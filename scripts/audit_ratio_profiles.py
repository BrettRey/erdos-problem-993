#!/usr/bin/env python3
"""Audit ratio-profile separation under repeated hard-component products."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from indpoly import is_log_concave, is_unimodal  # noqa: E402
from scripts.analyze_prufer_corpus import (  # noqa: E402
    lc_defects,
    mean_independent_set_size,
    mode_interval,
)
from scripts.audit_forest_products import default_components, poly_mul  # noqa: E402


def log_ratio(num: int, den: int) -> float:
    return math.log(num) - math.log(den) if num else float("-inf")


def ratio_to_float(log_value: float | None) -> float | None:
    if log_value is None:
        return None
    if log_value == float("-inf"):
        return 0.0
    return math.exp(log_value)


def first_descent(poly: list[int]) -> int | None:
    for k in range(1, len(poly)):
        if poly[k] < poly[k - 1]:
            return k
    return None


def max_transition_ratio(poly: list[int], start: int) -> tuple[int | None, float | None]:
    """Return (transition_j, log_ratio) for max a[j+1]/a[j]."""
    best_j: int | None = None
    best_log: float | None = None
    for j in range(start, len(poly) - 1):
        den = poly[j]
        num = poly[j + 1]
        if den == 0:
            continue
        current_log = log_ratio(num, den)
        if best_log is None or current_log > best_log:
            best_j = j
            best_log = current_log
    return best_j, best_log


def ratio_profile(poly: list[int], n: int, *, unsafe_threshold: float) -> dict[str, Any]:
    mode_first, mode_last = mode_interval(poly)
    descent = first_descent(poly)
    if descent is None:
        tail_j, tail_log = None, None
    else:
        tail_j, tail_log = max_transition_ratio(poly, descent)

    defects = lc_defects(poly)
    best_lc_k: int | None = None
    best_lc_log: float | None = None
    for defect in defects:
        k = defect["k"]
        if k + 1 >= len(poly) or poly[k] == 0:
            continue
        current_log = log_ratio(poly[k + 1], poly[k])
        if best_lc_log is None or current_log > best_lc_log:
            best_lc_k = k
            best_lc_log = current_log

    max_tail_ratio = ratio_to_float(tail_log)
    max_lc_bump_right_ratio = ratio_to_float(best_lc_log)
    if tail_log is not None and best_lc_log is not None:
        lc_right_below_tail = tail_log > best_lc_log
        separation_ratio = ratio_to_float(best_lc_log - tail_log)
    else:
        lc_right_below_tail = None
        separation_ratio = None

    mu = mean_independent_set_size(poly)
    return {
        "n": n,
        "alpha": len(poly) - 1,
        "mode_first": mode_first,
        "mode_last": mode_last,
        "first_descent": descent,
        "mu": mu,
        "low_mode_surplus": n // 3 + 1 - mode_first,
        "mean_gap_n_over_3": n / 3 - mu,
        "mean_ratio_n_over_3": mu / (n / 3) if n else None,
        "unimodal": is_unimodal(poly),
        "log_concave": is_log_concave(poly),
        "lc_defect_count": len(defects),
        "max_tail_transition": tail_j,
        "max_tail_ratio": max_tail_ratio,
        "max_tail_reserve": 1.0 - max_tail_ratio if max_tail_ratio is not None else None,
        "max_tail_distance_from_mode": tail_j - mode_first if tail_j is not None else None,
        "max_lc_bump_k": best_lc_k,
        "max_lc_bump_right_ratio": max_lc_bump_right_ratio,
        "max_lc_bump_right_reserve": 1.0 - max_lc_bump_right_ratio
        if max_lc_bump_right_ratio is not None
        else None,
        "max_lc_bump_distance_from_mode": best_lc_k - mode_first if best_lc_k is not None else None,
        "max_lc_bump_distance_from_tail_max": best_lc_k - tail_j
        if best_lc_k is not None and tail_j is not None
        else None,
        "lc_bump_right_below_tail_max": lc_right_below_tail,
        "lc_bump_to_tail_ratio": separation_ratio,
        "unsafe_lc_bump_right_ratio": (
            best_lc_log >= math.log(unsafe_threshold)
            if best_lc_log is not None
            else False
        ),
    }


def build_base_components() -> list[dict[str, Any]]:
    components = {item["label"]: item for item in default_components()}

    def base(labels: list[str]) -> dict[str, Any]:
        poly = [1]
        n = 0
        d_leaf_le1 = True
        for label in labels:
            item = components[label]
            poly = poly_mul(poly, item["poly"])
            n += item["n"]
            d_leaf_le1 = d_leaf_le1 and item["d_leaf_le1"]
        return {"labels": labels, "poly": poly, "n": n, "d_leaf_le1": d_leaf_le1}

    labels = [
        "galvin_T_21_11",
        "bautista_ramos_TG_8_8",
        "li_T_3_50_50",
        "li_Tstar_3_50_50",
        "ramos_epoch11_top_nm_rank6829",
        "ramos_epoch11_top_lc_rank3309",
        "n28_top_near_miss",
        "n28_worst_lc",
    ]
    bases = [base([label]) for label in labels]
    bases.extend(
        [
            base(["galvin_T_21_11", "bautista_ramos_TG_8_8"]),
            base(["galvin_T_21_11", "li_Tstar_3_50_50"]),
            base(["bautista_ramos_TG_8_8", "li_Tstar_3_50_50"]),
            base(["galvin_T_21_11", "bautista_ramos_TG_8_8", "li_Tstar_3_50_50"]),
            base(["n28_top_near_miss", "galvin_T_21_11"]),
            base(["n28_top_near_miss", "bautista_ramos_TG_8_8"]),
        ]
    )
    return bases


def product_power_rows(
    bases: list[dict[str, Any]],
    *,
    max_power: int,
    mixed_max_power: int,
    unsafe_threshold: float,
) -> list[dict[str, Any]]:
    rows = []
    rank = 1
    for base in bases:
        poly = [1]
        n = 0
        base_max_power = mixed_max_power if len(base["labels"]) > 1 else max_power
        for power in range(1, base_max_power + 1):
            poly = poly_mul(poly, base["poly"])
            n += base["n"]
            row = {
                "rank": rank,
                "base_labels": base["labels"],
                "power": power,
                "d_leaf_le1": base["d_leaf_le1"],
                **ratio_profile(poly, n, unsafe_threshold=unsafe_threshold),
            }
            rows.append(row)
            rank += 1
    return rows


def summarize_numeric(rows: list[dict[str, Any]], key: str) -> dict[str, float | int | None]:
    values = [row[key] for row in rows if row[key] is not None]
    if not values:
        return {"min": None, "max": None, "mean": None, "median": None}
    values = sorted(values)
    return {
        "min": values[0],
        "max": values[-1],
        "mean": sum(values) / len(values),
        "median": values[len(values) // 2],
    }


def compact(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "base_labels",
        "power",
        "n",
        "alpha",
        "mode_first",
        "low_mode_surplus",
        "mean_gap_n_over_3",
        "mean_ratio_n_over_3",
        "d_leaf_le1",
        "unimodal",
        "log_concave",
        "lc_defect_count",
        "max_tail_transition",
        "max_tail_ratio",
        "max_tail_reserve",
        "max_tail_distance_from_mode",
        "max_lc_bump_k",
        "max_lc_bump_right_ratio",
        "max_lc_bump_right_reserve",
        "max_lc_bump_distance_from_mode",
        "max_lc_bump_distance_from_tail_max",
        "lc_bump_right_below_tail_max",
        "lc_bump_to_tail_ratio",
        "unsafe_lc_bump_right_ratio",
    ]
    return {key: row[key] for key in keys}


def build_summary(
    rows: list[dict[str, Any]],
    *,
    max_power: int,
    mixed_max_power: int,
    unsafe_threshold: float,
    top: int,
) -> dict[str, Any]:
    lc_rows = [row for row in rows if row["lc_defect_count"] > 0]
    dleaf_rows = [row for row in rows if row["d_leaf_le1"]]
    return {
        "source": {
            "kind": "hard_component_product_powers",
            "max_power": max_power,
            "mixed_max_power": mixed_max_power,
            "unsafe_threshold": unsafe_threshold,
        },
        "processed": len(rows),
        "counts": {
            "non_unimodal": sum(1 for row in rows if not row["unimodal"]),
            "non_log_concave": sum(1 for row in rows if not row["log_concave"]),
            "d_leaf_le1": len(dleaf_rows),
            "d_leaf_le1_low_mode": sum(1 for row in dleaf_rows if row["low_mode_surplus"] >= 0),
            "d_leaf_le1_mean_lt_n_over_3": sum(1 for row in dleaf_rows if row["mean_gap_n_over_3"] > 0),
            "lc_rows": len(lc_rows),
            "lc_bump_right_below_tail_max": sum(1 for row in lc_rows if row["lc_bump_right_below_tail_max"]),
            "unsafe_lc_bump_right_ratio": sum(1 for row in rows if row["unsafe_lc_bump_right_ratio"]),
        },
        "numeric": {
            key: summarize_numeric(rows, key)
            for key in [
                "max_tail_ratio",
                "max_tail_reserve",
                "max_lc_bump_right_ratio",
                "max_lc_bump_right_reserve",
                "lc_bump_to_tail_ratio",
                "max_lc_bump_distance_from_mode",
                "max_lc_bump_distance_from_tail_max",
                "low_mode_surplus",
                "mean_gap_n_over_3",
                "mean_ratio_n_over_3",
            ]
        },
        "top_by_max_tail_ratio": [
            compact(row) for row in sorted(rows, key=lambda row: row["max_tail_ratio"] or -1, reverse=True)[:top]
        ],
        "top_by_lc_bump_right_ratio": [
            compact(row)
            for row in sorted(rows, key=lambda row: row["max_lc_bump_right_ratio"] or -1, reverse=True)[:top]
        ],
        "top_by_lc_bump_to_tail_ratio": [
            compact(row) for row in sorted(rows, key=lambda row: row["lc_bump_to_tail_ratio"] or -1, reverse=True)[:top]
        ],
        "unsafe_lc_bump_rows": [compact(row) for row in rows if row["unsafe_lc_bump_right_ratio"]][:top],
        "non_unimodal_rows": [compact(row) for row in rows if not row["unimodal"]][:top],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-power", type=int, default=12, help="Maximum repeated product power")
    parser.add_argument(
        "--mixed-max-power",
        type=int,
        default=4,
        help="Maximum repeated product power for mixed hard-component bases",
    )
    parser.add_argument(
        "--unsafe-threshold",
        type=float,
        default=0.95,
        help="Flag LC bump right-ratios at or above this threshold",
    )
    parser.add_argument("--top", type=int, default=10, help="Rows retained per top table")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/ratio_profile_power_audit_2026-07-03.json"),
        help="Output JSON path",
    )
    args = parser.parse_args()
    if args.max_power < 1:
        parser.error("--max-power must be positive")
    if args.mixed_max_power < 1:
        parser.error("--mixed-max-power must be positive")
    if not 0.0 < args.unsafe_threshold < 1.0:
        parser.error("--unsafe-threshold must be between 0 and 1")

    rows = product_power_rows(
        build_base_components(),
        max_power=args.max_power,
        mixed_max_power=args.mixed_max_power,
        unsafe_threshold=args.unsafe_threshold,
    )
    summary = build_summary(
        rows,
        max_power=args.max_power,
        mixed_max_power=args.mixed_max_power,
        unsafe_threshold=args.unsafe_threshold,
        top=args.top,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Analyzed {summary['processed']} product powers; "
        f"non-LC={summary['counts']['non_log_concave']}, "
        f"non-unimodal={summary['counts']['non_unimodal']}, "
        f"unsafe-LC-bumps={summary['counts']['unsafe_lc_bump_right_ratio']}, "
        f"wrote {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
