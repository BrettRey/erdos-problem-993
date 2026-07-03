#!/usr/bin/env python3
"""Audit forest/product stresses built from known hard tree polynomials."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from graph6 import parse_graph6  # noqa: E402
from indpoly import (  # noqa: E402
    independence_poly,
    is_log_concave,
    is_unimodal,
    log_concavity_ratio,
    near_miss_ratio,
)
from scripts.analyze_prufer_corpus import (  # noqa: E402
    lc_defects,
    make_bautista_ramos_tree,
    make_galvin_tree,
    make_li_tree,
    mean_independent_set_size,
    mode_interval,
    structural_metrics,
)


def poly_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj:
                out[i + j] += ai * bj
    return out


def graph6_adj(graph6: str) -> list[list[int]]:
    return parse_graph6(graph6.encode("ascii"))[1]


def component(label: str, adj: list[list[int]]) -> dict[str, Any]:
    poly = independence_poly(len(adj), adj)
    return {
        "label": label,
        "n": len(adj),
        "alpha": len(poly) - 1,
        "poly": poly,
        **structural_metrics(adj),
    }


def default_components() -> list[dict[str, Any]]:
    ramos = json.loads(Path("results/ramos_sun_60_epoch11_summary.json").read_text())
    n28 = json.loads(Path("results/analysis_n28_modal_lc_nm.json").read_text())
    return [
        component("ramos_epoch11_top_nm_rank6829", graph6_adj(ramos["top_by_near_miss_ratio"][0]["graph6"])),
        component("ramos_epoch11_top_lc_rank3309", graph6_adj(ramos["top_by_lc_ratio"][0]["graph6"])),
        component("galvin_T_21_11", make_galvin_tree(21, 11)),
        component("galvin_T_1_50", make_galvin_tree(1, 50)),
        component("bautista_ramos_TG_8_8", make_bautista_ramos_tree(8, 8)),
        component("li_T_3_50_50", make_li_tree(50, 50, starred=False)),
        component("li_Tstar_3_50_50", make_li_tree(50, 50, starred=True)),
        component("n28_top_near_miss", graph6_adj(n28["top_near_misses"][0]["graph6"])),
        component("n28_worst_lc", graph6_adj(n28["top_lc_failures"][0]["graph6"])),
    ]


def sequence_metrics(poly: list[int], n: int) -> dict[str, Any]:
    mode_first, mode_last = mode_interval(poly)
    mu = mean_independent_set_size(poly)
    nm_ratio, nm_pos = near_miss_ratio(poly)
    lc_ratio, lc_pos = log_concavity_ratio(poly)
    defects = lc_defects(poly)
    defect_distances_from_mode = [defect["k"] - mode_first for defect in defects]
    defect_distances_from_near_miss = [
        defect["k"] - nm_pos for defect in defects if nm_pos != -1
    ]
    return {
        "n": n,
        "alpha": len(poly) - 1,
        "mode_first": mode_first,
        "mode_last": mode_last,
        "mu": mu,
        "ceil_mu_minus_mode": math.ceil(mu) - mode_first,
        "low_mode_surplus": n // 3 + 1 - mode_first,
        "mean_gap_n_over_3": n / 3 - mu,
        "mean_ratio_n_over_3": mu / (n / 3) if n else None,
        "unimodal": is_unimodal(poly),
        "log_concave": is_log_concave(poly),
        "near_miss_ratio": nm_ratio,
        "near_miss_pos": nm_pos,
        "near_miss_reserve": 1.0 - nm_ratio if nm_pos != -1 else None,
        "lc_ratio": lc_ratio,
        "lc_pos": lc_pos,
        "lc_defect_count": len(defects),
        "lc_defect_min_distance_from_mode": min(defect_distances_from_mode)
        if defect_distances_from_mode
        else None,
        "lc_defect_min_distance_from_near_miss": min(defect_distances_from_near_miss)
        if defect_distances_from_near_miss
        else None,
    }


def compact_row(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "labels",
        "n",
        "alpha",
        "mode_first",
        "low_mode_surplus",
        "mu",
        "mean_gap_n_over_3",
        "mean_ratio_n_over_3",
        "ceil_mu_minus_mode",
        "unimodal",
        "log_concave",
        "near_miss_ratio",
        "near_miss_reserve",
        "near_miss_pos",
        "lc_ratio",
        "lc_pos",
        "lc_defect_count",
        "lc_defect_min_distance_from_mode",
        "lc_defect_min_distance_from_near_miss",
        "d_leaf_le1",
    ]
    return {key: row[key] for key in keys}


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


def build_rows(components: list[dict[str, Any]], max_factors: int) -> list[dict[str, Any]]:
    rows = []
    rank = 1
    for size in range(1, max_factors + 1):
        for combo in itertools.combinations_with_replacement(components, size):
            poly = [1]
            for item in combo:
                poly = poly_mul(poly, item["poly"])
            n = sum(item["n"] for item in combo)
            row = {
                "rank": rank,
                "labels": [item["label"] for item in combo],
                "component_count": size,
                "d_leaf_le1": all(item["d_leaf_le1"] for item in combo),
                **sequence_metrics(poly, n),
            }
            rows.append(row)
            rank += 1
    return rows


def build_summary(
    rows: list[dict[str, Any]],
    components: list[dict[str, Any]],
    *,
    max_factors: int,
    top: int,
) -> dict[str, Any]:
    dleaf_rows = [row for row in rows if row["d_leaf_le1"]]
    return {
        "source": {
            "component_labels": [item["label"] for item in components],
            "component_count": len(components),
            "max_factors": max_factors,
        },
        "processed": len(rows),
        "counts": {
            "non_unimodal": sum(1 for row in rows if not row["unimodal"]),
            "non_log_concave": sum(1 for row in rows if not row["log_concave"]),
            "low_mode": sum(1 for row in rows if row["low_mode_surplus"] >= 0),
            "mean_lt_n_over_3": sum(1 for row in rows if row["mean_gap_n_over_3"] > 0),
            "mode_le_ceil_mu": sum(1 for row in rows if row["ceil_mu_minus_mode"] >= 0),
            "d_leaf_le1": sum(1 for row in rows if row["d_leaf_le1"]),
            "d_leaf_le1_non_unimodal": sum(1 for row in dleaf_rows if not row["unimodal"]),
            "d_leaf_le1_non_log_concave": sum(1 for row in dleaf_rows if not row["log_concave"]),
            "d_leaf_le1_low_mode": sum(1 for row in dleaf_rows if row["low_mode_surplus"] >= 0),
            "d_leaf_le1_mean_lt_n_over_3": sum(1 for row in dleaf_rows if row["mean_gap_n_over_3"] > 0),
            "d_leaf_le1_mode_le_ceil_mu": sum(1 for row in dleaf_rows if row["ceil_mu_minus_mode"] >= 0),
            "lc_defects_all_after_mode": sum(
                1
                for row in rows
                if row["lc_defect_count"] > 0 and row["lc_defect_min_distance_from_mode"] > 0
            ),
            "lc_defects_all_after_near_miss": sum(
                1
                for row in rows
                if row["lc_defect_count"] > 0
                and row["lc_defect_min_distance_from_near_miss"] is not None
                and row["lc_defect_min_distance_from_near_miss"] > 0
            ),
        },
        "numeric": {
            key: summarize_numeric(rows, key)
            for key in [
                "low_mode_surplus",
                "mean_gap_n_over_3",
                "mean_ratio_n_over_3",
                "near_miss_ratio",
                "near_miss_reserve",
                "lc_ratio",
                "lc_defect_min_distance_from_mode",
                "lc_defect_min_distance_from_near_miss",
            ]
        },
        "top_by_near_miss_ratio": [
            compact_row(row) for row in sorted(rows, key=lambda r: r["near_miss_ratio"], reverse=True)[:top]
        ],
        "top_by_lc_ratio": [
            compact_row(row) for row in sorted(rows, key=lambda r: r["lc_ratio"], reverse=True)[:top]
        ],
        "bottom_by_mean_gap_n_over_3": [
            compact_row(row) for row in sorted(rows, key=lambda r: r["mean_gap_n_over_3"])[:top]
        ],
        "bottom_by_low_mode_surplus": [
            compact_row(row) for row in sorted(rows, key=lambda r: r["low_mode_surplus"])[:top]
        ],
        "non_unimodal_rows": [compact_row(row) for row in rows if not row["unimodal"]][:top],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-factors", type=int, default=3, help="Maximum forest components per product")
    parser.add_argument("--top", type=int, default=10, help="Rows retained per top/bottom table")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/forest_product_stress_audit_2026-07-03.json"),
        help="Output JSON path",
    )
    args = parser.parse_args()
    if args.max_factors < 1:
        parser.error("--max-factors must be positive")

    components = default_components()
    rows = build_rows(components, args.max_factors)
    summary = build_summary(rows, components, max_factors=args.max_factors, top=args.top)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Analyzed {summary['processed']} products; "
        f"non-LC={summary['counts']['non_log_concave']}, "
        f"non-unimodal={summary['counts']['non_unimodal']}, "
        f"wrote {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
