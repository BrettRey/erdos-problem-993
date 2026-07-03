#!/usr/bin/env python3
"""Beam-search hard-component products for ratio-profile stress."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.audit_forest_products import default_components, poly_mul  # noqa: E402
from scripts.audit_ratio_profiles import ratio_profile  # noqa: E402


def ratio_key(row: dict[str, Any], key: str) -> float:
    value = row[key]
    return float(value) if value is not None else -1.0


def component_counts(labels: list[str], counts: tuple[int, ...]) -> dict[str, int]:
    return {label: count for label, count in zip(labels, counts) if count}


def compact(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "component_counts",
        "component_count",
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


def make_row(
    *,
    labels: list[str],
    counts: tuple[int, ...],
    poly: list[int],
    n: int,
    d_leaf_le1: bool,
    unsafe_threshold: float,
) -> dict[str, Any]:
    return {
        "component_counts": component_counts(labels, counts),
        "component_count": sum(counts),
        "d_leaf_le1": d_leaf_le1,
        **ratio_profile(poly, n, unsafe_threshold=unsafe_threshold),
    }


def choose_frontier(states: list[dict[str, Any]], beam_width: int) -> list[dict[str, Any]]:
    selected: dict[tuple[int, ...], dict[str, Any]] = {}
    objective_keys = [
        "max_tail_ratio",
        "max_lc_bump_right_ratio",
        "lc_bump_to_tail_ratio",
        "mean_ratio_n_over_3",
    ]
    for key in objective_keys:
        for state in sorted(states, key=lambda item: ratio_key(item["row"], key), reverse=True)[:beam_width]:
            selected[state["counts"]] = state
    for state in states:
        row = state["row"]
        if not row["unimodal"] or row["unsafe_lc_bump_right_ratio"]:
            selected[state["counts"]] = state
    return list(selected.values())


def search_rows(
    components: list[dict[str, Any]],
    *,
    max_factors: int,
    beam_width: int,
    unsafe_threshold: float,
) -> tuple[list[dict[str, Any]], dict[int, int]]:
    labels = [component["label"] for component in components]
    zero_counts = (0,) * len(components)
    frontier = [
        {
            "counts": zero_counts,
            "last_index": 0,
            "poly": [1],
            "n": 0,
            "d_leaf_le1": True,
        }
    ]
    rows: list[dict[str, Any]] = []
    depth_counts: dict[int, int] = {}

    for depth in range(1, max_factors + 1):
        candidates: dict[tuple[int, ...], dict[str, Any]] = {}
        for state in frontier:
            for index in range(state["last_index"], len(components)):
                component = components[index]
                counts = list(state["counts"])
                counts[index] += 1
                counts_tuple = tuple(counts)
                if counts_tuple in candidates:
                    continue
                poly = poly_mul(state["poly"], component["poly"])
                n = state["n"] + component["n"]
                d_leaf_le1 = state["d_leaf_le1"] and component["d_leaf_le1"]
                row = make_row(
                    labels=labels,
                    counts=counts_tuple,
                    poly=poly,
                    n=n,
                    d_leaf_le1=d_leaf_le1,
                    unsafe_threshold=unsafe_threshold,
                )
                candidates[counts_tuple] = {
                    "counts": counts_tuple,
                    "last_index": index,
                    "poly": poly,
                    "n": n,
                    "d_leaf_le1": d_leaf_le1,
                    "row": row,
                }
        rows.extend(state["row"] for state in candidates.values())
        depth_counts[depth] = len(candidates)
        frontier = choose_frontier(list(candidates.values()), beam_width)

    return rows, depth_counts


def build_summary(
    rows: list[dict[str, Any]],
    components: list[dict[str, Any]],
    *,
    max_factors: int,
    beam_width: int,
    unsafe_threshold: float,
    depth_counts: dict[int, int],
    top: int,
) -> dict[str, Any]:
    dleaf_rows = [row for row in rows if row["d_leaf_le1"]]
    lc_rows = [row for row in rows if row["lc_defect_count"] > 0]
    return {
        "source": {
            "kind": "hard_component_product_beam_search",
            "component_labels": [component["label"] for component in components],
            "component_count": len(components),
            "max_factors": max_factors,
            "beam_width": beam_width,
            "unsafe_threshold": unsafe_threshold,
            "depth_candidates": depth_counts,
        },
        "processed": len(rows),
        "counts": {
            "non_unimodal": sum(1 for row in rows if not row["unimodal"]),
            "non_log_concave": sum(1 for row in rows if not row["log_concave"]),
            "unsafe_lc_bump_right_ratio": sum(1 for row in rows if row["unsafe_lc_bump_right_ratio"]),
            "d_leaf_le1": len(dleaf_rows),
            "d_leaf_le1_low_mode": sum(1 for row in dleaf_rows if row["low_mode_surplus"] >= 0),
            "d_leaf_le1_mean_lt_n_over_3": sum(1 for row in dleaf_rows if row["mean_gap_n_over_3"] > 0),
            "lc_rows": len(lc_rows),
            "lc_bump_right_below_tail_max": sum(1 for row in lc_rows if row["lc_bump_right_below_tail_max"]),
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
            compact(row) for row in sorted(rows, key=lambda row: ratio_key(row, "max_tail_ratio"), reverse=True)[:top]
        ],
        "top_by_lc_bump_right_ratio": [
            compact(row)
            for row in sorted(rows, key=lambda row: ratio_key(row, "max_lc_bump_right_ratio"), reverse=True)[
                :top
            ]
        ],
        "top_by_lc_bump_to_tail_ratio": [
            compact(row)
            for row in sorted(rows, key=lambda row: ratio_key(row, "lc_bump_to_tail_ratio"), reverse=True)[:top]
        ],
        "unsafe_lc_bump_rows": [compact(row) for row in rows if row["unsafe_lc_bump_right_ratio"]][:top],
        "non_unimodal_rows": [compact(row) for row in rows if not row["unimodal"]][:top],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-factors", type=int, default=8, help="Maximum product size")
    parser.add_argument("--beam-width", type=int, default=80, help="States retained per objective")
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
        default=Path("results/product_ratio_beam_search_2026-07-03.json"),
        help="Output JSON path",
    )
    args = parser.parse_args()
    if args.max_factors < 1:
        parser.error("--max-factors must be positive")
    if args.beam_width < 1:
        parser.error("--beam-width must be positive")
    if not 0.0 < args.unsafe_threshold < 1.0:
        parser.error("--unsafe-threshold must be between 0 and 1")

    components = default_components()
    rows, depth_counts = search_rows(
        components,
        max_factors=args.max_factors,
        beam_width=args.beam_width,
        unsafe_threshold=args.unsafe_threshold,
    )
    summary = build_summary(
        rows,
        components,
        max_factors=args.max_factors,
        beam_width=args.beam_width,
        unsafe_threshold=args.unsafe_threshold,
        depth_counts=depth_counts,
        top=args.top,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Beam-searched {summary['processed']} products; "
        f"non-LC={summary['counts']['non_log_concave']}, "
        f"non-unimodal={summary['counts']['non_unimodal']}, "
        f"unsafe-LC-bumps={summary['counts']['unsafe_lc_bump_right_ratio']}, "
        f"wrote {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
