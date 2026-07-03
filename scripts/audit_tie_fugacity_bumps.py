#!/usr/bin/env python3
"""Audit tie-fugacity margins at post-descent ratio bumps."""

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
from scripts.analyze_prufer_corpus import lc_defects, mode_interval  # noqa: E402
from scripts.audit_forest_products import default_components, poly_mul  # noqa: E402
from scripts.audit_ratio_profiles import build_base_components, ratio_profile  # noqa: E402
from scripts.search_product_ratio_stress import choose_frontier, make_row  # noqa: E402


def first_descent(poly: list[int]) -> int | None:
    for k in range(1, len(poly)):
        if poly[k] < poly[k - 1]:
            return k
    return None


def log_ratio(num: int, den: int) -> float:
    if num == 0:
        return float("-inf")
    return math.log(num) - math.log(den)


def ratio_to_float(log_value: float) -> float:
    if log_value == float("-inf"):
        return 0.0
    if log_value > 709:
        return float("inf")
    if log_value < -745:
        return 0.0
    return math.exp(log_value)


def mean_at_log_lambda(poly: list[int], log_lambda: float) -> float:
    log_weights = [math.log(value) + k * log_lambda for k, value in enumerate(poly)]
    offset = max(log_weights)
    total = 0.0
    weighted = 0.0
    for k, log_weight in enumerate(log_weights):
        weight = math.exp(log_weight - offset)
        total += weight
        weighted += k * weight
    return weighted / total


def tie_event(poly: list[int], *, k: int, event_type: str) -> dict[str, Any]:
    """Return the tie-fugacity margin for the transition k -> k+1."""
    right_log_ratio = log_ratio(poly[k + 1], poly[k])
    log_lambda = -right_log_ratio
    mu_lambda = mean_at_log_lambda(poly, log_lambda)
    return {
        "event_type": event_type,
        "k": k,
        "right_ratio": ratio_to_float(right_log_ratio),
        "log_lambda": log_lambda,
        "lambda": ratio_to_float(log_lambda),
        "mu_lambda": mu_lambda,
        "mu_minus_k": mu_lambda - k,
    }


def bump_events(poly: list[int]) -> list[dict[str, Any]]:
    descent = first_descent(poly)
    if descent is None:
        return []
    events = []
    for defect in lc_defects(poly):
        k = defect["k"]
        if k >= descent and k + 1 < len(poly):
            events.append(
                {
                    **tie_event(poly, k=k, event_type="post_descent_lc_bump"),
                    "lc_ratio": defect["ratio"],
                }
            )
    return events


def crossing_events(poly: list[int]) -> list[dict[str, Any]]:
    """Strict post-descent upward coefficient transitions."""
    descent = first_descent(poly)
    if descent is None:
        return []
    events = []
    for k in range(descent, len(poly) - 1):
        if poly[k + 1] > poly[k]:
            event = tie_event(poly, k=k, event_type="post_descent_upward_transition")
            if k >= 1:
                event["lc_bump_at_k"] = poly[k - 1] * poly[k + 1] > poly[k] * poly[k]
            else:
                event["lc_bump_at_k"] = None
            events.append(event)
    return events


def entry_metadata(
    *,
    source_kind: str,
    labels_key: str,
    labels: list[str],
    n: int,
    poly: list[int],
    d_leaf_le1: bool,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    mode_first, mode_last = mode_interval(poly)
    metadata = {
        "source_kind": source_kind,
        labels_key: labels,
        "n": n,
        "alpha": len(poly) - 1,
        "mode_first": mode_first,
        "mode_last": mode_last,
        "first_descent": first_descent(poly),
        "d_leaf_le1": d_leaf_le1,
        "unimodal": is_unimodal(poly),
        "log_concave": is_log_concave(poly),
    }
    if extra:
        metadata.update(extra)
    return metadata


def product_power_entries(
    *,
    max_power: int,
    mixed_max_power: int,
) -> list[dict[str, Any]]:
    entries = []
    for base in build_base_components():
        poly = [1]
        n = 0
        base_max_power = mixed_max_power if len(base["labels"]) > 1 else max_power
        for power in range(1, base_max_power + 1):
            poly = poly_mul(poly, base["poly"])
            n += base["n"]
            entries.append(
                {
                    "poly": poly[:],
                    "metadata": entry_metadata(
                        source_kind="product_power",
                        labels_key="base_labels",
                        labels=base["labels"],
                        n=n,
                        poly=poly,
                        d_leaf_le1=base["d_leaf_le1"],
                        extra={"power": power},
                    ),
                }
            )
    return entries


def beam_entries(
    *,
    max_factors: int,
    beam_width: int,
    unsafe_threshold: float,
) -> list[dict[str, Any]]:
    components = default_components()
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
    entries = []
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
                component_counts = row["component_counts"]
                candidates[counts_tuple] = {
                    "counts": counts_tuple,
                    "last_index": index,
                    "poly": poly,
                    "n": n,
                    "d_leaf_le1": d_leaf_le1,
                    "row": row,
                    "metadata": entry_metadata(
                        source_kind="product_beam",
                        labels_key="component_counts",
                        labels=[],
                        n=n,
                        poly=poly,
                        d_leaf_le1=d_leaf_le1,
                        extra={
                            "component_counts": component_counts,
                            "component_count": sum(counts_tuple),
                        },
                    ),
                }
        entries.extend(
            {"poly": state["poly"], "metadata": state["metadata"]} for state in candidates.values()
        )
        frontier = choose_frontier(list(candidates.values()), beam_width)
    return entries


def compact_event(entry: dict[str, Any], event: dict[str, Any]) -> dict[str, Any]:
    metadata = entry["metadata"]
    keys = [
        "source_kind",
        "base_labels",
        "component_counts",
        "component_count",
        "power",
        "n",
        "alpha",
        "mode_first",
        "first_descent",
        "d_leaf_le1",
        "unimodal",
        "log_concave",
    ]
    out = {key: metadata[key] for key in keys if key in metadata}
    out.update(event)
    return out


def summarize_numeric(events: list[dict[str, Any]], key: str) -> dict[str, float | None]:
    values = [event[key] for event in events if event.get(key) is not None]
    if not values:
        return {"min": None, "max": None, "mean": None, "median": None}
    values = sorted(values)
    return {
        "min": values[0],
        "max": values[-1],
        "mean": sum(values) / len(values),
        "median": values[len(values) // 2],
    }


def build_summary(
    entries: list[dict[str, Any]],
    *,
    max_power: int,
    mixed_max_power: int,
    beam_max_factors: int,
    beam_width: int,
    unsafe_threshold: float,
    top: int,
) -> dict[str, Any]:
    bump_rows = []
    crossing_rows = []
    for entry in entries:
        poly = entry["poly"]
        bump_rows.extend(compact_event(entry, event) for event in bump_events(poly))
        crossing_rows.extend(compact_event(entry, event) for event in crossing_events(poly))

    dleaf_entries = [entry for entry in entries if entry["metadata"]["d_leaf_le1"]]
    return {
        "source": {
            "kind": "tie_fugacity_bump_audit",
            "product_power": {
                "max_power": max_power,
                "mixed_max_power": mixed_max_power,
            },
            "product_beam": {
                "max_factors": beam_max_factors,
                "beam_width": beam_width,
            },
            "unsafe_threshold": unsafe_threshold,
        },
        "processed": len(entries),
        "counts": {
            "products": len(entries),
            "product_power_rows": sum(
                1 for entry in entries if entry["metadata"]["source_kind"] == "product_power"
            ),
            "product_beam_rows": sum(
                1 for entry in entries if entry["metadata"]["source_kind"] == "product_beam"
            ),
            "d_leaf_le1_products": len(dleaf_entries),
            "non_unimodal_products": sum(1 for entry in entries if not entry["metadata"]["unimodal"]),
            "non_log_concave_products": sum(1 for entry in entries if not entry["metadata"]["log_concave"]),
            "post_descent_lc_bump_events": len(bump_rows),
            "post_descent_upward_transition_events": len(crossing_rows),
            "negative_lc_bump_tie_margins": sum(1 for event in bump_rows if event["mu_minus_k"] < 0),
            "negative_crossing_tie_margins": sum(
                1 for event in crossing_rows if event["mu_minus_k"] < 0
            ),
        },
        "numeric": {
            "lc_bump_mu_minus_k": summarize_numeric(bump_rows, "mu_minus_k"),
            "lc_bump_right_ratio": summarize_numeric(bump_rows, "right_ratio"),
            "lc_bump_log_lambda": summarize_numeric(bump_rows, "log_lambda"),
            "crossing_mu_minus_k": summarize_numeric(crossing_rows, "mu_minus_k"),
            "crossing_right_ratio": summarize_numeric(crossing_rows, "right_ratio"),
        },
        "top_min_lc_bump_margin": sorted(bump_rows, key=lambda event: event["mu_minus_k"])[:top],
        "top_max_lc_bump_right_ratio": sorted(
            bump_rows, key=lambda event: event["right_ratio"], reverse=True
        )[:top],
        "post_descent_upward_transition_rows": crossing_rows[:top],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-power", type=int, default=20, help="Maximum single-base product power")
    parser.add_argument(
        "--mixed-max-power",
        type=int,
        default=4,
        help="Maximum repeated product power for mixed bases",
    )
    parser.add_argument("--beam-max-factors", type=int, default=8, help="Maximum beam product size")
    parser.add_argument("--beam-width", type=int, default=80, help="Beam states retained per objective")
    parser.add_argument(
        "--unsafe-threshold",
        type=float,
        default=0.95,
        help="Shared threshold used by imported beam frontier scoring",
    )
    parser.add_argument("--top", type=int, default=10, help="Rows retained per table")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/tie_fugacity_bump_audit_2026-07-03.json"),
        help="Output JSON path",
    )
    args = parser.parse_args()
    if args.max_power < 1:
        parser.error("--max-power must be positive")
    if args.mixed_max_power < 1:
        parser.error("--mixed-max-power must be positive")
    if args.beam_max_factors < 1:
        parser.error("--beam-max-factors must be positive")
    if args.beam_width < 1:
        parser.error("--beam-width must be positive")
    if not 0.0 < args.unsafe_threshold < 1.0:
        parser.error("--unsafe-threshold must be between 0 and 1")

    entries = product_power_entries(max_power=args.max_power, mixed_max_power=args.mixed_max_power)
    entries.extend(
        beam_entries(
            max_factors=args.beam_max_factors,
            beam_width=args.beam_width,
            unsafe_threshold=args.unsafe_threshold,
        )
    )
    summary = build_summary(
        entries,
        max_power=args.max_power,
        mixed_max_power=args.mixed_max_power,
        beam_max_factors=args.beam_max_factors,
        beam_width=args.beam_width,
        unsafe_threshold=args.unsafe_threshold,
        top=args.top,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Audited {summary['processed']} products; "
        f"LC-bump events={summary['counts']['post_descent_lc_bump_events']}, "
        f"upward-transition events={summary['counts']['post_descent_upward_transition_events']}, "
        f"negative bump margins={summary['counts']['negative_lc_bump_tie_margins']}, "
        f"wrote {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
