#!/usr/bin/env python3
"""Analyze conditional local ratios for signed PB reserve rows."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from probe_signed_pb_reserve import Block, first_descent, grouped_pmf


DEFAULT_INPUTS = [
    Path("results/signed_pb_reserve_optimizer_2026-07-03.json"),
    Path("results/signed_pb_reserve_optimizer_balanced_f05_2026-07-03.json"),
    Path("results/signed_pb_reserve_optimizer_balanced_f10_2026-07-03.json"),
    Path("results/signed_pb_reserve_optimizer_balanced_f25_2026-07-03.json"),
    Path("results/signed_pb_reserve_probe_2026-07-03.json"),
]


def as_blocks(raw: list[list[float]]) -> list[Block]:
    return [(int(count), float(p)) for count, p in raw]


def trim_zero_tail(pmf: np.ndarray) -> np.ndarray:
    active = np.flatnonzero(pmf > 0.0)
    if len(active) == 0:
        return pmf
    return pmf[: int(active[-1]) + 1]


def ratio(seq: np.ndarray, index: int) -> float:
    if index < 0:
        return 0.0
    if index >= len(seq) - 1:
        return 0.0
    denominator = float(seq[index])
    if denominator <= 0.0:
        return math.nan
    return float(seq[index + 1] / denominator)


def reverse_ratio(seq: np.ndarray, index: int) -> float:
    """Return seq[index - 1] / seq[index], with zero outside support."""
    if index <= 0:
        return 0.0
    if index >= len(seq):
        return 0.0
    denominator = float(seq[index])
    if denominator <= 0.0:
        return math.nan
    return float(seq[index - 1] / denominator)


def conditional_y(
    x_pmf: np.ndarray,
    y_pmf: np.ndarray,
    z: int,
) -> tuple[np.ndarray, np.ndarray]:
    lo = max(0, -z)
    hi = min(len(y_pmf) - 1, len(x_pmf) - 1 - z)
    if lo > hi:
        return np.array([], dtype=np.int64), np.array([], dtype=np.float64)
    y_values = np.arange(lo, hi + 1, dtype=np.int64)
    weights = x_pmf[z + y_values] * y_pmf[y_values]
    total = float(weights.sum())
    if total <= 0.0:
        return y_values, weights
    return y_values, weights / total


def conditional_stats(values: np.ndarray, weights: np.ndarray) -> dict[str, float]:
    if len(values) == 0 or float(weights.sum()) <= 0.0:
        return {"mean": math.nan, "variance": math.nan}
    mean = float(np.dot(values, weights))
    variance = float(np.dot((values - mean) ** 2, weights))
    return {"mean": mean, "variance": variance}


def total_variation(
    values_a: np.ndarray,
    weights_a: np.ndarray,
    values_b: np.ndarray,
    weights_b: np.ndarray,
) -> float:
    keys = sorted(set(int(v) for v in values_a) | set(int(v) for v in values_b))
    map_a = {int(v): float(w) for v, w in zip(values_a, weights_a)}
    map_b = {int(v): float(w) for v, w in zip(values_b, weights_b)}
    return 0.5 * sum(abs(map_a.get(key, 0.0) - map_b.get(key, 0.0)) for key in keys)


def expectation(values: np.ndarray, weights: np.ndarray, fn) -> float:
    if len(values) == 0 or float(weights.sum()) <= 0.0:
        return math.nan
    total = 0.0
    mass = 0.0
    for value, weight_raw in zip(values, weights):
        weight = float(weight_raw)
        if weight <= 0.0:
            continue
        term = fn(int(value))
        if not math.isfinite(term):
            continue
        total += weight * term
        mass += weight
    if mass <= 0.0:
        return math.nan
    return total / mass


def lower_x_boundary_term(x_pmf: np.ndarray, y_pmf: np.ndarray, z: int) -> float:
    """Return the c_{z+1} term with X=0 missing from the X-side ratio core."""
    y_index = -z - 1
    if 0 <= y_index < len(y_pmf):
        return float(x_pmf[0] * y_pmf[y_index])
    return 0.0


def upper_y_boundary_term(x_pmf: np.ndarray, y_pmf: np.ndarray, z: int) -> float:
    """Return the c_{z+1} term with Y=max(Y) missing from the Y-side ratio core."""
    y_index = len(y_pmf) - 1
    x_index = z + 1 + y_index
    if 0 <= x_index < len(x_pmf):
        return float(x_pmf[x_index] * y_pmf[y_index])
    return 0.0


def analyze_row(row: dict[str, Any], source: str) -> dict[str, Any] | None:
    x_blocks = as_blocks(row["x_blocks"])
    y_blocks = as_blocks(row["y_blocks"])
    x_pmf = trim_zero_tail(grouped_pmf(x_blocks))
    y_pmf = trim_zero_tail(grouped_pmf(y_blocks))
    signed = np.convolve(x_pmf, y_pmf[::-1])
    support_start = -(len(y_pmf) - 1)
    descent_index = first_descent(signed)
    if descent_index is None or descent_index >= len(signed) - 1:
        return None

    z = support_start + descent_index
    c_prev = float(signed[descent_index - 1])
    c_z = float(signed[descent_index])
    previous_ratio = float(signed[descent_index] / signed[descent_index - 1])
    pressure = float(signed[descent_index + 1] / signed[descent_index])
    effective_ratio = pressure / previous_ratio
    effective_drop = 1.0 - effective_ratio
    y_values, weights = conditional_y(x_pmf, y_pmf, z)
    y_prev_values, prev_weights = conditional_y(x_pmf, y_pmf, z - 1)
    y_stats = conditional_stats(y_values, weights)
    y_prev_stats = conditional_stats(y_prev_values, prev_weights)

    x_ratio_core = expectation(y_values, weights, lambda y: ratio(x_pmf, z + y))
    x_boundary = lower_x_boundary_term(x_pmf, y_pmf, z) / c_z
    x_ratio_identity = x_ratio_core + x_boundary
    x_previous_core = expectation(
        y_prev_values,
        prev_weights,
        lambda y: ratio(x_pmf, z - 1 + y),
    )
    x_previous_boundary = lower_x_boundary_term(x_pmf, y_pmf, z - 1) / c_prev
    x_previous_identity = x_previous_core + x_previous_boundary
    y_ratio_core = expectation(
        y_values,
        weights,
        lambda y: reverse_ratio(y_pmf, y),
    )
    y_boundary = upper_y_boundary_term(x_pmf, y_pmf, z) / c_z
    y_ratio_identity = y_ratio_core + y_boundary
    identity_errors = [
        abs(pressure - x_ratio_identity),
        abs(previous_ratio - x_previous_identity),
        abs(pressure - y_ratio_identity),
    ]
    finite_identity_errors = [error for error in identity_errors if math.isfinite(error)]

    x_index_mean = z + y_stats["mean"]
    x_prev_index_mean = z - 1 + y_prev_stats["mean"]
    x_var = float(row["x_variance"])
    y_var = float(row["y_variance"])
    variance = float(row["variance"])
    side_fraction = min(x_var, y_var) / variance if variance > 0.0 else math.nan

    inv_index = expectation(
        y_values,
        weights,
        lambda y: 1.0 / float(z + y + 1) if z + y + 1 > 0 else 0.0,
    )
    conditional_newton_factor = expectation(
        y_values,
        weights,
        lambda y: float(z + y) / float(z + y + 1) if z + y + 1 > 0 else 0.0,
    )

    return {
        "source": source,
        "x_n": row["x_n"],
        "y_n": row["y_n"],
        "x_variance": x_var,
        "y_variance": y_var,
        "variance": variance,
        "side_variance_fraction": side_fraction,
        "signed_mean": row["signed_mean"],
        "first_descent_value": z,
        "first_descent_minus_signed_mean": z - float(row["signed_mean"]),
        "previous_ratio": previous_ratio,
        "pressure": pressure,
        "reserve": 1.0 - pressure,
        "variance_times_reserve": variance * (1.0 - pressure),
        "effective_ratio_drop_factor": effective_ratio,
        "effective_ratio_drop": effective_drop,
        "variance_times_effective_ratio_drop": variance * effective_drop,
        "x_identity_pressure": x_ratio_identity,
        "x_identity_pressure_core": x_ratio_core,
        "x_identity_pressure_lower_boundary": x_boundary,
        "x_identity_previous_ratio": x_previous_identity,
        "x_identity_previous_ratio_core": x_previous_core,
        "x_identity_previous_lower_boundary": x_previous_boundary,
        "y_identity_pressure": y_ratio_identity,
        "y_identity_pressure_core": y_ratio_core,
        "y_identity_pressure_upper_boundary": y_boundary,
        "max_identity_error": max(finite_identity_errors, default=math.nan),
        "conditional_y_mean_at_descent": y_stats["mean"],
        "conditional_y_variance_at_descent": y_stats["variance"],
        "conditional_x_mean_at_descent": x_index_mean,
        "conditional_x_variance_at_descent": y_stats["variance"],
        "conditional_y_mean_previous": y_prev_stats["mean"],
        "conditional_x_mean_previous": x_prev_index_mean,
        "conditional_y_total_variation_from_previous": total_variation(
            y_values,
            weights,
            y_prev_values,
            prev_weights,
        ),
        "conditional_x_index_plus_one_over_variance": (x_index_mean + 1.0) / variance,
        "conditional_y_index_plus_one_over_variance": (y_stats["mean"] + 1.0)
        / variance,
        "conditional_expected_inverse_x_index_plus_one": inv_index,
        "variance_times_expected_inverse_x_index_plus_one": variance * inv_index,
        "conditional_newton_factor_average": conditional_newton_factor,
        "x_blocks": row["x_blocks"],
        "y_blocks": row["y_blocks"],
    }


def extract_rows(path: Path, *, top: int) -> list[tuple[str, dict[str, Any]]]:
    data = json.loads(path.read_text())
    rows: list[tuple[str, dict[str, Any]]] = []
    for item in data.get("best_by_variance_cutoff", []):
        best = item.get("best")
        if best is not None:
            rows.append((f"{path.name}:cutoff={item.get('variance_cutoff')}", best))
    for index, row in enumerate(data.get("best_by_variance_times_reserve", [])[:top]):
        rows.append((f"{path.name}:best_vres_rank={index + 1}", row))
    for index, row in enumerate(data.get("top_smallest_v_reserve", [])[:top]):
        rows.append((f"{path.name}:top_probe_rank={index + 1}", row))
    return rows


def row_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        tuple((int(count), round(float(p), 12)) for count, p in row["x_blocks"]),
        tuple((int(count), round(float(p), 12)) for count, p in row["y_blocks"]),
    )


def compact_analysis(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "source",
        "x_n",
        "y_n",
        "variance",
        "x_variance",
        "y_variance",
        "side_variance_fraction",
        "first_descent_value",
        "first_descent_minus_signed_mean",
        "previous_ratio",
        "pressure",
        "reserve",
        "variance_times_reserve",
        "effective_ratio_drop_factor",
        "effective_ratio_drop",
        "variance_times_effective_ratio_drop",
        "conditional_x_mean_at_descent",
        "conditional_y_mean_at_descent",
        "conditional_x_index_plus_one_over_variance",
        "conditional_y_index_plus_one_over_variance",
        "variance_times_expected_inverse_x_index_plus_one",
        "conditional_y_total_variation_from_previous",
        "max_identity_error",
    ]
    return {key: row[key] for key in keys}


def parse_paths(raw: list[str] | None) -> list[Path]:
    if not raw:
        return DEFAULT_INPUTS
    return [Path(item) for item in raw]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="*")
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/signed_pb_conditional_analysis_2026-07-03.json"),
    )
    args = parser.parse_args()

    seen = set()
    analyses: list[dict[str, Any]] = []
    for path in parse_paths(args.inputs):
        if not path.exists():
            raise SystemExit(f"missing input: {path}")
        for source, row in extract_rows(path, top=args.top):
            key = row_key(row)
            if key in seen:
                continue
            seen.add(key)
            analysis = analyze_row(row, source)
            if analysis is not None:
                analyses.append(analysis)

    analyses.sort(key=lambda row: row["variance_times_reserve"])
    summary = {
        "source": {
            "kind": "signed_pb_conditional_ratio_analysis",
            "inputs": [str(path) for path in parse_paths(args.inputs)],
            "top_per_section": args.top,
        },
        "processed": len(analyses),
        "max_identity_error": max(
            (
                row["max_identity_error"]
                for row in analyses
                if math.isfinite(row["max_identity_error"])
            ),
            default=math.nan,
        ),
        "best_by_variance_times_reserve": [
            compact_analysis(row) for row in analyses[: args.top]
        ],
        "largest_conditional_index_over_variance": [
            compact_analysis(row)
            for row in sorted(
                analyses,
                key=lambda row: row["conditional_x_index_plus_one_over_variance"],
                reverse=True,
            )[: args.top]
        ],
        "largest_conditional_tv_shift": [
            compact_analysis(row)
            for row in sorted(
                analyses,
                key=lambda row: row["conditional_y_total_variation_from_previous"],
                reverse=True,
            )[: args.top]
        ],
        "analyses": analyses,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Wrote {args.out}; processed={len(analyses)}, "
        f"max_identity_error={summary['max_identity_error']:.3e}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
