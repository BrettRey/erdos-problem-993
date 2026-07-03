#!/usr/bin/env python3
"""Floating-point probe for long broom-handle crossing pressure.

This is not a certificate script.  It studies the dominant product term
``(1+x)^s I(P_l)`` for brooms with many leaves and long handles.  The
hub-included term ``x I(P_{l-1})`` is exponentially small near the
first-descent region when ``s`` is large, so the probe is meant to test
which asymptotic denominator controls the reserve.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover - environment guard
    raise SystemExit("This exploratory probe requires numpy.") from exc

from scripts.scan_crossing_families import parse_int_values  # noqa: E402


def path_poly_shape(length: int) -> np.ndarray:
    """Return path independence coefficients up to an arbitrary scale."""
    if length < 0:
        raise ValueError("length must be nonnegative")
    if length == 0:
        return np.array([1.0], dtype=np.float64)
    if length == 1:
        return np.array([1.0, 1.0], dtype=np.float64)

    prev2 = np.array([1.0], dtype=np.float64)
    prev1 = np.array([1.0, 1.0], dtype=np.float64)
    for _ in range(2, length + 1):
        size = max(len(prev1), len(prev2) + 1)
        current = np.zeros(size, dtype=np.float64)
        current[: len(prev1)] += prev1
        current[1 : 1 + len(prev2)] += prev2

        # Global rescaling preserves all coefficient ratios and prevents overflow.
        scale = float(current.max())
        if scale > 1e100:
            current /= scale
            prev1 /= scale
            prev2 /= scale
        prev2, prev1 = prev1, current
    return prev1 / float(prev1.max())


def binomial_shape(s: int) -> np.ndarray:
    """Return binomial coefficients C(s,k) up to an arbitrary scale."""
    if s < 0:
        raise ValueError("s must be nonnegative")
    coeffs = np.zeros(s + 1, dtype=np.float64)
    mid = s // 2
    coeffs[mid] = 1.0
    for k in range(mid, s):
        coeffs[k + 1] = coeffs[k] * (s - k) / (k + 1)
    for k in range(mid, 0, -1):
        coeffs[k - 1] = coeffs[k] * k / (s - k + 1)
    return coeffs


def first_descent(seq: np.ndarray) -> int | None:
    for k in range(1, len(seq)):
        if seq[k] < seq[k - 1]:
            return k
    return None


def pressure(seq: np.ndarray) -> tuple[float, int | None, int | None]:
    descent = first_descent(seq)
    if descent is None:
        return 0.0, None, None

    left = seq[descent:-1]
    right = seq[descent + 1 :]
    mask = left > 0
    ratios = np.full_like(left, np.nan, dtype=np.float64)
    ratios[mask] = right[mask] / left[mask]
    if np.all(np.isnan(ratios)):
        return 0.0, descent, None
    rel = int(np.nanargmax(ratios))
    return float(ratios[rel]), descent, descent + rel


def distribution_stats(seq: np.ndarray) -> tuple[float, float]:
    total = float(seq.sum())
    positions = np.arange(len(seq), dtype=np.float64)
    mean = float((positions * seq).sum() / total)
    variance = float((((positions - mean) ** 2) * seq).sum() / total)
    return mean, variance


def evaluate(s: int, arm: int) -> dict[str, Any]:
    seq = np.convolve(binomial_shape(s), path_poly_shape(arm))
    crossing_pressure, descent, crossing_pos = pressure(seq)
    mean, variance = distribution_stats(seq)
    reserve = 1.0 - crossing_pressure
    n = s + arm + 1
    return {
        "s": s,
        "arm": arm,
        "n": n,
        "first_descent": descent,
        "crossing_pos": crossing_pos,
        "mean": mean,
        "variance": variance,
        "descent_minus_mean": None if descent is None else descent - mean,
        "crossing_pressure": crossing_pressure,
        "crossing_reserve": reserve,
        "s_times_reserve": s * reserve,
        "n_times_reserve": n * reserve,
        "variance_times_reserve": variance * reserve,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s", type=int, default=5000)
    parser.add_argument(
        "--arm-values",
        default="10,20,50,100,193,240,500,1000,2000,5000,10000,20000,50000",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/broom_variance_scaling_probe_2026-07-03.json"),
    )
    args = parser.parse_args()
    if args.s < 1:
        parser.error("--s must be positive")

    arms = parse_int_values(args.arm_values)
    if any(arm < 1 for arm in arms):
        parser.error("all arm lengths must be positive")

    rows = [evaluate(args.s, arm) for arm in arms]
    summary = {
        "source": {
            "kind": "broom_variance_scaling_probe",
            "s": args.s,
            "arm_values": arms,
            "approximation": "dominant product term (1+x)^s I(P_l); hub term omitted",
        },
        "processed": len(rows),
        "top_by_crossing_pressure": sorted(
            rows, key=lambda row: row["crossing_pressure"], reverse=True
        ),
        "rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {args.out}; processed={len(rows)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
