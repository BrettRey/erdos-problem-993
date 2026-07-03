#!/usr/bin/env python3
"""Probe the two-sided sparse Poisson limit for PB reserve.

Grouped Poisson-binomial optimizers push hard examples toward a nearly
deterministic shift plus many tiny success probabilities and a few tiny failure
probabilities.  After factoring the near-deterministic variables, the limiting
local law is the difference

    Z = Pois(lambda) - Pois(eta).

This script scans that Skellam-type limit for small values of
``variance * first-post-descent reserve``.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover - environment guard
    raise SystemExit("This exploratory probe requires numpy.") from exc


def poisson_shape(mu: float, *, tail_sigma: float, tail_extra: int) -> np.ndarray:
    if mu < 0:
        raise ValueError("Poisson mean must be nonnegative")
    if mu == 0:
        return np.array([1.0], dtype=np.float64)

    mode = int(math.floor(mu))
    k_max = max(
        mode + 1,
        int(math.ceil(mu + tail_sigma * math.sqrt(mu + 1.0) + tail_extra)),
    )
    pmf = np.zeros(k_max + 1, dtype=np.float64)
    pmf[mode] = 1.0
    for k in range(mode, k_max):
        pmf[k + 1] = pmf[k] * mu / float(k + 1)
    for k in range(mode, 0, -1):
        pmf[k - 1] = pmf[k] * float(k) / mu
    return pmf / float(pmf.sum())


def first_descent(seq: np.ndarray) -> int | None:
    for k in range(1, len(seq)):
        if seq[k] < seq[k - 1]:
            return k
    return None


def skellam_shape(
    lam: float,
    eta: float,
    *,
    tail_sigma: float,
    tail_extra: int,
    trim_relative: float,
) -> tuple[np.ndarray, int]:
    left = poisson_shape(lam, tail_sigma=tail_sigma, tail_extra=tail_extra)
    right = poisson_shape(eta, tail_sigma=tail_sigma, tail_extra=tail_extra)
    pmf = np.convolve(left, right[::-1])
    support_start = -(len(right) - 1)

    threshold = float(pmf.max()) * trim_relative
    active = np.flatnonzero(pmf >= threshold)
    if len(active) == 0:
        return pmf / float(pmf.sum()), support_start
    lo = max(0, int(active[0]) - 2)
    hi = min(len(pmf) - 1, int(active[-1]) + 2)
    trimmed = pmf[lo : hi + 1]
    return trimmed / float(trimmed.sum()), support_start + lo


def metric(
    lam: float,
    eta: float,
    *,
    tail_sigma: float,
    tail_extra: int,
    trim_relative: float,
) -> dict[str, Any] | None:
    pmf, support_start = skellam_shape(
        lam,
        eta,
        tail_sigma=tail_sigma,
        tail_extra=tail_extra,
        trim_relative=trim_relative,
    )
    descent_idx = first_descent(pmf)
    if descent_idx is None or descent_idx >= len(pmf) - 1 or pmf[descent_idx] <= 0:
        return None
    pressure = float(pmf[descent_idx + 1] / pmf[descent_idx])
    variance = lam + eta
    return {
        "lambda": lam,
        "eta": eta,
        "mean": lam - eta,
        "variance": variance,
        "first_descent": support_start + descent_idx,
        "post_descent_pressure": pressure,
        "post_descent_reserve": 1.0 - pressure,
        "variance_times_reserve": variance * (1.0 - pressure),
    }


def parse_floats(spec: str) -> list[float]:
    return [float(part.strip()) for part in spec.split(",") if part.strip()]


def eta_grid(variance: float, grid_size: int, boundary_exponents: int) -> list[float]:
    values = {0.0, variance}
    if grid_size > 1:
        for i in range(grid_size):
            values.add(variance * i / float(grid_size - 1))
    for exponent in range(1, boundary_exponents + 1):
        eps = variance * 10.0 ** (-exponent)
        values.add(eps)
        values.add(max(0.0, variance - eps))
    return sorted(values)


def compact(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "lambda": row["lambda"],
        "eta": row["eta"],
        "mean": row["mean"],
        "variance": row["variance"],
        "first_descent": row["first_descent"],
        "post_descent_pressure": row["post_descent_pressure"],
        "post_descent_reserve": row["post_descent_reserve"],
        "variance_times_reserve": row["variance_times_reserve"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variance-values", default="1,2,5,10,20,50,100")
    parser.add_argument("--grid-size", type=int, default=801)
    parser.add_argument("--boundary-exponents", type=int, default=8)
    parser.add_argument("--tail-sigma", type=float, default=12.0)
    parser.add_argument("--tail-extra", type=int, default=30)
    parser.add_argument("--trim-relative", type=float, default=1e-15)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/skellam_reserve_probe_2026-07-03.json"),
    )
    args = parser.parse_args()
    if args.grid_size < 2:
        parser.error("--grid-size must be at least 2")
    if args.boundary_exponents < 0:
        parser.error("--boundary-exponents must be nonnegative")

    rows: list[dict[str, Any]] = []
    best_by_variance = []
    for variance in parse_floats(args.variance_values):
        variance_rows = []
        for eta in eta_grid(variance, args.grid_size, args.boundary_exponents):
            lam = variance - eta
            row = metric(
                lam,
                eta,
                tail_sigma=args.tail_sigma,
                tail_extra=args.tail_extra,
                trim_relative=args.trim_relative,
            )
            if row is None:
                continue
            rows.append(row)
            variance_rows.append(row)
        best = min(variance_rows, key=lambda row: row["variance_times_reserve"])
        best_by_variance.append({"variance": variance, "best": compact(best)})

    near_sharp_v1 = []
    for exponent in range(1, args.boundary_exponents + 1):
        eta = 10.0 ** (-exponent)
        row = metric(
            1.0 - eta,
            eta,
            tail_sigma=args.tail_sigma,
            tail_extra=args.tail_extra,
            trim_relative=args.trim_relative,
        )
        if row is not None:
            near_sharp_v1.append(compact(row))

    summary = {
        "source": {
            "kind": "skellam_sparse_limit_reserve_probe",
            "variance_values": parse_floats(args.variance_values),
            "grid_size": args.grid_size,
            "boundary_exponents": args.boundary_exponents,
            "tail_sigma": args.tail_sigma,
            "tail_extra": args.tail_extra,
            "trim_relative": args.trim_relative,
        },
        "processed": len(rows),
        "best_by_variance": best_by_variance,
        "near_sharp_v1_boundary": near_sharp_v1,
        "top_smallest_variance_times_reserve": [
            compact(row)
            for row in sorted(rows, key=lambda row: row["variance_times_reserve"])[
                : args.top
            ]
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {args.out}; processed={len(rows)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
