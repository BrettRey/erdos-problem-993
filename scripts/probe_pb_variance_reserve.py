#!/usr/bin/env python3
"""Adversarial probe for the Poisson-binomial variance reserve.

This is an exploratory falsification script for the candidate lemma in issue #5:

    first post-descent reserve >= c / Var.

It searches binomial grids and deterministic random Poisson-binomial samples for
small values of ``Var * reserve``.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
from typing import Any

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover - environment guard
    raise SystemExit("This exploratory probe requires numpy.") from exc


def first_descent(seq: np.ndarray) -> int | None:
    for k in range(1, len(seq)):
        if seq[k] < seq[k - 1]:
            return k
    return None


def pmf_metric(pmf: np.ndarray, mean: float, variance: float) -> dict[str, Any] | None:
    descent = first_descent(pmf)
    if descent is None or descent >= len(pmf) - 1 or pmf[descent] <= 0:
        return None
    pressure = float(pmf[descent + 1] / pmf[descent])
    reserve = 1.0 - pressure
    return {
        "mean": mean,
        "variance": variance,
        "first_descent": descent,
        "crossing_pressure": pressure,
        "crossing_reserve": reserve,
        "variance_times_reserve": variance * reserve,
    }


def poisson_binomial_pmf(ps: list[float]) -> np.ndarray:
    pmf = np.array([1.0], dtype=np.float64)
    for p in ps:
        pmf = np.convolve(pmf, np.array([1.0 - p, p], dtype=np.float64))
        total = float(pmf.sum())
        if total == 0.0:
            raise FloatingPointError("probability mass underflowed to zero")
        if pmf.max() > 1e100 or total < 1e-100:
            pmf /= total
    return pmf / float(pmf.sum())


def binomial_pmf_shape(n: int, p: float) -> np.ndarray:
    mode_guess = min(max(int((n + 1) * p), 0), n)
    pmf = np.zeros(n + 1, dtype=np.float64)
    pmf[mode_guess] = 1.0
    if p == 0.0 or p == 1.0:
        return pmf
    odds = p / (1.0 - p)
    for k in range(mode_guess, n):
        pmf[k + 1] = pmf[k] * (n - k) / (k + 1) * odds
    for k in range(mode_guess, 0, -1):
        pmf[k - 1] = pmf[k] * k / (n - k + 1) / odds
    return pmf / float(pmf.sum())


def binomial_metric(n: int, p: float) -> dict[str, Any] | None:
    mean = n * p
    variance = n * p * (1.0 - p)
    metric = pmf_metric(binomial_pmf_shape(n, p), mean, variance)
    if metric is None:
        return None
    metric.update({"kind": "binomial", "n": n, "p": p})
    return metric


def best_by_variance_cutoff(rows: list[dict[str, Any]], cutoffs: list[float]) -> list[dict[str, Any]]:
    out = []
    for cutoff in cutoffs:
        candidates = [row for row in rows if row["variance"] >= cutoff]
        if not candidates:
            out.append({"variance_cutoff": cutoff, "best": None})
            continue
        best = min(candidates, key=lambda row: row["variance_times_reserve"])
        out.append({"variance_cutoff": cutoff, "best": compact(best)})
    return out


def compact(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "kind",
        "n",
        "p",
        "family",
        "mean",
        "variance",
        "first_descent",
        "crossing_pressure",
        "crossing_reserve",
        "variance_times_reserve",
        "min_p",
        "max_p",
        "mean_p",
    ]
    return {key: row[key] for key in keys if key in row}


def random_parameter_sets(rng: random.Random, n: int, samples: int) -> list[tuple[str, list[float]]]:
    out: list[tuple[str, list[float]]] = []
    for _ in range(samples):
        out.append(("uniform", [rng.random() for _ in range(n)]))
        out.append(("beta_0.2_0.2", [rng.betavariate(0.2, 0.2) for _ in range(n)]))
        out.append(("beta_0.5_5", [rng.betavariate(0.5, 5.0) for _ in range(n)]))
        out.append(("beta_5_0.5", [rng.betavariate(5.0, 0.5) for _ in range(n)]))
        a = rng.random()
        b = rng.random()
        mix = rng.random()
        out.append(("two_point", [a if rng.random() < mix else b for _ in range(n)]))
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--p-grid-size", type=int, default=999)
    parser.add_argument("--random-samples", type=int, default=100)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/pb_variance_reserve_probe_2026-07-03.json"),
    )
    args = parser.parse_args()
    rng = random.Random(args.seed)

    binomial_rows: list[dict[str, Any]] = []
    p_values = np.linspace(0.001, 0.999, args.p_grid_size)
    for n in [10, 20, 50, 100, 200, 500, 1000, 5000]:
        for p_raw in p_values:
            metric = binomial_metric(n, float(p_raw))
            if metric is not None:
                binomial_rows.append(metric)

    random_rows: list[dict[str, Any]] = []
    for n in [20, 50, 100, 200, 300]:
        for family, ps in random_parameter_sets(rng, n, args.random_samples):
            pmf = poisson_binomial_pmf(ps)
            metric = pmf_metric(
                pmf,
                mean=sum(ps),
                variance=sum(p * (1.0 - p) for p in ps),
            )
            if metric is None:
                continue
            metric.update(
                {
                    "kind": "random_poisson_binomial",
                    "family": family,
                    "n": n,
                    "min_p": min(ps),
                    "max_p": max(ps),
                    "mean_p": sum(ps) / len(ps),
                }
            )
            random_rows.append(metric)

    all_rows = binomial_rows + random_rows
    summary = {
        "source": {
            "kind": "pb_variance_reserve_probe",
            "seed": args.seed,
            "p_grid_size": args.p_grid_size,
            "random_samples": args.random_samples,
        },
        "processed": {
            "binomial_rows": len(binomial_rows),
            "random_rows": len(random_rows),
            "total_rows": len(all_rows),
        },
        "best_by_variance_cutoff": best_by_variance_cutoff(
            all_rows, [1, 2, 5, 10, 20, 50, 100, 500]
        ),
        "top_smallest_variance_times_reserve": [
            compact(row)
            for row in sorted(all_rows, key=lambda row: row["variance_times_reserve"])[: args.top]
        ],
        "top_random_smallest_variance_times_reserve": [
            compact(row)
            for row in sorted(random_rows, key=lambda row: row["variance_times_reserve"])[
                : args.top
            ]
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Wrote {args.out}; binomial={len(binomial_rows)}, random={len(random_rows)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
