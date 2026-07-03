#!/usr/bin/env python3
"""Sanity probe for the low-probability PB reserve lemma."""

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


def poisson_binomial_pmf(ps: list[float]) -> np.ndarray:
    pmf = np.array([1.0], dtype=np.float64)
    for p in ps:
        pmf = np.convolve(pmf, np.array([1.0 - p, p], dtype=np.float64))
        total = float(pmf.sum())
        if total == 0.0:
            raise FloatingPointError("probability mass underflowed to zero")
        pmf /= total
    return pmf / float(pmf.sum())


def first_descent(seq: np.ndarray) -> int | None:
    for k in range(1, len(seq)):
        if seq[k] < seq[k - 1]:
            return k
    return None


def metric(ps: list[float], kind: str) -> dict[str, Any] | None:
    pmf = poisson_binomial_pmf(ps)
    descent = first_descent(pmf)
    if descent is None or descent >= len(pmf) - 1 or pmf[descent] <= 0:
        return None
    mean = sum(ps)
    variance = sum(p * (1.0 - p) for p in ps)
    if variance <= 0:
        return None
    pressure = float(pmf[descent + 1] / pmf[descent])
    reserve = 1.0 - pressure
    return {
        "kind": kind,
        "n": len(ps),
        "mean": mean,
        "variance": variance,
        "first_descent": descent,
        "pressure": pressure,
        "reserve": reserve,
        "variance_times_reserve": variance * reserve,
        "newton_floor": 1.0 / (descent + 1),
        "mu_floor": 1.0 / (mean + 3.0),
        "variance_floor": 1.0 / (2.0 * variance + 3.0),
        "target_floor": 1.0 / (5.0 * variance),
        "target_ok": variance < 1.0 or reserve > 1.0 / (5.0 * variance),
        "min_p": min(ps),
        "max_p": max(ps),
    }


def random_parameter_sets(rng: random.Random, n: int, samples: int) -> list[tuple[str, list[float]]]:
    out: list[tuple[str, list[float]]] = []
    for _ in range(samples):
        out.append(("uniform_0_half", [0.5 * rng.random() for _ in range(n)]))
        out.append(("beta_sparse", [0.5 * rng.betavariate(0.35, 4.0) for _ in range(n)]))
        out.append(("beta_mid", [0.5 * rng.betavariate(2.0, 2.0) for _ in range(n)]))
        p = 0.5 * rng.random()
        out.append(("binomial_random_p", [p for _ in range(n)]))
    return out


def compact(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "kind",
        "n",
        "mean",
        "variance",
        "first_descent",
        "pressure",
        "reserve",
        "variance_times_reserve",
        "newton_floor",
        "mu_floor",
        "variance_floor",
        "target_floor",
        "min_p",
        "max_p",
    ]
    return {key: row[key] for key in keys}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--samples", type=int, default=200)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/low_probability_pb_reserve_probe_2026-07-03.json"),
    )
    args = parser.parse_args()
    rng = random.Random(args.seed)

    rows: list[dict[str, Any]] = []
    for n in [5, 10, 20, 50, 100, 200, 300]:
        for kind, ps in random_parameter_sets(rng, n, args.samples):
            row = metric(ps, kind)
            if row is not None:
                rows.append(row)

    failures = [row for row in rows if not row["target_ok"]]
    v_ge_1 = [row for row in rows if row["variance"] >= 1.0]
    summary = {
        "source": {
            "kind": "low_probability_pb_reserve_probe",
            "seed": args.seed,
            "samples": args.samples,
            "n_values": [5, 10, 20, 50, 100, 200, 300],
        },
        "processed": len(rows),
        "counts": {
            "variance_ge_1": len(v_ge_1),
            "target_failures": len(failures),
        },
        "top_smallest_v_reserve": [
            compact(row)
            for row in sorted(v_ge_1, key=lambda row: row["variance_times_reserve"])[
                : args.top
            ]
        ],
        "closest_to_target": [
            compact(row)
            for row in sorted(
                v_ge_1,
                key=lambda row: row["reserve"] - row["target_floor"],
            )[: args.top]
        ],
        "failures": [compact(row) for row in failures[: args.top]],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Wrote {args.out}; processed={len(rows)}, failures={len(failures)}",
        flush=True,
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
