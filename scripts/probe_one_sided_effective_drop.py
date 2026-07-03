#!/usr/bin/env python3
"""Probe one-sided effective ratio drop for low-probability PB laws."""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path
from typing import Any

import numpy as np

from probe_signed_pb_reserve import (
    Block,
    first_descent,
    grouped_pmf,
    random_grouped_blocks,
)


def metric(blocks: list[Block], kind: str, details: dict[str, Any] | None = None) -> dict[str, Any] | None:
    pmf = grouped_pmf(blocks)
    descent = first_descent(pmf)
    if descent is None or descent <= 0 or descent >= len(pmf) - 1:
        return None
    previous_ratio = float(pmf[descent] / pmf[descent - 1])
    pressure = float(pmf[descent + 1] / pmf[descent])
    effective_ratio = pressure / previous_ratio
    effective_drop = 1.0 - effective_ratio
    mean = sum(count * p for count, p in blocks)
    variance = sum(count * p * (1.0 - p) for count, p in blocks)
    if variance <= 0.0:
        return None
    row = {
        "kind": kind,
        "n": sum(count for count, _ in blocks),
        "blocks": [[count, p] for count, p in blocks],
        "mean": mean,
        "variance": variance,
        "first_descent": descent,
        "previous_ratio": previous_ratio,
        "pressure": pressure,
        "reserve": 1.0 - pressure,
        "variance_times_reserve": variance * (1.0 - pressure),
        "effective_ratio_drop_factor": effective_ratio,
        "effective_ratio_drop": effective_drop,
        "variance_times_effective_ratio_drop": variance * effective_drop,
        "mode_variance_ratio": (descent + 1.0) / variance,
        "quarter_effective_ok": variance < 1.0 or effective_drop >= 1.0 / (4.0 * variance),
        "mode_variance_4_ok": variance < 1.0 or descent + 1.0 <= 4.0 * variance,
    }
    if details:
        row.update(details)
    return row


def compact(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "kind",
        "n",
        "mean",
        "variance",
        "first_descent",
        "previous_ratio",
        "pressure",
        "variance_times_reserve",
        "effective_ratio_drop_factor",
        "variance_times_effective_ratio_drop",
        "mode_variance_ratio",
        "lambda",
        "p",
        "blocks",
    ]
    return {key: row[key] for key in keys if key in row}


def add_random_rows(rows: list[dict[str, Any]], *, rng: random.Random, samples: int, max_groups: int) -> None:
    n_values = [5, 10, 20, 50, 100, 200, 300]
    for _ in range(samples):
        n = rng.choice(n_values)
        row = metric(random_grouped_blocks(rng, n, max_groups), "random_grouped")
        if row is not None:
            rows.append(row)


def add_binomial_grid(rows: list[dict[str, Any]], *, grid_size: int) -> None:
    for n in [5, 10, 20, 50, 100, 200, 300, 500]:
        for p_raw in np.linspace(0.002, 0.5, grid_size):
            p = float(p_raw)
            row = metric([(n, p)], "binomial_grid", {"p": p})
            if row is not None:
                rows.append(row)


def add_finite_poisson(rows: list[dict[str, Any]]) -> None:
    lambdas = [1.0, 1.01, 1.05, 1.1, 1.25, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0]
    approximants = [50, 100, 200, 500, 1000, 5000]
    for lam in lambdas:
        for base in approximants:
            n = max(1, int(math.ceil(base * lam)))
            p = min(0.5, lam / n)
            row = metric([(n, p)], "finite_poisson", {"lambda": lam, "p": p})
            if row is not None:
                rows.append(row)


def best_by_cutoff(rows: list[dict[str, Any]], cutoffs: list[float], key: str, *, reverse: bool = False) -> list[dict[str, Any]]:
    out = []
    for cutoff in cutoffs:
        candidates = [row for row in rows if row["variance"] >= cutoff]
        best = sorted(candidates, key=lambda row: row[key], reverse=reverse)
        out.append({"variance_cutoff": cutoff, "best": compact(best[0]) if best else None})
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--random-samples", type=int, default=4000)
    parser.add_argument("--max-groups", type=int, default=6)
    parser.add_argument("--binomial-grid-size", type=int, default=501)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/one_sided_effective_drop_probe_2026-07-03.json"),
    )
    args = parser.parse_args()

    rng = random.Random(args.seed)
    rows: list[dict[str, Any]] = []
    add_random_rows(rows, rng=rng, samples=args.random_samples, max_groups=args.max_groups)
    add_binomial_grid(rows, grid_size=args.binomial_grid_size)
    add_finite_poisson(rows)

    v_ge_1 = [row for row in rows if row["variance"] >= 1.0]
    quarter_failures = [row for row in v_ge_1 if not row["quarter_effective_ok"]]
    mode_variance_failures = [row for row in v_ge_1 if not row["mode_variance_4_ok"]]
    summary = {
        "source": {
            "kind": "one_sided_low_probability_effective_drop_probe",
            "seed": args.seed,
            "random_samples": args.random_samples,
            "max_groups": args.max_groups,
            "binomial_grid_size": args.binomial_grid_size,
        },
        "processed": len(rows),
        "variance_ge_1": len(v_ge_1),
        "quarter_effective_failures": len(quarter_failures),
        "mode_variance_4_failures": len(mode_variance_failures),
        "best_effective_drop_by_cutoff": best_by_cutoff(
            v_ge_1,
            [1, 1.25, 1.5, 2, 5, 10, 20, 50],
            "variance_times_effective_ratio_drop",
        ),
        "largest_mode_variance_ratio_by_cutoff": best_by_cutoff(
            v_ge_1,
            [1, 1.25, 1.5, 2, 5, 10, 20, 50],
            "mode_variance_ratio",
            reverse=True,
        ),
        "smallest_effective_drop_rows": [
            compact(row)
            for row in sorted(v_ge_1, key=lambda row: row["variance_times_effective_ratio_drop"])[
                : args.top
            ]
        ],
        "largest_mode_variance_ratio_rows": [
            compact(row)
            for row in sorted(v_ge_1, key=lambda row: row["mode_variance_ratio"], reverse=True)[
                : args.top
            ]
        ],
        "quarter_effective_failures_rows": [
            compact(row)
            for row in sorted(quarter_failures, key=lambda row: row["variance_times_effective_ratio_drop"])[
                : args.top
            ]
        ],
        "mode_variance_4_failures_rows": [
            compact(row)
            for row in sorted(mode_variance_failures, key=lambda row: row["mode_variance_ratio"], reverse=True)[
                : args.top
            ]
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Wrote {args.out}; processed={len(rows)}, "
        f"quarter_failures={len(quarter_failures)}, "
        f"mode_variance_failures={len(mode_variance_failures)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
