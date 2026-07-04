#!/usr/bin/env python3
"""Probe the elementary D=4 variance bound.

For low-probability Bernoulli sums, write odds w_i = p_i/(1-p_i) <= 1.
The remaining one-sided localization gap reduces to the candidate inequality

    e_1 >= e_0, e_2 >= e_1, e_3 >= e_2
        =>  sum_i w_i/(1+w_i)^2 >= 5/4,

or, more narrowly,

    first strict descent of e_k(w) occurs at D=4
        =>  sum_i w_i/(1+w_i)^2 >= 5/4.

It also tracks the still stronger-looking nondegenerate shortcut
e_2 > 0 and e_3 >= e_2, because the low-support equality cases are the
only known obstruction to the naive e_3 >= e_2 formulation.

This script is a falsification probe, not a proof.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
from typing import Any


def elementary_prefix(weights: list[float], degree: int = 5) -> list[float]:
    coeffs = [1.0] + [0.0] * degree
    for w in weights:
        for k in range(degree, 0, -1):
            coeffs[k] += coeffs[k - 1] * w
    return coeffs


def variance_from_odds(weights: list[float]) -> float:
    return sum(w / ((1.0 + w) ** 2) for w in weights)


def mean_from_odds(weights: list[float]) -> float:
    return sum(w / (1.0 + w) for w in weights)


def metric(weights: list[float], kind: str, details: dict[str, Any] | None = None) -> dict[str, Any]:
    coeffs = elementary_prefix(weights)
    mean = mean_from_odds(weights)
    variance = variance_from_odds(weights)
    first_descent_prefix = None
    for k in range(1, len(coeffs)):
        if coeffs[k] < coeffs[k - 1]:
            first_descent_prefix = k
            break
    pre_d4_guard = (
        coeffs[1] >= coeffs[0]
        and coeffs[2] >= coeffs[1]
        and coeffs[3] >= coeffs[2]
    )
    d4_case = (
        len([w for w in weights if w > 0.0]) >= 5
        and pre_d4_guard
        and coeffs[4] < coeffs[3]
    )
    weak_guard = coeffs[3] >= coeffs[2]
    weak_positive_guard = coeffs[2] > 1e-12 and weak_guard
    row = {
        "kind": kind,
        "n": len(weights),
        "support_size": len([w for w in weights if w > 0.0]),
        "e1": coeffs[1],
        "e2": coeffs[2],
        "e3": coeffs[3],
        "e4": coeffs[4],
        "e5": coeffs[5],
        "e3_minus_e2": coeffs[3] - coeffs[2],
        "e4_minus_e3": coeffs[4] - coeffs[3],
        "mean": mean,
        "mean_minus_5_over_2": mean - 2.5,
        "variance": variance,
        "variance_minus_5_over_4": variance - 1.25,
        "first_descent_prefix": first_descent_prefix,
        "weak_guard_e3_ge_e2": weak_guard,
        "weak_bound_ok": not weak_guard or variance + 1e-12 >= 1.25,
        "weak_positive_guard": weak_positive_guard,
        "weak_positive_bound_ok": not weak_positive_guard or variance + 1e-12 >= 1.25,
        "pre_d4_guard": pre_d4_guard,
        "pre_d4_bound_ok": not pre_d4_guard or variance + 1e-12 >= 1.25,
        "d4_case": d4_case,
        "d4_bound_ok": not d4_case or variance + 1e-12 >= 1.25,
        "weights": weights,
    }
    if details:
        row.update(details)
    return row


def random_weights(rng: random.Random, n: int, family: str) -> list[float]:
    if family == "uniform":
        return [rng.random() for _ in range(n)]
    if family == "edge":
        return [rng.betavariate(0.25, 0.25) for _ in range(n)]
    if family == "near_one":
        return [1.0 - rng.random() ** 5 for _ in range(n)]
    if family == "sparse":
        return [rng.betavariate(0.3, 5.0) for _ in range(n)]
    return [rng.betavariate(1.0, 2.0) for _ in range(n)]


def compact(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "kind",
        "n",
        "support_size",
        "mean",
        "mean_minus_5_over_2",
        "variance",
        "variance_minus_5_over_4",
        "first_descent_prefix",
        "e1",
        "e2",
        "e3",
        "e4",
        "e5",
        "e3_minus_e2",
        "e4_minus_e3",
        "m",
        "r",
        "family",
        "weights",
    ]
    return {key: row[key] for key in keys if key in row}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--samples-per-n", type=int, default=1000)
    parser.add_argument("--max-n", type=int, default=80)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/d4_variance_bound_probe_2026-07-04.json"),
    )
    args = parser.parse_args()

    rng = random.Random(args.seed)
    rows: list[dict[str, Any]] = []
    families = ["uniform", "edge", "near_one", "sparse", "mixed"]
    for n in range(3, args.max_n + 1):
        for _ in range(args.samples_per_n):
            family = rng.choice(families)
            rows.append(metric(random_weights(rng, n, family), "random", {"family": family}))

    for m in range(0, args.max_n + 1):
        for i in range(1001):
            r = i / 1000.0
            weights = [1.0] * m
            if r > 0.0:
                weights.append(r)
            if weights:
                rows.append(metric(weights, "ones_plus_remainder", {"m": m, "r": r}))

    weak_feasible = [row for row in rows if row["weak_guard_e3_ge_e2"]]
    weak_failures = [row for row in rows if not row["weak_bound_ok"]]
    weak_positive_cases = [row for row in rows if row["weak_positive_guard"]]
    weak_positive_failures = [row for row in rows if not row["weak_positive_bound_ok"]]
    pre_d4_cases = [row for row in rows if row["pre_d4_guard"]]
    pre_d4_failures = [row for row in rows if not row["pre_d4_bound_ok"]]
    d4_cases = [row for row in rows if row["d4_case"]]
    d4_failures = [row for row in rows if not row["d4_bound_ok"]]
    summary = {
        "source": {
            "kind": "d4_variance_bound_probe",
            "seed": args.seed,
            "samples_per_n": args.samples_per_n,
            "max_n": args.max_n,
        },
        "processed": len(rows),
        "weak_guard_e3_ge_e2_rows": len(weak_feasible),
        "weak_guard_failures": len(weak_failures),
        "weak_positive_guard_cases": len(weak_positive_cases),
        "weak_positive_guard_failures": len(weak_positive_failures),
        "pre_d4_guard_cases": len(pre_d4_cases),
        "pre_d4_guard_failures": len(pre_d4_failures),
        "d4_cases": len(d4_cases),
        "d4_failures": len(d4_failures),
        "best_weak_positive_by_variance": [
            compact(row)
            for row in sorted(weak_positive_cases, key=lambda row: row["variance"])[: args.top]
        ],
        "best_weak_positive_by_mean": [
            compact(row)
            for row in sorted(weak_positive_cases, key=lambda row: row["mean"])[: args.top]
        ],
        "best_pre_d4_by_variance": [
            compact(row)
            for row in sorted(pre_d4_cases, key=lambda row: row["variance"])[: args.top]
        ],
        "best_pre_d4_by_mean": [
            compact(row)
            for row in sorted(pre_d4_cases, key=lambda row: row["mean"])[: args.top]
        ],
        "best_d4_by_variance": [
            compact(row)
            for row in sorted(d4_cases, key=lambda row: row["variance"])[: args.top]
        ],
        "best_d4_by_mean": [
            compact(row)
            for row in sorted(d4_cases, key=lambda row: row["mean"])[: args.top]
        ],
        "closest_d4_to_boundary": [
            compact(row)
            for row in sorted(
                d4_cases,
                key=lambda row: (abs(row["variance_minus_5_over_4"]), row["variance"]),
            )[: args.top]
        ],
        "d4_failures_rows": [
            compact(row)
            for row in sorted(d4_failures, key=lambda row: row["variance_minus_5_over_4"])[
                : args.top
            ]
        ],
        "pre_d4_guard_failures_rows": [
            compact(row)
            for row in sorted(pre_d4_failures, key=lambda row: row["variance_minus_5_over_4"])[
                : args.top
            ]
        ],
        "weak_positive_guard_failures_rows": [
            compact(row)
            for row in sorted(weak_positive_failures, key=lambda row: row["variance_minus_5_over_4"])[
                : args.top
            ]
        ],
        "weak_guard_failures_rows": [
            compact(row)
            for row in sorted(weak_failures, key=lambda row: row["variance_minus_5_over_4"])[
                : args.top
            ]
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Wrote {args.out}; processed={len(rows)}, "
        f"weak_positive_cases={len(weak_positive_cases)}, "
        f"weak_positive_failures={len(weak_positive_failures)}, "
        f"pre_d4_cases={len(pre_d4_cases)}, "
        f"pre_d4_failures={len(pre_d4_failures)}, "
        f"d4_cases={len(d4_cases)}, d4_failures={len(d4_failures)}, "
        f"weak_failures={len(weak_failures)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
