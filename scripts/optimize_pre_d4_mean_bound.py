#!/usr/bin/env python3
"""Adversarially probe the pre-D4 elementary mean bound.

Target lemma:

    0 <= w_i <= 1,
    e_1 >= 1, e_2 >= e_1, e_3 >= e_2
        =>  sum_i w_i/(1+w_i) >= 5/2.

This script minimizes the mean under a quadratic penalty for violating the
three pre-D4 inequalities.  It is a falsification/diagnostic harness, not a
proof.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from probe_d4_variance_bound import elementary_prefix, mean_from_odds, variance_from_odds


def as_float_list(values: Any) -> list[float]:
    return [float(value) for value in values]


def pre_d4_violations(weights: list[float]) -> list[float]:
    coeffs = elementary_prefix(weights)
    return [
        1.0 - coeffs[1],
        coeffs[1] - coeffs[2],
        coeffs[2] - coeffs[3],
    ]


def feasible(weights: list[float], tolerance: float) -> bool:
    return max(pre_d4_violations(weights)) <= tolerance


def objective(weights: Any, *, penalty: float) -> float:
    w = as_float_list(weights)
    violation_penalty = sum(max(0.0, violation) ** 2 for violation in pre_d4_violations(w))
    return mean_from_odds(w) + penalty * violation_penalty


def compact(
    weights: list[float],
    *,
    n: int,
    seed: int,
    fun: float,
    tolerance: float,
    row_kind: str,
) -> dict[str, Any]:
    coeffs = elementary_prefix(weights)
    violations = pre_d4_violations(weights)
    mean = mean_from_odds(weights)
    variance = variance_from_odds(weights)
    return {
        "n": n,
        "seed": seed,
        "row_kind": row_kind,
        "fun": float(fun),
        "feasible": feasible(weights, tolerance),
        "max_violation": max(violations),
        "violations": violations,
        "mean": mean,
        "mean_minus_5_over_2": mean - 2.5,
        "variance": variance,
        "variance_minus_5_over_4": variance - 1.25,
        "e1": coeffs[1],
        "e2": coeffs[2],
        "e3": coeffs[3],
        "e4": coeffs[4],
        "weights_desc": sorted(weights, reverse=True),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=12)
    parser.add_argument("--restarts", type=int, default=4)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--maxiter", type=int, default=600)
    parser.add_argument("--popsize", type=int, default=12)
    parser.add_argument("--penalty", type=float, default=1e6)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/pre_d4_mean_optimizer_2026-07-04.json"),
    )
    args = parser.parse_args()

    try:
        from scipy.optimize import differential_evolution
    except ImportError as exc:
        raise SystemExit("scipy is required for this optimizer") from exc

    rows: list[dict[str, Any]] = []
    for n in range(args.min_n, args.max_n + 1):
        for restart in range(args.restarts):
            seed = args.seed + 1009 * n + restart
            result = differential_evolution(
                lambda x: objective(x, penalty=args.penalty),
                [(0.0, 1.0)] * n,
                seed=seed,
                maxiter=args.maxiter,
                popsize=args.popsize,
                polish=True,
                tol=1e-9,
                workers=1,
            )
            weights = as_float_list(result.x)
            rows.append(
                compact(
                    weights,
                    n=n,
                    seed=seed,
                    fun=float(result.fun),
                    tolerance=args.tolerance,
                    row_kind="optimizer",
                )
            )
        if n >= 5:
            weights = [1.0] * 5 + [0.0] * (n - 5)
            rows.append(
                compact(
                    weights,
                    n=n,
                    seed=-1,
                    fun=mean_from_odds(weights),
                    tolerance=args.tolerance,
                    row_kind="structured_boundary",
                )
            )

    feasible_rows = [row for row in rows if row["feasible"]]
    failures = [
        row
        for row in feasible_rows
        if row["mean"] + args.tolerance < 2.5 or row["variance"] + args.tolerance < 1.25
    ]
    summary = {
        "source": {
            "kind": "pre_d4_mean_optimizer",
            "min_n": args.min_n,
            "max_n": args.max_n,
            "restarts": args.restarts,
            "seed": args.seed,
            "maxiter": args.maxiter,
            "popsize": args.popsize,
            "penalty": args.penalty,
            "tolerance": args.tolerance,
        },
        "processed": len(rows),
        "feasible": len(feasible_rows),
        "failures": len(failures),
        "best_feasible_by_mean": sorted(feasible_rows, key=lambda row: row["mean"])[:20],
        "best_feasible_by_variance": sorted(feasible_rows, key=lambda row: row["variance"])[:20],
        "best_by_n": [
            min([row for row in feasible_rows if row["n"] == n], key=lambda row: row["mean"])
            for n in range(args.min_n, args.max_n + 1)
            if any(row["n"] == n for row in feasible_rows)
        ],
        "failures_rows": sorted(failures, key=lambda row: row["mean"])[:20],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Wrote {args.out}; processed={len(rows)}, "
        f"feasible={len(feasible_rows)}, failures={len(failures)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
