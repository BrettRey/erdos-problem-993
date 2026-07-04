#!/usr/bin/env python3
"""Breaker harness for the local-mode mean lemma.

Target under attack:

    0 <= w_i <= 1, e_2(w) > 0, e_3(w) >= e_2(w)
        =>  sum_i w_i/(1+w_i) >= 5/2.

This is a falsification script, not a proof.  It combines:

* grouped constrained optimization over multiplicity patterns;
* ungrouped constrained SLSQP restarts;
* structured boundary families around five fair weights;
* coarse two-group grids.

If a feasible row with mean < 5/2 is found, the JSON output contains it as a
counterexample candidate.  Otherwise the output is only search evidence.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable


TOL = 1e-9


def elementary_prefix(weights: Iterable[float], degree: int = 4) -> list[float]:
    coeffs = [1.0] + [0.0] * degree
    for w in weights:
        for k in range(degree, 0, -1):
            coeffs[k] += coeffs[k - 1] * w
    return coeffs


def mean_from_odds(weights: Iterable[float]) -> float:
    return sum(w / (1.0 + w) for w in weights)


def variance_from_odds(weights: Iterable[float]) -> float:
    return sum(w / ((1.0 + w) ** 2) for w in weights)


def expand_groups(counts: list[int], values: list[float]) -> list[float]:
    weights: list[float] = []
    for count, value in zip(counts, values):
        weights.extend([float(value)] * count)
    return weights


def rational(value: float, max_denominator: int = 1000) -> str:
    return str(Fraction(value).limit_denominator(max_denominator))


def row_from_weights(
    *,
    weights: list[float],
    source: str,
    details: dict[str, Any] | None = None,
) -> dict[str, Any]:
    coeffs = elementary_prefix(weights)
    mean = mean_from_odds(weights)
    variance = variance_from_odds(weights)
    e2 = coeffs[2]
    e3 = coeffs[3]
    feasible = e2 > TOL and e3 + TOL >= e2
    row: dict[str, Any] = {
        "source": source,
        "support_size": sum(1 for w in weights if w > TOL),
        "n": len(weights),
        "feasible": feasible,
        "constraint_margin_e3_minus_e2": e3 - e2,
        "e1": coeffs[1],
        "e2": e2,
        "e3": e3,
        "e4": coeffs[4],
        "mean": mean,
        "mean_minus_5_over_2": mean - 2.5,
        "variance": variance,
        "variance_minus_5_over_4": variance - 1.25,
        "weights_desc": sorted(weights, reverse=True),
        "weights_rational_approx": [rational(w) for w in sorted(weights, reverse=True)],
    }
    if details:
        row.update(details)
    return row


def compact(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "source",
        "n",
        "support_size",
        "feasible",
        "constraint_margin_e3_minus_e2",
        "mean",
        "mean_minus_5_over_2",
        "variance",
        "variance_minus_5_over_4",
        "e1",
        "e2",
        "e3",
        "e4",
        "counts",
        "values",
        "values_rational_approx",
        "grid_a",
        "grid_b",
        "weights_desc",
        "weights_rational_approx",
        "success",
        "message",
        "fun",
        "seed",
        "start_kind",
    ]
    return {key: row[key] for key in keys if key in row}


def integer_partitions_at_most(total: int, parts: int) -> Iterable[list[int]]:
    """Nonincreasing positive partitions of total into at most parts parts."""
    if parts <= 0:
        return

    def rec(remaining: int, max_next: int, prefix: list[int]) -> Iterable[list[int]]:
        if remaining == 0:
            yield prefix
            return
        if len(prefix) == parts:
            return
        for value in range(min(max_next, remaining), 0, -1):
            yield from rec(remaining - value, value, prefix + [value])

    yield from rec(total, total, [])


def seed_vectors(counts: list[int], rng: random.Random) -> list[list[float]]:
    g = len(counts)
    seeds: list[list[float]] = []
    seeds.append([1.0] * min(g, 1) + [0.2] * max(0, g - 1))
    seeds.append([1.0] * g)
    seeds.append([min(1.0, 3.0 / max(1, sum(counts) - 2))] * g)
    seeds.append([0.95] * g)
    seeds.append([0.5] * g)
    for _ in range(2):
        seeds.append([rng.random() for _ in range(g)])
    return [seed[:g] for seed in seeds]


def optimize_grouped(
    *,
    counts: list[int],
    rng: random.Random,
    restarts: int,
    maxiter: int,
    scipy_minimize: Any,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    def objective(values: Any) -> float:
        return mean_from_odds(expand_groups(counts, [float(v) for v in values]))

    def constraint(values: Any) -> float:
        coeffs = elementary_prefix(expand_groups(counts, [float(v) for v in values]))
        return coeffs[3] - coeffs[2]

    starts = seed_vectors(counts, rng)
    while len(starts) < restarts:
        starts.append([rng.random() for _ in counts])
    for start_index, start in enumerate(starts[:restarts]):
        result = scipy_minimize(
            objective,
            start,
            method="SLSQP",
            bounds=[(0.0, 1.0)] * len(counts),
            constraints=[{"type": "ineq", "fun": constraint}],
            options={"maxiter": maxiter, "ftol": 1e-12, "disp": False},
        )
        values = [float(max(0.0, min(1.0, v))) for v in result.x]
        weights = expand_groups(counts, values)
        rows.append(
            row_from_weights(
                weights=weights,
                source="grouped_slsqp",
                details={
                    "counts": counts,
                    "values": values,
                    "values_rational_approx": [rational(value) for value in values],
                    "success": bool(result.success),
                    "message": str(result.message),
                    "fun": float(result.fun),
                    "start_kind": f"seed_{start_index}",
                },
            )
        )
    return rows


def optimize_ungrouped(
    *,
    n: int,
    rng: random.Random,
    restarts: int,
    maxiter: int,
    scipy_minimize: Any,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    def objective(values: Any) -> float:
        return mean_from_odds([float(v) for v in values])

    def constraint(values: Any) -> float:
        coeffs = elementary_prefix([float(v) for v in values])
        return coeffs[3] - coeffs[2]

    starts: list[list[float]] = []
    starts.append([1.0] * min(n, 5) + [0.0] * max(0, n - 5))
    starts.append([min(1.0, 3.0 / max(1, n - 2))] * n)
    starts.append([0.95] * n)
    starts.append([0.5] * n)
    while len(starts) < restarts:
        family = rng.choice(["uniform", "near_one", "edge"])
        if family == "near_one":
            starts.append([1.0 - rng.random() ** 4 for _ in range(n)])
        elif family == "edge":
            starts.append([rng.betavariate(0.25, 0.25) for _ in range(n)])
        else:
            starts.append([rng.random() for _ in range(n)])
    for start_index, start in enumerate(starts[:restarts]):
        result = scipy_minimize(
            objective,
            start,
            method="SLSQP",
            bounds=[(0.0, 1.0)] * n,
            constraints=[{"type": "ineq", "fun": constraint}],
            options={"maxiter": maxiter, "ftol": 1e-12, "disp": False},
        )
        weights = [float(max(0.0, min(1.0, v))) for v in result.x]
        rows.append(
            row_from_weights(
                weights=weights,
                source="ungrouped_slsqp",
                details={
                    "success": bool(result.success),
                    "message": str(result.message),
                    "fun": float(result.fun),
                    "seed": start_index,
                    "start_kind": f"ungrouped_seed_{start_index}",
                },
            )
        )
    return rows


def structured_rows(max_n: int, grid_step: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    grid_points = [i * grid_step for i in range(int(round(1.0 / grid_step)) + 1)]
    for n in range(3, max_n + 1):
        if n >= 5:
            rows.append(
                row_from_weights(
                    weights=[1.0] * 5 + [0.0] * (n - 5),
                    source="five_ones_plus_zero",
                )
            )
        for m in range(0, n + 1):
            for r in grid_points:
                weights = [1.0] * m + [r] * (n - m)
                rows.append(
                    row_from_weights(
                        weights=weights,
                        source="ones_plus_equal_remainder_grid",
                        details={"counts": [m, n - m], "values": [1.0, r]},
                    )
                )

    for a_count in range(1, max_n + 1):
        for b_count in range(1, max_n - a_count + 1):
            n = a_count + b_count
            if n > max_n:
                continue
            for a in grid_points:
                for b in grid_points:
                    weights = [a] * a_count + [b] * b_count
                    rows.append(
                        row_from_weights(
                            weights=weights,
                            source="two_group_grid",
                            details={
                                "counts": [a_count, b_count],
                                "values": [a, b],
                                "grid_a": a,
                                "grid_b": b,
                            },
                        )
                    )
    return rows


def differential_evolution_rows(
    *,
    max_groups: int,
    max_support: int,
    restarts: int,
    maxiter: int,
    rng: random.Random,
    differential_evolution: Any,
    penalty: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    candidate_patterns = [
        [5],
        [4, 1],
        [3, 2],
        [3, 1, 1],
        [2, 2, 1],
        [2, 1, 1, 1],
        [1, 1, 1, 1, 1],
    ]
    for n in range(6, max_support + 1):
        candidate_patterns.append([n])
        candidate_patterns.append([n - 1, 1])
        candidate_patterns.append([n - 2, 2])
        if max_groups >= 3:
            candidate_patterns.append([n - 2, 1, 1])
        if max_groups >= 4 and n >= 7:
            candidate_patterns.append([n - 3, 1, 1, 1])

    seen: set[tuple[int, ...]] = set()
    patterns: list[list[int]] = []
    for pattern in candidate_patterns:
        if sum(pattern) <= max_support and len(pattern) <= max_groups:
            key = tuple(sorted(pattern, reverse=True))
            if key not in seen:
                seen.add(key)
                patterns.append(list(key))

    for pattern in patterns:
        for restart in range(restarts):
            seed = rng.randrange(1, 2**31 - 1)

            def objective(values: Any) -> float:
                weights = expand_groups(pattern, [float(v) for v in values])
                coeffs = elementary_prefix(weights)
                violation = max(0.0, coeffs[2] - coeffs[3])
                return mean_from_odds(weights) + penalty * violation * violation

            result = differential_evolution(
                objective,
                [(0.0, 1.0)] * len(pattern),
                seed=seed,
                maxiter=maxiter,
                popsize=10,
                polish=True,
                tol=1e-9,
                workers=1,
            )
            values = [float(max(0.0, min(1.0, v))) for v in result.x]
            rows.append(
                row_from_weights(
                    weights=expand_groups(pattern, values),
                    source="grouped_differential_evolution_penalty",
                    details={
                        "counts": pattern,
                        "values": values,
                        "values_rational_approx": [rational(value) for value in values],
                        "success": bool(result.success),
                        "message": str(result.message),
                        "fun": float(result.fun),
                        "seed": seed,
                    },
                )
            )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--max-support", type=int, default=16)
    parser.add_argument("--max-groups", type=int, default=4)
    parser.add_argument("--group-restarts", type=int, default=5)
    parser.add_argument("--ungrouped-restarts", type=int, default=5)
    parser.add_argument("--ungrouped-max-n", type=int, default=12)
    parser.add_argument("--slsqp-maxiter", type=int, default=400)
    parser.add_argument("--grid-step", type=float, default=0.05)
    parser.add_argument("--de-restarts", type=int, default=1)
    parser.add_argument("--de-maxiter", type=int, default=160)
    parser.add_argument("--de-penalty", type=float, default=1e7)
    parser.add_argument("--top", type=int, default=25)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/local_mode_mean_breaker_2026-07-04.json"),
    )
    args = parser.parse_args()

    if args.max_support < 3:
        parser.error("--max-support must be at least 3")
    if args.max_groups < 1:
        parser.error("--max-groups must be positive")
    if not (0.0 < args.grid_step <= 0.5):
        parser.error("--grid-step must be in (0, 0.5]")

    try:
        from scipy.optimize import differential_evolution, minimize
    except ImportError as exc:
        raise SystemExit("scipy is required for local_mode_mean_breaker.py") from exc

    rng = random.Random(args.seed)
    rows: list[dict[str, Any]] = []
    rows.extend(structured_rows(args.max_support, args.grid_step))

    grouped_patterns = []
    for total in range(3, args.max_support + 1):
        grouped_patterns.extend(integer_partitions_at_most(total, args.max_groups))
    for counts in grouped_patterns:
        rows.extend(
            optimize_grouped(
                counts=counts,
                rng=rng,
                restarts=args.group_restarts,
                maxiter=args.slsqp_maxiter,
                scipy_minimize=minimize,
            )
        )

    for n in range(3, min(args.ungrouped_max_n, args.max_support) + 1):
        rows.extend(
            optimize_ungrouped(
                n=n,
                rng=rng,
                restarts=args.ungrouped_restarts,
                maxiter=args.slsqp_maxiter,
                scipy_minimize=minimize,
            )
        )

    rows.extend(
        differential_evolution_rows(
            max_groups=args.max_groups,
            max_support=args.max_support,
            restarts=args.de_restarts,
            maxiter=args.de_maxiter,
            rng=rng,
            differential_evolution=differential_evolution,
            penalty=args.de_penalty,
        )
    )

    feasible_rows = [row for row in rows if row["feasible"]]
    counterexamples = [
        row for row in feasible_rows if row["mean"] < 2.5 - 1e-7
    ]
    near_counterexamples = [
        row for row in feasible_rows if row["mean"] < 2.5 - 1e-9
    ]
    best_by_source = []
    for source in sorted({row["source"] for row in rows}):
        source_feasible = [row for row in feasible_rows if row["source"] == source]
        if source_feasible:
            best_by_source.append(compact(min(source_feasible, key=lambda row: row["mean"])))

    summary = {
        "source": {
            "kind": "local_mode_mean_breaker",
            "seed": args.seed,
            "max_support": args.max_support,
            "max_groups": args.max_groups,
            "group_restarts": args.group_restarts,
            "ungrouped_restarts": args.ungrouped_restarts,
            "ungrouped_max_n": args.ungrouped_max_n,
            "slsqp_maxiter": args.slsqp_maxiter,
            "grid_step": args.grid_step,
            "de_restarts": args.de_restarts,
            "de_maxiter": args.de_maxiter,
            "de_penalty": args.de_penalty,
        },
        "processed": len(rows),
        "feasible": len(feasible_rows),
        "counterexamples": len(counterexamples),
        "near_counterexamples_strict_tolerance": len(near_counterexamples),
        "best_feasible_by_mean": [
            compact(row)
            for row in sorted(feasible_rows, key=lambda row: row["mean"])[: args.top]
        ],
        "best_feasible_by_constraint_margin": [
            compact(row)
            for row in sorted(
                feasible_rows,
                key=lambda row: (abs(row["constraint_margin_e3_minus_e2"]), row["mean"]),
            )[: args.top]
        ],
        "best_by_source": best_by_source,
        "counterexample_candidates": [
            compact(row)
            for row in sorted(counterexamples, key=lambda row: row["mean"])[: args.top]
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    best_mean = min((row["mean"] for row in feasible_rows), default=math.nan)
    print(
        f"Wrote {args.out}; processed={len(rows)}, feasible={len(feasible_rows)}, "
        f"counterexamples={len(counterexamples)}, best_mean={best_mean:.12g}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
