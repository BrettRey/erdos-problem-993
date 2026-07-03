#!/usr/bin/env python3
"""Heuristic optimizer for effective signed ratio drop."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
from typing import Any

from analyze_signed_conditionals import analyze_row, compact_analysis
from optimize_signed_pb_reserve import (
    State,
    is_feasible,
    mutate,
    parse_floats,
    parse_pairs,
    seed_states,
    state_from_eval,
    state_key,
)
from probe_signed_pb_reserve import signed_metric


def evaluate(state: State) -> dict[str, Any] | None:
    x_blocks, y_blocks = state
    metric = signed_metric(
        x_blocks=x_blocks,
        y_blocks=y_blocks,
        kind="signed_ratio_drop_optimizer",
    )
    if metric is None:
        return None
    return analyze_row(metric, "optimizer_candidate")


def objective(ev: dict[str, Any]) -> float:
    return ev["variance_times_effective_ratio_drop"]


def score(
    ev: dict[str, Any],
    *,
    variance_cutoff: float,
    min_side_variance: float,
    min_side_variance_fraction: float,
) -> float:
    penalty = 0.0
    if ev["variance"] < variance_cutoff:
        penalty += 100.0 + 10.0 * (variance_cutoff - ev["variance"])
    if not is_feasible(
        ev,
        variance_cutoff=min(ev["variance"], variance_cutoff),
        min_side_variance=min_side_variance,
        min_side_variance_fraction=min_side_variance_fraction,
    ):
        required = max(min_side_variance, min_side_variance_fraction * ev["variance"])
        shortfall = max(0.0, required - ev["x_variance"]) + max(
            0.0,
            required - ev["y_variance"],
        )
        penalty += 100.0 + 20.0 * shortfall
    return objective(ev) + penalty


def compact_optimizer(row: dict[str, Any] | None) -> dict[str, Any] | None:
    if row is None:
        return None
    out = compact_analysis(row)
    out["x_blocks"] = row["x_blocks"]
    out["y_blocks"] = row["y_blocks"]
    return out


def run_one(
    *,
    x_n: int,
    y_n: int,
    variance_cutoff: float,
    max_groups: int,
    population_size: int,
    generations: int,
    seed: int,
    min_side_variance: float,
    min_side_variance_fraction: float,
) -> dict[str, Any]:
    rng = random.Random(seed)
    population: list[dict[str, Any]] = []
    for state in seed_states(x_n, y_n, max_groups, rng):
        ev = evaluate(state)
        if ev is None:
            continue
        ev["score"] = score(
            ev,
            variance_cutoff=variance_cutoff,
            min_side_variance=min_side_variance,
            min_side_variance_fraction=min_side_variance_fraction,
        )
        population.append(ev)
    population.sort(key=lambda ev: ev["score"])
    population = population[:population_size]

    history = []
    scale = 1.0
    for generation in range(generations + 1):
        best_feasible = min(
            (
                ev
                for ev in population
                if is_feasible(
                    ev,
                    variance_cutoff=variance_cutoff,
                    min_side_variance=min_side_variance,
                    min_side_variance_fraction=min_side_variance_fraction,
                )
            ),
            key=objective,
            default=None,
        )
        if generation in {0, 1, generations} or generation % 50 == 0:
            history.append(
                {
                    "generation": generation,
                    "best_score": population[0]["score"],
                    "best_feasible": compact_optimizer(best_feasible),
                }
            )
        if generation == generations:
            break

        children: list[dict[str, Any]] = []
        elite = population[: max(4, population_size // 5)]
        for parent in elite:
            parent_state = state_from_eval(parent)
            for _ in range(4):
                ev = evaluate(mutate(parent_state, max_groups, rng, scale))
                if ev is None:
                    continue
                ev["score"] = score(
                    ev,
                    variance_cutoff=variance_cutoff,
                    min_side_variance=min_side_variance,
                    min_side_variance_fraction=min_side_variance_fraction,
                )
                children.append(ev)

        immigrant_count = max(1, population_size // 10)
        immigrant_seeds = seed_states(x_n, y_n, max_groups, rng)
        for state in rng.sample(immigrant_seeds, min(immigrant_count, len(immigrant_seeds))):
            ev = evaluate(state)
            if ev is not None:
                ev["score"] = score(
                    ev,
                    variance_cutoff=variance_cutoff,
                    min_side_variance=min_side_variance,
                    min_side_variance_fraction=min_side_variance_fraction,
                )
                children.append(ev)

        population.extend(children)
        best_by_key: dict[
            tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]],
            dict[str, Any],
        ] = {}
        for ev in population:
            key = state_key(ev)
            incumbent = best_by_key.get(key)
            if incumbent is None or ev["score"] < incumbent["score"]:
                best_by_key[key] = ev
        population = sorted(best_by_key.values(), key=lambda ev: ev["score"])[:population_size]
        scale = max(0.05, scale * 0.995)

    feasible = [
        ev
        for ev in population
        if is_feasible(
            ev,
            variance_cutoff=variance_cutoff,
            min_side_variance=min_side_variance,
            min_side_variance_fraction=min_side_variance_fraction,
        )
    ]
    best = min(feasible, key=objective, default=None)
    return {
        "x_n": x_n,
        "y_n": y_n,
        "variance_cutoff": variance_cutoff,
        "min_side_variance": min_side_variance,
        "min_side_variance_fraction": min_side_variance_fraction,
        "max_groups": max_groups,
        "seed": seed,
        "best": compact_optimizer(best),
        "history": history,
    }


def best_by_cutoff(rows: list[dict[str, Any]], cutoffs: list[float]) -> list[dict[str, Any]]:
    out = []
    for cutoff in cutoffs:
        candidates = [row for row in rows if row["variance"] >= cutoff]
        best = min(candidates, key=lambda row: row["variance_times_effective_ratio_drop"], default=None)
        out.append({"variance_cutoff": cutoff, "best": compact_optimizer(best)})
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-pairs", default="20:20,50:50,100:100,200:200,50:200,200:50")
    parser.add_argument("--variance-cutoffs", default="1,2,5,10,20,50")
    parser.add_argument("--max-groups", type=int, default=6)
    parser.add_argument("--population-size", type=int, default=70)
    parser.add_argument("--generations", type=int, default=180)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--candidate-constant", type=float, default=0.25)
    parser.add_argument("--min-side-variance", type=float, default=0.0)
    parser.add_argument("--min-side-variance-fraction", type=float, default=0.0)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/signed_pb_ratio_drop_optimizer_2026-07-03.json"),
    )
    args = parser.parse_args()
    if args.max_groups < 1:
        parser.error("--max-groups must be positive")
    if args.min_side_variance < 0.0:
        parser.error("--min-side-variance must be nonnegative")
    if not 0.0 <= args.min_side_variance_fraction <= 0.5:
        parser.error("--min-side-variance-fraction must be in [0, 0.5]")

    n_pairs = parse_pairs(args.n_pairs)
    cutoffs = parse_floats(args.variance_cutoffs)
    results = []
    run_index = 0
    for x_n, y_n in n_pairs:
        for cutoff in cutoffs:
            if cutoff >= (x_n + y_n) / 4:
                continue
            results.append(
                run_one(
                    x_n=x_n,
                    y_n=y_n,
                    variance_cutoff=cutoff,
                    max_groups=args.max_groups,
                    population_size=args.population_size,
                    generations=args.generations,
                    seed=args.seed + run_index,
                    min_side_variance=args.min_side_variance,
                    min_side_variance_fraction=args.min_side_variance_fraction,
                )
            )
            run_index += 1

    best_rows = [result["best"] for result in results if result["best"] is not None]
    failures = [
        row
        for row in best_rows
        if row["variance"] >= 1.0
        and row["variance_times_effective_ratio_drop"] < args.candidate_constant
    ]
    summary = {
        "source": {
            "kind": "signed_pb_effective_ratio_drop_optimizer",
            "n_pairs": n_pairs,
            "variance_cutoffs": cutoffs,
            "max_groups": args.max_groups,
            "population_size": args.population_size,
            "generations": args.generations,
            "seed": args.seed,
            "candidate_constant": args.candidate_constant,
            "min_side_variance": args.min_side_variance,
            "min_side_variance_fraction": args.min_side_variance_fraction,
        },
        "processed": len(results),
        "candidate_failures": len(failures),
        "best_by_variance_cutoff": best_by_cutoff(best_rows, cutoffs),
        "best_by_effective_ratio_drop": sorted(
            best_rows,
            key=lambda row: row["variance_times_effective_ratio_drop"],
        )[: args.top],
        "candidate_failure_rows": sorted(
            failures,
            key=lambda row: row["variance_times_effective_ratio_drop"],
        )[: args.top],
        "results": results,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Wrote {args.out}; runs={len(results)}, "
        f"candidate_failures={len(failures)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
