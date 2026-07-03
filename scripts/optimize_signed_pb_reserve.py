#!/usr/bin/env python3
"""Heuristic optimizer for signed low-probability PB reserve.

The search space is a pair of grouped low-probability Bernoulli sums

    X = sum_g Bin(x_g, p_g),    Y = sum_h Bin(y_h, q_h),

with all probabilities at most 1/2.  The objective is to minimize
``V * reserve`` for the signed law ``X - Y`` at its first descent.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from copy import deepcopy
from pathlib import Path
from typing import Any

from probe_signed_pb_reserve import Block, compact, signed_metric


State = tuple[list[Block], list[Block]]


def normalize_blocks(blocks: list[Block]) -> list[Block]:
    clean = [(count, min(max(p, 1e-6), 0.5)) for count, p in blocks if count > 0]
    clean.sort(key=lambda item: item[1])
    merged: list[Block] = []
    for count, p in clean:
        if merged and abs(merged[-1][1] - p) < 1e-9:
            old_count, old_p = merged[-1]
            merged[-1] = (old_count + count, old_p)
        else:
            merged.append((count, p))
    return merged


def evaluate(state: State) -> dict[str, Any] | None:
    x_blocks, y_blocks = state
    return signed_metric(
        x_blocks=normalize_blocks(x_blocks),
        y_blocks=normalize_blocks(y_blocks),
        kind="signed_grouped_optimizer",
    )


def side_variance_requirement(
    ev: dict[str, Any],
    *,
    min_side_variance: float,
    min_side_variance_fraction: float,
) -> float:
    return max(min_side_variance, min_side_variance_fraction * ev["variance"])


def side_variance_shortfall(
    ev: dict[str, Any],
    *,
    min_side_variance: float,
    min_side_variance_fraction: float,
) -> float:
    required = side_variance_requirement(
        ev,
        min_side_variance=min_side_variance,
        min_side_variance_fraction=min_side_variance_fraction,
    )
    return max(0.0, required - ev["x_variance"]) + max(0.0, required - ev["y_variance"])


def is_feasible(
    ev: dict[str, Any],
    *,
    variance_cutoff: float,
    min_side_variance: float,
    min_side_variance_fraction: float,
) -> bool:
    return (
        ev["variance"] >= variance_cutoff
        and side_variance_shortfall(
            ev,
            min_side_variance=min_side_variance,
            min_side_variance_fraction=min_side_variance_fraction,
        )
        <= 0.0
    )


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
    shortfall = side_variance_shortfall(
        ev,
        min_side_variance=min_side_variance,
        min_side_variance_fraction=min_side_variance_fraction,
    )
    if shortfall > 0.0:
        penalty += 100.0 + 20.0 * shortfall
    return ev["variance_times_reserve"] + penalty


def split_counts(n: int, groups: int, rng: random.Random) -> list[int]:
    groups = min(groups, n)
    if groups <= 1:
        return [n]
    cuts = sorted(rng.sample(range(1, n), groups - 1))
    return [cuts[0], *[cuts[i] - cuts[i - 1] for i in range(1, len(cuts))], n - cuts[-1]]


def random_blocks(n: int, max_groups: int, rng: random.Random, family: str) -> list[Block]:
    counts = split_counts(n, rng.randint(1, min(max_groups, n)), rng)
    probs = []
    for _ in counts:
        if family == "sparse":
            probs.append(0.5 * rng.betavariate(0.35, 4.0))
        elif family == "middle":
            probs.append(0.5 * rng.betavariate(2.0, 2.0))
        elif family == "edge":
            probs.append(0.5 * rng.betavariate(0.4, 0.4))
        else:
            probs.append(0.5 * rng.random())
    return normalize_blocks(list(zip(counts, probs)))


def seed_states(x_n: int, y_n: int, max_groups: int, rng: random.Random) -> list[State]:
    seeds: list[State] = []
    base_probs = [0.005, 0.02, 0.05, 0.1, 0.2, 0.35, 0.49, 0.5]
    for p in base_probs:
        for q in base_probs:
            seeds.append(([(x_n, p)], [(y_n, q)]))

    families = ["uniform", "sparse", "middle", "edge"]
    for _ in range(40):
        seeds.append(
            (
                random_blocks(x_n, max_groups, rng, rng.choice(families)),
                random_blocks(y_n, max_groups, rng, rng.choice(families)),
            )
        )
    for _ in range(20):
        seeds.append(
            (
                random_blocks(x_n, max_groups, rng, "sparse"),
                random_blocks(y_n, max_groups, rng, "sparse"),
            )
        )
    return seeds


def mutate_counts(blocks: list[Block], rng: random.Random) -> list[Block]:
    blocks = deepcopy(blocks)
    if len(blocks) < 2:
        return blocks
    i, j = rng.sample(range(len(blocks)), 2)
    if blocks[i][0] <= 1:
        return blocks
    move = rng.randint(1, min(blocks[i][0] - 1, max(1, blocks[i][0] // 5)))
    blocks[i] = (blocks[i][0] - move, blocks[i][1])
    blocks[j] = (blocks[j][0] + move, blocks[j][1])
    return normalize_blocks(blocks)


def mutate_probs(blocks: list[Block], rng: random.Random, scale: float) -> list[Block]:
    blocks = deepcopy(blocks)
    i = rng.randrange(len(blocks))
    count, p = blocks[i]
    p = min(max(p, 1e-6), 0.5 - 1e-9)
    logit = math.log(p / (0.5 - p))
    logit += rng.gauss(0.0, scale)
    p2 = 0.5 / (1.0 + math.exp(-logit))
    blocks[i] = (count, p2)
    return normalize_blocks(blocks)


def split_or_merge(blocks: list[Block], max_groups: int, rng: random.Random) -> list[Block]:
    blocks = deepcopy(blocks)
    if len(blocks) < max_groups and (len(blocks) == 1 or rng.random() < 0.7):
        choices = [i for i, (count, _) in enumerate(blocks) if count >= 2]
        if not choices:
            return blocks
        i = rng.choice(choices)
        count, p = blocks.pop(i)
        left = rng.randint(1, count - 1)
        p = min(max(p, 1e-6), 0.5 - 1e-9)
        logit = math.log(p / (0.5 - p))
        delta = rng.gauss(0.0, 0.7)
        p1 = 0.5 / (1.0 + math.exp(-(logit - delta)))
        p2 = 0.5 / (1.0 + math.exp(-(logit + delta)))
        blocks.extend([(left, p1), (count - left, p2)])
        return normalize_blocks(blocks)

    if len(blocks) >= 2:
        i = rng.randrange(len(blocks) - 1)
        c1, p1 = blocks.pop(i)
        c2, p2 = blocks.pop(i)
        blocks.append((c1 + c2, (c1 * p1 + c2 * p2) / (c1 + c2)))
    return normalize_blocks(blocks)


def mutate_side(blocks: list[Block], max_groups: int, rng: random.Random, scale: float) -> list[Block]:
    roll = rng.random()
    if roll < 0.65:
        return mutate_probs(blocks, rng, scale)
    if roll < 0.85:
        return mutate_counts(blocks, rng)
    return split_or_merge(blocks, max_groups, rng)


def mutate(state: State, max_groups: int, rng: random.Random, scale: float) -> State:
    x_blocks, y_blocks = state
    if rng.random() < 0.5:
        return mutate_side(x_blocks, max_groups, rng, scale), y_blocks
    return x_blocks, mutate_side(y_blocks, max_groups, rng, scale)


def state_from_eval(ev: dict[str, Any]) -> State:
    return (
        [(int(count), float(p)) for count, p in ev["x_blocks"]],
        [(int(count), float(p)) for count, p in ev["y_blocks"]],
    )


def state_key(ev: dict[str, Any]) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]:
    return (
        tuple((count, round(p * 10**9)) for count, p in ev["x_blocks"]),
        tuple((count, round(p * 10**9)) for count, p in ev["y_blocks"]),
    )


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
            key=lambda ev: ev["variance_times_reserve"],
            default=None,
        )
        if generation in {0, 1, generations} or generation % 50 == 0:
            history.append(
                {
                    "generation": generation,
                    "best_score": population[0]["score"],
                    "best_feasible": compact(best_feasible) if best_feasible else None,
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
    best = min(feasible, key=lambda ev: ev["variance_times_reserve"], default=None)
    return {
        "x_n": x_n,
        "y_n": y_n,
        "variance_cutoff": variance_cutoff,
        "min_side_variance": min_side_variance,
        "min_side_variance_fraction": min_side_variance_fraction,
        "max_groups": max_groups,
        "seed": seed,
        "best": compact(best) if best else None,
        "history": history,
    }


def parse_pairs(spec: str) -> list[tuple[int, int]]:
    out = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" not in part:
            value = int(part)
            out.append((value, value))
            continue
        left, right = part.split(":", 1)
        out.append((int(left), int(right)))
    return out


def parse_floats(spec: str) -> list[float]:
    return [float(part.strip()) for part in spec.split(",") if part.strip()]


def best_by_cutoff(rows: list[dict[str, Any]], cutoffs: list[float]) -> list[dict[str, Any]]:
    out = []
    for cutoff in cutoffs:
        candidates = [row for row in rows if row["variance"] >= cutoff]
        best = min(candidates, key=lambda row: row["variance_times_reserve"], default=None)
        out.append({"variance_cutoff": cutoff, "best": compact(best) if best else None})
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
    parser.add_argument("--min-side-variance", type=float, default=0.0)
    parser.add_argument("--min-side-variance-fraction", type=float, default=0.0)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/signed_pb_reserve_optimizer_2026-07-03.json"),
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
    summary = {
        "source": {
            "kind": "signed_pb_variance_reserve_grouped_optimizer",
            "n_pairs": n_pairs,
            "variance_cutoffs": cutoffs,
            "max_groups": args.max_groups,
            "population_size": args.population_size,
            "generations": args.generations,
            "seed": args.seed,
            "min_side_variance": args.min_side_variance,
            "min_side_variance_fraction": args.min_side_variance_fraction,
        },
        "processed": len(results),
        "best_by_variance_cutoff": best_by_cutoff(best_rows, cutoffs),
        "best_by_variance_times_reserve": sorted(
            best_rows, key=lambda row: row["variance_times_reserve"]
        )[: args.top],
        "results": results,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {args.out}; runs={len(results)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
