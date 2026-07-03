#!/usr/bin/env python3
"""Heuristic optimizer for small Poisson-binomial variance reserve.

The search space is grouped Bernoulli sums: a state is a list of blocks
``(count, p)`` representing ``count`` independent Bernoulli variables with the
same success probability.  This covers binomial laws and many heterogeneous
Poisson-binomial laws while keeping convolution cheap.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from copy import deepcopy
from pathlib import Path
from typing import Any

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover - environment guard
    raise SystemExit("This optimizer requires numpy.") from exc


Block = tuple[int, float]


def first_descent(seq: np.ndarray) -> int | None:
    for k in range(1, len(seq)):
        if seq[k] < seq[k - 1]:
            return k
    return None


def binomial_shape(n: int, p: float) -> np.ndarray:
    if n == 0:
        return np.array([1.0], dtype=np.float64)
    if p <= 0.0:
        out = np.zeros(n + 1, dtype=np.float64)
        out[0] = 1.0
        return out
    if p >= 1.0:
        out = np.zeros(n + 1, dtype=np.float64)
        out[n] = 1.0
        return out

    mode_guess = min(max(int((n + 1) * p), 0), n)
    pmf = np.zeros(n + 1, dtype=np.float64)
    pmf[mode_guess] = 1.0
    odds = p / (1.0 - p)
    for k in range(mode_guess, n):
        pmf[k + 1] = pmf[k] * (n - k) / (k + 1) * odds
    for k in range(mode_guess, 0, -1):
        pmf[k - 1] = pmf[k] * k / (n - k + 1) / odds
    return pmf / float(pmf.sum())


def grouped_pmf(blocks: list[Block]) -> np.ndarray:
    pmf = np.array([1.0], dtype=np.float64)
    for count, p in blocks:
        if count <= 0:
            continue
        pmf = np.convolve(pmf, binomial_shape(count, p))
        pmf /= float(pmf.sum())
    return pmf


def evaluate(blocks: list[Block]) -> dict[str, Any] | None:
    n = sum(count for count, _ in blocks)
    mean = sum(count * p for count, p in blocks)
    variance = sum(count * p * (1.0 - p) for count, p in blocks)
    if n < 2 or variance <= 0:
        return None
    pmf = grouped_pmf(blocks)
    descent = first_descent(pmf)
    if descent is None or descent >= len(pmf) - 1 or pmf[descent] <= 0:
        return None
    pressure = float(pmf[descent + 1] / pmf[descent])
    reserve = 1.0 - pressure
    return {
        "n": n,
        "blocks": normalize_blocks(blocks),
        "mean": mean,
        "variance": variance,
        "first_descent": descent,
        "crossing_pressure": pressure,
        "crossing_reserve": reserve,
        "variance_times_reserve": variance * reserve,
    }


def normalize_blocks(blocks: list[Block]) -> list[Block]:
    clean = [(count, min(max(p, 1e-6), 1.0 - 1e-6)) for count, p in blocks if count > 0]
    clean.sort(key=lambda item: item[1])
    merged: list[Block] = []
    for count, p in clean:
        if merged and abs(merged[-1][1] - p) < 1e-9:
            old_count, old_p = merged[-1]
            merged[-1] = (old_count + count, old_p)
        else:
            merged.append((count, p))
    return merged


def compact(ev: dict[str, Any]) -> dict[str, Any]:
    return {
        "n": ev["n"],
        "blocks": [[count, p] for count, p in ev["blocks"]],
        "mean": ev["mean"],
        "variance": ev["variance"],
        "first_descent": ev["first_descent"],
        "crossing_pressure": ev["crossing_pressure"],
        "crossing_reserve": ev["crossing_reserve"],
        "variance_times_reserve": ev["variance_times_reserve"],
    }


def seed_blocks(n: int, groups: int, rng: random.Random) -> list[list[Block]]:
    seeds: list[list[Block]] = []
    for p in [0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.95, 0.98]:
        seeds.append([(n, p)])
    for _ in range(30):
        cuts = sorted(rng.sample(range(1, n), min(groups - 1, n - 1)))
        counts = [cuts[0], *[cuts[i] - cuts[i - 1] for i in range(1, len(cuts))], n - cuts[-1]]
        probs = sorted(rng.random() for _ in counts)
        seeds.append(list(zip(counts, probs)))
    for _ in range(30):
        cuts = sorted(rng.sample(range(1, n), min(groups - 1, n - 1)))
        counts = [cuts[0], *[cuts[i] - cuts[i - 1] for i in range(1, len(cuts))], n - cuts[-1]]
        probs = sorted(rng.betavariate(0.4, 0.4) for _ in counts)
        seeds.append(list(zip(counts, probs)))
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
    logit = math.log(p / (1.0 - p))
    logit += rng.gauss(0.0, scale)
    p2 = 1.0 / (1.0 + math.exp(-logit))
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
        delta = rng.gauss(0.0, 0.5)
        logit = math.log(p / (1.0 - p))
        p1 = 1.0 / (1.0 + math.exp(-(logit - delta)))
        p2 = 1.0 / (1.0 + math.exp(-(logit + delta)))
        blocks.extend([(left, p1), (count - left, p2)])
        return normalize_blocks(blocks)

    if len(blocks) >= 2:
        i = rng.randrange(len(blocks) - 1)
        c1, p1 = blocks.pop(i)
        c2, p2 = blocks.pop(i)
        blocks.append((c1 + c2, (c1 * p1 + c2 * p2) / (c1 + c2)))
    return normalize_blocks(blocks)


def mutate(blocks: list[Block], max_groups: int, rng: random.Random, scale: float) -> list[Block]:
    roll = rng.random()
    if roll < 0.65:
        return mutate_probs(blocks, rng, scale)
    if roll < 0.85:
        return mutate_counts(blocks, rng)
    return split_or_merge(blocks, max_groups, rng)


def score(ev: dict[str, Any], variance_cutoff: float) -> float:
    penalty = 0.0
    if ev["variance"] < variance_cutoff:
        penalty = (variance_cutoff - ev["variance"]) * 10.0 + 100.0
    return ev["variance_times_reserve"] + penalty


def run_one(
    *,
    n: int,
    variance_cutoff: float,
    max_groups: int,
    population_size: int,
    generations: int,
    seed: int,
) -> dict[str, Any]:
    rng = random.Random(seed)
    population: list[dict[str, Any]] = []
    for blocks in seed_blocks(n, max_groups, rng):
        ev = evaluate(blocks)
        if ev is None:
            continue
        ev["score"] = score(ev, variance_cutoff)
        population.append(ev)
    population.sort(key=lambda ev: ev["score"])
    population = population[:population_size]

    history = []
    scale = 1.0
    for generation in range(generations + 1):
        best_feasible = min(
            (ev for ev in population if ev["variance"] >= variance_cutoff),
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
            for _ in range(4):
                child_blocks = mutate(parent["blocks"], max_groups, rng, scale)
                ev = evaluate(child_blocks)
                if ev is None:
                    continue
                ev["score"] = score(ev, variance_cutoff)
                children.append(ev)
        immigrant_count = max(1, population_size // 10)
        immigrant_seeds = seed_blocks(n, max_groups, rng)
        for blocks in rng.sample(immigrant_seeds, min(immigrant_count, len(immigrant_seeds))):
            ev = evaluate(blocks)
            if ev is not None:
                ev["score"] = score(ev, variance_cutoff)
                children.append(ev)
        population.extend(children)

        best_by_key: dict[tuple[tuple[int, int], ...], dict[str, Any]] = {}
        for ev in population:
            key = tuple((count, round(p * 10**9)) for count, p in ev["blocks"])
            incumbent = best_by_key.get(key)
            if incumbent is None or ev["score"] < incumbent["score"]:
                best_by_key[key] = ev
        population = sorted(best_by_key.values(), key=lambda ev: ev["score"])[:population_size]
        scale = max(0.05, scale * 0.995)

    feasible = [ev for ev in population if ev["variance"] >= variance_cutoff]
    best = min(feasible, key=lambda ev: ev["variance_times_reserve"], default=None)
    return {
        "n": n,
        "variance_cutoff": variance_cutoff,
        "max_groups": max_groups,
        "seed": seed,
        "best": compact(best) if best else None,
        "history": history,
    }


def parse_ints(spec: str) -> list[int]:
    out = set()
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a_s, b_s = part.split("-", 1)
            a = int(a_s)
            b = int(b_s)
            out.update(range(a, b + 1))
        else:
            out.add(int(part))
    return sorted(out)


def parse_floats(spec: str) -> list[float]:
    return [float(part.strip()) for part in spec.split(",") if part.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-values", default="20,50,100,200,300")
    parser.add_argument("--variance-cutoffs", default="1,2,5,10,20,50")
    parser.add_argument("--max-groups", type=int, default=6)
    parser.add_argument("--population-size", type=int, default=80)
    parser.add_argument("--generations", type=int, default=250)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/pb_variance_reserve_optimizer_2026-07-03.json"),
    )
    args = parser.parse_args()
    if args.max_groups < 1:
        parser.error("--max-groups must be positive")

    results = []
    for i, n in enumerate(parse_ints(args.n_values)):
        for j, cutoff in enumerate(parse_floats(args.variance_cutoffs)):
            if cutoff >= n / 4:
                continue
            results.append(
                run_one(
                    n=n,
                    variance_cutoff=cutoff,
                    max_groups=args.max_groups,
                    population_size=args.population_size,
                    generations=args.generations,
                    seed=args.seed + 1000 * i + j,
                )
            )

    best_rows = [result["best"] for result in results if result["best"] is not None]
    summary = {
        "source": {
            "kind": "pb_variance_reserve_grouped_optimizer",
            "n_values": parse_ints(args.n_values),
            "variance_cutoffs": parse_floats(args.variance_cutoffs),
            "max_groups": args.max_groups,
            "population_size": args.population_size,
            "generations": args.generations,
            "seed": args.seed,
        },
        "processed": len(results),
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
