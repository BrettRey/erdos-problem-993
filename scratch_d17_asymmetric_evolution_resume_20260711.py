#!/usr/bin/env python3
"""Crash-safe replay of the interrupted D17 asymmetric evolutionary search.

This preserves the interrupted run's convention: branching tuples are read
from the joined root outward.  The lower-level spherical recurrence reads
leaves toward the root, so every tuple is reversed exactly once at that API
boundary.
"""

from __future__ import annotations

import argparse
import ast
import json
import math
import random
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np

from scratch_d17_checkpointed_search_20260711 import (
    RankedCandidate,
    append_record,
    exact_joined_polynomial,
    exact_valley,
    strict_post_descent_pressure,
    treehood_certificate,
    utc_now,
)
from scratch_d17_spherical_search_20260711 import (
    add,
    multiply,
    spherical_states,
)


Branching = tuple[int, ...]
Pair = tuple[Branching, Branching]


def side_order(root_outward: Branching) -> int:
    level = 1
    total = 1
    for branch in root_outward:
        level *= branch
        total += level
    return total


def pair_order(pair: Pair) -> int:
    return side_order(pair[0]) + side_order(pair[1])


def normalize_pair(pair: Pair) -> Pair:
    return tuple(sorted((tuple(pair[0]), tuple(pair[1]))))  # type: ignore[return-value]


class FloatEvaluator:
    def __init__(self, *, relative_floor: float, descent_tolerance: float, min_separation: int):
        self.relative_floor = relative_floor
        self.descent_tolerance = descent_tolerance
        self.min_separation = min_separation
        self.state_cache: dict[Branching, Any] = {}
        self.pair_cache: dict[Pair, RankedCandidate | None] = {}

    def states(self, root_outward: Branching):
        if root_outward not in self.state_cache:
            self.state_cache[root_outward] = spherical_states(tuple(reversed(root_outward)))[:2]
        return self.state_cache[root_outward]

    def evaluate(self, pair: Pair) -> RankedCandidate | None:
        key = normalize_pair(pair)
        if key in self.pair_cache:
            return self.pair_cache[key]
        left_excluded, left_selected = self.states(key[0])
        right_excluded, right_selected = self.states(key[1])
        full = add(
            multiply(left_excluded, right_excluded),
            add(
                multiply(left_selected, right_excluded),
                multiply(left_excluded, right_selected),
            ),
        )
        probability = np.maximum(full.values, 0.0)
        probability /= probability.sum()
        ranked = strict_post_descent_pressure(
            probability,
            relative_floor=self.relative_floor,
            descent_tolerance=self.descent_tolerance,
            min_separation=self.min_separation,
        )
        if ranked is None:
            self.pair_cache[key] = None
            return None
        candidate = RankedCandidate(
            score=float(ranked["best_later_ratio"]),
            left=key[0],
            right=key[1],
            order=pair_order(key),
            first_descent_at=int(ranked["first_descent_at"]),
            first_descent_ratio=float(ranked["first_descent_ratio"]),
            best_later_at=int(ranked["best_later_at"]),
            best_later_ratio=float(ranked["best_later_ratio"]),
            rebound=float(ranked["rebound"]),
        )
        self.pair_cache[key] = candidate
        return candidate


def seed_population(max_order: int) -> list[Pair]:
    rows: list[Branching] = [
        *((2,) * depth for depth in range(4, 12)),
        (4, 4, 2, 1, 2, 2, 1, 2, 4, 4),
        (2, 3, 2, 1, 2, 1, 1),
        (3, 2, 2, 1, 1, 2, 2),
    ]
    return [
        (left, right)
        for left in rows
        for right in rows
        if side_order(left) + side_order(right) <= max_order
    ]


def random_injection(rng: random.Random, max_order: int) -> Pair | None:
    sides = []
    for _ in range(2):
        depth = rng.randint(4, 13)
        sides.append(tuple(rng.choice((1, 1, 2, 2, 2, 3, 4)) for _ in range(depth)))
    pair: Pair = (sides[0], sides[1])
    return pair if pair_order(pair) <= max_order else None


def mutate(parent: Pair, rng: random.Random, max_order: int) -> Pair | None:
    rows = [list(parent[0]), list(parent[1])]
    side = rng.randrange(2)
    row = rows[side]
    operation = rng.randrange(5)
    if operation == 0:
        row[rng.randrange(len(row))] = rng.randint(1, 6)
    elif operation == 1 and len(row) < 15:
        row.insert(rng.randrange(len(row) + 1), rng.randint(1, 5))
    elif operation == 2 and len(row) > 3:
        row.pop(rng.randrange(len(row)))
    elif operation == 3:
        rng.shuffle(row)
    else:
        j = rng.randrange(len(row))
        row[j] = max(1, min(8, row[j] + rng.choice((-1, 1))))
    child: Pair = (tuple(rows[0]), tuple(rows[1]))
    return child if pair_order(child) <= max_order else None


def load_checkpoint(path: Path, run_id: str) -> dict[str, Any] | None:
    if not path.exists():
        return None
    last = None
    for line in path.read_text(encoding="utf-8").splitlines():
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            continue
        if record.get("run_id") == run_id and record.get("kind") == "evolution_checkpoint":
            last = record
    return last


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run-id", default="asym-evolution-recovered-20260711")
    parser.add_argument("--iterations", type=int, default=5_000)
    parser.add_argument("--max-order", type=int, default=14_000)
    parser.add_argument("--seed", type=int, default=20_260_711)
    parser.add_argument("--checkpoint-every", type=int, default=25)
    parser.add_argument("--population", type=int, default=40)
    parser.add_argument("--parent-pool", type=int, default=20)
    parser.add_argument("--relative-floor", type=float, default=1e-12)
    parser.add_argument("--descent-tolerance", type=float, default=1e-9)
    parser.add_argument("--ascent-tolerance", type=float, default=1e-7)
    parser.add_argument("--min-separation", type=int, default=1)
    args = parser.parse_args()

    evaluator = FloatEvaluator(
        relative_floor=args.relative_floor,
        descent_tolerance=args.descent_tolerance,
        min_separation=args.min_separation,
    )
    checkpoint = load_checkpoint(args.output, args.run_id)
    if checkpoint is None:
        rng = random.Random(args.seed)
        population = seed_population(args.max_order)
        start_iteration = 0
        best: RankedCandidate | None = None
        append_record(
            args.output,
            {
                "kind": "evolution_start",
                "run_id": args.run_id,
                "at": utc_now(),
                "parameters": vars(args) | {"output": str(args.output)},
                "branching_convention": "root outward; reversed once for leaves-to-root recurrence",
                "recovered_objective": "seed 20260711, 5000 iterations, total order <= 14000",
            },
        )
    else:
        rng = random.Random()
        rng.setstate(ast.literal_eval(checkpoint["rng_state_repr"]))
        population = [
            (tuple(pair[0]), tuple(pair[1])) for pair in checkpoint["population"]
        ]
        start_iteration = int(checkpoint["next_iteration"])
        raw_best = checkpoint.get("best")
        best = RankedCandidate(**raw_best) if raw_best is not None else None
        append_record(
            args.output,
            {
                "kind": "evolution_resume",
                "run_id": args.run_id,
                "at": utc_now(),
                "start_iteration": start_iteration,
            },
        )

    started = time.monotonic()
    for iteration in range(start_iteration, args.iterations):
        unique_population = set(map(normalize_pair, population))
        ranked = [
            candidate
            for pair in unique_population
            if (candidate := evaluator.evaluate(pair)) is not None
        ]
        ranked.sort(key=lambda item: (item.score, item.rebound), reverse=True)
        if not ranked:
            raise RuntimeError("evolution population has no rankable candidate")
        population = [(item.left, item.right) for item in ranked[: args.population]]
        champion = ranked[0]
        if best is None or (champion.score, champion.rebound) > (best.score, best.rebound):
            best = champion
            append_record(
                args.output,
                {
                    "kind": "evolution_best",
                    "run_id": args.run_id,
                    "at": utc_now(),
                    "iteration": iteration,
                    "candidate": asdict(best),
                },
            )
            print(json.dumps({"iteration": iteration, "best": asdict(best)}, sort_keys=True), flush=True)

        if champion.score > 1.0 + args.ascent_tolerance:
            # Convert root-outward tuples to the leaves-to-root exact API.
            left = tuple(reversed(champion.left))
            right = tuple(reversed(champion.right))
            exact_started = time.monotonic()
            coefficients = exact_joined_polynomial(left, right)
            valley = exact_valley(coefficients)
            exact_record: dict[str, Any] = {
                "kind": "exact_replay",
                "run_id": args.run_id,
                "at": utc_now(),
                "iteration": iteration,
                "candidate_root_outward": asdict(champion),
                "exact_seconds": time.monotonic() - exact_started,
                "exact_degree": len(coefficients) - 1,
                "exact_valley": valley,
            }
            if valley is not None:
                exact_record["treehood"] = treehood_certificate(left, right)
            append_record(args.output, exact_record)
            if valley is not None:
                print(json.dumps(exact_record, sort_keys=True), flush=True)
                return

        parent = rng.choice(population[: min(args.parent_pool, len(population))])
        child = mutate(parent, rng, args.max_order)
        if child is not None:
            population.append(child)
        if iteration % 30 == 0:
            for _ in range(12):
                injected = random_injection(rng, args.max_order)
                if injected is not None:
                    population.append(injected)

        if (iteration + 1) % args.checkpoint_every == 0 or iteration + 1 == args.iterations:
            record = {
                "kind": "evolution_checkpoint",
                "run_id": args.run_id,
                "at": utc_now(),
                "next_iteration": iteration + 1,
                "elapsed_seconds_this_invocation": time.monotonic() - started,
                "population": population,
                "rng_state_repr": repr(rng.getstate()),
                "best": asdict(best) if best is not None else None,
                "evaluated_pairs": len(evaluator.pair_cache),
                "cached_states": len(evaluator.state_cache),
            }
            append_record(args.output, record)
            print(
                json.dumps(
                    {
                        "checkpoint": iteration + 1,
                        "best": asdict(best) if best is not None else None,
                        "evaluated_pairs": len(evaluator.pair_cache),
                    },
                    sort_keys=True,
                ),
                flush=True,
            )

    append_record(
        args.output,
        {
            "kind": "evolution_complete",
            "run_id": args.run_id,
            "at": utc_now(),
            "iterations": args.iterations,
            "best": asdict(best) if best is not None else None,
        },
    )


if __name__ == "__main__":
    main()
