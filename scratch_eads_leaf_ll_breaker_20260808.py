#!/usr/bin/env python3
"""Evolutionary kill test for an EADS orientation-transition lemma.

For a vertex ``v`` put

    A_v = I(T-v),   B_v = x I(T-N[v]).

A split is bad when the modal intervals of ``A_v`` and ``B_v`` are at
distance at least two.  It is left-oriented when the modes of ``A_v`` lie
to the right of those of ``B_v``.  Exact enumeration through order 18 found
no tree satisfying both

1. every leaf split is bad; and
2. some edge has two left-oriented bad endpoints.

This program is a falsification search for that statement, not a proof.  It
keeps exact integer polynomials and records a complete certificate if it
finds a counterexample.  The search is seeded on both sides of the apparent
barrier: the exceptional order-18 tree with every leaf split bad, and double
stars, whose centers can both be left-bad while their leaves remain good.
"""

from __future__ import annotations

import argparse
import json
import random
import time
from copy import deepcopy
from pathlib import Path

from nm_optimizer import _random_tree, generate_seeds, mutate
from scratch_eads_local_constraint_breaker_20260808 import (
    encode_split,
    fingerprint,
    signed_orientation,
)
from scratch_existential_split_mode_20260807 import tree_vertex_splits


ORDER_18_ALL_LEAVES_BAD_EDGES = (
    (0, 14), (0, 15), (1, 14), (1, 16), (2, 15), (2, 17),
    (3, 16), (4, 16), (5, 16), (6, 16), (7, 16),
    (8, 17), (9, 17), (10, 17), (11, 17), (12, 17), (13, 17),
)

# Exact all-vertex-bad order-60 obstruction found on 2026-08-08.  It lies
# directly on the conjectured alternating-orientation boundary and is a much
# stronger seed than the older order-18 all-leaves-bad example.
ORDER_60_ALL_VERTICES_BAD_EDGES = (
    (0, 1),
    *((0, v) for v in range(14, 19)),
    (1, 2), (1, 3), (1, 4), (1, 5), (1, 7), (1, 9), (1, 10),
    *((2, v) for v in range(19, 24)),
    *((3, v) for v in range(24, 28)),
    (4, 12),
    *((4, v) for v in range(28, 31)),
    (5, 6),
    *((5, v) for v in range(31, 34)),
    (6, 8), (6, 11), (6, 13),
    *((7, v) for v in range(34, 39)),
    *((8, v) for v in range(39, 43)),
    *((9, v) for v in range(43, 47)),
    *((10, v) for v in range(47, 51)),
    *((11, v) for v in range(51, 55)),
    *((13, v) for v in range(55, 60)),
)


def from_edges(n: int, edges: tuple[tuple[int, int], ...]) -> list[list[int]]:
    adj: list[list[int]] = [[] for _ in range(n)]
    for left, right in edges:
        adj[left].append(right)
        adj[right].append(left)
    return adj


def double_star(n: int, left_leaves: int) -> list[list[int]]:
    """Return the double star with adjacent centers 0 and 1."""
    if n < 4 or not 1 <= left_leaves <= n - 3:
        raise ValueError("a double star needs at least one leaf at each center")
    adj: list[list[int]] = [[] for _ in range(n)]
    adj[0].append(1)
    adj[1].append(0)
    for vertex in range(2, n):
        center = 0 if vertex < 2 + left_leaves else 1
        adj[center].append(vertex)
        adj[vertex].append(center)
    return adj


def evaluate(adj: list[list[int]]) -> dict[str, object]:
    full, splits = tree_vertex_splits(adj)
    leaves = [v for v, row in enumerate(adj) if len(row) == 1]
    leaf_floor = min(splits[v].distance for v in leaves)

    # This signed margin is a smoother search coordinate than the categorical
    # orientation: it reaches 2 exactly when a vertex is left-bad.
    left_margins = [
        split.a_modes[0] - split.b_modes[1]
        for split in splits
    ]
    left_edge_pressure, pressure_edge = max(
        (
            (min(left_margins[u], left_margins[v]), (u, v))
            for u, row in enumerate(adj)
            for v in row
            if u < v
        ),
        default=(-10**9, None),
    )
    target = leaf_floor >= 2 and left_edge_pressure >= 2
    return {
        "target": target,
        "leaf_floor": leaf_floor,
        "left_edge_pressure": left_edge_pressure,
        "pressure_edge": pressure_edge,
        "polynomial": list(full),
        "splits": splits,
        "left_margins": left_margins,
    }


def encode_candidate(
    adj: list[list[int]], evaluation: dict[str, object], generation: int,
) -> dict[str, object]:
    splits = evaluation["splits"]
    assert isinstance(splits, list)
    return {
        "n": len(adj),
        "generation": generation,
        "target": evaluation["target"],
        "leaf_floor": evaluation["leaf_floor"],
        "left_edge_pressure": evaluation["left_edge_pressure"],
        "pressure_edge": list(evaluation["pressure_edge"])
        if evaluation["pressure_edge"] is not None else None,
        "edges": [list(edge) for edge in fingerprint(adj)],
        "degrees": [len(row) for row in adj],
        "polynomial": evaluation["polynomial"],
        "splits": [encode_split(split) for split in splits],
    }


def rank_key(item: tuple[list[list[int]], dict[str, object]]) -> tuple[int, ...]:
    _, ev = item
    leaf = int(ev["leaf_floor"])
    edge = int(ev["left_edge_pressure"])
    return min(leaf, edge), leaf + edge, leaf, edge


def select_diverse(
    candidates: list[tuple[list[list[int]], dict[str, object]]],
    population_size: int,
    rng: random.Random,
) -> list[tuple[list[list[int]], dict[str, object]]]:
    """Retain balanced, leaf-side, and left-edge-side search niches."""
    unique: dict[tuple[tuple[int, int], ...], tuple[list[list[int]], dict[str, object]]] = {}
    for candidate in candidates:
        unique.setdefault(fingerprint(candidate[0]), candidate)
    pool = list(unique.values())
    selected: list[tuple[list[list[int]], dict[str, object]]] = []
    selected_ids: set[tuple[tuple[int, int], ...]] = set()

    orderings = [
        sorted(pool, key=rank_key, reverse=True),
        sorted(
            pool,
            key=lambda item: (
                item[1]["leaf_floor"], item[1]["left_edge_pressure"],
            ),
            reverse=True,
        ),
        sorted(
            pool,
            key=lambda item: (
                item[1]["left_edge_pressure"], item[1]["leaf_floor"],
            ),
            reverse=True,
        ),
    ]
    cursor = [0, 0, 0]
    lane = 0
    while len(selected) < min(population_size, len(pool)):
        ordering = orderings[lane]
        while cursor[lane] < len(ordering):
            candidate = ordering[cursor[lane]]
            cursor[lane] += 1
            fp = fingerprint(candidate[0])
            if fp not in selected_ids:
                selected_ids.add(fp)
                selected.append(candidate)
                break
        lane = (lane + 1) % len(orderings)
        if all(cursor[j] >= len(orderings[j]) for j in range(3)):
            break
    rng.shuffle(selected)
    return selected


def initial_population(
    args: argparse.Namespace, rng: random.Random,
) -> list[tuple[list[list[int]], dict[str, object]]]:
    trees: list[list[list[int]]] = [
        from_edges(18, ORDER_18_ALL_LEAVES_BAD_EDGES),
        from_edges(60, ORDER_60_ALL_VERTICES_BAD_EDGES),
    ]
    while len(trees) < args.population:
        n = rng.randint(args.min_n, args.max_n)
        lane = len(trees) % 3
        if lane == 0 and n >= 4:
            trees.append(double_star(n, rng.randint(1, n - 3)))
        elif lane == 1:
            _, adj, _ = generate_seeds(n, 1, rng)[0]
            trees.append(adj if len(adj) == n else _random_tree(n, rng))
        else:
            trees.append(_random_tree(n, rng))
    return [(adj, evaluate(adj)) for adj in trees]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=14)
    parser.add_argument("--max-n", type=int, default=220)
    parser.add_argument("--population", type=int, default=120)
    parser.add_argument("--generations", type=int, default=1000)
    parser.add_argument("--offspring", type=int, default=7)
    parser.add_argument("--mutations", type=int, default=4)
    parser.add_argument("--seed", type=int, default=202608085)
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/eads_leaf_ll_breaker_20260808.json"),
    )
    args = parser.parse_args()
    rng = random.Random(args.seed)
    started = time.time()
    population = initial_population(args, rng)
    evaluated = len(population)
    archive: dict[
        tuple[tuple[int, int], ...],
        tuple[list[list[int]], dict[str, object], int],
    ] = {}
    certificate = None

    for generation in range(args.generations + 1):
        for adj, ev in population:
            fp = fingerprint(adj)
            archive.setdefault(fp, (deepcopy(adj), ev, generation))
            if ev["target"]:
                certificate = encode_candidate(adj, ev, generation)
                print(json.dumps({"event": "counterexample", **certificate}), flush=True)
                break
        if certificate is not None:
            break

        if generation % 10 == 0:
            balanced = max(population, key=rank_key)
            best_leaf = max(
                population,
                key=lambda item: (
                    item[1]["leaf_floor"], item[1]["left_edge_pressure"],
                ),
            )
            best_edge = max(
                population,
                key=lambda item: (
                    item[1]["left_edge_pressure"], item[1]["leaf_floor"],
                ),
            )
            print(json.dumps({
                "event": "generation",
                "generation": generation,
                "evaluated": evaluated,
                "balanced": rank_key(balanced),
                "best_leaf": [
                    best_leaf[1]["leaf_floor"],
                    best_leaf[1]["left_edge_pressure"],
                    len(best_leaf[0]),
                ],
                "best_edge": [
                    best_edge[1]["leaf_floor"],
                    best_edge[1]["left_edge_pressure"],
                    len(best_edge[0]),
                ],
            }), flush=True)
        if generation == args.generations:
            break

        parents = select_diverse(population, max(6, args.population // 4), rng)
        candidates = list(parents)
        attempts = 0
        target_size = args.population * args.offspring
        while len(candidates) < target_size and attempts < target_size * 12:
            attempts += 1
            parent_adj, _ = rng.choice(parents)
            child = deepcopy(parent_adj)
            for _ in range(rng.randint(1, args.mutations)):
                _, child = mutate(len(child), child, rng)
            candidates.append((child, evaluate(child)))
            evaluated += 1
        population = select_diverse(candidates, args.population, rng)

    best_rows = sorted(archive.values(), key=lambda row: rank_key((row[0], row[1])), reverse=True)[:10]
    report = {
        "claim_scope": "evolutionary exact falsification search; not a proof",
        "statement_tested": (
            "if every leaf split is bad, no edge has two left-oriented "
            "bad endpoints"
        ),
        "seed": args.seed,
        "parameters": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "population": args.population,
            "generations": args.generations,
            "offspring": args.offspring,
            "mutations": args.mutations,
        },
        "evaluated": evaluated,
        "counterexample": certificate,
        "best": [encode_candidate(adj, ev, gen) for adj, ev, gen in best_rows],
        "elapsed_seconds": time.time() - started,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "evaluated": evaluated,
        "counterexample": certificate is not None,
        "elapsed_seconds": report["elapsed_seconds"],
        "out": str(args.out),
    }), flush=True)
    raise SystemExit(1 if certificate is not None else 0)


if __name__ == "__main__":
    main()
