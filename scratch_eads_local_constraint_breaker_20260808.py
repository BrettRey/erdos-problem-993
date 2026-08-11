#!/usr/bin/env python3
"""Evolutionary falsification search for two candidate EADS local lemmas.

For a vertex ``v`` put

    A_v = I(T-v),   B_v = x I(T-N[v]).

A bad split has modal distance at least two.  Call it right-oriented when
the modes of ``B_v`` lie to the right of those of ``A_v``, and left-oriented
otherwise.  Exact small-tree data suggest, but do not prove, that

1. no edge has two right-oriented bad endpoints; and
2. every left-oriented bad vertex has degree at least four.

The second assertion is closely related to the open vertex-deletion mode
stability problem, so both assertions must be treated as conjectural.  This
script searches directly for either kind of exact counterexample.  It uses
integer independence-polynomial arithmetic throughout and records the full
tree and split certificate when it succeeds.
"""

from __future__ import annotations

import argparse
import json
import random
import time
from copy import deepcopy
from pathlib import Path

from nm_optimizer import _random_tree, generate_seeds, mutate
from scratch_existential_split_mode_20260807 import (
    VertexSplit,
    tree_vertex_splits,
)


def signed_orientation(split: VertexSplit) -> tuple[str, int]:
    """Return orientation and signed modal separation."""
    if split.b_modes[0] > split.a_modes[1]:
        return "right", split.b_modes[0] - split.a_modes[1]
    if split.a_modes[0] > split.b_modes[1]:
        return "left", split.a_modes[0] - split.b_modes[1]
    return "overlap", 0


def evaluate(adj: list[list[int]]) -> dict[str, object]:
    full, splits = tree_vertex_splits(adj)
    oriented = [signed_orientation(split) for split in splits]

    left_margins = [
        split.a_modes[0] - split.b_modes[1]
        for split in splits
    ]
    right_margins = [
        split.b_modes[0] - split.a_modes[1]
        for split in splits
    ]
    low_left = max(
        (
            left_margins[v]
            for v in range(len(adj))
            if len(adj[v]) <= 3
        ),
        default=-1,
    )
    right_edge = max(
        (
            min(right_margins[u], right_margins[v])
            for u, neighbors in enumerate(adj)
            for v in neighbors
            if u < v
        ),
        default=-1,
    )
    objective = max(low_left, right_edge)
    return {
        "objective": objective,
        "low_degree_left_separation": low_left,
        "right_edge_separation": right_edge,
        "polynomial": list(full),
        "splits": splits,
        "oriented": oriented,
    }


def fingerprint(adj: list[list[int]]) -> tuple[tuple[int, int], ...]:
    return tuple(
        (u, v)
        for u, neighbors in enumerate(adj)
        for v in neighbors
        if u < v
    )


def encode_split(split: VertexSplit) -> dict[str, object]:
    orientation, separation = signed_orientation(split)
    return {
        "vertex": split.vertex,
        "orientation": orientation,
        "separation": separation,
        "a_modes": list(split.a_modes),
        "b_modes": list(split.b_modes),
        "a_poly": list(split.a_poly),
        "b_poly": list(split.b_poly),
    }


def encode_candidate(
    adj: list[list[int]], evaluation: dict[str, object], generation: int,
) -> dict[str, object]:
    splits = evaluation["splits"]
    assert isinstance(splits, list)
    return {
        "n": len(adj),
        "generation": generation,
        "objective": evaluation["objective"],
        "low_degree_left_separation": evaluation[
            "low_degree_left_separation"
        ],
        "right_edge_separation": evaluation["right_edge_separation"],
        "edges": [list(edge) for edge in fingerprint(adj)],
        "degrees": [len(neighbors) for neighbors in adj],
        "polynomial": evaluation["polynomial"],
        "splits": [encode_split(split) for split in splits],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=28)
    parser.add_argument("--max-n", type=int, default=180)
    parser.add_argument("--population", type=int, default=80)
    parser.add_argument("--generations", type=int, default=600)
    parser.add_argument("--offspring", type=int, default=5)
    parser.add_argument("--mutations", type=int, default=3)
    parser.add_argument("--seed", type=int, default=202608084)
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/eads_local_constraint_breaker_20260808.json"),
    )
    args = parser.parse_args()
    rng = random.Random(args.seed)
    started = time.time()

    population: list[tuple[list[list[int]], dict[str, object]]] = []
    certificate_path = Path("results/eads_counterexample_n60_20260808.json")
    if certificate_path.exists():
        row = json.loads(certificate_path.read_text())
        seeded = [[] for _ in range(int(row["order"]))]
        for left, right in row["edges"]:
            seeded[left].append(right)
            seeded[right].append(left)
        population.append((seeded, evaluate(seeded)))
    sizes = [
        rng.randint(args.min_n, args.max_n)
        for _ in range(args.population - len(population))
    ]
    for n in sizes:
        # Half structured seeds and half uniform labelled trees.
        if len(population) % 2 == 0:
            _, adj, _ = generate_seeds(n, 1, rng)[0]
            if len(adj) != n:
                adj = _random_tree(n, rng)
        else:
            adj = _random_tree(n, rng)
        population.append((adj, evaluate(adj)))

    archive: dict[tuple[tuple[int, int], ...], tuple[list[list[int]], dict[str, object], int]] = {}
    certificate = None
    evaluated = len(population)

    for generation in range(args.generations + 1):
        for adj, ev in population:
            fp = fingerprint(adj)
            incumbent = archive.get(fp)
            if incumbent is None or ev["objective"] > incumbent[1]["objective"]:
                archive[fp] = (deepcopy(adj), ev, generation)
            if ev["objective"] >= 2:
                certificate = encode_candidate(adj, ev, generation)
                break
        if certificate is not None:
            print(json.dumps({"event": "counterexample", **certificate}), flush=True)
            break

        population.sort(
            key=lambda item: (
                item[1]["objective"],
                item[1]["low_degree_left_separation"],
                item[1]["right_edge_separation"],
            ),
            reverse=True,
        )
        if generation % 10 == 0:
            best = population[0][1]
            print(json.dumps({
                "event": "generation",
                "generation": generation,
                "evaluated": evaluated,
                "best_objective": best["objective"],
                "best_low_degree_left": best["low_degree_left_separation"],
                "best_right_edge": best["right_edge_separation"],
                "best_n": len(population[0][0]),
            }), flush=True)
        if generation == args.generations:
            break

        elite_count = max(4, args.population // 5)
        parents = population[:elite_count]
        candidates = list(parents)
        seen = {fingerprint(adj) for adj, _ in candidates}
        attempts = 0
        target = args.population * args.offspring
        while len(candidates) < target and attempts < target * 10:
            attempts += 1
            parent_adj, _ = rng.choice(parents)
            child = deepcopy(parent_adj)
            for _ in range(rng.randint(1, args.mutations)):
                _, child = mutate(len(child), child, rng)
            fp = fingerprint(child)
            if fp in seen:
                continue
            seen.add(fp)
            candidates.append((child, evaluate(child)))
            evaluated += 1
        candidates.sort(
            key=lambda item: (
                item[1]["objective"],
                item[1]["low_degree_left_separation"],
                item[1]["right_edge_separation"],
                rng.random(),
            ),
            reverse=True,
        )
        population = candidates[:args.population]

    best_rows = sorted(
        archive.values(),
        key=lambda row: row[1]["objective"],
        reverse=True,
    )[:10]
    report = {
        "claim_scope": "evolutionary exact falsification search; not a proof",
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
