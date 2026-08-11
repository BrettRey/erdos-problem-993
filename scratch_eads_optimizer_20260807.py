#!/usr/bin/env python3
"""Evolutionary kill test for the existential adjacent-mode split (EADS).

For every vertex v of a tree T, compare the modal intervals of

    A_v = I(T-v),             B_v = x I(T-N[v]).

EADS asserts that some vertex has interval distance at most one.  This
search directly minimizes the number of such vertices, using exact integer
polynomials throughout.  It is a counterexample search, not a proof.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from nm_optimizer import _adj_fingerprint, _random_tree, mutate_leaf_move, mutate_spr
from scratch_existential_split_mode_20260807 import tree_vertex_splits


def caterpillar_seed(loads: tuple[int, ...]) -> list[list[int]]:
    """Build a caterpillar used to seed the low-good basin."""
    n = len(loads) + sum(loads)
    adj: list[list[int]] = [[] for _ in range(n)]
    for u in range(len(loads) - 1):
        adj[u].append(u + 1)
        adj[u + 1].append(u)
    leaf = len(loads)
    for support, count in enumerate(loads):
        for _ in range(count):
            adj[support].append(leaf)
            adj[leaf].append(support)
            leaf += 1
    return adj


def evaluate(adj: list[list[int]]) -> dict[str, object]:
    poly, splits = tree_vertex_splits(adj)
    distances = [split.distance for split in splits]
    good = [split.vertex for split in splits if split.distance <= 1]
    # Lexicographic pressure: first eliminate good vertices, then push their
    # distances upward, then push every split upward.
    ordered = sorted(distances)
    return {
        "good_count": len(good),
        "good_vertices": good,
        "distances": distances,
        "ordered_distances": ordered,
        "polynomial": list(poly),
    }


def key(ev: dict[str, object]) -> tuple[object, ...]:
    ordered = ev["ordered_distances"]
    assert isinstance(ordered, list)
    # A larger distance vector is better after minimizing the number <= 1.
    return (ev["good_count"], tuple(-value for value in ordered))


def encode(adj: list[list[int]], ev: dict[str, object], restart: int, step: int) -> dict[str, object]:
    return {
        "restart": restart,
        "step": step,
        "n": len(adj),
        "edges": [
            [u, v]
            for u, neighbors in enumerate(adj)
            for v in neighbors
            if u < v
        ],
        **ev,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=24)
    parser.add_argument("--max-n", type=int, default=100)
    parser.add_argument("--restarts", type=int, default=60)
    parser.add_argument("--steps", type=int, default=500)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument(
        "--seed-loads", default="3,0,4,1,5,6",
        help="comma-separated caterpillar loads; empty disables the seed",
    )
    parser.add_argument(
        "--best-restarts", action="store_true",
        help="start later restarts from the best tree instead of a random tree",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/eads_optimizer_20260807.json"),
    )
    args = parser.parse_args()
    rng = random.Random(args.seed)

    best_adj: list[list[int]] | None = None
    best_ev: dict[str, object] | None = None
    best_restart = -1
    best_step = -1
    evaluated = 0
    seed_adj = None
    if args.seed_loads.strip():
        seed_adj = caterpillar_seed(tuple(
            int(value) for value in args.seed_loads.split(",")
        ))

    for restart in range(args.restarts):
        if restart == 0 and seed_adj is not None:
            adj = [neighbors.copy() for neighbors in seed_adj]
            n = len(adj)
        elif args.best_restarts and best_adj is not None:
            adj = [neighbors.copy() for neighbors in best_adj]
            n = len(adj)
            for _ in range(1 + restart % 7):
                if rng.random() < 0.55:
                    _, adj = mutate_leaf_move(n, adj, rng)
                else:
                    _, adj = mutate_spr(n, adj, rng)
        else:
            n = rng.randint(args.min_n, args.max_n)
            adj = _random_tree(n, rng)
        ev = evaluate(adj)
        evaluated += 1
        for step in range(args.steps):
            if rng.random() < 0.55:
                _, candidate = mutate_leaf_move(n, adj, rng)
            else:
                _, candidate = mutate_spr(n, adj, rng)
            candidate_ev = evaluate(candidate)
            evaluated += 1
            if key(candidate_ev) <= key(ev) or rng.random() < 0.015:
                adj, ev = candidate, candidate_ev
            if best_ev is None or key(ev) < key(best_ev):
                best_adj = [neighbors.copy() for neighbors in adj]
                best_ev = ev.copy()
                best_restart = restart
                best_step = step
            if ev["good_count"] == 0:
                break
        if best_ev is not None:
            print(json.dumps({
                "restart": restart,
                "evaluated": evaluated,
                "best_good_count": best_ev["good_count"],
                "best_distance_prefix": best_ev["ordered_distances"][:12],
            }), flush=True)
        if best_ev is not None and best_ev["good_count"] == 0:
            break

    assert best_adj is not None and best_ev is not None
    result = {
        "claim_scope": "exact evolutionary computation; evidence only",
        "statement_tested": (
            "Every tree has a vertex v whose I(T-v) and xI(T-N[v]) "
            "modal intervals have distance at most one."
        ),
        "seed": args.seed,
        "evaluated": evaluated,
        "counterexample_found": best_ev["good_count"] == 0,
        "best": encode(best_adj, best_ev, best_restart, best_step),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "evaluated": evaluated,
        "counterexample_found": result["counterexample_found"],
        "best_good_count": best_ev["good_count"],
        "out": str(args.out),
    }), flush=True)
    raise SystemExit(1 if result["counterexample_found"] else 0)


if __name__ == "__main__":
    main()
