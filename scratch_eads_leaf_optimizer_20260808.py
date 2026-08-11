#!/usr/bin/env python3
"""Adversarial test of the existential-leaf strengthening of EADS.

The statement under test is that every nontrivial tree has a leaf v for which

    A_v = I(T-v),              B_v = x I(T-N[v])

have modal intervals at distance at most one.  This implies EADS and has the
special leaf form A_v=I(H), B_v=xI(H-u), where H=T-v and u is its support.
All split polynomials and all objective comparisons are computed exactly.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from nm_optimizer import _random_tree, mutate_leaf_move, mutate_spr
from scratch_existential_split_mode_20260807 import tree_vertex_splits


def evaluate(adj: list[list[int]]) -> dict[str, object]:
    polynomial, splits = tree_vertex_splits(adj)
    leaves = [vertex for vertex, neighbors in enumerate(adj) if len(neighbors) == 1]
    distance = {split.vertex: split.distance for split in splits}
    leaf_profile = sorted(distance[leaf] for leaf in leaves)
    good_leaves = sorted(leaf for leaf in leaves if distance[leaf] <= 1)
    good_vertices = sorted(split.vertex for split in splits if split.distance <= 1)
    return {
        "leaves": leaves,
        "leaf_profile": leaf_profile,
        "good_leaves": good_leaves,
        "good_vertices": good_vertices,
        "distances": [split.distance for split in splits],
        "polynomial": list(polynomial),
    }


def key(evaluation: dict[str, object]) -> tuple[object, ...]:
    profile = evaluation["leaf_profile"]
    assert isinstance(profile, list)
    return (
        len(evaluation["good_leaves"]),  # type: ignore[arg-type]
        tuple(-int(value) for value in profile),
        len(evaluation["good_vertices"]),  # type: ignore[arg-type]
    )


def load_seed(path: Path) -> list[list[int]] | None:
    if not path.exists():
        return None
    data = json.loads(path.read_text())
    record = data.get("best")
    if not isinstance(record, dict) or "edges" not in record:
        return None
    adj: list[list[int]] = [[] for _ in range(int(record["n"]))]
    for left, right in record["edges"]:
        adj[left].append(right)
        adj[right].append(left)
    return adj


def encode(
    restart: int,
    step: int,
    adj: list[list[int]],
    evaluation: dict[str, object],
) -> dict[str, object]:
    return {
        "restart": restart,
        "step": step,
        "n": len(adj),
        "edges": [
            [left, right]
            for left, neighbors in enumerate(adj)
            for right in neighbors
            if left < right
        ],
        **evaluation,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=20)
    parser.add_argument("--max-n", type=int, default=120)
    parser.add_argument("--restarts", type=int, default=100)
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260808)
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/eads_leaf_optimizer_20260808.json"),
    )
    args = parser.parse_args()
    rng = random.Random(args.seed)
    seed_paths = [
        Path("results/eads_seeded_optimizer_20260807.json"),
        Path("results/eads_seeded_n60_optimizer_20260807.json"),
        Path("results/eads_center_ball_r3_optimizer_smoke_20260808.json"),
    ]
    seeds = [seed for path in seed_paths if (seed := load_seed(path))]

    best_adj: list[list[int]] | None = None
    best_evaluation: dict[str, object] | None = None
    best_restart = -1
    best_step = -1
    evaluated = 0

    for restart in range(args.restarts):
        if restart < len(seeds):
            adj = [neighbors.copy() for neighbors in seeds[restart]]
        elif best_adj is not None and restart % 3 == 0:
            adj = [neighbors.copy() for neighbors in best_adj]
            for _ in range(1 + restart % 11):
                if rng.random() < 0.55:
                    _, adj = mutate_leaf_move(len(adj), adj, rng)
                else:
                    _, adj = mutate_spr(len(adj), adj, rng)
        else:
            adj = _random_tree(rng.randint(args.min_n, args.max_n), rng)
        n = len(adj)
        evaluation = evaluate(adj)
        evaluated += 1

        for step in range(args.steps):
            if rng.random() < 0.55:
                _, candidate = mutate_leaf_move(n, adj, rng)
            else:
                _, candidate = mutate_spr(n, adj, rng)
            candidate_evaluation = evaluate(candidate)
            evaluated += 1
            if key(candidate_evaluation) <= key(evaluation) or rng.random() < 0.01:
                adj, evaluation = candidate, candidate_evaluation
            if best_evaluation is None or key(evaluation) < key(best_evaluation):
                best_adj = [neighbors.copy() for neighbors in adj]
                best_evaluation = evaluation.copy()
                best_restart = restart
                best_step = step
            if not evaluation["good_leaves"]:
                break

        assert best_evaluation is not None
        print(json.dumps({
            "event": "restart",
            "restart": restart,
            "evaluated": evaluated,
            "best_good_leaves": len(best_evaluation["good_leaves"]),
            "best_leaf_profile": best_evaluation["leaf_profile"],
            "best_total_good": len(best_evaluation["good_vertices"]),
        }), flush=True)
        if not best_evaluation["good_leaves"]:
            break

    assert best_adj is not None and best_evaluation is not None
    result = {
        "claim_scope": "exact evolutionary computation; evidence only",
        "statement_tested": (
            "Every nontrivial tree has a leaf whose I(T-v) and "
            "xI(T-N[v]) modal intervals have distance at most one."
        ),
        "seed": args.seed,
        "evaluated": evaluated,
        "counterexample_found": not best_evaluation["good_leaves"],
        "best": encode(
            best_restart, best_step, best_adj, best_evaluation,
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "evaluated": evaluated,
        "counterexample_found": result["counterexample_found"],
        "out": str(args.out),
    }), flush=True)
    raise SystemExit(1 if result["counterexample_found"] else 0)


if __name__ == "__main__":
    main()
