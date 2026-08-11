#!/usr/bin/env python3
"""Adversarial test of a canonical center-ball version of EADS.

For a tree T let C(T) be its one- or two-vertex graph-theoretic center.  For a
chosen radius r, this script searches for a tree on which every vertex at
distance at most r from C(T) has modal-split distance at least two.  Such a
tree refutes the proposed strengthening

    some v at distance at most r from C(T) is an EADS certificate.

The calculation is exact.  A failure of the strengthening is not necessarily
a failure of EADS itself; the latter permits a certificate anywhere in T.
"""

from __future__ import annotations

import argparse
import json
import random
from collections import deque
from pathlib import Path

from nm_optimizer import _random_tree, mutate_leaf_move, mutate_spr
from scratch_existential_split_mode_20260807 import tree_vertex_splits


def centers(adj: list[list[int]]) -> set[int]:
    """Return the graph-theoretic center(s) of a tree."""
    degree = [len(neighbors) for neighbors in adj]
    queue = deque(v for v, value in enumerate(degree) if value <= 1)
    remaining = len(adj)
    while remaining > 2:
        for _ in range(len(queue)):
            vertex = queue.popleft()
            remaining -= 1
            for neighbor in adj[vertex]:
                degree[neighbor] -= 1
                if degree[neighbor] == 1:
                    queue.append(neighbor)
    return set(queue)


def evaluate(adj: list[list[int]], radius: int) -> dict[str, object]:
    polynomial, splits = tree_vertex_splits(adj)
    center = centers(adj)
    ball = set(center)
    frontier = set(center)
    for _ in range(radius):
        frontier = {
            neighbor
            for vertex in frontier
            for neighbor in adj[vertex]
        } - ball
        ball.update(frontier)
    distance = {split.vertex: split.distance for split in splits}
    ball_profile = sorted(distance[vertex] for vertex in ball)
    good = sorted(split.vertex for split in splits if split.distance <= 1)
    ball_good = sorted(vertex for vertex in ball if distance[vertex] <= 1)
    return {
        "centers": sorted(center),
        "center_ball": sorted(ball),
        "center_ball_profile": ball_profile,
        "center_ball_good": ball_good,
        "good_vertices": good,
        "distances": [split.distance for split in splits],
        "polynomial": list(polynomial),
    }


def key(evaluation: dict[str, object]) -> tuple[object, ...]:
    profile = evaluation["center_ball_profile"]
    assert isinstance(profile, list)
    # First eliminate certificates from the center ball.  At a tie, maximize
    # its sorted distance profile and then reduce certificates globally.
    return (
        len(evaluation["center_ball_good"]),  # type: ignore[arg-type]
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
    parser.add_argument("--radius", type=int, default=1)
    parser.add_argument("--seed", type=int, default=20260808)
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/eads_center_ball_optimizer_20260808.json"),
    )
    args = parser.parse_args()
    rng = random.Random(args.seed)
    seed_paths = [
        Path("results/eads_seeded_optimizer_20260807.json"),
        Path("results/eads_seeded_n60_optimizer_20260807.json"),
        Path("results/eads_center_ball_optimizer_smoke_20260808.json"),
        Path("results/eads_center_ball_r2_optimizer_smoke_20260808.json"),
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
            for _ in range(1 + restart % 9):
                if rng.random() < 0.55:
                    _, adj = mutate_leaf_move(len(adj), adj, rng)
                else:
                    _, adj = mutate_spr(len(adj), adj, rng)
        else:
            adj = _random_tree(rng.randint(args.min_n, args.max_n), rng)
        n = len(adj)
        evaluation = evaluate(adj, args.radius)
        evaluated += 1

        for step in range(args.steps):
            if rng.random() < 0.55:
                _, candidate = mutate_leaf_move(n, adj, rng)
            else:
                _, candidate = mutate_spr(n, adj, rng)
            candidate_evaluation = evaluate(candidate, args.radius)
            evaluated += 1
            if key(candidate_evaluation) <= key(evaluation) or rng.random() < 0.01:
                adj, evaluation = candidate, candidate_evaluation
            if best_evaluation is None or key(evaluation) < key(best_evaluation):
                best_adj = [neighbors.copy() for neighbors in adj]
                best_evaluation = evaluation.copy()
                best_restart = restart
                best_step = step
            if not evaluation["center_ball_good"]:
                break

        assert best_evaluation is not None
        print(json.dumps({
            "event": "restart",
            "restart": restart,
            "evaluated": evaluated,
            "best_center_ball_good": len(best_evaluation["center_ball_good"]),
            "best_profile": best_evaluation["center_ball_profile"],
            "best_total_good": len(best_evaluation["good_vertices"]),
        }), flush=True)
        if not best_evaluation["center_ball_good"]:
            break

    assert best_adj is not None and best_evaluation is not None
    result = {
        "claim_scope": "exact evolutionary computation; evidence only",
        "statement_tested": (
            f"Some vertex at distance at most {args.radius} from the center "
            "of every tree has EADS modal-interval distance at most one."
        ),
        "seed": args.seed,
        "radius": args.radius,
        "evaluated": evaluated,
        "counterexample_found": not best_evaluation["center_ball_good"],
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
