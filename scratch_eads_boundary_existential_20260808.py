#!/usr/bin/env python3
"""Exact scan for an EADS-good vertex near the boundary of every tree.

For a vertex v, its deletion split is

    I(T) = I(T-v) + x I(T-N[v]).

Call v good when the modal intervals of the two summands have distance at
most one.  This script tests two structural statements:

1. some leaf or support of a leaf is good;
2. among the endpoints and supports of one deterministic diameter, some
   vertex is good.

Both statements are evidence-only.  A minimal-counterexample proof would
also need to justify the forest polynomials produced by the deletions.
"""

from __future__ import annotations

import argparse
import json
import time
from collections import Counter, deque
from pathlib import Path

from scratch_existential_split_mode_20260807 import (
    encode_split,
    tree_vertex_splits,
)
from trees import trees_geng_raw


def farthest(adj: list[list[int]], start: int) -> int:
    """Return the least-labelled farthest vertex from start."""
    distance = [-1] * len(adj)
    distance[start] = 0
    queue = deque([start])
    while queue:
        vertex = queue.popleft()
        for neighbor in adj[vertex]:
            if distance[neighbor] >= 0:
                continue
            distance[neighbor] = distance[vertex] + 1
            queue.append(neighbor)
    maximum = max(distance)
    return next(v for v, value in enumerate(distance) if value == maximum)


def diameter_boundary(adj: list[list[int]]) -> set[int]:
    if len(adj) <= 1:
        return {0} if adj else set()
    left = farthest(adj, 0)
    right = farthest(adj, left)
    return {left, adj[left][0], right, adj[right][0]}


def encoded_tree(
    order: int,
    tree_index: int,
    graph6: bytes,
    adj: list[list[int]],
    polynomial: tuple[int, ...] | list[int],
    splits,
    candidates: set[int],
) -> dict[str, object]:
    return {
        "n": order,
        "tree_index": tree_index,
        "graph6": graph6.decode("ascii"),
        "edges": [
            [left, right]
            for left, row in enumerate(adj)
            for right in row
            if left < right
        ],
        "independence_polynomial": list(polynomial),
        "candidates": sorted(candidates),
        "candidate_splits": {
            str(vertex): encode_split(splits[vertex], True)
            for vertex in sorted(candidates)
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=2)
    parser.add_argument("--max-n", type=int, default=20)
    parser.add_argument("--stop-on-first", action="store_true")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/eads_boundary_existential_n20_20260808.json"),
    )
    args = parser.parse_args()

    started = time.time()
    orders: list[dict[str, object]] = []
    boundary_failures: list[dict[str, object]] = []
    diameter_failures: list[dict[str, object]] = []
    boundary_good_histogram: Counter[int] = Counter()
    total_trees = 0

    for order in range(args.min_n, args.max_n + 1):
        order_trees = 0
        order_boundary_failures = 0
        order_diameter_failures = 0
        order_histogram: Counter[int] = Counter()
        minimum_good: dict[str, object] | None = None
        minimum_good_count: int | None = None

        for tree_index, (_, adj, graph6) in enumerate(trees_geng_raw(order)):
            polynomial, splits = tree_vertex_splits(adj)
            leaves = {v for v, row in enumerate(adj) if len(row) == 1}
            boundary = leaves | {adj[leaf][0] for leaf in leaves}
            good = {v for v in boundary if splits[v].distance <= 1}
            good_count = len(good)
            order_histogram[good_count] += 1
            boundary_good_histogram[good_count] += 1

            if minimum_good_count is None or good_count < minimum_good_count:
                minimum_good_count = good_count
                minimum_good = {
                    **encoded_tree(
                        order, tree_index, graph6, adj, polynomial, splits,
                        boundary,
                    ),
                    "good_boundary_vertices": sorted(good),
                }

            if not good:
                order_boundary_failures += 1
                if len(boundary_failures) < 20:
                    boundary_failures.append(minimum_good)

            diameter = diameter_boundary(adj)
            diameter_good = {
                v for v in diameter if splits[v].distance <= 1
            }
            if not diameter_good:
                order_diameter_failures += 1
                if len(diameter_failures) < 20:
                    diameter_failures.append({
                        **encoded_tree(
                            order, tree_index, graph6, adj, polynomial,
                            splits, diameter,
                        ),
                        "good_diameter_boundary_vertices": [],
                    })

            order_trees += 1
            if args.stop_on_first and (
                order_boundary_failures or order_diameter_failures
            ):
                break

        total_trees += order_trees
        row = {
            "n": order,
            "trees": order_trees,
            "boundary_failures": order_boundary_failures,
            "diameter_boundary_failures": order_diameter_failures,
            "boundary_good_histogram": {
                str(count): frequency
                for count, frequency in sorted(order_histogram.items())
            },
            "minimum_good_count": minimum_good_count,
            "minimum_good_witness": minimum_good,
            "elapsed_seconds": time.time() - started,
        }
        orders.append(row)
        print(json.dumps({"event": "order", **row}), flush=True)
        if args.stop_on_first and (
            order_boundary_failures or order_diameter_failures
        ):
            break

    result = {
        "claim_scope": "exhaustive exact computation; evidence only",
        "min_n": args.min_n,
        "max_n": args.max_n,
        "total_trees": total_trees,
        "boundary_failure_count": sum(
            row["boundary_failures"] for row in orders
        ),
        "diameter_boundary_failure_count": sum(
            row["diameter_boundary_failures"] for row in orders
        ),
        "boundary_good_histogram": {
            str(count): frequency
            for count, frequency in sorted(boundary_good_histogram.items())
        },
        "boundary_failures": boundary_failures,
        "diameter_boundary_failures": diameter_failures,
        "orders": orders,
        "elapsed_seconds": time.time() - started,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "trees": total_trees,
        "boundary_failures": result["boundary_failure_count"],
        "diameter_boundary_failures": result[
            "diameter_boundary_failure_count"
        ],
        "out": str(args.out),
    }), flush=True)
    raise SystemExit(
        1 if result["boundary_failure_count"]
        or result["diameter_boundary_failure_count"] else 0
    )


if __name__ == "__main__":
    main()
