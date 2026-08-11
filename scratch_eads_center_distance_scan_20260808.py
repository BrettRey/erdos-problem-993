#!/usr/bin/env python3
"""Exhaustively locate the nearest EADS certificate to a tree center.

This is a diagnostic for canonical-center strengthenings of EADS, not a test
of EADS itself.  For every non-isomorphic tree in the requested order range,
the script computes all exact deletion splits, finds every EADS-good vertex,
and records its minimum graph distance from the center of the tree.
"""

from __future__ import annotations

import argparse
import json
import time
from collections import Counter, deque
from pathlib import Path

from scratch_eads_center_ball_optimizer_20260808 import centers
from scratch_existential_split_mode_20260807 import (
    encode_split,
    tree_vertex_splits,
)
from trees import trees_geng_raw


def distances_from(adj: list[list[int]], sources: set[int]) -> list[int]:
    distance = [-1] * len(adj)
    queue = deque(sources)
    for source in sources:
        distance[source] = 0
    while queue:
        vertex = queue.popleft()
        for neighbor in adj[vertex]:
            if distance[neighbor] < 0:
                distance[neighbor] = distance[vertex] + 1
                queue.append(neighbor)
    return distance


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=2)
    parser.add_argument("--max-n", type=int, default=18)
    parser.add_argument("--record-distance", type=int, default=2)
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/eads_center_distance_n18_20260808.json"),
    )
    args = parser.parse_args()
    started = time.time()
    orders: list[dict[str, object]] = []
    witnesses: list[dict[str, object]] = []
    total = 0
    global_histogram: Counter[int] = Counter()

    for n in range(args.min_n, args.max_n + 1):
        histogram: Counter[int] = Counter()
        order_trees = 0
        for tree_index, (_, adj, graph6) in enumerate(trees_geng_raw(n)):
            polynomial, splits = tree_vertex_splits(adj)
            center = centers(adj)
            distance = distances_from(adj, center)
            good = [split.vertex for split in splits if split.distance <= 1]
            nearest = min(distance[vertex] for vertex in good)
            histogram[nearest] += 1
            global_histogram[nearest] += 1
            order_trees += 1
            total += 1
            if nearest >= args.record_distance:
                witnesses.append({
                    "n": n,
                    "tree_index": tree_index,
                    "graph6": graph6.decode("ascii"),
                    "centers": sorted(center),
                    "nearest_good_distance": nearest,
                    "good_vertices": good,
                    "edges": [
                        [left, right]
                        for left, neighbors in enumerate(adj)
                        for right in neighbors
                        if left < right
                    ],
                    "independence_polynomial": list(polynomial),
                    "splits": [encode_split(split, True) for split in splits],
                })
        row = {
            "n": n,
            "trees": order_trees,
            "nearest_good_distance_histogram": {
                str(key): value for key, value in sorted(histogram.items())
            },
            "elapsed_seconds": time.time() - started,
        }
        orders.append(row)
        print(json.dumps({"event": "order", **row}), flush=True)

    result = {
        "claim_scope": "exhaustive exact computation; diagnostic only",
        "statement_tested": (
            "Distance from a graph-theoretic center to the nearest vertex "
            "whose EADS modal-interval distance is at most one."
        ),
        "min_n": args.min_n,
        "max_n": args.max_n,
        "total_trees": total,
        "nearest_good_distance_histogram": {
            str(key): value for key, value in sorted(global_histogram.items())
        },
        "record_distance": args.record_distance,
        "witness_count": len(witnesses),
        "witnesses": witnesses,
        "orders": orders,
        "elapsed_seconds": time.time() - started,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "trees": total,
        "witnesses": len(witnesses),
        "out": str(args.out),
    }), flush=True)


if __name__ == "__main__":
    main()
