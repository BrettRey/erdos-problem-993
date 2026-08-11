#!/usr/bin/env python3
"""Exact kill-test for leaf/support complementarity in EADS.

For every leaf ell with support u, test whether ell or u has an adjacent-mode
deletion split.  This local statement is a structural kill-test only: even if
it held, using it in an induction would still require control of the forest
polynomials created by vertex deletion.  This script provides finite exact
evidence only.
"""

from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path

from scratch_existential_split_mode_20260807 import (
    encode_split,
    tree_vertex_splits,
)
from trees import trees_geng_raw


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=2)
    parser.add_argument("--max-n", type=int, default=20)
    parser.add_argument("--stop-on-first", action="store_true")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/eads_pendant_pair_n20_20260808.json"),
    )
    args = parser.parse_args()

    started = time.time()
    total_trees = 0
    total_pairs = 0
    failures: list[dict[str, object]] = []
    orders: list[dict[str, object]] = []
    distance_pairs: Counter[tuple[int, int]] = Counter()

    for order in range(args.min_n, args.max_n + 1):
        order_trees = 0
        order_pairs = 0
        order_failures = 0
        order_distances: Counter[tuple[int, int]] = Counter()
        for tree_index, (_, adj, graph6) in enumerate(trees_geng_raw(order)):
            polynomial, splits = tree_vertex_splits(adj)
            for leaf, neighbors in enumerate(adj):
                if len(neighbors) != 1:
                    continue
                support = neighbors[0]
                pair = (splits[leaf].distance, splits[support].distance)
                order_distances[pair] += 1
                distance_pairs[pair] += 1
                order_pairs += 1
                if min(pair) <= 1:
                    continue
                order_failures += 1
                failures.append({
                    "n": order,
                    "tree_index": tree_index,
                    "graph6": graph6.decode("ascii"),
                    "leaf": leaf,
                    "support": support,
                    "distance_pair": list(pair),
                    "edges": [
                        [left, right]
                        for left, row in enumerate(adj)
                        for right in row
                        if left < right
                    ],
                    "independence_polynomial": list(polynomial),
                    "leaf_split": encode_split(splits[leaf], True),
                    "support_split": encode_split(splits[support], True),
                })
                if args.stop_on_first:
                    break
            order_trees += 1
            if args.stop_on_first and order_failures:
                break

        total_trees += order_trees
        total_pairs += order_pairs
        row = {
            "n": order,
            "trees": order_trees,
            "pendant_pairs": order_pairs,
            "failures": order_failures,
            "distance_pairs": {
                f"{left},{right}": count
                for (left, right), count in sorted(order_distances.items())
            },
            "elapsed_seconds": time.time() - started,
        }
        orders.append(row)
        print(json.dumps({"event": "order", **row}), flush=True)
        if args.stop_on_first and order_failures:
            break

    result = {
        "claim_scope": "exhaustive exact computation; evidence only",
        "statement_tested": (
            "For every leaf ell with support u, at least one of ell and u "
            "has deletion-split modal distance at most one."
        ),
        "min_n": args.min_n,
        "max_n": args.max_n,
        "total_trees": total_trees,
        "total_pendant_pairs": total_pairs,
        "failure_count": len(failures),
        "distance_pairs": {
            f"{left},{right}": count
            for (left, right), count in sorted(distance_pairs.items())
        },
        "failures": failures,
        "orders": orders,
        "elapsed_seconds": time.time() - started,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "trees": total_trees,
        "pendant_pairs": total_pairs,
        "failures": len(failures),
        "out": str(args.out),
    }), flush=True)
    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()
