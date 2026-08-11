#!/usr/bin/env python3
"""Exhaustively test whether every nontrivial tree has an EADS-good leaf."""

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
    parser.add_argument("--max-n", type=int, default=18)
    parser.add_argument("--stop-on-first", action="store_true")
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/eads_leaf_n18_20260808.json"),
    )
    args = parser.parse_args()
    started = time.time()
    orders: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    total_trees = 0
    total_leaf_splits = 0
    global_min_good_histogram: Counter[int] = Counter()
    global_total_good_histogram: Counter[int] = Counter()
    global_min_total_good: int | None = None
    global_min_total_records: list[dict[str, object]] = []

    for n in range(args.min_n, args.max_n + 1):
        order_trees = 0
        leaf_splits = 0
        failure_count = 0
        min_good_histogram: Counter[int] = Counter()
        total_good_histogram: Counter[int] = Counter()
        order_min_total_good: int | None = None
        for tree_index, (_, adj, graph6) in enumerate(trees_geng_raw(n)):
            polynomial, splits = tree_vertex_splits(adj)
            leaves = [
                vertex for vertex, neighbors in enumerate(adj)
                if len(neighbors) == 1
            ]
            leaf_split_records = [splits[leaf] for leaf in leaves]
            best = min(split.distance for split in leaf_split_records)
            total_good = sum(split.distance <= 1 for split in splits)
            min_good_histogram[best] += 1
            global_min_good_histogram[best] += 1
            total_good_histogram[total_good] += 1
            global_total_good_histogram[total_good] += 1
            order_min_total_good = (
                total_good if order_min_total_good is None
                else min(order_min_total_good, total_good)
            )
            if global_min_total_good is None or total_good < global_min_total_good:
                global_min_total_good = total_good
                global_min_total_records = []
            if total_good == global_min_total_good and len(global_min_total_records) < 10:
                global_min_total_records.append({
                    "n": n,
                    "tree_index": tree_index,
                    "graph6": graph6.decode("ascii"),
                    "total_good": total_good,
                    "good_vertices": [
                        split.vertex for split in splits if split.distance <= 1
                    ],
                    "edges": [
                        [left, right]
                        for left, neighbors in enumerate(adj)
                        for right in neighbors
                        if left < right
                    ],
                })
            order_trees += 1
            leaf_splits += len(leaves)
            if best > 1:
                failure_count += 1
                failures.append({
                    "n": n,
                    "tree_index": tree_index,
                    "graph6": graph6.decode("ascii"),
                    "leaves": leaves,
                    "best_leaf_distance": best,
                    "edges": [
                        [left, right]
                        for left, neighbors in enumerate(adj)
                        for right in neighbors
                        if left < right
                    ],
                    "independence_polynomial": list(polynomial),
                    "leaf_splits": [
                        encode_split(split, True)
                        for split in leaf_split_records
                    ],
                })
                if args.stop_on_first:
                    break
        total_trees += order_trees
        total_leaf_splits += leaf_splits
        row = {
            "n": n,
            "trees": order_trees,
            "leaf_splits": leaf_splits,
            "failures": failure_count,
            "best_leaf_distance_histogram": {
                str(key): value
                for key, value in sorted(min_good_histogram.items())
            },
            "total_good_histogram": {
                str(key): value
                for key, value in sorted(total_good_histogram.items())
            },
            "min_total_good": order_min_total_good,
            "elapsed_seconds": time.time() - started,
        }
        orders.append(row)
        print(json.dumps({"event": "order", **row}), flush=True)
        if args.stop_on_first and failure_count:
            break

    result = {
        "claim_scope": "exhaustive exact computation; evidence only",
        "statement_tested": (
            "Every nontrivial tree has a leaf whose I(T-v) and "
            "xI(T-N[v]) modal intervals have distance at most one."
        ),
        "min_n": args.min_n,
        "max_n": args.max_n,
        "total_trees": total_trees,
        "total_leaf_splits": total_leaf_splits,
        "best_leaf_distance_histogram": {
            str(key): value
            for key, value in sorted(global_min_good_histogram.items())
        },
        "total_good_histogram": {
            str(key): value
            for key, value in sorted(global_total_good_histogram.items())
        },
        "min_total_good": global_min_total_good,
        "min_total_good_records": global_min_total_records,
        "failure_count": len(failures),
        "failures": failures,
        "orders": orders,
        "elapsed_seconds": time.time() - started,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "trees": total_trees,
        "leaf_splits": total_leaf_splits,
        "failures": len(failures),
        "out": str(args.out),
    }), flush=True)
    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()
