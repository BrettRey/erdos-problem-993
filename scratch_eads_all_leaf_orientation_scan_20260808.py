#!/usr/bin/env python3
"""Exact scan of bad-edge orientations conditional on every leaf being bad.

For every non-isomorphic tree in the requested order range, compute all
vertex-deletion splits exactly.  Among trees whose leaf splits all have modal
distance at least two, count edges with two left-bad or two right-bad
endpoints.  A same-orientation edge is a counterexample to the candidate
alternation lemma.  The output is an independently replayable JSON artifact;
absence of a witness in a finite range is computational evidence only.
"""

from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path

from scratch_eads_local_constraint_breaker_20260808 import (
    encode_split,
    signed_orientation,
)
from scratch_existential_split_mode_20260807 import tree_vertex_splits
from trees import trees_geng_raw


def encode_witness(
    order: int,
    tree_index: int,
    graph6: bytes,
    adj: list[list[int]],
    polynomial: tuple[int, ...],
    splits,
    oriented: list[tuple[str, int]],
    same_edges: list[tuple[int, int, str]],
) -> dict[str, object]:
    return {
        "n": order,
        "tree_index": tree_index,
        "graph6": graph6.decode("ascii"),
        "edges": [
            [u, v]
            for u, row in enumerate(adj)
            for v in row
            if u < v
        ],
        "degrees": [len(row) for row in adj],
        "polynomial": list(polynomial),
        "same_orientation_bad_edges": [list(edge) for edge in same_edges],
        "splits": [encode_split(split) for split in splits],
        "orientations": [
            {"vertex": v, "orientation": direction, "separation": separation}
            for v, (direction, separation) in enumerate(oriented)
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=19)
    parser.add_argument("--max-n", type=int, default=20)
    parser.add_argument("--stop-on-first", action="store_true")
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/eads_all_leaf_orientation_n20_20260808.json"),
    )
    args = parser.parse_args()
    started = time.time()
    total_trees = 0
    total_all_leaf_bad = 0
    total_ll = 0
    total_rr = 0
    witnesses: list[dict[str, object]] = []
    orders: list[dict[str, object]] = []

    for order in range(args.min_n, args.max_n + 1):
        order_started = time.time()
        order_trees = 0
        order_all_leaf_bad = 0
        order_ll = 0
        order_rr = 0
        orientation_histogram: Counter[str] = Counter()
        for tree_index, (_, adj, graph6) in enumerate(trees_geng_raw(order)):
            polynomial, splits = tree_vertex_splits(adj)
            leaves = [v for v, row in enumerate(adj) if len(row) == 1]
            if min(splits[v].distance for v in leaves) >= 2:
                order_all_leaf_bad += 1
                oriented = [signed_orientation(split) for split in splits]
                for direction, _ in oriented:
                    orientation_histogram[direction] += 1
                same_edges: list[tuple[int, int, str]] = []
                for u, row in enumerate(adj):
                    for v in row:
                        if u >= v:
                            continue
                        if (
                            oriented[u][0] == oriented[v][0]
                            and oriented[u][0] in {"left", "right"}
                            and oriented[u][1] >= 2
                            and oriented[v][1] >= 2
                        ):
                            direction = oriented[u][0]
                            same_edges.append((u, v, direction))
                            if direction == "left":
                                order_ll += 1
                            else:
                                order_rr += 1
                if same_edges and len(witnesses) < 20:
                    witnesses.append(encode_witness(
                        order, tree_index, graph6, adj, polynomial,
                        splits, oriented, same_edges,
                    ))
                if same_edges and args.stop_on_first:
                    order_trees += 1
                    break
            order_trees += 1

        total_trees += order_trees
        total_all_leaf_bad += order_all_leaf_bad
        total_ll += order_ll
        total_rr += order_rr
        row = {
            "n": order,
            "trees": order_trees,
            "all_leaf_bad_trees": order_all_leaf_bad,
            "left_left_bad_edges": order_ll,
            "right_right_bad_edges": order_rr,
            "orientation_histogram_on_all_leaf_bad_trees": dict(
                sorted(orientation_histogram.items())
            ),
            "elapsed_seconds": time.time() - order_started,
        }
        orders.append(row)
        print(json.dumps({"event": "order", **row}), flush=True)
        if args.stop_on_first and (order_ll or order_rr):
            break

    report = {
        "claim_scope": "exhaustive exact computation; finite-range evidence only",
        "statement_tested": (
            "if every leaf split is bad, no edge has two bad endpoints "
            "of the same orientation"
        ),
        "min_n": args.min_n,
        "max_n": args.max_n,
        "total_trees": total_trees,
        "all_leaf_bad_trees": total_all_leaf_bad,
        "left_left_bad_edges": total_ll,
        "right_right_bad_edges": total_rr,
        "witnesses": witnesses,
        "orders": orders,
        "elapsed_seconds": time.time() - started,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "trees": total_trees,
        "all_leaf_bad_trees": total_all_leaf_bad,
        "left_left_bad_edges": total_ll,
        "right_right_bad_edges": total_rr,
        "elapsed_seconds": report["elapsed_seconds"],
        "out": str(args.out),
    }), flush=True)
    raise SystemExit(1 if total_ll or total_rr else 0)


if __name__ == "__main__":
    main()
