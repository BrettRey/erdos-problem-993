#!/usr/bin/env python3
"""Exact EADS search by joining two adversarial tree blocks.

Single-edge evolutionary moves cannot duplicate a coordinated modal phase.
This search takes the current two-certificate tree and exhausts every pair of
attachment vertices under vertex coalescence and bridges of prescribed length.
It tests both the EADS proof route and Erdos Problem 993 with exact integers.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from scratch_eads_pareto_optimizer_20260808 import (
    adjacency,
    encode,
    make_individual,
    total_key,
)


def load_block(path: Path, leader: str) -> list[list[int]]:
    data = json.loads(path.read_text())
    record = data["leaders"][leader]
    return adjacency(int(record["n"]), record["edges"])


def bridge(
    left: list[list[int]], right: list[list[int]], u: int, v: int,
    internal_vertices: int,
) -> list[list[int]]:
    """Join two trees by a path with ``internal_vertices`` new vertices."""
    n_left = len(left)
    n_right = len(right)
    out = [[] for _ in range(n_left + n_right + internal_vertices)]
    for a, neighbors in enumerate(left):
        for b in neighbors:
            if a < b:
                out[a].append(b)
                out[b].append(a)
    for a, neighbors in enumerate(right):
        for b in neighbors:
            if a < b:
                aa, bb = n_left + a, n_left + b
                out[aa].append(bb)
                out[bb].append(aa)
    path = [u]
    path.extend(range(n_left + n_right, len(out)))
    path.append(n_left + v)
    for a, b in zip(path, path[1:]):
        out[a].append(b)
        out[b].append(a)
    return out


def coalesce(
    left: list[list[int]], right: list[list[int]], u: int, v: int,
) -> list[list[int]]:
    """Identify vertex ``u`` of the first tree with ``v`` of the second."""
    n_left = len(left)
    mapping = {v: u}
    next_vertex = n_left
    for old in range(len(right)):
        if old != v:
            mapping[old] = next_vertex
            next_vertex += 1
    out = [neighbors.copy() for neighbors in left]
    out.extend([] for _ in range(next_vertex - n_left))
    for a, neighbors in enumerate(right):
        for b in neighbors:
            if a < b:
                aa, bb = mapping[a], mapping[b]
                out[aa].append(bb)
                out[bb].append(aa)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", type=Path,
        default=Path("results/eads_pareto_optimizer_20260808.json"),
    )
    parser.add_argument("--leader", default="total_good")
    parser.add_argument("--right-input", type=Path)
    parser.add_argument("--right-leader")
    parser.add_argument(
        "--path-internals", type=int, nargs="*", default=[0, 1, 2, 3],
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/eads_block_gluing_20260808.json"),
    )
    args = parser.parse_args()
    left_block = load_block(args.input, args.leader)
    right_input = args.right_input or args.input
    right_leader = args.right_leader or args.leader
    right_block = load_block(right_input, right_leader)
    operations = ["coalesce"] + [f"bridge_{r}" for r in args.path_internals]
    best: dict[str, dict[str, object]] = {}
    evaluated = 0
    eads_counterexample = None
    main_counterexample = None

    for operation in operations:
        for u in range(len(left_block)):
            for v in range(len(right_block)):
                if operation == "coalesce":
                    adj = coalesce(left_block, right_block, u, v)
                else:
                    internal = int(operation.removeprefix("bridge_"))
                    adj = bridge(left_block, right_block, u, v, internal)
                item = make_individual(adj, f"{operation}:{u}:{v}")
                item["operation"] = operation
                item["attachment_vertices"] = [u, v]
                evaluated += 1
                incumbent = best.get(operation)
                if incumbent is None or total_key(item) < total_key(incumbent):
                    best[operation] = item
                ev = item["evaluation"]
                assert isinstance(ev, dict)
                if not ev["unimodal"]:
                    main_counterexample = item
                    break
                if ev["good_count"] == 0:
                    eads_counterexample = item
                    break
            if eads_counterexample or main_counterexample:
                break
        leader = best[operation]
        print(json.dumps({
            "event": "operation",
            "operation": operation,
            "evaluated": evaluated,
            "best_good": leader["evaluation"]["good_count"],
            "best_label": leader["label"],
        }), flush=True)
        if eads_counterexample or main_counterexample:
            break

    def encoded(item: dict[str, object] | None) -> dict[str, object] | None:
        if item is None:
            return None
        return {
            "operation": item["operation"],
            "attachment_vertices": item["attachment_vertices"],
            **encode(item),
        }

    result = {
        "claim_scope": "exact finite block-gluing computation; evidence only",
        "input": str(args.input),
        "leader": args.leader,
        "right_input": str(right_input),
        "right_leader": right_leader,
        "block_orders": [len(left_block), len(right_block)],
        "operations": operations,
        "evaluated": evaluated,
        "erdos_993_counterexample_found": main_counterexample is not None,
        "eads_counterexample_found": eads_counterexample is not None,
        "erdos_993_counterexample": encoded(main_counterexample),
        "eads_counterexample": encoded(eads_counterexample),
        "best_by_operation": {
            operation: encoded(item) for operation, item in best.items()
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "evaluated": evaluated,
        "erdos_993_counterexample_found": main_counterexample is not None,
        "eads_counterexample_found": eads_counterexample is not None,
        "out": str(args.out),
    }), flush=True)
    raise SystemExit(1 if eads_counterexample or main_counterexample else 0)


if __name__ == "__main__":
    main()
