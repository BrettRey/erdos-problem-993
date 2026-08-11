#!/usr/bin/env python3
"""Replace the final EADS phase anchor by every small rooted-tree gadget.

The current extremal tree has exactly two good vertices: a leaf and its
degree-two support.  Removing that pendant pair leaves a fixed 21-vertex base.
This scan attaches every rooted unlabeled tree in the requested order range at
the former parent, testing whether a nontrivial finite gadget can absorb both
remaining certificates.  All polynomial computations are exact.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from scratch_eads_pareto_optimizer_20260808 import (
    adjacency,
    encode,
    make_individual,
    scaled_total_key,
    total_key,
)
from trees import trees_geng_raw


def load_extremal(path: Path) -> tuple[list[list[int]], int, list[int]]:
    data = json.loads(path.read_text())
    record = data["leaders"]["total_good"]
    adj = adjacency(int(record["n"]), record["edges"])
    good = list(record["good_vertices"])
    leaf = next(v for v in good if len(adj[v]) == 1)
    support = adj[leaf][0]
    if support not in good or len(adj[support]) != 2:
        raise ValueError("expected a good pendant leaf/support pair")
    parent = next(v for v in adj[support] if v != leaf)
    removed = {leaf, support}
    mapping: dict[int, int] = {}
    for old in range(len(adj)):
        if old not in removed:
            mapping[old] = len(mapping)
    edges = [
        [mapping[u], mapping[v]]
        for u, neighbors in enumerate(adj)
        for v in neighbors
        if u < v and u not in removed and v not in removed
    ]
    return adjacency(len(mapping), edges), mapping[parent], [leaf, support]


def attach(
    base: list[list[int]], parent: int, gadget: list[list[int]], root: int,
) -> list[list[int]]:
    offset = len(base)
    out = [neighbors.copy() for neighbors in base]
    out.extend([] for _ in gadget)
    for u, neighbors in enumerate(gadget):
        for v in neighbors:
            if u < v:
                out[offset + u].append(offset + v)
                out[offset + v].append(offset + u)
    out[parent].append(offset + root)
    out[offset + root].append(parent)
    return out


def gadgets(n: int):
    if n == 1:
        yield "@", [[]]
    else:
        for _, adj, graph6 in trees_geng_raw(n):
            yield graph6.decode("ascii"), adj


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", type=Path,
        default=Path("results/eads_potential_optimizer_20260808.json"),
    )
    parser.add_argument("--min-gadget-n", type=int, default=1)
    parser.add_argument("--max-gadget-n", type=int, default=13)
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/eads_anchor_gadget_scan_20260808.json"),
    )
    args = parser.parse_args()
    base, parent, removed = load_extremal(args.input)
    best = None
    scaled_best = None
    eads_counterexample = None
    main_counterexample = None
    evaluated = 0
    rows = []

    for n in range(args.min_gadget_n, args.max_gadget_n + 1):
        order_evaluated = 0
        order_min_good = None
        order_hist: dict[int, int] = {}
        for graph6, gadget in gadgets(n):
            for root in range(n):
                item = make_individual(
                    attach(base, parent, gadget, root),
                    f"gadget_n{n}:{graph6}:root{root}",
                )
                item["gadget_graph6"] = graph6
                item["gadget_root"] = root
                evaluated += 1
                order_evaluated += 1
                ev = item["evaluation"]
                assert isinstance(ev, dict)
                good = int(ev["good_count"])
                order_hist[good] = order_hist.get(good, 0) + 1
                order_min_good = good if order_min_good is None else min(order_min_good, good)
                if best is None or total_key(item) < total_key(best):
                    best = item
                if scaled_best is None or scaled_total_key(item) < scaled_total_key(scaled_best):
                    scaled_best = item
                if not ev["unimodal"]:
                    main_counterexample = item
                    break
                if good == 0:
                    eads_counterexample = item
                    break
            if eads_counterexample or main_counterexample:
                break
        row = {
            "gadget_order": n,
            "evaluated": order_evaluated,
            "min_good": order_min_good,
            "good_count_histogram": {str(k): v for k, v in sorted(order_hist.items())},
        }
        rows.append(row)
        print(json.dumps({"event": "order", **row}), flush=True)
        if eads_counterexample or main_counterexample:
            break
    assert best is not None and scaled_best is not None

    def record(item: dict[str, object] | None):
        if item is None:
            return None
        return {
            "gadget_graph6": item["gadget_graph6"],
            "gadget_root": item["gadget_root"],
            **encode(item),
        }

    result = {
        "claim_scope": "exhaustive rooted-gadget computation in a fixed base; evidence only",
        "input": str(args.input),
        "base_order": len(base),
        "attachment_parent": parent,
        "removed_original_vertices": removed,
        "gadget_order_range": [args.min_gadget_n, args.max_gadget_n],
        "evaluated": evaluated,
        "orders": rows,
        "erdos_993_counterexample_found": main_counterexample is not None,
        "eads_counterexample_found": eads_counterexample is not None,
        "erdos_993_counterexample": record(main_counterexample),
        "eads_counterexample": record(eads_counterexample),
        "best": record(best),
        "scaled_best": record(scaled_best),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done", "evaluated": evaluated,
        "best_good": best["evaluation"]["good_count"],
        "eads_counterexample_found": eads_counterexample is not None,
        "erdos_993_counterexample_found": main_counterexample is not None,
        "out": str(args.out),
    }), flush=True)
    raise SystemExit(1 if eads_counterexample or main_counterexample else 0)


if __name__ == "__main__":
    main()
