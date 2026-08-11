#!/usr/bin/env python3
"""Replace an interior two-certificate anchor by all small two-terminal trees."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from scratch_eads_anchor_gadget_scan_20260808 import gadgets
from scratch_eads_decorated_spider_optimizer_20260808 import build, seed_specs
from scratch_eads_pareto_optimizer_20260808 import (
    encode,
    evaluate,
    make_individual,
    scaled_total_key,
    total_key,
)


def remove_anchor(adj: list[list[int]]) -> tuple[list[list[int]], list[int], list[int]]:
    ev = evaluate(adj)
    good = list(ev["good_vertices"])
    leaves = [v for v in good if len(adj[v]) == 1]
    if len(good) != 2 or len(leaves) != 1:
        raise ValueError((good, leaves))
    leaf = leaves[0]
    support = adj[leaf][0]
    boundaries = [v for v in adj[support] if v != leaf]
    if len(boundaries) != 2:
        raise ValueError((support, boundaries))
    removed = {leaf, support}
    mapping: dict[int, int] = {}
    for old in range(len(adj)):
        if old not in removed:
            mapping[old] = len(mapping)
    out = [[] for _ in mapping]
    for u, neighbors in enumerate(adj):
        for v in neighbors:
            if u < v and u not in removed and v not in removed:
                uu, vv = mapping[u], mapping[v]
                out[uu].append(vv)
                out[vv].append(uu)
    return out, [mapping[v] for v in boundaries], [leaf, support]


def attach(
    base: list[list[int]], boundaries: list[int], gadget: list[list[int]],
    terminals: tuple[int, int],
) -> list[list[int]]:
    offset = len(base)
    out = [neighbors.copy() for neighbors in base]
    out.extend([] for _ in gadget)
    for u, neighbors in enumerate(gadget):
        for v in neighbors:
            if u < v:
                out[offset + u].append(offset + v)
                out[offset + v].append(offset + u)
    for boundary, terminal in zip(boundaries, terminals):
        out[boundary].append(offset + terminal)
        out[offset + terminal].append(boundary)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-gadget-n", type=int, default=10)
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/eads_two_terminal_gadget_scan_20260808.json"),
    )
    args = parser.parse_args()
    source_spec = seed_specs()[-1]
    source = build(source_spec)
    base, boundaries, removed = remove_anchor(source)
    best = None
    scaled_best = None
    eads_counterexample = None
    main_counterexample = None
    evaluated = 0
    rows = []
    for n in range(1, args.max_gadget_n + 1):
        order_count = 0
        order_min = None
        hist: dict[int, int] = {}
        for graph6, gadget in gadgets(n):
            for left in range(n):
                for right in range(n):
                    item = make_individual(
                        attach(base, boundaries, gadget, (left, right)),
                        f"gadget_n{n}:{graph6}:term{left},{right}",
                    )
                    item["gadget_graph6"] = graph6
                    item["gadget_terminals"] = [left, right]
                    evaluated += 1
                    order_count += 1
                    ev = item["evaluation"]
                    assert isinstance(ev, dict)
                    good = int(ev["good_count"])
                    hist[good] = hist.get(good, 0) + 1
                    order_min = good if order_min is None else min(order_min, good)
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
            if eads_counterexample or main_counterexample:
                break
        row = {
            "gadget_order": n,
            "evaluated": order_count,
            "min_good": order_min,
            "good_count_histogram": {str(k): v for k, v in sorted(hist.items())},
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
            "gadget_terminals": item["gadget_terminals"],
            **encode(item),
        }

    result = {
        "claim_scope": "exhaustive two-terminal gadget computation in a fixed base; evidence only",
        "source_spec": source_spec,
        "source_order": len(source),
        "base_order": len(base),
        "boundaries": boundaries,
        "removed_source_vertices": removed,
        "max_gadget_order": args.max_gadget_n,
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
