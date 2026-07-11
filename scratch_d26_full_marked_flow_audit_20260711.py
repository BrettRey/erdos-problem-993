#!/usr/bin/env python3
"""Exact audit of the full nonlocal D26 fallback/marked flow.

For rank r, a source (A,C) in I_r x I_(r+2) has demand r+2.  For every
deletion D subset C of size r+1 it reaches

* a fallback target F_(A,D), of capacity one; and
* a marked target G_(D,j), of capacity N_j, for every matching block j
  empty in A.

Here N_j is the number of marked pairs (B,b) with B in I_(r+1) and
b in B intersect Q_j.  Since an independent set occupies at most one vertex
of a matching block, N_j is simply the number of rank-(r+1) sets occupying
block j.  The total target capacity is

    i_r i_(r+1) + (r+1) i_(r+1)^2,

so the global Hall cut is exactly GSB.  This script checks every cut at once
by exact integer maximum flow.  It is deliberately independent of the saved
support-weight enumerator whose old adjacency indexing was defective.
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import networkx as nx
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_flow

from scratch_d17_checkpointed_search_20260711 import append_record, utc_now
from scratch_d26_marked_pair_obstruction_certificate_20260711 import (
    alpha_tree,
    graph_data,
    independent_masks,
    prefix_last,
)
from scratch_d26_support_exchange_stress_20260711 import maximum_matching_blocks


def deletion_masks(mask: int) -> tuple[int, ...]:
    return tuple(mask ^ (1 << vertex) for vertex in range(mask.bit_length()) if mask >> vertex & 1)


def audit_case(tree: nx.Graph, rank: int) -> dict[str, object]:
    adjacency = graph_data(tree)  # Explicitly indexed by labels 0,...,n-1.
    blocks = maximum_matching_blocks(tree)
    alpha = len(blocks)
    assert alpha == alpha_tree(tree)

    block_of = [-1] * len(tree)
    for block, vertices in enumerate(blocks):
        for vertex in vertices:
            block_of[vertex] = block
    assert all(block >= 0 for block in block_of)

    level_r = independent_masks(adjacency, rank)
    level_r1 = independent_masks(adjacency, rank + 1)
    level_r2 = independent_masks(adjacency, rank + 2)
    i_r, i_r1, i_r2 = len(level_r), len(level_r1), len(level_r2)

    n_marked = [0] * alpha
    for chosen in level_r1:
        occupied = 0
        for vertex in range(len(tree)):
            if chosen >> vertex & 1:
                occupied |= 1 << block_of[vertex]
        assert occupied.bit_count() == rank + 1
        for block in range(alpha):
            if occupied >> block & 1:
                n_marked[block] += 1
    assert sum(n_marked) == (rank + 1) * i_r1

    support_r: dict[int, int] = {}
    for chosen in level_r:
        occupied = 0
        for vertex in range(len(tree)):
            if chosen >> vertex & 1:
                occupied |= 1 << block_of[vertex]
        assert occupied.bit_count() == rank
        support_r[chosen] = occupied

    demand = (rank + 2) * i_r * i_r2
    target_capacity = i_r * i_r1 + i_r1 * sum(n_marked)
    global_slack = target_capacity - demand
    assert global_slack >= 0

    # Node layout: source; all pair sources; fallback targets; marked targets;
    # sink.  Dictionaries keep the construction transparent and auditable.
    next_node = 1
    pair_node: dict[tuple[int, int], int] = {}
    for chosen_r in level_r:
        for chosen_r2 in level_r2:
            pair_node[(chosen_r, chosen_r2)] = next_node
            next_node += 1

    fallback_node: dict[tuple[int, int], int] = {}
    for chosen_r in level_r:
        for deleted in level_r1:
            fallback_node[(chosen_r, deleted)] = next_node
            next_node += 1

    marked_node: dict[tuple[int, int], int] = {}
    for deleted in level_r1:
        for block in range(alpha):
            marked_node[(deleted, block)] = next_node
            next_node += 1

    sink = next_node
    node_count = sink + 1
    infinite = demand + 1

    rows: list[int] = []
    columns: list[int] = []
    capacities: list[int] = []

    def add_edge(source: int, target: int, capacity: int) -> None:
        assert capacity >= 0
        if capacity:
            rows.append(source)
            columns.append(target)
            capacities.append(capacity)

    for (chosen_r, chosen_r2), node in pair_node.items():
        add_edge(0, node, rank + 2)
        empty_blocks = [
            block for block in range(alpha) if not (support_r[chosen_r] >> block) & 1
        ]
        for deleted in deletion_masks(chosen_r2):
            add_edge(node, fallback_node[(chosen_r, deleted)], infinite)
            for block in empty_blocks:
                add_edge(node, marked_node[(deleted, block)], infinite)

    for node in fallback_node.values():
        add_edge(node, sink, 1)
    for (deleted, block), node in marked_node.items():
        add_edge(node, sink, n_marked[block])

    network = coo_matrix(
        (capacities, (rows, columns)), shape=(node_count, node_count), dtype="int64"
    ).tocsr()
    result = maximum_flow(network, 0, sink, method="dinic")
    flow_value = int(result.flow_value)
    assert flow_value == demand, (
        nx.to_graph6_bytes(tree, header=False).decode().strip(),
        rank,
        flow_value,
        demand,
    )

    return {
        "alpha": alpha,
        "rank": rank,
        "i_r": i_r,
        "i_r1": i_r1,
        "i_r2": i_r2,
        "demand": demand,
        "target_capacity": target_capacity,
        "global_slack": global_slack,
        "flow_value": flow_value,
        "nodes": node_count,
        "arcs": len(capacities),
        "blocks": [list(block) for block in blocks],
        "n_marked": n_marked,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-order", type=int, default=2)
    parser.add_argument("--max-order", type=int, default=11)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/d26_full_marked_flow_audit_20260711.jsonl"),
    )
    args = parser.parse_args()

    cases = 0
    trees = 0
    sharp: dict[str, object] | None = None
    largest: dict[str, object] | None = None
    for n in range(args.min_order, args.max_order + 1):
        for tree_index, tree in enumerate(nx.nonisomorphic_trees(n)):
            trees += 1
            alpha = alpha_tree(tree)
            for rank in range(max(0, prefix_last(alpha) + 1)):
                row = audit_case(tree, rank)
                cases += 1
                tagged = {
                    "n": n,
                    "tree_index": tree_index,
                    "graph6": nx.to_graph6_bytes(tree, header=False).decode().strip(),
                    **row,
                }
                if sharp is None or row["global_slack"] < sharp["global_slack"]:
                    sharp = tagged
                if largest is None or row["arcs"] > largest["arcs"]:
                    largest = tagged

    assert sharp is not None and largest is not None
    record = {
        "at": utc_now(),
        "kind": "d26_full_marked_flow_audit",
        "certificate": "passed",
        "order_range": [args.min_order, args.max_order],
        "trees": trees,
        "rank_cases": cases,
        "sharp_global_cut": sharp,
        "largest_network": largest,
        "statement": "Every full fallback/marked network saturated exactly. Min-cut surplus itself is zero because the source cut has capacity equal to demand; global_slack records the nontrivial all-target cut surplus.",
    }
    append_record(args.output, record)
    print(json.dumps(record, sort_keys=True))


if __name__ == "__main__":
    main()
