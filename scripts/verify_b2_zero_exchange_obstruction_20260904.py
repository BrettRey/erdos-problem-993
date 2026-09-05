#!/usr/bin/env python3
"""Verify the first obstruction to the naive five-repair injection."""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import networkx as nx


GRAPH6 = "Hs`?GGC"
DEFAULT_OUTPUT = Path("results/b2_zero_exchange_obstruction_20260904.json")


def independent_masks(tree: nx.Graph) -> list[int]:
    adjacency = [0] * tree.number_of_nodes()
    for left, right in tree.edges():
        adjacency[left] |= 1 << right
        adjacency[right] |= 1 << left

    masks: list[int] = []
    for mask in range(1 << tree.number_of_nodes()):
        remaining = mask
        while remaining:
            bit = remaining & -remaining
            vertex = bit.bit_length() - 1
            remaining -= bit
            if adjacency[vertex] & remaining:
                break
        else:
            masks.append(mask)
    return masks


def run_verification() -> dict[str, object]:
    tree = nx.from_graph6_bytes(GRAPH6.encode())
    independent = independent_masks(tree)
    alpha = max(mask.bit_count() for mask in independent)
    maximum = [mask for mask in independent if mask.bit_count() == alpha]

    profiles: dict[int, tuple[int, int]] = {}
    blocked_by_defect: dict[int, list[int]] = {}
    for defect in range(5):
        size = alpha - defect
        layer = [mask for mask in independent if mask.bit_count() == size]
        extendable = [
            mask
            for mask in layer
            if any(mask & ~maximum_mask == 0 for maximum_mask in maximum)
        ]
        blocked = [mask for mask in layer if mask not in set(extendable)]
        profiles[defect] = (len(extendable), len(blocked))
        blocked_by_defect[defect] = blocked

    target_blocked = blocked_by_defect[4]
    repair_neighborhoods: list[set[int]] = []
    for blocked in target_blocked:
        repairs: set[int] = set()
        for maximum_mask in maximum:
            outside = blocked & ~maximum_mask
            available = maximum_mask & ~blocked
            common = blocked & maximum_mask
            replacement_count = outside.bit_count()
            vertices = [
                vertex
                for vertex in range(tree.number_of_nodes())
                if available & (1 << vertex)
            ]
            for replacements in combinations(vertices, replacement_count):
                repaired = common
                for vertex in replacements:
                    repaired |= 1 << vertex
                repairs.add(repaired)
        repair_neighborhoods.append(repairs)

    source = ("source",)
    sink = ("sink",)
    flow_graph = nx.DiGraph()
    for index, repairs in enumerate(repair_neighborhoods):
        blocked_node = ("blocked", index)
        flow_graph.add_edge(source, blocked_node, capacity=5)
        for repaired in repairs:
            extendable_node = ("extendable", repaired)
            flow_graph.add_edge(blocked_node, extendable_node, capacity=1)
            flow_graph.add_edge(extendable_node, sink, capacity=1)
    flow_value = nx.maximum_flow_value(flow_graph, source, sink)

    extendable = [profiles[defect][0] for defect in range(5)]
    blocked = [profiles[defect][1] for defect in range(5)]
    result = {
        "certificate_date": "2026-09-04",
        "kind": "naive_b2_zero_exchange_injection_obstruction",
        "claim_status": "exact counterexample to this injection, not to 5*b_4<=e_4",
        "graph6": GRAPH6,
        "n": tree.number_of_nodes(),
        "alpha": alpha,
        "edges": [list(edge) for edge in sorted(tree.edges())],
        "e_0_to_4": extendable,
        "b_0_to_4": blocked,
        "repair_rule": (
            "For each blocked defect-4 set S and every maximum set M, replace "
            "S\\M by an equally large subset of M\\S, retaining S intersection M."
        ),
        "repair_neighborhood_sizes": [
            len(repairs) for repairs in repair_neighborhoods
        ],
        "required_distinct_repair_slots": 5 * len(target_blocked),
        "maximum_distinct_repair_slots": flow_value,
        "global_five_bound_lhs": 5 * blocked[4],
        "global_five_bound_rhs": extendable[4],
    }

    expected = {
        "n": 9,
        "alpha": 7,
        "e_0_to_4": [1, 7, 21, 35, 35],
        "b_0_to_4": [0, 0, 0, 2, 6],
        "repair_neighborhood_sizes": [5, 5, 5, 5, 5, 5],
        "required_distinct_repair_slots": 30,
        "maximum_distinct_repair_slots": 26,
        "global_five_bound_lhs": 30,
        "global_five_bound_rhs": 35,
    }
    for key, value in expected.items():
        if result[key] != value:
            raise AssertionError(f"{key}: {result[key]!r} != {value!r}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run_verification()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "PASS: naive five-repair injection max flow "
        f"{result['maximum_distinct_repair_slots']} < "
        f"{result['required_distinct_repair_slots']}; global bound survives"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
