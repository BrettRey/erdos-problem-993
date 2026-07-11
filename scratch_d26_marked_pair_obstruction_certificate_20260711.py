#!/usr/bin/env python3
"""Exact D26 obstruction to the one-spare-block marked-pair injection.

The proposed marked source has size

    (r+2) i_r i_(r+2),

and a no-subset-sum-one source must be exported from its symmetric-difference
fiber.  The first natural marked repair gives, for fixed (A,D=C-c), one
canonical alternating-closure target for each matching block empty in A,
with that closure label used as the codomain mark, plus the single unmarked
fallback (A,D).  There are alpha-r+1 such target labels.

This script certifies an exact prefix witness with six hard marked extensions
but only alpha-r+1=5 labels.  Thus the mark does not repair the old Hall
obstruction if it merely records one canonical spare-block exchange.  Any
further marked rule must create multiple c-dependent, invertible target
states inside a single spare block; block labels plus fallback are impossible.
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from math import ceil
from pathlib import Path

import networkx as nx

from scratch_d17_checkpointed_search_20260711 import append_record, utc_now


WITNESS_G6 = "JpD?GCA?_C?"
WITNESS_A = frozenset({0, 3, 6})
WITNESS_D = frozenset({7, 8, 9, 10})


def independent(adjacency: list[int], chosen: int) -> bool:
    return all(not (adjacency[v] & chosen) for v in range(len(adjacency)) if chosen >> v & 1)


def independent_masks(adjacency: list[int], rank: int) -> list[int]:
    return [
        sum(1 << v for v in chosen)
        for chosen in combinations(range(len(adjacency)), rank)
        if independent(adjacency, sum(1 << v for v in chosen))
    ]


def addable_mask(adjacency: list[int], chosen: int) -> int:
    blocked = chosen
    for v in range(len(adjacency)):
        if chosen >> v & 1:
            blocked |= adjacency[v]
    return ((1 << len(adjacency)) - 1) ^ blocked


def component_imbalances(adjacency: list[int], a_mask: int, c_mask: int) -> tuple[int, ...]:
    difference = a_mask ^ c_mask
    seen = 0
    imbalances = []
    while difference & ~seen:
        seed = (difference & ~seen) & -(difference & ~seen)
        stack = seed
        seen |= seed
        a_count = 0
        c_count = 0
        while stack:
            bit = stack & -stack
            stack -= bit
            vertex = bit.bit_length() - 1
            a_count += bool(a_mask & bit)
            c_count += bool(c_mask & bit)
            new = adjacency[vertex] & difference & ~seen
            seen |= new
            stack |= new
        imbalances.append(c_count - a_count)
    return tuple(sorted(imbalances))


def subset_sums(values: tuple[int, ...]) -> set[int]:
    sums = {0}
    for value in values:
        sums |= {old + value for old in tuple(sums)}
    return sums


def is_hard(adjacency: list[int], a_mask: int, c_mask: int) -> bool:
    return 1 not in subset_sums(component_imbalances(adjacency, a_mask, c_mask))


def alpha_tree(graph: nx.Graph) -> int:
    return len(graph) - len(nx.max_weight_matching(graph, maxcardinality=True))


def prefix_last(alpha: int) -> int:
    return ceil((2 * alpha - 1) / 3) - 2


def graph_data(graph: nx.Graph) -> list[int]:
    return [sum(1 << u for u in graph.neighbors(v)) for v in range(len(graph))]


def hard_extensions(adjacency: list[int], a_mask: int, d_mask: int) -> list[int]:
    result = []
    extensions = addable_mask(adjacency, d_mask)
    while extensions:
        bit = extensions & -extensions
        extensions -= bit
        if is_hard(adjacency, a_mask, d_mask | bit):
            result.append(bit.bit_length() - 1)
    return result


def alternating_closure(
    graph: nx.Graph,
    chosen: frozenset[int],
    maximum_set: frozenset[int],
    matching_map: dict[int, int],
    label: int,
) -> frozenset[int]:
    outside = chosen - maximum_set
    w = {label}
    while True:
        p = {x for x in outside if any(graph.has_edge(x, y) for y in w)}
        new_w = {label} | {matching_map[x] for x in p}
        if new_w == w:
            break
        w = new_w
    assert len(w) == len(p) + 1
    target = (chosen - p) | w
    assert len(target) == len(chosen) + 1
    assert all(not graph.has_edge(u, v) for u, v in combinations(target, 2))
    return frozenset(target)


def star_audit(maximum_leaves: int = 32) -> list[dict[str, int]]:
    rows = []
    for leaves in range(4, maximum_leaves + 1):
        graph = nx.star_graph(leaves)
        adjacency = graph_data(graph)
        a_mask = 1  # centre
        d_mask = (1 << 1) | (1 << 2)  # two leaves
        hard = hard_extensions(adjacency, a_mask, d_mask)
        # The only hard extensions are the leaves outside D.
        assert hard == list(range(3, leaves + 1))
        alpha = leaves
        rank = 1
        assert rank <= prefix_last(alpha)
        capacity = alpha - rank + 1
        assert len(hard) == leaves - 2 <= capacity
        rows.append(
            {
                "leaves": leaves,
                "hard_demand": len(hard),
                "spare_block_plus_fallback_capacity": capacity,
            }
        )
    return rows


def exhaustive_first_failure(max_order: int) -> dict[str, object] | None:
    tested_pairs = 0
    for n in range(2, max_order + 1):
        for tree_index, graph in enumerate(nx.generators.nonisomorphic_trees(n)):
            adjacency = graph_data(graph)
            alpha = alpha_tree(graph)
            for rank in range(max(0, prefix_last(alpha) + 1)):
                for a_mask in independent_masks(adjacency, rank):
                    for d_mask in independent_masks(adjacency, rank + 1):
                        hard = hard_extensions(adjacency, a_mask, d_mask)
                        tested_pairs += 1
                        capacity = alpha - rank + 1
                        if len(hard) > capacity:
                            return {
                                "order": n,
                                "tree_index": tree_index,
                                "graph6": nx.to_graph6_bytes(graph, header=False).decode().strip(),
                                "alpha": alpha,
                                "rank": rank,
                                "a_mask": a_mask,
                                "d_mask": d_mask,
                                "hard_extensions": hard,
                                "demand": len(hard),
                                "capacity": capacity,
                                "tested_pairs_through_first_failure": tested_pairs,
                            }
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exhaustive-max-order", type=int, default=11)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/d26_marked_pair_obstruction_20260711.jsonl"),
    )
    args = parser.parse_args()

    graph = nx.from_graph6_bytes(WITNESS_G6.encode())
    adjacency = graph_data(graph)
    alpha = alpha_tree(graph)
    rank = len(WITNESS_A)
    a_mask = sum(1 << v for v in WITNESS_A)
    d_mask = sum(1 << v for v in WITNESS_D)
    assert alpha == 7 and rank == 3 and rank == prefix_last(alpha)
    assert independent(adjacency, a_mask) and independent(adjacency, d_mask)

    hard = hard_extensions(adjacency, a_mask, d_mask)
    assert hard == [0, 1, 2, 3, 4, 5]
    profiles = {}
    for marked in hard:
        c_mask = d_mask | (1 << marked)
        imbalances = component_imbalances(adjacency, a_mask, c_mask)
        assert sum(imbalances) == 2
        assert 1 not in subset_sums(imbalances)
        profiles[str(marked)] = list(imbalances)

    maximum_set = frozenset({0, 3, 5, 7, 8, 9, 10})
    complement = set(graph) - maximum_set
    matching_map = {1: 5, 2: 0, 4: 3, 6: 10}
    assert set(matching_map) == complement
    assert len(set(matching_map.values())) == len(matching_map)
    assert all(graph.has_edge(x, y) for x, y in matching_map.items())
    blocks = tuple(sorted((x, y)) for x, y in matching_map.items()) + tuple(
        (vertex,) for vertex in sorted(maximum_set - set(matching_map.values()))
    )
    assert len(blocks) == alpha
    assert sum(bool(WITNESS_A.intersection(block)) for block in blocks) == rank

    outside = WITNESS_A - maximum_set
    labels = sorted(maximum_set - (WITNESS_A & maximum_set) - {matching_map[x] for x in outside})
    assert labels == [5, 7, 8, 9]
    marked_targets = {
        label: sorted(alternating_closure(graph, WITNESS_A, maximum_set, matching_map, label))
        for label in labels
    }
    assert len({tuple(value) for value in marked_targets.values()}) == len(labels)
    demand = len(hard)
    capacity = len(labels) + 1  # four marked closures plus (A,D)
    assert demand == 6 > 5 == capacity

    stars = star_audit()
    first_failure = exhaustive_first_failure(args.exhaustive_max_order)
    assert first_failure is not None
    assert first_failure["graph6"] == WITNESS_G6
    assert first_failure["alpha"] == alpha
    assert first_failure["rank"] == rank
    assert first_failure["a_mask"] == a_mask
    assert first_failure["d_mask"] == d_mask
    assert first_failure["hard_extensions"] == hard

    result = {
        "at": utc_now(),
        "kind": "d26_marked_pair_obstruction",
        "certificate": "passed",
        "network_model": "component switches for subset-sum-one fibers; for each hard deletion base (A,D), one canonical marked alternating-closure target per A-empty matching block plus one unmarked fallback",
        "tree": WITNESS_G6,
        "edges": sorted([list(edge) for edge in graph.edges()]),
        "alpha": alpha,
        "rank": rank,
        "prefix_last": prefix_last(alpha),
        "A": sorted(WITNESS_A),
        "D": sorted(WITNESS_D),
        "hard_marked_extensions": hard,
        "component_imbalance_profiles": profiles,
        "matching_blocks": [list(block) for block in blocks],
        "spare_labels": labels,
        "marked_closure_targets": marked_targets,
        "demand": demand,
        "marked_plus_fallback_capacity": capacity,
        "hall_deficit": demand - capacity,
        "star_checks": len(stars),
        "star_minimum_slack": min(
            row["spare_block_plus_fallback_capacity"] - row["hard_demand"]
            for row in stars
        ),
        "exhaustive_first_failure": first_failure,
        "route_conclusion": "A globally invertible marked rule cannot use only the spare-block/closure label and the one unmarked fallback. It must encode c through multiple distinct targets within at least one spare block, or use a nonlocal fractional flow.",
    }
    append_record(args.output, result)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
