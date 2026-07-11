#!/usr/bin/env python3
"""Exact prefix obstruction to the unmarked D26 support exchange.

Take seven matched double blocks D_i=(2i,2i+1), joined in a path through six
unmatched singleton blocks q_i=14+i.  Both edges at an internal double block
hit the same endpoint, so the displayed matching has no augmenting path and
is maximum.  At alpha=13, rank r=7 is the last required prefix rank.

For the support R consisting of all seven double blocks,

    w(R)=128,
    sum_{j in R,q not in R} w(R-j+q)=864 < 7*128=896.

Thus both the bare support-exchange statement and its deletion-cavity
sufficient condition are false in the exact prefix.  Adding one fallback
copy of w(R) repairs this witness:

    864 + 128 >= 896.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import networkx as nx

from scratch_d17_checkpointed_search_20260711 import append_record, utc_now
from scratch_d26_support_exchange_stress_20260711 import (
    exchange_row,
    make_csp,
    prefix_last,
    support_weight,
)


RANK = 7
SINGLETONS = 6


def witness() -> tuple[nx.Graph, tuple[tuple[int, ...], ...], int]:
    tree = nx.Graph()
    blocks: list[tuple[int, ...]] = [
        (2 * block, 2 * block + 1) for block in range(RANK)
    ]
    tree.add_nodes_from(range(2 * RANK + SINGLETONS))
    tree.add_edges_from(blocks)
    for index in range(SINGLETONS):
        singleton = 2 * RANK + index
        tree.add_edge(singleton, 2 * index)
        tree.add_edge(singleton, 2 * (index + 1))
        blocks.append((singleton,))
    support = (1 << RANK) - 1
    return tree, tuple(blocks), support


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/d26_support_exchange_obstruction_20260711.jsonl"),
    )
    args = parser.parse_args()

    tree, blocks, support = witness()
    matching = {(2 * block, 2 * block + 1) for block in range(RANK)}
    maximum = nx.max_weight_matching(tree, maxcardinality=True)
    alpha = len(tree) - len(maximum)

    assert len(tree) == 20 and tree.number_of_edges() == 19 and nx.is_tree(tree)
    assert all(tree.has_edge(*edge) for edge in matching)
    assert len(maximum) == len(matching) == RANK
    assert alpha == len(blocks) == 13
    assert prefix_last(alpha) == RANK

    csp = make_csp(tree, blocks)
    row = exchange_row(csp, support)
    assert row["rank"] == RANK
    assert row["weight"] == 2**RANK == 128

    # Each q_i has two incident deleted blocks, giving weight 2^(r-2), and
    # five nonincident deleted blocks, giving weight 2^(r-3).
    swap_weights: dict[str, int] = {}
    for deleted in range(RANK):
        for singleton in range(RANK, RANK + SINGLETONS):
            swapped = (support ^ (1 << deleted)) | (1 << singleton)
            value = support_weight(csp, swapped)
            path_index = singleton - RANK
            expected = 2 ** (RANK - 2) if deleted in (path_index, path_index + 1) else 2 ** (RANK - 3)
            assert value == expected
            swap_weights[f"{deleted},{singleton}"] = value

    expected_neighbour = SINGLETONS * (
        2 * 2 ** (RANK - 2) + (RANK - 2) * 2 ** (RANK - 3)
    )
    assert expected_neighbour == 864
    assert row["neighbour_weight"] == expected_neighbour
    assert row["exchange_gap"] == -32
    assert Fraction(row["neighbour_weight"], row["weight"]) == Fraction(27, 4)

    # Every deletion cavity has weight 2^(r-1); hence CE is exactly twice
    # the normalized neighbour weight on this witness.
    assert all(denominator == 64 for _, _, denominator in row["cavity_rows"])
    assert row["cavity_sum"] == Fraction(27, 2)
    assert row["cavity_gap"] == Fraction(-1, 2)

    # The orientation blocker sufficient inequality fails by 64: each q_i
    # sees two independent fair endpoint choices, so its expected charge is
    # r/4 + 1/2 = 9/4 per orientation.
    orientation_blocker_charge = SINGLETONS * row["weight"] * Fraction(RANK + 2, 4)
    assert orientation_blocker_charge == 1728
    assert orientation_blocker_charge - 2 * RANK * row["weight"] == -64

    fallback_gap = row["neighbour_weight"] + row["weight"] - RANK * row["weight"]
    assert fallback_gap == 96 > 0

    result = {
        "at": utc_now(),
        "kind": "d26_support_exchange_obstruction",
        "certificate": "passed",
        "tree_order": len(tree),
        "tree_edges": sorted([sorted(edge) for edge in tree.edges()]),
        "graph6": nx.to_graph6_bytes(tree, header=False).decode().strip(),
        "matching": sorted([list(edge) for edge in matching]),
        "matching_size": len(matching),
        "alpha": alpha,
        "rank": RANK,
        "prefix_last": prefix_last(alpha),
        "blocks": [list(block) for block in blocks],
        "support": support,
        "support_blocks": list(range(RANK)),
        "weight": row["weight"],
        "swap_weights": swap_weights,
        "neighbour_weight": row["neighbour_weight"],
        "bare_target": RANK * row["weight"],
        "bare_exchange_gap": row["exchange_gap"],
        "normalized_neighbour_weight": str(Fraction(27, 4)),
        "cavity_sum": str(row["cavity_sum"]),
        "cavity_target": 2 * RANK,
        "cavity_gap": str(row["cavity_gap"]),
        "orientation_blocker_charge": int(orientation_blocker_charge),
        "orientation_blocker_gap": -64,
        "fallback_corrected_gap": fallback_gap,
        "route_conclusion": "Bare support exchange SE, deletion-cavity CE, and the orientation blocker sufficient inequality are false in the required prefix. The minimally marked continuation is neighbour_weight + w(R) >= r w(R), using one fallback copy.",
    }
    append_record(args.output, result)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
