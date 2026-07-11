#!/usr/bin/env python3
"""Conditioned-DP stress test for the D26 block-support exchange.

Fix a maximum matching of a tree.  Its matched edges and unmatched vertices
are blocks Q_0,...,Q_(alpha-1).  For R a set of blocks, w(R) is the number of
independent transversals which occupy exactly the blocks in R.  The live D26
candidate is

    sum_{j in R, q not in R} w(R-j+q) >= |R| w(R).            (SE)

The exhaustive helper enumerates every independent set and is exponential in
the order.  This script instead evaluates a requested w(R) by a linear tree
DP on the contracted matching blocks.  It also tests the stronger sufficient
cavity statement

    sum_{j in R} sum_{q not in R} w(R-j+q)/w(R-j) >= 2 |R|,  (CE)

using exact rational arithmetic.  The deterministic corpus includes the P36
safe-pair obstruction, random labelled trees, and adversarial all-double
block trees (paths, combs, stars, and random contracted trees) with randomized
endpoint labels.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import dataclass
from fractions import Fraction
from math import ceil
from pathlib import Path

import networkx as nx

from scratch_d17_checkpointed_search_20260711 import append_record, utc_now


@dataclass(frozen=True)
class BlockCSP:
    blocks: tuple[tuple[int, ...], ...]
    # Each item is (neighbour block, endpoint here, endpoint there).
    adjacency: tuple[tuple[tuple[int, int, int], ...], ...]


def prefix_last(alpha: int) -> int:
    return ceil((2 * alpha - 1) / 3) - 2


def maximum_matching_blocks(tree: nx.Graph) -> tuple[tuple[int, ...], ...]:
    matching = nx.max_weight_matching(tree, maxcardinality=True)
    blocks = [tuple(sorted(edge)) for edge in matching]
    covered = {vertex for edge in matching for vertex in edge}
    blocks.extend((vertex,) for vertex in range(len(tree)) if vertex not in covered)
    blocks.sort()
    return tuple(blocks)


def make_csp(tree: nx.Graph, blocks: tuple[tuple[int, ...], ...]) -> BlockCSP:
    block_of = [-1] * len(tree)
    for index, block in enumerate(blocks):
        for vertex in block:
            assert block_of[vertex] == -1
            block_of[vertex] = index
    assert all(index >= 0 for index in block_of)

    adjacency: list[list[tuple[int, int, int]]] = [[] for _ in blocks]
    internal_edges = 0
    for u, v in tree.edges():
        bu, bv = block_of[u], block_of[v]
        if bu == bv:
            internal_edges += 1
            assert len(blocks[bu]) == 2 and {u, v} == set(blocks[bu])
        else:
            adjacency[bu].append((bv, u, v))
            adjacency[bv].append((bu, v, u))
    assert internal_edges == sum(len(block) == 2 for block in blocks)

    contracted = nx.Graph()
    contracted.add_nodes_from(range(len(blocks)))
    for here, entries in enumerate(adjacency):
        for there, _, _ in entries:
            if here < there:
                contracted.add_edge(here, there)
    assert nx.is_tree(contracted)
    return BlockCSP(blocks, tuple(tuple(entries) for entries in adjacency))


def support_weight(csp: BlockCSP, support: int) -> int:
    """Count independent transversals with exactly ``support`` occupied."""

    def states(block: int) -> tuple[int, ...]:
        if not (support >> block) & 1:
            return (-1,)  # Empty block.
        return csp.blocks[block]

    def visit(block: int, parent: int) -> dict[int, int]:
        values = {state: 1 for state in states(block)}
        for child, endpoint_here, endpoint_child in csp.adjacency[block]:
            if child == parent:
                continue
            child_values = visit(child, block)
            child_total = sum(child_values.values())
            for state in tuple(values):
                allowed = child_total
                if state == endpoint_here:
                    allowed -= child_values.get(endpoint_child, 0)
                values[state] *= allowed
        return values

    return sum(visit(0, -1).values())


def exchange_row(csp: BlockCSP, support: int) -> dict[str, object]:
    alpha = len(csp.blocks)
    rank = support.bit_count()
    weight = support_weight(csp, support)
    assert weight > 0
    selected = [j for j in range(alpha) if support >> j & 1]
    empty = [q for q in range(alpha) if not (support >> q & 1)]

    neighbour_weight = 0
    cavity = Fraction(0)
    cavity_rows = []
    for j in selected:
        deleted = support ^ (1 << j)
        deleted_weight = support_weight(csp, deleted)
        extension_weight = 0
        for q in empty:
            swapped_weight = support_weight(csp, deleted | (1 << q))
            neighbour_weight += swapped_weight
            extension_weight += swapped_weight
        cavity += Fraction(extension_weight, deleted_weight)
        cavity_rows.append((j, extension_weight, deleted_weight))

    return {
        "alpha": alpha,
        "rank": rank,
        "support": support,
        "weight": weight,
        "neighbour_weight": neighbour_weight,
        "exchange_gap": neighbour_weight - rank * weight,
        "cavity_sum": cavity,
        "cavity_gap": cavity - 2 * rank,
        "cavity_rows": cavity_rows,
    }


def random_tree(n: int, rng: random.Random) -> nx.Graph:
    return nx.from_prufer_sequence([rng.randrange(n) for _ in range(n - 2)])


def doubled_block_tree(
    contracted: nx.Graph, rng: random.Random
) -> tuple[nx.Graph, tuple[tuple[int, int], ...]]:
    """Expand every contracted vertex to a matched edge with random ports."""
    alpha = len(contracted)
    tree = nx.Graph()
    blocks = tuple((2 * block, 2 * block + 1) for block in range(alpha))
    tree.add_nodes_from(range(2 * alpha))
    tree.add_edges_from(blocks)
    for u, v in contracted.edges():
        tree.add_edge(2 * u + rng.randrange(2), 2 * v + rng.randrange(2))
    assert nx.is_tree(tree)
    assert len(nx.max_weight_matching(tree, maxcardinality=True)) == alpha
    return tree, blocks


def comb_tree(alpha: int) -> nx.Graph:
    """A contracted comb, truncated to exactly ``alpha`` vertices."""
    spine = (alpha + 1) // 2
    graph = nx.path_graph(spine)
    next_vertex = spine
    for vertex in range(spine):
        if next_vertex == alpha:
            break
        graph.add_edge(vertex, next_vertex)
        next_vertex += 1
    return graph


def sampled_support(alpha: int, rng: random.Random, last_only: bool) -> int:
    last = prefix_last(alpha)
    assert last >= 1
    rank = last if last_only else rng.randint(max(1, last - 3), last)
    return sum(1 << block for block in rng.sample(range(alpha), rank))


def p36_row() -> dict[str, object]:
    tree = nx.path_graph(36)
    blocks = tuple((2 * block, 2 * block + 1) for block in range(18))
    support_blocks = list(range(0, 17, 2)) + [17]
    support = sum(1 << block for block in support_blocks)
    row = exchange_row(make_csp(tree, blocks), support)
    assert row["rank"] == 10 and row["weight"] == 768
    assert row["neighbour_weight"] == 35392
    assert row["exchange_gap"] == 27712
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--random-trials", type=int, default=1200)
    parser.add_argument("--block-trials", type=int, default=1200)
    parser.add_argument("--min-n", type=int, default=15)
    parser.add_argument("--max-n", type=int, default=30)
    parser.add_argument("--seed", type=int, default=260993)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/d26_support_exchange_stress_20260711.jsonl"),
    )
    args = parser.parse_args()
    rng = random.Random(args.seed)

    sharp_exchange: tuple[Fraction, dict[str, object]] | None = None
    sharp_cavity: tuple[Fraction, dict[str, object]] | None = None
    checks = 0

    def audit(tree: nx.Graph, blocks: tuple[tuple[int, ...], ...], support: int, label: str) -> None:
        nonlocal sharp_exchange, sharp_cavity, checks
        csp = make_csp(tree, blocks)
        row = exchange_row(csp, support)
        checks += 1
        assert row["exchange_gap"] >= 0, (label, row, blocks, sorted(tree.edges()))
        assert row["cavity_gap"] >= 0, (label, row, blocks, sorted(tree.edges()))
        exchange_ratio = Fraction(row["exchange_gap"], row["weight"])
        if sharp_exchange is None or exchange_ratio < sharp_exchange[0]:
            sharp_exchange = (exchange_ratio, {"label": label, **row})
        cavity_gap = row["cavity_gap"]
        assert isinstance(cavity_gap, Fraction)
        if sharp_cavity is None or cavity_gap < sharp_cavity[0]:
            sharp_cavity = (cavity_gap, {"label": label, **row})

    p36 = p36_row()

    for trial in range(args.random_trials):
        n = rng.randint(args.min_n, args.max_n)
        tree = random_tree(n, rng)
        blocks = maximum_matching_blocks(tree)
        if prefix_last(len(blocks)) < 1:
            continue
        support = sampled_support(len(blocks), rng, trial % 2 == 0)
        audit(tree, blocks, support, f"random-tree-{trial}")

    families = ("path", "comb", "star", "random")
    for trial in range(args.block_trials):
        alpha = rng.randint(7, 24)
        family = families[trial % len(families)]
        if family == "path":
            contracted = nx.path_graph(alpha)
        elif family == "comb":
            contracted = comb_tree(alpha)
        elif family == "star":
            contracted = nx.star_graph(alpha - 1)
        else:
            contracted = random_tree(alpha, rng)
        tree, blocks = doubled_block_tree(contracted, rng)
        support = sampled_support(alpha, rng, trial % 3 != 0)
        audit(tree, blocks, support, f"block-{family}-{trial}")

    assert sharp_exchange is not None and sharp_cavity is not None
    result = {
        "at": utc_now(),
        "kind": "d26_support_exchange_stress",
        "certificate": "passed",
        "seed": args.seed,
        "checks": checks,
        "random_tree_trials": args.random_trials,
        "all_double_block_trials": args.block_trials,
        "order_range": [args.min_n, args.max_n],
        "p36": {
            key: str(value) if isinstance(value, Fraction) else value
            for key, value in p36.items()
        },
        "sharp_exchange": {
            key: str(value) if isinstance(value, Fraction) else value
            for key, value in sharp_exchange[1].items()
            if key != "cavity_rows"
        },
        "sharp_cavity": {
            key: str(value) if isinstance(value, Fraction) else value
            for key, value in sharp_cavity[1].items()
            if key != "cavity_rows"
        },
        "statement": "SE and the sufficient deletion-cavity inequality CE passed this deterministic stress corpus; this is falsification evidence, not a proof.",
    }
    append_record(args.output, result)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
