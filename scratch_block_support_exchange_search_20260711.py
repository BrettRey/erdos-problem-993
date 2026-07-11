#!/usr/bin/env python3
"""Stress-test the maximum-matching block support exchange inequality.

For a maximum matching of a tree, its edges and unmatched vertices give
alpha(T) blocks.  Let w(R) count independent sets whose occupied-block
support is exactly R.  The candidate inequality is

    sum_{j in R, q not in R} w(R-j+q) >= |R| w(R).

The script uses exact integer arithmetic and samples random labelled trees.
"""

from __future__ import annotations

import argparse
import random
from collections import defaultdict
from math import ceil

import networkx as nx


def matching_blocks(tree: nx.Graph) -> tuple[list[tuple[int, ...]], list[int]]:
    matching = nx.max_weight_matching(tree, maxcardinality=True)
    blocks = [tuple(sorted(edge)) for edge in matching]
    covered = {v for edge in matching for v in edge}
    blocks.extend((v,) for v in tree if v not in covered)
    blocks.sort()
    block_of = [0] * len(tree)
    for index, block in enumerate(blocks):
        for vertex in block:
            block_of[vertex] = index
    return blocks, block_of


def support_weights(tree: nx.Graph, block_of: list[int]) -> dict[int, int]:
    # ``nonisomorphic_trees`` does not promise node-iteration order
    # ``0,1,...,n-1``.  This list is indexed by the vertex label below, so it
    # must be constructed in label order rather than graph iteration order.
    adjacency = [
        sum(1 << u for u in tree.neighbors(v)) for v in range(len(tree))
    ]
    weights: dict[int, int] = defaultdict(int)

    def visit(available: int, chosen_support: int) -> None:
        if not available:
            weights[chosen_support] += 1
            return
        bit = available & -available
        vertex = bit.bit_length() - 1
        visit(available ^ bit, chosen_support)
        visit(
            available & ~bit & ~adjacency[vertex],
            chosen_support | (1 << block_of[vertex]),
        )

    visit((1 << len(tree)) - 1, 0)
    return weights


def check_tree(tree: nx.Graph) -> dict[str, object] | None:
    blocks, block_of = matching_blocks(tree)
    alpha = len(blocks)
    weights = support_weights(tree, block_of)
    last = ceil((2 * alpha - 1) / 3) - 2
    all_blocks = (1 << alpha) - 1
    for support, weight in weights.items():
        rank = support.bit_count()
        if rank < 1 or rank > last:
            continue
        outside = all_blocks ^ support
        neighbor_weight = 0
        selected = support
        while selected:
            j = selected & -selected
            selected ^= j
            empty = outside
            while empty:
                q = empty & -empty
                empty ^= q
                neighbor_weight += weights.get((support ^ j) | q, 0)
        if neighbor_weight < rank * weight:
            return {
                "n": len(tree),
                "graph6": nx.to_graph6_bytes(tree, header=False).decode().strip(),
                "blocks": blocks,
                "alpha": alpha,
                "rank": rank,
                "support": support,
                "weight": weight,
                "neighbor_weight": neighbor_weight,
                "gap": neighbor_weight - rank * weight,
            }
    return None


def random_tree(n: int, rng: random.Random) -> nx.Graph:
    prufer = [rng.randrange(n) for _ in range(n - 2)]
    return nx.from_prufer_sequence(prufer)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=5000)
    parser.add_argument("--min-n", type=int, default=15)
    parser.add_argument("--max-n", type=int, default=28)
    parser.add_argument("--seed", type=int, default=993)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    sharpest: tuple[float, dict[str, object]] | None = None
    for trial in range(1, args.trials + 1):
        n = rng.randint(args.min_n, args.max_n)
        tree = random_tree(n, rng)
        failure = check_tree(tree)
        if failure is not None:
            print({"status": "failure", "trial": trial, **failure}, flush=True)
            return
        if trial % 100 == 0:
            print({"status": "progress", "trial": trial}, flush=True)
    print({"status": "passed", "trials": args.trials}, flush=True)


if __name__ == "__main__":
    main()
