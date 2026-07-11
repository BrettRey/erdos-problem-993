#!/usr/bin/env python3
"""Search one-edge augmentations of trees for a nonunimodal base graph.

If ``G=T+uv`` with ``uv`` a nonedge of the tree, then

    I(G;x) = I(T;x) - x^2 I(T-N_T[u]-N_T[v];x).

This makes it cheap to scan unicyclic graphs exactly.  A hit would be a
candidate input for the Bencs stable-path-tree transformation; it is not by
itself a counterexample to the tree conjecture.
"""

from __future__ import annotations

import argparse
import random

import networkx as nx

from indpoly import independence_poly, is_unimodal


def adjacency_list(graph: nx.Graph) -> list[list[int]]:
    return [sorted(graph.neighbors(v)) for v in range(len(graph))]


def forest_after_closed_neighborhoods(
    tree: nx.Graph, u: int, v: int
) -> list[list[int]]:
    removed = {u, v, *tree.neighbors(u), *tree.neighbors(v)}
    kept = [x for x in tree if x not in removed]
    relabel = {x: i for i, x in enumerate(kept)}
    return [
        sorted(relabel[y] for y in tree.neighbors(x) if y in relabel)
        for x in kept
    ]


def add_edge_polynomial(
    tree: nx.Graph, tree_poly: list[int], u: int, v: int
) -> list[int]:
    remainder_adj = forest_after_closed_neighborhoods(tree, u, v)
    remainder = independence_poly(len(remainder_adj), remainder_adj)
    out = list(tree_poly)
    for rank, value in enumerate(remainder):
        out[rank + 2] -= value
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    assert all(value > 0 for value in out)
    return out


def rebound_score(poly: list[int]) -> tuple[int, int, int, int]:
    first = next(
        (k for k in range(len(poly) - 1) if poly[k + 1] < poly[k]), -1
    )
    if first < 0:
        return (0, 1, -1, -1)
    best_num, best_den, best_k = 0, 1, -1
    for k in range(first + 1, len(poly) - 1):
        if poly[k + 1] * best_den > best_num * poly[k]:
            best_num, best_den, best_k = poly[k + 1], poly[k], k
    return best_num, best_den, best_k, first


def random_tree(n: int, rng: random.Random, hub_bias: float) -> nx.Graph:
    if rng.random() >= hub_bias:
        return nx.from_prufer_sequence([rng.randrange(n) for _ in range(n - 2)])
    hubs = max(2, int(n**0.5))
    return nx.from_prufer_sequence(
        [rng.randrange(hubs) if rng.random() < 0.85 else rng.randrange(n)
         for _ in range(n - 2)]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=2000)
    parser.add_argument("--min-n", type=int, default=15)
    parser.add_argument("--max-n", type=int, default=60)
    parser.add_argument("--edges-per-tree", type=int, default=200)
    parser.add_argument("--hub-bias", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=993)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    best: tuple[float, dict[str, object]] | None = None
    for trial in range(1, args.trials + 1):
        n = rng.randint(args.min_n, args.max_n)
        tree = random_tree(n, rng, args.hub_bias)
        tree_adj = adjacency_list(tree)
        tree_poly = independence_poly(n, tree_adj)
        nonedges = list(nx.non_edges(tree))
        rng.shuffle(nonedges)
        for u, v in nonedges[: args.edges_per_tree]:
            poly = add_edge_polynomial(tree, tree_poly, u, v)
            score = rebound_score(poly)
            ratio = score[0] / score[1]
            row = {
                "trial": trial,
                "tree_graph6": nx.to_graph6_bytes(tree, header=False).decode().strip(),
                "added_edge": [u, v],
                "order": n,
                "score": score,
                "ratio": ratio,
            }
            if best is None or ratio > best[0]:
                best = (ratio, row)
            if not is_unimodal(poly):
                print({"status": "nonunimodal_unicyclic", **row,
                       "coefficients": poly}, flush=True)
                return
        if trial % 100 == 0:
            print({"status": "progress", "trial": trial, "best": best}, flush=True)
    print({"status": "passed", "trials": args.trials, "best": best}, flush=True)


if __name__ == "__main__":
    main()
