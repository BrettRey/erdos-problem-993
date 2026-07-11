#!/usr/bin/env python3
"""Scan one-edge closures of the exact order-139 connector near miss."""

from __future__ import annotations

import argparse

import networkx as nx

from indpoly import independence_poly, is_unimodal
from scratch_connector_path_search_20260711 import build_tree, state_catalog
from scratch_unicyclic_stable_path_probe_20260711 import (
    add_edge_polynomial,
    rebound_score,
)


LEFT = [
    ("J??????oD~?", 10, 11),
    ("J????A?oB}?", 2, 11),
    ("K??????_F?B|", 11, 12),
    ("K???????C?^}", 11, 12),
    ("K???????FoO}", 0, 12),
    ("J??????wC~?", 1, 11),
    ("I????A?~o", 9, 10),
    ("K???????C?^}", 1, 12),
    ("J??????{?^_", 10, 11),
]
RIGHT = [
    ("J????B_FCe?", 6, 11),
    ("K??????_E?N{", 9, 12),
    ("J???????F~_", 0, 11),
]


def build() -> nx.Graph:
    states = state_catalog()
    lookup = {(state.graph6, state.root, state.order): i for i, state in enumerate(states)}
    left = [lookup[row] for row in LEFT]
    right = [lookup[row] for row in RIGHT]
    adjacency = build_tree(states, left, right, 2)
    graph = nx.Graph()
    graph.add_nodes_from(range(len(adjacency)))
    graph.add_edges_from((u, v) for u, row in enumerate(adjacency) for v in row if u < v)
    assert nx.is_tree(graph) and len(graph) == 139
    return graph


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    tree = build()
    adjacency = [sorted(tree.neighbors(v)) for v in tree]
    tree_poly = independence_poly(len(tree), adjacency)
    candidates = list(nx.non_edges(tree))
    if not args.all:
        core = {v for v in tree if tree.degree(v) >= 3} | {0, 1, 2}
        candidates = [(u, v) for u, v in candidates if u in core or v in core]
    best = None
    for index, (u, v) in enumerate(candidates, 1):
        poly = add_edge_polynomial(tree, tree_poly, u, v)
        score = rebound_score(poly)
        ratio = score[0] / score[1]
        row = {"edge": (u, v), "distance": nx.shortest_path_length(tree, u, v),
               "score": score, "ratio": ratio}
        if best is None or ratio > best[0]:
            best = ratio, row
        if not is_unimodal(poly):
            print({"status": "hit", **row, "coefficients": poly}, flush=True)
            return
        if index % 500 == 0:
            print({"status": "progress", "tested": index, "total": len(candidates),
                   "best": best}, flush=True)
    print({"status": "passed", "tested": len(candidates), "best": best}, flush=True)


if __name__ == "__main__":
    main()
