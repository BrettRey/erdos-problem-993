#!/usr/bin/env python3
"""Exact search for nonunimodal Bencs stable-path trees.

For an ordered graph G and root u, Bencs' stable-path tree recursively attaches
the stable-path trees rooted at the ordered neighbours u_i inside
G-{u,u_1,...,u_{i-1}}.  This driver computes the rooted-tree DP directly from
that recursion, without materialising the (potentially enormous) tree, and
checks the resulting exact coefficient sequence for a valley.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import permutations

import networkx as nx

from indpoly import _polyadd, _polymul, is_unimodal


def product(polys: list[list[int]]) -> list[int]:
    out = [1]
    for poly in polys:
        out = _polymul(out, poly)
    return out


def stable_path_states(
    adjacency: tuple[int, ...], order: tuple[int, ...], root: int
) -> tuple[list[int], list[int], int]:
    """Return exclude-root, select-root, and order of T^<_{G,root}."""

    position = [0] * len(order)
    for index, vertex in enumerate(order):
        position[vertex] = index

    @lru_cache(maxsize=None)
    def recurse(mask: int, vertex: int) -> tuple[tuple[int, ...], tuple[int, ...], int]:
        neighbors_mask = adjacency[vertex] & mask
        neighbors = sorted(
            (v for v in range(len(adjacency)) if neighbors_mask >> v & 1),
            key=position.__getitem__,
        )
        removed = 1 << vertex
        children = []
        for neighbor in neighbors:
            child_mask = mask & ~removed
            children.append(recurse(child_mask, neighbor))
            removed |= 1 << neighbor
        excluded = product([list(total) for total, _, _ in children])
        selected = [0] + product([list(exc) for _, exc, _ in children])
        total = _polyadd(excluded, selected)
        tree_order = 1 + sum(size for _, _, size in children)
        return tuple(total), tuple(excluded), tree_order

    total, excluded, tree_order = recurse((1 << len(adjacency)) - 1, root)
    selected = [
        value - (excluded[index] if index < len(excluded) else 0)
        for index, value in enumerate(total)
    ]
    return list(excluded), selected, tree_order


def graph_masks(graph: nx.Graph) -> tuple[int, ...]:
    vertices = sorted(graph)
    relabel = {vertex: index for index, vertex in enumerate(vertices)}
    return tuple(
        sum(1 << relabel[w] for w in graph.neighbors(vertex))
        for vertex in vertices
    )


def valleys(coefficients: list[int]) -> list[tuple[int, int, int]]:
    seen_descent = False
    result = []
    for index in range(len(coefficients) - 1):
        if coefficients[index + 1] < coefficients[index]:
            seen_descent = True
        elif seen_descent and coefficients[index + 1] > coefficients[index]:
            result.append((index, coefficients[index], coefficients[index + 1]))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas-max-order", type=int, default=7)
    parser.add_argument("--all-orders-through", type=int, default=7)
    parser.add_argument("--max-orderings", type=int, default=5040)
    args = parser.parse_args()

    tested_graphs = tested_instances = 0
    best = None
    for graph in nx.graph_atlas_g():
        n = len(graph)
        if n < 2 or n > args.atlas_max_order or not nx.is_connected(graph):
            continue
        adjacency = graph_masks(graph)
        graph_poly = None
        # Nonunimodal bases are prioritized, but all connected atlas graphs are
        # valid because induced factors can themselves introduce a valley.
        tested_graphs += 1
        order_iter = permutations(range(n))
        for order_index, order in enumerate(order_iter):
            if order_index >= args.max_orderings:
                break
            if n > args.all_orders_through and order_index:
                break
            for root in range(n):
                excluded, selected, tree_order = stable_path_states(
                    adjacency, order, root
                )
                coefficients = _polyadd(excluded, selected)
                tested_instances += 1
                pressure = 0.0
                seen_descent = False
                for index in range(len(coefficients) - 1):
                    if coefficients[index + 1] < coefficients[index]:
                        seen_descent = True
                    elif seen_descent and coefficients[index]:
                        pressure = max(
                            pressure, coefficients[index + 1] / coefficients[index]
                        )
                item = (pressure, n, order, root, tree_order, coefficients)
                if best is None or item[0] > best[0]:
                    best = item
                if not is_unimodal(coefficients):
                    print(
                        {
                            "counterexample": True,
                            "base_graph6": nx.to_graph6_bytes(
                                graph, header=False
                            ).decode().strip(),
                            "base_order": n,
                            "vertex_order": order,
                            "root": root,
                            "tree_order": tree_order,
                            "coefficients": coefficients,
                            "valleys": valleys(coefficients),
                        },
                        flush=True,
                    )
                    return
        if tested_graphs % 100 == 0:
            print(
                {
                    "tested_graphs": tested_graphs,
                    "tested_instances": tested_instances,
                    "best": best,
                },
                flush=True,
            )
    print(
        {
            "counterexample": False,
            "tested_graphs": tested_graphs,
            "tested_instances": tested_instances,
            "best": best,
        },
        flush=True,
    )


if __name__ == "__main__":
    main()
