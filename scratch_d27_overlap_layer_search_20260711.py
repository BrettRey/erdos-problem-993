#!/usr/bin/env python3
"""Exact common-set-layer audit for the D27 far covariance regrouping.

For ordered rank-r sets A,B, group the far quadratic by C=A intersection B.
The direct pair loop is avoided by an upper zeta transform over the Boolean
lattice, followed by Moebius inversion.  At n=16 all arrays have only 65,536
rows, so exhaustive tree scans are practical.
"""

from __future__ import annotations

import argparse
from itertools import combinations
from math import ceil

import numpy as np

from indpoly import independence_poly
from scratch_d18_covariance_search_20260711 import covariance_profile
from trees import trees


def independent_rank_masks(adj: list[list[int]], rank: int) -> list[int]:
    edge_masks = [sum(1 << neighbor for neighbor in neighbors) for neighbors in adj]
    out = []
    for vertices in combinations(range(len(adj)), rank):
        mask = sum(1 << vertex for vertex in vertices)
        if all(not (edge_masks[vertex] & mask) for vertex in vertices):
            out.append(mask)
    return out


def far_matrix(adj: list[list[int]]) -> np.ndarray:
    order = len(adj)
    matrix = np.zeros((order, order), dtype=np.int64)
    for source in range(order):
        distance = [-1] * order
        distance[source] = 0
        queue = [source]
        for vertex in queue:
            for neighbor in adj[vertex]:
                if distance[neighbor] < 0:
                    distance[neighbor] = distance[vertex] + 1
                    queue.append(neighbor)
        for target in range(source + 1, order):
            if distance[target] >= 3:
                matrix[source, target] = matrix[target, source] = 1
    return matrix


def upper_zeta(array: np.ndarray, order: int) -> None:
    """In place: f(C) <- sum_(A superset C) f(A)."""
    tail_shape = array.shape[1:]
    for bit in range(order):
        block = 1 << bit
        view = array.reshape((-1, 2, block, *tail_shape))
        view[:, 0] += view[:, 1]


def upper_mobius(array: np.ndarray, order: int) -> None:
    """Invert ``upper_zeta`` in place."""
    tail_shape = array.shape[1:]
    for bit in range(order):
        block = 1 << bit
        view = array.reshape((-1, 2, block, *tail_shape))
        view[:, 0] -= view[:, 1]


def exact_layers(adj: list[list[int]], rank: int) -> np.ndarray:
    order = len(adj)
    size = 1 << order
    far = far_matrix(adj)
    count = np.zeros(size, dtype=np.int64)
    marked_far = np.zeros(size, dtype=np.int64)
    marked_vertex = np.zeros((size, order), dtype=np.int64)

    closed = [
        (1 << vertex) | sum(1 << neighbor for neighbor in adj[vertex])
        for vertex in range(order)
    ]
    for mask in independent_rank_masks(adj, rank):
        addable = np.asarray(
            [int(not (mask & neighborhood)) for neighborhood in closed],
            dtype=np.int64,
        )
        count[mask] = 1
        marked_vertex[mask] = addable
        marked_far[mask] = int(addable @ far @ addable // 2)

    upper_zeta(count, order)
    upper_zeta(marked_far, order)
    upper_zeta(marked_vertex, order)
    cross = np.sum((marked_vertex @ far) * marked_vertex, axis=1) // 2
    cumulative = count * marked_far - cross
    layers = cumulative.copy()
    upper_mobius(layers, order)
    return layers


def audit_tree(adj: list[list[int]], rank: int) -> tuple[int, int | None]:
    layers = exact_layers(adj, rank)
    positive = np.flatnonzero(layers > 0)
    profile_qfar = covariance_profile(adj)["rows"][rank]["qfar"]
    if int(np.sum(layers)) != profile_qfar:
        raise AssertionError((int(np.sum(layers)), profile_qfar))
    if len(positive):
        index = int(positive[np.argmax(layers[positive])])
        return int(layers[index]), index
    return int(np.max(layers)), None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--min-n', type=int, default=9)
    parser.add_argument('--max-n', type=int, default=16)
    args = parser.parse_args()
    best = None
    for order in range(args.min_n, args.max_n + 1):
        tree_count = 0
        for graph6, adj in trees(order):
            poly = independence_poly(order, adj)
            alpha = len(poly) - 1
            last = ceil((2 * alpha - 1) / 3) - 2
            for rank in range(4, last + 1):
                value, common = audit_tree(adj, rank)
                row = (value, graph6, alpha, rank, common)
                if best is None or row[0] > best[0]:
                    best = row
                if common is not None:
                    print(
                        {
                            'failure': row,
                            'adj': adj,
                            'poly': poly,
                        },
                        flush=True,
                    )
                    raise SystemExit(1)
            tree_count += 1
        print({'n': order, 'trees': tree_count, 'best_layer': best}, flush=True)
    print({'failure': None, 'best_layer': best}, flush=True)


if __name__ == '__main__':
    main()
