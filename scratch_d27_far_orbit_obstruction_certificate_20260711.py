#!/usr/bin/env python3
"""Exact certificate for two insufficient D27 switching regroupings.

Neither a complete symmetric-difference orbit nor an exact-common-set layer
has a fixed nonpositive far quadratic, even in the required prefix.  The
global qfar remains negative on both witnesses.
"""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations

from indpoly import independence_poly
from scratch_d18_covariance_search_20260711 import covariance_profile
from scratch_d27_overlap_layer_search_20260711 import exact_layers


ORBIT_TREE = [
    [9, 10], [9], [9], [10], [10], [10], [10], [10], [10],
    [0, 1, 2], [0, 3, 4, 5, 6, 7, 8],
]
ORBIT_U = sum(1 << vertex for vertex in (0, 2, 4, 5, 6, 7, 8, 9))

COMMON_TREE = [
    [10, 11, 12], [10], [10], [10], [11], [11], [11], [12], [12], [12],
    [0, 1, 2, 3], [0, 4, 5, 6], [0, 7, 8, 9],
]


def independent_sets(adj: list[list[int]], rank: int) -> list[int]:
    out = []
    for vertices in combinations(range(len(adj)), rank):
        mask = sum(1 << vertex for vertex in vertices)
        if all(
            not (mask & (1 << u) and mask & (1 << v))
            for u, neighbors in enumerate(adj)
            for v in neighbors
            if u < v
        ):
            out.append(mask)
    return out


def addable_mask(adj: list[list[int]], chosen: int) -> int:
    return sum(
        1 << vertex
        for vertex, neighbors in enumerate(adj)
        if not (chosen & (1 << vertex))
        and all(not (chosen & (1 << neighbor)) for neighbor in neighbors)
    )


def far_pairs(adj: list[list[int]]) -> list[tuple[int, int]]:
    out = []
    for source in range(len(adj)):
        distance = [-1] * len(adj)
        distance[source] = 0
        queue = [source]
        for vertex in queue:
            for neighbor in adj[vertex]:
                if distance[neighbor] < 0:
                    distance[neighbor] = distance[vertex] + 1
                    queue.append(neighbor)
        out.extend(
            (source, target)
            for target in range(source + 1, len(adj))
            if distance[target] >= 3
        )
    return out


def orbit_certificate() -> dict:
    adj = ORBIT_TREE
    rank = 4
    sets = independent_sets(adj, rank)
    addable = {chosen: addable_mask(adj, chosen) for chosen in sets}
    far = far_pairs(adj)
    orbit_twice: defaultdict[tuple[int, int], int] = defaultdict(int)
    orbit_size: defaultdict[tuple[int, int], int] = defaultdict(int)
    for left in sets:
        for right in sets:
            difference = [
                ((addable[left] >> vertex) & 1)
                - ((addable[right] >> vertex) & 1)
                for vertex in range(len(adj))
            ]
            key = (left & right, left ^ right)
            orbit_twice[key] += sum(
                difference[u] * difference[v] for u, v in far
            )
            orbit_size[key] += 1

    key = (0, ORBIT_U)
    assert orbit_size[key] == 20
    assert orbit_twice[key] == 20
    profile = covariance_profile(adj)
    assert profile['alpha'] == 9 and profile['limit'] - 2 == rank
    assert profile['rows'][rank]['qfar'] == -8700

    # The induced symmetric-difference components have absolute color
    # imbalances 1,1,1,1,1,1 (one P3 and five singletons).
    assert independence_poly(len(adj), adj) == [
        1, 11, 45, 100, 146, 141, 90, 37, 9, 1,
    ]
    return {
        'n': len(adj),
        'alpha': profile['alpha'],
        'rank': rank,
        'balanced_orientations': orbit_size[key],
        'orbit_contribution': orbit_twice[key] // 2,
        'global_qfar': profile['rows'][rank]['qfar'],
    }


def common_layer_certificate() -> dict:
    adj = COMMON_TREE
    rank = 5
    poly = independence_poly(len(adj), adj)
    assert poly == [1, 13, 66, 175, 279, 300, 228, 123, 45, 10, 1]
    layers = exact_layers(adj, rank)
    assert int(layers[0]) == 144
    by_overlap = [0] * (rank + 1)
    for common, value in enumerate(layers):
        if value:
            by_overlap[common.bit_count()] += int(value)
    assert by_overlap == [144, -4743, -21717, -17946, -2718, 0]
    profile = covariance_profile(adj)
    assert profile['alpha'] == 10 and profile['limit'] - 2 == rank
    assert int(layers.sum()) == profile['rows'][rank]['qfar'] == -46980
    return {
        'n': len(adj),
        'alpha': profile['alpha'],
        'rank': rank,
        'empty_common_layer': int(layers[0]),
        'by_overlap': by_overlap,
        'global_qfar': profile['rows'][rank]['qfar'],
    }


def main() -> None:
    print(
        {
            'complete_orbit': orbit_certificate(),
            'exact_common_layer': common_layer_certificate(),
            'certificate': 'passed',
        },
        flush=True,
    )


if __name__ == '__main__':
    main()
