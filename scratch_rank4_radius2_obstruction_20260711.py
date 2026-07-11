#!/usr/bin/env python3
"""Exact obstruction to a radius-two-only extrapolation of rank-three GSB.

The two trees below have independence number nine and hence rank r=4 lies in
the prefix-GSB window.  They agree on every aggregate used by the audited
rank-three certificate,

    n, S2, S3, S4, P, P21, L,

but have different i_6 and therefore different rank-four GSB gaps.  Thus an
arbitrary-r proof cannot simply reuse that fixed radius-two moment tuple.
"""

from __future__ import annotations

from graph6 import parse_graph6
from indpoly import independence_poly


TREE_0 = [
    [9, 12],
    [10, 11],
    [10],
    [10],
    [11],
    [11],
    [12],
    [12],
    [12],
    [0],
    [1, 2, 3, 12],
    [1, 4, 5],
    [0, 6, 7, 8, 10],
]

TREE_1 = [
    [9],
    [9],
    [10, 11],
    [10],
    [11],
    [12],
    [12],
    [12],
    [12],
    [0, 1, 11, 12],
    [2, 3],
    [2, 4, 9],
    [5, 6, 7, 8, 9],
]

# NetworkX/nauty graph6 encodings of the displayed labelled trees.
GRAPH6 = ("L??????_B_H__y", "L??????o@_DA@{")


def assert_tree(adj: list[list[int]]) -> None:
    n = len(adj)
    assert all(v in adj[u] for u in range(n) for v in adj[u])
    edges = sum(map(len, adj)) // 2
    assert edges == n - 1
    seen = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                stack.append(v)
    assert len(seen) == n


def rank_three_moment_tuple(adj: list[list[int]]) -> tuple[int, ...]:
    degree = [len(neighbors) for neighbors in adj]
    s2 = sum(d**2 for d in degree)
    s3 = sum(d**3 for d in degree)
    s4 = sum(d**4 for d in degree)
    p = 0
    p21 = 0
    for u, neighbors in enumerate(adj):
        for v in neighbors:
            if u < v:
                p += degree[u] * degree[v]
                p21 += degree[u] * degree[v] * (degree[u] + degree[v])
    ell = sum(
        sum(degree[u] for u in neighbors) ** 2 for neighbors in adj
    )
    return len(adj), s2, s3, s4, p, p21, ell


def rank_four_gsb_gap(poly: list[int]) -> int:
    # r=4: 6 i_4 i_6 <= 5 i_5^2 + i_4 i_5.
    return 5 * poly[5] ** 2 + poly[4] * poly[5] - 6 * poly[4] * poly[6]


def main() -> None:
    for encoding, tree in zip(GRAPH6, (TREE_0, TREE_1), strict=True):
        assert_tree(tree)
        n, replay = parse_graph6(encoding.encode("ascii"))
        assert n == 13 and replay == tree

    moment_0 = rank_three_moment_tuple(TREE_0)
    moment_1 = rank_three_moment_tuple(TREE_1)
    assert moment_0 == moment_1 == (13, 66, 240, 1002, 75, 488, 392)

    poly_0 = independence_poly(13, TREE_0)
    poly_1 = independence_poly(13, TREE_1)
    assert poly_0 == [1, 13, 66, 175, 274, 269, 170, 69, 17, 2]
    assert poly_1 == [1, 13, 66, 175, 274, 269, 168, 64, 13, 1]
    assert len(poly_0) - 1 == len(poly_1) - 1 == 9

    gap_0 = rank_four_gsb_gap(poly_0)
    gap_1 = rank_four_gsb_gap(poly_1)
    assert (gap_0, gap_1) == (156_031, 159_319)
    assert gap_0 > 0 and gap_1 > 0

    print(
        {
            "rank_three_moment_tuple": moment_0,
            "independence_polynomials": [poly_0, poly_1],
            "rank_four_gsb_gaps": [gap_0, gap_1],
            "graph6": GRAPH6,
            "status": "passed",
        }
    )


if __name__ == "__main__":
    main()
