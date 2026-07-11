#!/usr/bin/env python3
"""Exact local-moment formula for rank-three extension counts in a tree.

For an independent triple S, e(S) is the number of vertices that can be
adjoined to S.  The function below computes

    N  = # independent triples,
    d1 = sum_S e(S),
    d2 = sum_S e(S)(e(S)-1)

using only degree and radius-two sums.  It is a replay/falsification artifact
for the July 2026 variance route, not a general-purpose implementation.
"""

from __future__ import annotations

from itertools import combinations
from math import comb

from scratch_extension_variance_dp_20260711 import extension_moment_jets
from trees import trees


def rank3_local_moments(
    adj: list[list[int]],
) -> tuple[int, int, int, dict[str, int]]:
    n = len(adj)
    t = n - 1
    degree = [len(neighbors) for neighbors in adj]
    neighbor_degree = [
        sum(degree[u] for u in neighbors) for neighbors in adj
    ]

    s2 = sum(d * d for d in degree)
    s3 = sum(d * d * d for d in degree)
    edge_product = sum(
        degree[u] * degree[v]
        for u, neighbors in enumerate(adj)
        for v in neighbors
        if u < v
    )
    edge_product_sum = sum(
        degree[u] * degree[v] * (degree[u] + degree[v])
        for u, neighbors in enumerate(adj)
        for v in neighbors
        if u < v
    )

    # j_v counts independent pairs in T-N[v], hence independent triples
    # containing v.
    j = [
        comb(t - degree[v], 2) - (t - neighbor_degree[v])
        for v in range(n)
    ]
    n3 = sum(j) // 3

    a1 = sum(degree[v] * j[v] for v in range(n))
    p1 = (
        t * sum(comb(d, 2) for d in degree)
        - 2 * edge_product
        + s2
    )
    q1 = sum(comb(d, 3) for d in degree)
    x1 = a1 - p1 + q1

    a2_single = sum(degree[v] ** 2 * j[v] for v in range(n))
    distance2_product = 0
    ap = 0
    for center, neighbors in enumerate(adj):
        for u, v in combinations(neighbors, 2):
            distance2_product += degree[u] * degree[v]
            extensions = t - degree[u] - degree[v]
            extension_degree_sum = (
                2 * t
                - degree[u]
                - degree[v]
                - neighbor_degree[u]
                - neighbor_degree[v]
                + degree[center]
            )
            ap += (
                (degree[u] + degree[v]) * extensions
                + extension_degree_sum
            )

    nonedge_degree_product = 2 * t * t - s2 // 2 - edge_product
    nonedge_degree_product_sum = (
        2 * t * s2 - s3 - edge_product_sum
    )
    pair_weighted_extensions = (
        (t - 1) * nonedge_degree_product
        - nonedge_degree_product_sum
        + distance2_product
    )
    a2 = a2_single + 2 * pair_weighted_extensions

    distance2_degree = [
        neighbor_degree[v] - degree[v] for v in range(n)
    ]
    p2 = p1 + 2 * sum(comb(z, 2) for z in distance2_degree)
    aq = sum(
        comb(degree[v] - 1, 2) * neighbor_degree[v]
        for v in range(n)
    )
    x2 = a2 + p2 - 2 * ap + 2 * aq - 5 * q1

    base = t - 2
    d1 = n3 * base - x1
    d2 = n3 * base * (base - 1) - (2 * base - 1) * x1 + x2
    pieces = {
        "A1": a1,
        "P1": p1,
        "Q1": q1,
        "X1": x1,
        "A2": a2,
        "P2": p2,
        "AP": ap,
        "AQ": aq,
        "X2": x2,
    }
    return n3, d1, d2, pieces


def main() -> None:
    for n in range(3, 15):
        count = 0
        for _, adj in trees(n, backend="geng"):
            got = rank3_local_moments(adj)[:3]
            jet = extension_moment_jets(adj)
            expected = tuple(
                polynomial[3] if len(polynomial) > 3 else 0
                for polynomial in (jet.value, jet.first, jet.second)
            )
            if got != expected:
                raise AssertionError((n, got, expected, adj))
            count += 1
        print({"n": n, "trees": count, "formula_mismatches": 0}, flush=True)


if __name__ == "__main__":
    main()
