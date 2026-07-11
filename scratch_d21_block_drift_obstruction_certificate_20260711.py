#!/usr/bin/env python3
"""Exact certificate for the D21 matching-block drift obstruction.

The requested averaged drift is valid on this example, but the smallest
prefix star K_{1,4} rules out a pointwise containment-edge/involution proof
with loss at most three.  It also checks the exact size-biased identity.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import ceil

from indpoly import independence_poly


def star(leaves: int) -> list[list[int]]:
    adj = [[] for _ in range(leaves + 1)]
    for leaf in range(1, leaves + 1):
        adj[0].append(leaf)
        adj[leaf].append(0)
    return adj


def independent(adj: list[list[int]], chosen: frozenset[int]) -> bool:
    return all(v not in chosen for u in chosen for v in adj[u])


def independent_sets(
    adj: list[list[int]], rank: int
) -> list[frozenset[int]]:
    return [
        frozenset(chosen)
        for chosen in combinations(range(len(adj)), rank)
        if independent(adj, frozenset(chosen))
    ]


def addable(adj: list[list[int]], chosen: frozenset[int]) -> frozenset[int]:
    blocked = set(chosen)
    for u in chosen:
        blocked.update(adj[u])
    return frozenset(set(range(len(adj))) - blocked)


def residual_edges(
    adj: list[list[int]], chosen: frozenset[int]
) -> tuple[tuple[int, int], ...]:
    available = addable(adj, chosen)
    return tuple(
        (u, v)
        for u in available
        for v in adj[u]
        if u < v and v in available
    )


def mean(values: list[int]) -> Fraction:
    return Fraction(sum(values), len(values))


def main() -> None:
    adj = star(4)
    polynomial = independence_poly(len(adj), adj)
    alpha = len(polynomial) - 1
    rank = 1
    prefix_last = ceil((2 * alpha - 1) / 3) - 2

    assert polynomial == [1, 5, 6, 4, 1]
    assert alpha == 4 and rank == prefix_last

    # One maximum-matching block and three unmatched singleton blocks.
    blocks = ((0, 1), (2,), (3,), (4,))
    assert len(blocks) == alpha
    assert {v for block in blocks for v in block} == set(range(5))

    levels = {
        k: independent_sets(adj, k)
        for k in (rank, rank + 1)
    }

    def e(chosen: frozenset[int]) -> int:
        return len(addable(adj, chosen))

    def b(chosen: frozenset[int]) -> int:
        return 2 * (alpha - len(chosen)) - e(chosen)

    # Every independent set occupies exactly |A| displayed blocks.
    for level in levels.values():
        for chosen in level:
            occupied = sum(bool(chosen.intersection(block)) for block in blocks)
            assert occupied == len(chosen)

    rank_one = levels[rank]
    rank_two = levels[rank + 1]
    center = frozenset({0})
    assert sorted(e(chosen) for chosen in rank_one) == [0, 3, 3, 3, 3]
    assert sorted(b(chosen) for chosen in rank_one) == [3, 3, 3, 3, 6]
    assert all(e(chosen) == 2 and b(chosen) == 2 for chosen in rank_two)

    containment_edges = [
        (chosen, chosen | {vertex})
        for chosen in rank_one
        for vertex in addable(adj, chosen)
    ]
    assert len(containment_edges) == 12
    assert all(parent != center for parent, _ in containment_edges)

    # The deterministic local identity holds on every containment edge.
    for parent, child in containment_edges:
        vertex = next(iter(child - parent))
        residual_degree = sum(
            neighbor in addable(adj, parent) for neighbor in adj[vertex]
        )
        assert b(child) - b(parent) == residual_degree - 1 == -1

    # Nevertheless, no map from the center set to any rank-two set can have
    # pointwise loss <= 3:  b({center})=6 while every child-level b is 2.
    assert b(center) == 6
    assert max(map(b, rank_two)) == 2 < b(center) - 3

    mu = mean([e(chosen) for chosen in rank_one])
    mean_b_one = mean([b(chosen) for chosen in rank_one])
    mean_b_two = mean([b(chosen) for chosen in rank_two])
    mean_b_size_biased = Fraction(
        sum(e(chosen) * b(chosen) for chosen in rank_one),
        sum(e(chosen) for chosen in rank_one),
    )
    second_e = mean([e(chosen) ** 2 for chosen in rank_one])
    variance = second_e - mu * mu
    eta = mean([len(residual_edges(adj, chosen)) for chosen in rank_one])

    assert mu == Fraction(12, 5)
    assert mean_b_one == Fraction(18, 5)
    assert mean_b_two == 2
    assert mean_b_size_biased == 3
    assert variance == Fraction(36, 25)
    assert eta == 0

    # Exact edge-size-biased identity and requested averaged drift.
    assert mean_b_size_biased == mean_b_one - variance / mu
    assert mean_b_two == mean_b_one - 1 + (2 * eta - variance) / mu
    assert mean_b_two - mean_b_one == Fraction(-8, 5) >= -3

    print(
        {
            "certificate": "passed",
            "tree": "K_1,4",
            "n": len(adj),
            "alpha": alpha,
            "rank": rank,
            "prefix_last": prefix_last,
            "independence_polynomial": polynomial,
            "rank_one_e": sorted(e(chosen) for chosen in rank_one),
            "rank_one_b": sorted(b(chosen) for chosen in rank_one),
            "rank_two_b": sorted(b(chosen) for chosen in rank_two),
            "mean_e": str(mu),
            "variance_e": str(variance),
            "mean_h": str(eta),
            "mean_b_rank_one": str(mean_b_one),
            "mean_b_size_biased": str(mean_b_size_biased),
            "mean_b_rank_two": str(mean_b_two),
            "mean_drift": str(mean_b_two - mean_b_one),
            "pointwise_loss_from_center": b(center) - max(map(b, rank_two)),
        }
    )


if __name__ == "__main__":
    main()
