#!/usr/bin/env python3
"""Exact obstruction to fixed-alpha star-corona extremality for D13.

Let D_{p,q} be a double star with adjacent hubs carrying p and q leaves.
Replace every quotient vertex by a matching edge.  Attach every quotient-leaf
edge to the hub endpoint of its hub block, but attach the hub-hub bridge to the
*leaf* endpoint of the right hub block.  This is the misaligned lift M_{p,q}.

At (p,q)=(20,21), alpha=43 and prefix rank r=10, its prefix-GSB ratio is
strictly larger than that of the proposed fixed-alpha extremizer
C_alpha=S(2^(alpha-1),1).  Both ratios are strictly below one, so this
obstructs only the compression/extremal lemma, not prefix GSB.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

from indpoly import _polyadd, _polymul, independence_poly
from scratch_d13_block_extremal_20260711 import (
    corona_star_poly,
    gsb_ratio,
    prefix_last,
)


def binomial_poly(n: int, scale: int = 1) -> list[int]:
    """Coefficients of (1+scale*x)^n."""
    return [comb(n, j) * scale**j for j in range(n + 1)]


def shift(poly: list[int], amount: int = 1) -> list[int]:
    return [0] * amount + poly


def misaligned_double_star_poly(p: int, q: int) -> list[int]:
    """Closed independence polynomial of M_{p,q}."""
    one_p = binomial_poly(p)
    one_q = binomial_poly(q)
    one_pq = binomial_poly(p + q)
    two_p = binomial_poly(p, 2)
    two_q = binomial_poly(q, 2)
    two_pq = binomial_poly(p + q, 2)

    # (1+x)^2(1+2x)^(p+q)
    term0 = _polymul(binomial_poly(2), two_pq)
    # x(1+x)^(q+1)(1+2x)^p
    term1 = shift(_polymul(binomial_poly(q + 1), two_p))
    # x(1+x)^p(1+2x)^q
    term2 = shift(_polymul(one_p, two_q))
    # x^2(1+x)^(p+q)
    term3 = shift(one_pq, 2)
    return _polyadd(_polyadd(term0, term1), _polyadd(term2, term3))


def misaligned_double_star_adj(p: int, q: int) -> list[list[int]]:
    """Build M_{p,q}, with its displayed perfect matching."""
    # Hub blocks are (0,1) and (2,3).  The bridge is 0--3, while all
    # quotient-leaf branches at the right hub attach to vertex 2.
    edges: list[tuple[int, int]] = [(0, 1), (2, 3), (0, 3)]
    next_vertex = 4
    for hub, count in ((0, p), (2, q)):
        for _ in range(count):
            base, leaf = next_vertex, next_vertex + 1
            next_vertex += 2
            edges.extend(((base, leaf), (hub, base)))
    adj = [[] for _ in range(next_vertex)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return adj


def main() -> None:
    p, q, rank = 20, 21, 10
    alpha = p + q + 2
    adj = misaligned_double_star_adj(p, q)
    assert len(adj) == 2 * alpha == 86
    assert sum(map(len, adj)) // 2 == len(adj) - 1

    # The displayed alpha internal block edges form a perfect matching, so a
    # bipartite tree has matching number alpha and independence number alpha.
    poly = independence_poly(len(adj), adj)
    formula = misaligned_double_star_poly(p, q)
    assert poly == formula
    assert len(poly) - 1 == alpha
    assert rank <= prefix_last(alpha) == 27

    extremal = corona_star_poly(alpha)
    got = gsb_ratio(poly, rank)
    proposed_bound = gsb_ratio(extremal, rank)

    assert poly[rank:rank + 3] == [
        1557817164957,
        9071748758487,
        46964708916216,
    ]
    assert extremal[rank:rank + 3] == [
        1735500102882,
        10274818745373,
        54064917032672,
    ]
    assert got == Fraction(
        155836993649129940784,
        163193680184864456073,
    )
    assert proposed_bound == Fraction(8989217536, 9413681835)
    assert got > proposed_bound
    assert got < 1

    comparison_cross = (
        got.numerator * proposed_bound.denominator
        - proposed_bound.numerator * got.denominator
    )
    gsb_gap = (
        (rank + 1) * poly[rank + 1] ** 2
        + poly[rank] * poly[rank + 1]
        - (rank + 2) * poly[rank] * poly[rank + 2]
    )
    assert comparison_cross == 16384653665596798453162512
    assert gsb_gap == 41445850477678897261648374

    print({
        "certificate": "passed",
        "p": p,
        "q": q,
        "n": len(adj),
        "alpha": alpha,
        "rank": rank,
        "prefix_last": prefix_last(alpha),
        "coefficients": poly[rank:rank + 3],
        "star_corona_coefficients": extremal[rank:rank + 3],
        "ratio": [got.numerator, got.denominator],
        "star_corona_ratio": [
            proposed_bound.numerator, proposed_bound.denominator,
        ],
        "comparison_cross": comparison_cross,
        "gsb_gap": gsb_gap,
    })


if __name__ == "__main__":
    main()
