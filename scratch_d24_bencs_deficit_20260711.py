#!/usr/bin/env python3
"""Exact D24 certificate for the Bencs--GSB coefficient decomposition.

For ``F=I(T;x)=sum_k a_k x^k``, put

    D_F(x,y) = x F'(x) F(y) - y F(x) F'(y).

Bencs' bivariate Christoffel--Darboux identity gives, for a tree,

    D_F(x,y) = sum_H d_H (x^p y^q-x^q y^p) G_H(x)G_H(y),

where H runs over nonempty connected induced subtrees, p >= q are its
bipartition sizes, d_H=p-q, and G_H=I(T-N[H];x).  This script verifies the
coefficient identity exhaustively on small trees and replays two exact
obstructions to termwise-positive extractions of prefix GSB.
"""

from __future__ import annotations

from fractions import Fraction

from indpoly import independence_poly
from trees import trees


Adj = list[list[int]]


def coefficient(poly: list[int], rank: int) -> int:
    return poly[rank] if 0 <= rank < len(poly) else 0


def adjacency(order: int, edges: list[tuple[int, int]]) -> Adj:
    adj = [[] for _ in range(order)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    for row in adj:
        row.sort()
    return adj


def connected_bipartition(adj: Adj, mask: int) -> tuple[int, int] | None:
    """Return ordered bipartition sizes p >= q, or None if disconnected."""
    root = (mask & -mask).bit_length() - 1
    colors = {root: 0}
    seen = 1 << root
    stack = [root]
    while stack:
        vertex = stack.pop()
        for neighbor in adj[vertex]:
            if not (mask >> neighbor) & 1:
                continue
            if (seen >> neighbor) & 1:
                assert colors[neighbor] != colors[vertex]
                continue
            seen |= 1 << neighbor
            colors[neighbor] = 1 - colors[vertex]
            stack.append(neighbor)
    if seen != mask:
        return None
    first = sum(color == 0 for color in colors.values())
    second = len(colors) - first
    return max(first, second), min(first, second)


def residual_mask(adj: Adj, mask: int) -> int:
    closed = mask
    for vertex in range(len(adj)):
        if (mask >> vertex) & 1:
            for neighbor in adj[vertex]:
                closed |= 1 << neighbor
    return ((1 << len(adj)) - 1) ^ closed


def induced_forest(adj: Adj, mask: int) -> Adj:
    vertices = [v for v in range(len(adj)) if (mask >> v) & 1]
    relabel = {vertex: index for index, vertex in enumerate(vertices)}
    return [
        [relabel[u] for u in adj[v] if u in relabel]
        for v in vertices
    ]


def subtree_records(adj: Adj) -> list[tuple[int, int, int, list[int], int]]:
    """Return (p,q,d,residual polynomial,mask) for connected induced H."""
    records = []
    for mask in range(1, 1 << len(adj)):
        parts = connected_bipartition(adj, mask)
        if parts is None:
            continue
        p, q = parts
        d = p - q
        residual = induced_forest(adj, residual_mask(adj, mask))
        records.append(
            (p, q, d, independence_poly(len(residual), residual), mask)
        )
    return records


def bencs_term(
    record: tuple[int, int, int, list[int], int], i: int, j: int
) -> int:
    p, q, d, residual, _ = record
    return d * (
        coefficient(residual, i - p) * coefficient(residual, j - q)
        - coefficient(residual, i - q) * coefficient(residual, j - p)
    )


def twice_local_gsb_term(
    record: tuple[int, int, int, list[int], int], rank: int
) -> int:
    """Twice the H contribution after polarizing the GSB deficit."""
    return (
        2 * bencs_term(record, rank + 1, rank)
        - (rank + 2) * bencs_term(record, rank + 2, rank)
    )


def gsb_deficit(poly: list[int], rank: int) -> int:
    return (
        (rank + 1) * poly[rank + 1] ** 2
        + poly[rank] * poly[rank + 1]
        - (rank + 2) * poly[rank] * poly[rank + 2]
    )


def convolution(poly: list[int], total: int) -> int:
    return sum(
        poly[i] * coefficient(poly, total - i)
        for i in range(len(poly))
    )


def energy(poly: list[int], total: int) -> Fraction:
    return sum(
        Fraction(
            (2 * i - total) ** 2
            * poly[i]
            * coefficient(poly, total - i),
            2,
        )
        for i in range(len(poly))
    )


def verify_bivariate_identity(max_order: int = 9) -> int:
    """Verify every coefficient of Bencs' identity through max_order."""
    checked = 0
    for order in range(1, max_order + 1):
        for _, adj in trees(order, "networkx"):
            poly = independence_poly(order, adj)
            records = subtree_records(adj)
            for i in range(len(poly)):
                for j in range(len(poly)):
                    rhs = sum(bencs_term(record, i, j) for record in records)
                    assert rhs == (i - j) * poly[i] * poly[j]
                    checked += 1
    return checked


def verify_no_charged_branch_failure(max_order: int = 10) -> int:
    """Exhaust the canonical singleton-charge claim through max_order.

    The square term has the canonical derivative decomposition

        (r+1)a_(r+1)^2 = a_(r+1) sum_v [x^r] I(T-N[v]).

    We assign that charge to the singleton H={v}.  The returned count is the
    number of connected-subtree/rank summands checked at the still-open ranks
    r>=4 in the prefix window.
    """
    checked = 0
    for order in range(1, max_order + 1):
        for _, adj in trees(order, "networkx"):
            poly = independence_poly(order, adj)
            alpha = len(poly) - 1
            tail_start = (2 * alpha + 1) // 3
            records = subtree_records(adj)
            for rank in range(4, max(4, tail_start - 1)):
                local_sum = 0
                for record in records:
                    p, q, _, residual, mask = record
                    del p, q
                    term = twice_local_gsb_term(record, rank)
                    if mask & (mask - 1) == 0:
                        term += (
                            2
                            * poly[rank + 1]
                            * coefficient(residual, rank)
                        )
                    assert term >= 0
                    local_sum += term
                    checked += 1
                assert local_sum == 2 * gsb_deficit(poly, rank)
    return checked


def star_obstruction() -> dict[str, object]:
    """Replay the first possible rank-four energy-core obstruction."""
    rank = 4
    adj = adjacency(10, [(0, leaf) for leaf in range(1, 10)])
    poly = independence_poly(10, adj)
    assert poly == [1, 10, 36, 84, 126, 126, 84, 36, 9, 1]
    alpha = len(poly) - 1
    tail_start = (2 * alpha + 1) // 3
    assert alpha == 9 and tail_start == 6 and rank == tail_start - 2

    # H is one leaf.  T-N[H] consists of eight isolated vertices.
    records = subtree_records(adj)
    leaf_mask = 1 << 1
    record = next(record for record in records if record[-1] == leaf_mask)
    p, q, d, residual, _ = record
    assert (p, q, d) == (1, 0, 1)
    assert residual == [1, 8, 28, 56, 70, 56, 28, 8, 1]
    first = bencs_term(record, rank + 1, rank)
    second = bencs_term(record, rank + 2, rank)
    local_twice = twice_local_gsb_term(record, rank)
    assert (first, second, local_twice) == (1764, 2352, -10584)

    # The diagonal-energy extraction with only nonnegative outer corrections:
    #
    # Delta = K + sum_{t>=1} A_t[(3r+4)(t+1)^2-2(r+1)]
    #             + sum_{t>=2} B_t[((2t+1)^2-9)/8].
    even_total = 2 * rank + 2
    odd_total = 2 * rank + 1
    core = (
        (rank + 1) * convolution(poly, even_total)
        + Fraction(9, 16) * convolution(poly, odd_total)
        - Fraction(3 * rank + 4, 4) * energy(poly, even_total)
        - Fraction(1, 8) * energy(poly, odd_total)
    )
    deficit = gsb_deficit(poly, rank)
    assert core == -180471
    assert deficit == 31752

    even_correction = sum(
        (
            (3 * rank + 4) * (t + 1) ** 2 - 2 * (rank + 1)
        )
        * poly[rank - t]
        * poly[rank + 2 + t]
        for t in range(1, rank + 1)
        if rank + 2 + t < len(poly)
    )
    odd_correction = sum(
        Fraction((2 * t + 1) ** 2 - 9, 8)
        * poly[rank - t]
        * poly[rank + 1 + t]
        for t in range(2, rank + 1)
        if rank + 1 + t < len(poly)
    )
    assert even_correction == 209172
    assert odd_correction == 3051
    assert core + even_correction + odd_correction == deficit

    return {
        "tree": "K_1,9",
        "rank": rank,
        "alpha": alpha,
        "tail_start": tail_start,
        "leaf_local_twice": local_twice,
        "energy_core": int(core),
        "positive_outer_correction": int(even_correction + odd_correction),
        "GSB_deficit": deficit,
    }


def charged_broom_obstruction() -> dict[str, object]:
    """Replay the order-11 failure of the canonical singleton charge."""
    rank = 4
    edges = (
        [(0, 1), (0, 9), (9, 10)]
        + [(1, leaf) for leaf in range(2, 9)]
    )
    adj = adjacency(11, edges)
    poly = independence_poly(11, adj)
    assert poly == [1, 11, 45, 105, 161, 161, 105, 43, 10, 1]
    alpha = len(poly) - 1
    tail_start = (2 * alpha + 1) // 3
    assert alpha == 9 and tail_start == 6 and rank == tail_start - 2

    # The induced path H=0-9-10 has bipartition 2+1.  Its closed
    # neighborhood also removes vertex 1, leaving seven isolated vertices.
    mask = (1 << 0) | (1 << 9) | (1 << 10)
    record = next(record for record in subtree_records(adj) if record[-1] == mask)
    p, q, d, residual, _ = record
    assert (p, q, d) == (2, 1, 1)
    assert residual == [1, 7, 21, 35, 35, 21, 7, 1]
    first = bencs_term(record, rank + 1, rank)
    second = bencs_term(record, rank + 2, rank)
    local_twice = twice_local_gsb_term(record, rank)
    assert (first, second, local_twice) == (490, 784, -3724)

    # H is not a singleton, so the derivative allocation assigns it no share
    # of (r+1)a_(r+1)^2.  This is therefore also the charged local term.
    assert mask & (mask - 1)
    assert gsb_deficit(poly, rank) == 54096
    return {
        "tree": "seven-leaf star with one arm extended by two vertices",
        "edges": edges,
        "rank": rank,
        "alpha": alpha,
        "tail_start": tail_start,
        "H": [0, 9, 10],
        "H_bipartition": [p, q],
        "residual_polynomial": residual,
        "Bencs_adjacent_term": first,
        "Bencs_two_step_term": second,
        "charged_local_twice": local_twice,
        "GSB_deficit": gsb_deficit(poly, rank),
    }


def main() -> None:
    coefficient_checks = verify_bivariate_identity()
    charged_checks = verify_no_charged_branch_failure()
    star = star_obstruction()
    broom = charged_broom_obstruction()
    print(
        {
            "bivariate_coefficient_checks_n_le_9": coefficient_checks,
            "nonnegative_charged_summands_n_le_10": charged_checks,
            "star_obstruction": star,
            "charged_broom_obstruction": broom,
            "certificate": "passed",
        }
    )


if __name__ == "__main__":
    main()
