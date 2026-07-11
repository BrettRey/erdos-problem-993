#!/usr/bin/env python3
"""Exact low-rank search for a prefix obstruction to aggregate far covariance.

The full D18 profile deliberately computes every coefficient.  This driver
truncates every cavity polynomial at ``rank + 2`` and is therefore suitable
for evolving large leaf blow-ups at the sharp rank-four boundary.
"""

from __future__ import annotations

from functools import lru_cache
from math import comb


def add(left: list[int], right: list[int], cap: int) -> list[int]:
    out = [0] * min(cap + 1, max(len(left), len(right)))
    for index, value in enumerate(left[: cap + 1]):
        out[index] += value
    for index, value in enumerate(right[: cap + 1]):
        out[index] += value
    return out


def mul(left: list[int], right: list[int], cap: int) -> list[int]:
    out = [0] * min(cap + 1, len(left) + len(right) - 1)
    for i, a_i in enumerate(left):
        if not a_i:
            continue
        for j, b_j in enumerate(right[: cap + 1 - i]):
            out[i + j] += a_i * b_j
    return out


def product(polys: list[list[int]], cap: int) -> list[int]:
    out = [1]
    for poly in polys:
        out = mul(out, poly, cap)
    return out


def coefficient(poly: list[int], rank: int) -> int:
    return poly[rank] if rank < len(poly) else 0


def cavity_messages(
    adj: list[list[int]], cap: int
) -> tuple[dict[tuple[int, int], list[int]], dict[tuple[int, int], list[int]]]:
    @lru_cache(maxsize=None)
    def message(vertex: int, parent: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
        children = [message(w, vertex) for w in adj[vertex] if w != parent]
        excluded = product([list(total) for total, _ in children], cap)
        selected_tail = product([list(exc) for _, exc in children], cap - 1)
        total = add(excluded, [0, *selected_tail], cap)
        return tuple(total), tuple(excluded)

    total: dict[tuple[int, int], list[int]] = {}
    excluded: dict[tuple[int, int], list[int]] = {}
    for vertex, neighbors in enumerate(adj):
        for parent in neighbors:
            p, e = message(vertex, parent)
            total[(vertex, parent)] = list(p)
            excluded[(vertex, parent)] = list(e)
    return total, excluded


def tree_polynomial(adj: list[list[int]], cap: int) -> list[int]:
    if not adj:
        return [1]
    total, excluded = cavity_messages(adj, cap)
    root = 0
    children = [w for w in adj[root]]
    root_excluded = product([total[(w, root)] for w in children], cap)
    root_selected = product([excluded[(w, root)] for w in children], cap - 1)
    return add(root_excluded, [0, *root_selected], cap)


def independence_number(adj: list[list[int]]) -> int:
    def visit(vertex: int, parent: int) -> tuple[int, int]:
        states = [visit(w, vertex) for w in adj[vertex] if w != parent]
        excluded = sum(max(a, b) for a, b in states)
        selected = 1 + sum(a for a, _ in states)
        return excluded, selected

    return max(visit(0, -1)) if adj else 0


def far_numerator(adj: list[list[int]], rank: int) -> dict[str, int]:
    """Return the exact unnormalised qfar at one rank."""
    cap = rank + 2
    poly = tree_polynomial(adj, cap)
    total, excluded = cavity_messages(adj, rank)
    one_polys = [
        product([excluded[(w, v)] for w in adj[v]], rank)
        for v in range(len(adj))
    ]
    one = [coefficient(marked, rank) for marked in one_polys]

    joint_two = 0
    for center, neighbors in enumerate(adj):
        # Sum R_u R_v product_(w != u,v) P_w without quadratic pair loops.
        # Here R_u is the u-side polynomial after deleting N[u] at the path
        # endpoint, and P_u is the unrestricted u-side cavity message.
        zero, marked_once, marked_twice = [1], [0], [0]
        for u in neighbors:
            unrestricted = total[(u, center)]
            endpoint = product(
                [excluded[(w, u)] for w in adj[u] if w != center], rank
            )
            marked_twice = add(
                mul(marked_twice, unrestricted, rank),
                mul(marked_once, endpoint, rank),
                rank,
            )
            marked_once = add(
                mul(marked_once, unrestricted, rank),
                mul(zero, endpoint, rank),
                rank,
            )
            zero = mul(zero, unrestricted, rank)
        joint_two += coefficient(marked_twice, rank)

    # Every unordered vertex pair is either equal, adjacent, distance two,
    # or far.  In a tree, a distance-two pair has one unique center.
    all_pairs = (sum(one) ** 2 - sum(value * value for value in one)) // 2
    edge_pairs = sum(
        one[u] * one[v]
        for u, neighbors in enumerate(adj)
        for v in neighbors
        if u < v
    )
    two_pairs = 0
    for neighbors in adj:
        values = [one[v] for v in neighbors]
        two_pairs += (sum(values) ** 2 - sum(value * value for value in values)) // 2
    product_far = all_pairs - edge_pairs - two_pairs

    count = coefficient(poly, rank)
    joint_all = comb(rank + 2, 2) * coefficient(poly, rank + 2)
    joint_far = joint_all - joint_two
    return {
        "rank": rank,
        "N": count,
        "d1": sum(one),
        "joint_far": joint_far,
        "product_far": product_far,
        "qfar": count * joint_far - product_far,
    }


if __name__ == "__main__":
    # Regression: the known deep-tail T(3,4,1) obstruction.
    from targeted import make_T_m_t_1

    _, tree = make_T_m_t_1(3, 4)
    row = far_numerator(tree, 13)
    assert row["qfar"] == 196002
    print(row)
