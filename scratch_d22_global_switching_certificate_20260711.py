#!/usr/bin/env python3
"""Exact certificates for the D22 global-switching lane.

The script records three falsifications that any Bencs-style covariance
proof has to survive:

* the zero-capacity shortcut ``2 q2 + qfar <= 0`` fails first at order 16;
* the stronger nonedge sign ``q2 + qfar <= 0`` fails in the prefix on the
  Galvin tree T(2,23,1);
* assigning the diagonal/edge collision capacity independently to each
  connected symmetric-difference component is false already at order 11.

It also audits the connected-component residual upper bound used in the
accompanying proof packet on every tree through order eight.
"""

from __future__ import annotations

from itertools import combinations
from math import ceil, comb

from graph6 import parse_graph6
from indpoly import independence_poly
from scratch_d18_covariance_search_20260711 import (
    cavity_messages,
    closed_neighborhood,
    coeff,
    covariance_profile,
    distance_classes,
    induced_poly,
    product,
)
from targeted import make_T_m_t_1
from trees import trees


Poly = list[int]


def poly_add(left: Poly, right: Poly, limit: int) -> Poly:
    return [
        (left[i] if i < len(left) else 0)
        + (right[i] if i < len(right) else 0)
        for i in range(min(limit + 1, max(len(left), len(right))))
    ]


def poly_mul(left: Poly, right: Poly, limit: int) -> Poly:
    out = [0] * min(limit + 1, len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right[: limit + 1 - i]):
            out[i + j] += x * y
    return out


def poly_pow(base: Poly, exponent: int, limit: int) -> Poly:
    out = [1]
    while exponent:
        if exponent & 1:
            out = poly_mul(out, base, limit)
        exponent //= 2
        if exponent:
            base = poly_mul(base, base, limit)
    return out


def galvin_orbit_data(m: int, t: int, rank: int) -> dict[str, int]:
    """Return exact marginal/covariance data using the four vertex orbits."""
    limit = rank + 2
    edge = [comb(t, k) * 2**k for k in range(min(t, limit) + 1)]
    leaves = [comb(t, k) for k in range(min(t, limit) + 1)]
    branch = poly_add(edge, [0] + leaves, limit)
    tree_poly = poly_add(
        poly_pow(branch, m, limit),
        [0] + poly_pow(edge, m, limit),
        limit,
    )

    edge_minus = [
        comb(t - 1, k) * 2**k
        for k in range(min(t - 1, limit) + 1)
    ]
    leaves_minus = [
        comb(t - 1, k) for k in range(min(t - 1, limit) + 1)
    ]
    branch_minus = poly_add(edge_minus, [0] + leaves_minus, limit)

    # A_v=I(T-N[v]) for root, level-one, level-two, and leaf vertices.
    a_root = comb(m * t, rank) * 2**rank
    a_w_poly = poly_mul(leaves, poly_pow(branch, m - 1, limit), limit)
    remaining_root = poly_add(
        poly_pow(branch, m - 1, limit),
        [0] + poly_pow(edge, m - 1, limit),
        limit,
    )
    a_x_poly = poly_mul(edge_minus, remaining_root, limit)
    a_y_poly = poly_add(
        poly_mul(branch_minus, poly_pow(branch, m - 1, limit), limit),
        [0]
        + poly_mul(edge_minus, poly_pow(edge, m - 1, limit), limit),
        limit,
    )
    a_w = coeff(a_w_poly, rank)
    a_x = coeff(a_x_poly, rank)
    a_y = coeff(a_y_poly, rank)

    n_rank = coeff(tree_poly, rank)
    n_rank_plus_two = coeff(tree_poly, rank + 2)
    d1 = a_root + m * a_w + m * t * (a_x + a_y)
    squares = a_root**2 + m * a_w**2 + m * t * (a_x**2 + a_y**2)
    edge_products = (
        m * a_root * a_w
        + m * t * a_w * a_x
        + m * t * a_x * a_y
    )
    nonedge_products = (d1**2 - squares) // 2 - edge_products
    q_nonedge = (
        n_rank * comb(rank + 2, 2) * n_rank_plus_two
        - nonedge_products
    )
    return {
        "N": n_rank,
        "i_rank_plus_two": n_rank_plus_two,
        "a_root": a_root,
        "a_w": a_w,
        "a_x": a_x,
        "a_y": a_y,
        "d1": d1,
        "squares": squares,
        "edge_products": edge_products,
        "q_nonedge": q_nonedge,
    }


def zero_capacity_obstruction() -> None:
    adj = [
        [13, 15], [13], [13], [13], [14, 15], [14], [14], [14],
        [14], [15], [15], [15], [15], [0, 1, 2, 3],
        [4, 5, 6, 7, 8], [0, 4, 9, 10, 11, 12],
    ]
    profile = covariance_profile(adj)
    row = profile["rows"][4]
    assert profile["n"] == 16 and profile["alpha"] == 13
    assert profile["limit"] == 9 and row["prefix"]
    assert row["N"] == 905 and row["d1"] == 7635
    assert row["q2"] == 316868 and row["qfar"] == -629342
    assert 2 * row["q2"] + row["qfar"] == 4394 > 0


def galvin_nonedge_obstruction() -> None:
    m, t, rank = 2, 23, 4
    _, adj = make_T_m_t_1(m, t)
    profile = covariance_profile(adj)
    row = profile["rows"][rank]
    orbit = galvin_orbit_data(m, t, rank)

    assert len(adj) == 95 and profile["alpha"] == 48
    assert profile["limit"] == 32 and row["prefix"]
    assert row["N"] == orbit["N"] == 2835141
    assert row["d1"] == orbit["d1"] == 240590350
    assert row["q2"] == 62005580307264
    assert row["qfar"] == -61454066298966
    assert row["q2"] + row["qfar"] == orbit["q_nonedge"]
    assert orbit["q_nonedge"] == 551514008298 > 0
    assert orbit["squares"] == 615288884289486
    assert orbit["edge_products"] == 404840482368080
    ordered_gap = (
        orbit["squares"]
        + 2 * orbit["edge_products"]
        - 2 * orbit["q_nonedge"]
    )
    assert ordered_gap == 1423866821009050 > 0
    assert ordered_gap + row["N"] * row["d1"] == 2105974386498400


def componentwise_collision_obstruction() -> None:
    n, adj = parse_graph6(b"J??????wC~?")
    assert n == 11 and independence_poly(n, adj)[-1] > 0
    assert len(independence_poly(n, adj)) - 1 == 9
    rank = 4
    assert rank == ceil((2 * 9 - 1) / 3) - 2

    h_vertices = {0, 1, 2, 9}
    removed = set(h_vertices)
    for vertex in h_vertices:
        removed.update(adj[vertex])
    remainder = induced_poly(adj, removed)
    assert remainder == [1, 6, 15, 20, 15, 6, 1]

    # H is K_(1,3), with endpoint color sizes p=3 and q=1.  If the
    # diagonal/edge target capacity is assigned component by component, its
    # residual is positive, so common-set/cross-component compensation is
    # indispensable.
    p, q = 3, 1
    local_residual = (
        comb(p, 2) * coeff(remainder, rank - q)
        * coeff(remainder, rank + 2 - p)
        + comb(q, 2) * coeff(remainder, rank - p)
        * coeff(remainder, rank + 2 - q)
        - p * q * coeff(remainder, rank + 1 - p)
        * coeff(remainder, rank + 1 - q)
    )
    assert local_residual == 525 > 0


def connected_subsets(adj: list[list[int]]):
    n = len(adj)
    for mask in range(1, 1 << n):
        vertices = [v for v in range(n) if mask >> v & 1]
        vertex_set = set(vertices)
        seen = {vertices[0]}
        queue = [vertices[0]]
        for vertex in queue:
            for neighbor in adj[vertex]:
                if neighbor in vertex_set and neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        if len(seen) == len(vertices):
            yield vertices, vertex_set


def bencs_residual_upper(
    adj: list[list[int]], rank: int, far_only: bool
) -> int:
    """Connected-component residual after the standard color switches.

    Pairs whose two marks lie in distinct symmetric-difference components
    cancel.  Pairs with both marks common to the two target sets contribute
    only negatively and are omitted here, hence this is an upper bound.
    """
    distances = distance_classes(adj)
    total = 0
    for vertices, vertex_set in connected_subsets(adj):
        color = {vertices[0]: 0}
        queue = [vertices[0]]
        for vertex in queue:
            for neighbor in adj[vertex]:
                if neighbor in vertex_set and neighbor not in color:
                    color[neighbor] = 1 - color[vertex]
                    queue.append(neighbor)
        p_side = [v for v in vertices if color[v] == 0]
        q_side = [v for v in vertices if color[v] == 1]

        removed = set(vertices)
        for vertex in vertices:
            removed.update(adj[vertex])
        remainder = induced_poly(adj, removed)

        def same_color_pairs(side: list[int]) -> int:
            return sum(
                not far_only or distances[u][v] >= 3
                for u, v in combinations(side, 2)
            )

        cross_pairs = sum(
            v not in adj[u] for u in p_side for v in q_side
        )
        p = len(p_side)
        q = len(q_side)
        total += (
            same_color_pairs(p_side)
            * coeff(remainder, rank - q)
            * coeff(remainder, rank + 2 - p)
            + same_color_pairs(q_side)
            * coeff(remainder, rank - p)
            * coeff(remainder, rank + 2 - q)
            - cross_pairs
            * coeff(remainder, rank + 1 - p)
            * coeff(remainder, rank + 1 - q)
        )
    return total


def exact_covariance_sum(
    adj: list[list[int]], rank: int, far_only: bool
) -> int:
    poly = independence_poly(len(adj), adj)
    n_rank = coeff(poly, rank)
    distances = distance_classes(adj)
    one = [
        induced_poly(adj, closed_neighborhood(adj, vertex))
        for vertex in range(len(adj))
    ]
    total = 0
    for u, v in combinations(range(len(adj)), 2):
        if distances[u][v] < (3 if far_only else 2):
            continue
        joint = induced_poly(
            adj,
            closed_neighborhood(adj, u) | closed_neighborhood(adj, v),
        )
        total += (
            n_rank * coeff(joint, rank)
            - coeff(one[u], rank) * coeff(one[v], rank)
        )
    return total


def audit_bencs_residual(max_n: int = 8) -> None:
    for n in range(1, max_n + 1):
        for _, adj in trees(n):
            poly = independence_poly(n, adj)
            for rank in range(len(poly)):
                for far_only in (False, True):
                    exact = exact_covariance_sum(adj, rank, far_only)
                    residual = bencs_residual_upper(adj, rank, far_only)
                    assert exact <= residual


def main() -> None:
    zero_capacity_obstruction()
    galvin_nonedge_obstruction()
    componentwise_collision_obstruction()
    audit_bencs_residual()
    print(
        {
            "zero_capacity_n16": "passed",
            "galvin_T_2_23_nonedge": "passed",
            "componentwise_n11": "passed",
            "bencs_residual_through_n8": "passed",
        },
        flush=True,
    )


if __name__ == "__main__":
    main()
