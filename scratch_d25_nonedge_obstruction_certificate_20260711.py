#!/usr/bin/env python3
"""Compact exact certificate refuting the D25 q_nonedge <= 0 lemma.

The witness is Galvin's tree T_(40,20,1), at rank r=20.  Polynomial
arithmetic is truncated at degree r+2 and uses arbitrary-precision integers.
Representative closed-neighborhood deletions are independently replayed by a
generic forest DP on the explicit 1641-vertex adjacency list.
"""

from __future__ import annotations

from math import ceil, comb

from targeted import make_T_m_t_1


M = 40
T = 20
R = 20
CAP = R + 2
EXPECTED_Q = int(
    "478465879328267257053095625419441855903392361962734248496390523870975697502462729069752063100"
)


def add(left: list[int], right: list[int]) -> list[int]:
    out = [0] * min(CAP + 1, max(len(left), len(right)))
    for i, value in enumerate(left[: len(out)]):
        out[i] += value
    for i, value in enumerate(right[: len(out)]):
        out[i] += value
    return out


def mul(left: list[int], right: list[int]) -> list[int]:
    out = [0] * min(CAP + 1, len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right[: len(out) - i]):
            out[i + j] += a * b
    return out


def power(base: list[int], exponent: int) -> list[int]:
    out = [1]
    while exponent:
        if exponent & 1:
            out = mul(out, base)
        exponent //= 2
        if exponent:
            base = mul(base, base)
    return out


def shift(poly: list[int]) -> list[int]:
    return [0] + poly[:CAP]


def coeff(poly: list[int], rank: int) -> int:
    return poly[rank] if rank < len(poly) else 0


def induced_forest_poly(
    adj: list[list[int]], removed: set[int], cap: int
) -> list[int]:
    """Generic truncated independence-polynomial DP for an induced forest."""
    global CAP
    assert cap == CAP
    n = len(adj)
    seen = [False] * n
    result = [1]
    for root in range(n):
        if root in removed or seen[root]:
            continue
        seen[root] = True
        parent = {root: -1}
        order = [root]
        for u in order:
            for v in adj[u]:
                if v not in removed and not seen[v]:
                    seen[v] = True
                    parent[v] = u
                    order.append(v)
        excluded: dict[int, list[int]] = {}
        selected: dict[int, list[int]] = {}
        for u in reversed(order):
            e, s = [1], [0, 1]
            for v in adj[u]:
                if parent.get(v) != u:
                    continue
                e = mul(e, add(excluded[v], selected[v]))
                s = mul(s, excluded[v])
            excluded[u], selected[u] = e, s
        result = mul(result, add(excluded[root], selected[root]))
    return result


def main() -> None:
    n, adj = make_T_m_t_1(M, T)
    alpha = M * (T + 1)
    limit = ceil((2 * alpha - 1) / 3)
    assert n == 1641 and alpha == 840 and limit == 560
    assert R >= 4 and R <= limit - 2
    assert sum(map(len, adj)) == 2 * (n - 1)

    p2 = [1, 2]
    leaf = [1, 1]
    e_t = power(p2, T)
    s_t = shift(power(leaf, T))
    branch = add(e_t, s_t)
    e_tm = power(p2, M * T)
    whole = add(power(branch, M), shift(e_tm))

    # Closed-neighborhood deletion polynomials for the four vertex orbits:
    # root, level-one w, level-two x, and leaf y.
    a_root_poly = e_tm
    a_w_poly = mul(power(leaf, T), power(branch, M - 1))
    root_other = add(
        power(branch, M - 1),
        shift(power(p2, (M - 1) * T)),
    )
    a_x_poly = mul(power(p2, T - 1), root_other)
    branch_less = add(power(p2, T - 1), shift(power(leaf, T - 1)))
    a_y_poly = add(
        mul(branch_less, power(branch, M - 1)),
        shift(power(p2, M * T - 1)),
    )

    # Independent generic replay on the explicit adjacency list.
    root, w, x, y = 0, 1, 1 + M, 1 + M + M * T
    generic_whole = induced_forest_poly(adj, set(), CAP)
    assert generic_whole == whole
    for vertex, formula in (
        (root, a_root_poly),
        (w, a_w_poly),
        (x, a_x_poly),
        (y, a_y_poly),
    ):
        removed = {vertex, *adj[vertex]}
        assert induced_forest_poly(adj, removed, CAP) == formula

    n_r = coeff(whole, R)
    i_next = coeff(whole, R + 1)
    i_two = coeff(whole, R + 2)
    a0, aw, ax, ay = (
        coeff(a_root_poly, R),
        coeff(a_w_poly, R),
        coeff(a_x_poly, R),
        coeff(a_y_poly, R),
    )
    d1 = a0 + M * aw + M * T * ax + M * T * ay
    assert d1 == (R + 1) * i_next

    sum_squares = a0 * a0 + M * aw * aw + M * T * ax * ax + M * T * ay * ay
    edge_products = M * a0 * aw + M * T * aw * ax + M * T * ax * ay
    all_pair_products = (d1 * d1 - sum_squares) // 2
    nonedge_products = all_pair_products - edge_products
    joint_nonedge = comb(R + 2, 2) * i_two
    q_nonedge = n_r * joint_nonedge - nonedge_products
    assert q_nonedge == EXPECTED_Q and q_nonedge > 0

    # The pivoted, full local budgets remain comfortably positive here.
    ordered_lc_gap = sum_squares + 2 * edge_products - 2 * q_nonedge
    gsb_gap = ordered_lc_gap + n_r * d1
    assert ordered_lc_gap > 0 and gsb_gap > 0
    assert ordered_lc_gap == d1 * d1 - 2 * n_r * joint_nonedge
    assert gsb_gap == d1 * d1 + n_r * d1 - 2 * n_r * joint_nonedge

    print(
        {
            "family": "T_(40,20,1)",
            "n": n,
            "alpha": alpha,
            "prefix_limit": limit,
            "rank": R,
            "N": n_r,
            "d1": d1,
            "a_orbits": {"root": a0, "w": aw, "x": ax, "y": ay},
            "q_nonedge": q_nonedge,
            "sum_squares": sum_squares,
            "edge_products": edge_products,
            "ordered_lc_gap": ordered_lc_gap,
            "gsb_gap": gsb_gap,
            "certificate": "passed",
        },
        flush=True,
    )


if __name__ == "__main__":
    main()
