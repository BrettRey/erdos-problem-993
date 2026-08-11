#!/usr/bin/env python3
"""Exact replay for once-subdivided double-star log-concavity."""

from __future__ import annotations

import argparse
from math import comb

from indpoly import independence_poly, is_log_concave, is_unimodal
from verify_double_star_log_concavity_20260808 import (
    binomial_poly,
    poly_add,
    poly_mul,
    shift,
)


def once_subdivided_formula(p: int, q: int) -> list[int]:
    """Return the polynomial after subdividing the double-star hub edge once."""
    base = binomial_poly(p + q + 1)
    hub_terms = poly_add(shift(binomial_poly(p)), shift(binomial_poly(q)))
    return poly_add(poly_add(base, hub_terms), [0, 0, 1])


def once_subdivided_factorization(p: int, q: int) -> list[int]:
    """Expand the log-concave core plus x^2 used in the proof."""
    if p > q:
        p, q = q, p
    d = q - p
    bracket = poly_add(
        binomial_poly(q + 1),
        poly_add(shift(binomial_poly(d)), [0, 1]),
    )
    core = poly_mul(binomial_poly(p), bracket)
    return poly_add(core, [0, 0, 1])


def once_subdivided_adj(p: int, q: int) -> list[list[int]]:
    """Build the tree with hubs 0,1 and their unique connector vertex 2."""
    adj = [[] for _ in range(p + q + 3)]

    def add_edge(left: int, right: int) -> None:
        adj[left].append(right)
        adj[right].append(left)

    add_edge(0, 2)
    add_edge(2, 1)
    next_vertex = 3
    for _ in range(p):
        add_edge(0, next_vertex)
        next_vertex += 1
    for _ in range(q):
        add_edge(1, next_vertex)
        next_vertex += 1
    return adj


def index_three_gap(p: int, q: int) -> int:
    """Return f_3^2-f_2 f_4 for the once-subdivided family."""
    s = p + q
    f2 = comb(s + 2, 2)
    f3 = comb(s + 1, 3) + comb(p, 2) + comb(q, 2)
    f4 = comb(s + 1, 4) + comb(p, 3) + comb(q, 3)
    return f3 * f3 - f2 * f4


def symmetric_scaled_gap(p: int, q: int) -> int:
    """Return the symmetric polynomial proved equal to 144 times the gap."""
    s = p + q
    z = (p - q) ** 2
    base = s * (s**5 + 6 * s**4 + s**3 - 9 * s**2 + 16 * s - 60)
    z_coefficient = 3 * (s * (s - 1) * (s + 4) + 12)
    return base + z_coefficient * z + 9 * z * z


def verify(max_leaves: int, dp_max_leaves: int) -> dict[str, int]:
    dp_cases = 0
    for p in range(dp_max_leaves + 1):
        for q in range(dp_max_leaves + 1):
            adj = once_subdivided_adj(p, q)
            assert independence_poly(len(adj), adj) == once_subdivided_formula(p, q)
            dp_cases += 1

    parameter_cases = 0
    minimum_gap: int | None = None
    for p in range(max_leaves + 1):
        for q in range(max_leaves + 1):
            formula = once_subdivided_formula(p, q)
            assert formula == once_subdivided_factorization(p, q)
            assert is_log_concave(formula)
            assert is_unimodal(formula)
            gap = index_three_gap(p, q)
            assert 144 * gap == symmetric_scaled_gap(p, q)
            minimum_gap = gap if minimum_gap is None else min(minimum_gap, gap)
            parameter_cases += 1

    return {
        "dp_cases": dp_cases,
        "parameter_cases": parameter_cases,
        "minimum_index_three_gap": minimum_gap or 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-leaves", type=int, default=100)
    parser.add_argument("--dp-max-leaves", type=int, default=8)
    args = parser.parse_args()
    result = verify(args.max_leaves, args.dp_max_leaves)
    print("once-subdivided formula DP cases:", result["dp_cases"])
    print("exact parameter cases:", result["parameter_cases"])
    print("minimum index-three Turan gap:", result["minimum_index_three_gap"])
    print("all checked polynomials log-concave and unimodal: true")


if __name__ == "__main__":
    main()
