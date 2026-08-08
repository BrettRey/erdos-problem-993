#!/usr/bin/env python3
"""Exact replay checks for the double-star log-concavity lemma.

The proof is algebraic; this script independently checks the graph formula,
the factorization used in the proof, and a finite parameter box using only
integer arithmetic.
"""

from __future__ import annotations

import argparse
from math import comb

from indpoly import independence_poly, is_log_concave, is_unimodal


def poly_add(left: list[int], right: list[int]) -> list[int]:
    out = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        out[index] += value
    for index, value in enumerate(right):
        out[index] += value
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_mul(left: list[int], right: list[int]) -> list[int]:
    out = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            out[i + j] += left_value * right_value
    return out


def binomial_poly(power: int) -> list[int]:
    return [comb(power, index) for index in range(power + 1)]


def shift(poly: list[int]) -> list[int]:
    return [0] + poly


def double_star_formula(p: int, q: int) -> list[int]:
    """Return I(D_{p,q};x), where the two hubs have p and q leaves."""
    return poly_add(
        binomial_poly(p + q),
        poly_add(shift(binomial_poly(p)), shift(binomial_poly(q))),
    )


def double_star_factorization(p: int, q: int) -> list[int]:
    """Return the factorization used in the proof, expanded exactly."""
    if p > q:
        p, q = q, p
    d = q - p
    star = poly_add(binomial_poly(p), [0, 1])
    g_poly = poly_mul(binomial_poly(d), star)
    h_poly = poly_add(g_poly, [0, 1])
    return poly_mul(binomial_poly(p), h_poly)


def double_star_adj(p: int, q: int) -> list[list[int]]:
    """Build the double star with adjacent hubs 0 and 1."""
    adj = [[] for _ in range(p + q + 2)]

    def add_edge(left: int, right: int) -> None:
        adj[left].append(right)
        adj[right].append(left)

    add_edge(0, 1)
    next_vertex = 2
    for _ in range(p):
        add_edge(0, next_vertex)
        next_vertex += 1
    for _ in range(q):
        add_edge(1, next_vertex)
        next_vertex += 1
    return adj


def exceptional_gap(q: int, d: int) -> int:
    """The only Turan gap made harder by adding the final x term."""
    return (comb(q, 2) + d) ** 2 - (q + 2) * (
        comb(q, 3) + comb(d, 2)
    )


def verify(max_leaves: int, dp_max_leaves: int) -> dict[str, int]:
    dp_cases = 0
    for p in range(dp_max_leaves + 1):
        for q in range(dp_max_leaves + 1):
            adj = double_star_adj(p, q)
            assert independence_poly(len(adj), adj) == double_star_formula(p, q)
            dp_cases += 1

    parameter_cases = 0
    minimum_gap: int | None = None
    for p in range(max_leaves + 1):
        for q in range(max_leaves + 1):
            formula = double_star_formula(p, q)
            assert formula == double_star_factorization(p, q)
            assert is_log_concave(formula)
            assert is_unimodal(formula)
            small, large = sorted((p, q))
            gap = exceptional_gap(large, large - small)
            minimum_gap = gap if minimum_gap is None else min(minimum_gap, gap)
            parameter_cases += 1

    # Exact checks of the two endpoint factorizations in the proof.
    for q in range(max_leaves + 1):
        assert 12 * exceptional_gap(q, 0) == (
            q * (q - 1) * (q * q - 3 * q + 8)
        )
        assert 12 * exceptional_gap(q, q) == (
            q * (q + 1) * (q * q + q + 4)
        )

    return {
        "dp_cases": dp_cases,
        "parameter_cases": parameter_cases,
        "minimum_exceptional_gap": minimum_gap or 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-leaves", type=int, default=100)
    parser.add_argument("--dp-max-leaves", type=int, default=8)
    args = parser.parse_args()
    result = verify(args.max_leaves, args.dp_max_leaves)
    print("double-star formula DP cases:", result["dp_cases"])
    print("exact parameter cases:", result["parameter_cases"])
    print("minimum exceptional Turan gap:", result["minimum_exceptional_gap"])
    print("all checked polynomials log-concave and unimodal: true")


if __name__ == "__main__":
    main()
