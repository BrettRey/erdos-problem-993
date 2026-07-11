#!/usr/bin/env python3
"""Exact certificate for prefix GSB at rank three.

The target is

    4 i_4(T)^2 + i_3(T)i_4(T) - 5 i_3(T)i_5(T) >= 0,

equivalently ``mu_4 <= mu_3 + 1``.  The script checks the universal
symbolic certificates for trees with at most three nonleaf vertices and for
all trees of order at least 15 with at least four nonleaf vertices.  It then
performs the finite exact base check for orders 10 through 14.

All symbolic comparisons are coefficientwise over polynomials in
nonnegative shifted variables.  No floating-point comparisons occur.
"""

from __future__ import annotations

from math import comb

import sympy as sp

from indpoly import independence_poly
from trees import trees_geng_raw


def bernstein_coefficients(
    expression: sp.Expr,
    variable: sp.Symbol,
    lower: sp.Expr,
    upper: sp.Expr,
    degree: int,
) -> list[sp.Expr]:
    y = sp.symbols("y")
    power = sp.Poly(
        sp.expand(expression.subs(variable, lower + y * (upper - lower))),
        y,
    )
    coefficients = [power.coeff_monomial(y**j) for j in range(degree + 1)]
    return [
        sp.factor(
            sum(
                coefficients[j]
                * sp.binomial(k, j)
                / sp.binomial(degree, j)
                for j in range(k + 1)
            )
        )
        for k in range(degree + 1)
    ]


def assert_shift_nonnegative(
    expression: sp.Expr,
    substitutions: dict[sp.Symbol, sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> None:
    polynomial = sp.Poly(sp.expand(expression.subs(substitutions)), *variables)
    assert all(coefficient >= 0 for _, coefficient in polynomial.terms())
    assert polynomial.as_expr() != 0


def verify_two_nonleaf_core() -> None:
    """Adjacent double stars, including every prefix case alpha>=7."""
    m, p, u, y = sp.symbols("m p u y")
    i3 = (m**3 - m - 6 * p) / 6
    i4 = (m - 2) * (m**3 - m - 12 * p) / 24
    i5 = (
        m**5
        - 5 * m**4
        + 5 * m**3
        - 20 * m**2 * p
        + 5 * m**2
        + 90 * m * p
        - 6 * m
        + 10 * p**2
        - 110 * p
    ) / 120
    deficit = sp.factor(4 * i4**2 + i3 * i4 - 5 * i3 * i5)
    # If a,b>=1 and a+b=m, then m-1<=ab<=m^2/4.  Prefix rank three
    # means alpha=m>=7.
    for coefficient in bernstein_coefficients(
        72 * deficit, p, m - 1, m**2 / 4, 3
    ):
        assert_shift_nonnegative(coefficient, {m: u + 7}, (u,))


def h3_deficit() -> tuple[sp.Expr, sp.Symbol, sp.Symbol, sp.Symbol]:
    """Return 72 times the deficit for a weighted three-vertex core path."""
    a, b, c, m, p = sp.symbols("a b c m p")
    t = a + b + c + 1
    leaf_count = a + b + c - 1
    s2 = (a + 1) ** 2 + (b + 1) ** 2 + (c + 1) ** 2 + leaf_count
    s3 = (a + 1) ** 3 + (b + 1) ** 3 + (c + 1) ** 3 + leaf_count
    s4 = (a + 1) ** 4 + (b + 1) ** 4 + (c + 1) ** 4 + leaf_count
    edge_product = (
        (a + 1) * (c + 1)
        + (b + 1) * (c + 1)
        + a * (a + 1)
        + b * (b + 1)
        + (c - 1) * (c + 1)
    )
    edge_product_sum = (
        (a + 1) * (c + 1) * (a + c + 2)
        + (b + 1) * (c + 1) * (b + c + 2)
        + a * (a + 1) * (a + 2)
        + b * (b + 1) * (b + 2)
        + (c - 1) * (c + 1) * (c + 2)
    )
    neighbor_square = (
        (a + c + 1) ** 2
        + (b + c + 1) ** 2
        + t**2
        + a * (a + 1) ** 2
        + b * (b + 1) ** 2
        + (c - 1) * (c + 1) ** 2
    )
    i3 = (3 * s2 + t**3 - 6 * t**2 - t) / 6
    i4 = (
        -24 * edge_product
        + 12 * t * s2
        - 4 * s3
        + t**4
        - 14 * t**3
        + 23 * t**2
        - 2 * t
    ) / 24
    i5 = (
        60 * neighbor_square
        - 120 * t * edge_product
        + 60 * edge_product_sum
        + 30 * t**2 * s2
        - 90 * t * s2
        + 55 * s2
        - 20 * t * s3
        - 30 * s3
        + 5 * s4
        + t**5
        - 25 * t**4
        + 125 * t**3
        - 115 * t**2
        - 6 * t
    ) / 120
    raw = sp.expand(72 * (4 * i4**2 + i3 * i4 - 5 * i3 * i5))
    symmetric, remainder, mapping = sp.symmetrize(raw, [a, b], formal=True)
    assert remainder == 0
    s1, s2_symbol = (entry[0] for entry in mapping)
    return sp.factor(symmetric.subs({s1: m, s2_symbol: p})), m, p, c


def verify_three_nonleaf_core() -> None:
    """Every three-nonleaf core path in the rank-three prefix."""
    deficit, m, p, c = h3_deficit()
    u, v = sp.symbols("u v")
    coefficients = bernstein_coefficients(deficit, p, m - 1, m**2 / 4, 3)
    # c=1 has matching number two, so prefix alpha>=7 gives m>=6.
    # c>=2 has matching number three and m+c>=8.  The five orthants
    # below cover those integer cases.
    orthants = ((6, 1), (5, 2), (4, 4), (3, 5), (2, 6))
    for coefficient in coefficients:
        for m0, c0 in orthants:
            assert_shift_nonnegative(
                coefficient, {m: u + m0, c: v + c0}, (u, v)
            )


def rank3_degree_formulas() -> tuple[sp.Expr, ...]:
    t, s2, s3, s4, edge_product, edge_product_sum, neighbor_square = sp.symbols(
        "t S2 S3 S4 P P21 L"
    )
    i3 = (3 * s2 + t**3 - 6 * t**2 - t) / 6
    i4 = (
        -24 * edge_product
        + 12 * t * s2
        - 4 * s3
        + t**4
        - 14 * t**3
        + 23 * t**2
        - 2 * t
    ) / 24
    i5 = (
        60 * neighbor_square
        - 120 * t * edge_product
        + 60 * edge_product_sum
        + 30 * t**2 * s2
        - 90 * t * s2
        + 55 * s2
        - 20 * t * s3
        - 30 * s3
        + 5 * s4
        + t**5
        - 25 * t**4
        + 125 * t**3
        - 115 * t**2
        - 6 * t
    ) / 120
    return (
        t,
        s2,
        s3,
        s4,
        edge_product,
        edge_product_sum,
        neighbor_square,
        i3,
        i4,
        i5,
    )


def verify_degree_formulas_against_dp() -> None:
    """Replay the three closed formulas against the tree DP through n=14."""
    for n in range(5, 15):
        for _, adj, _ in trees_geng_raw(n):
            t = n - 1
            degree = [len(neighbors) for neighbors in adj]
            s2 = sum(d**2 for d in degree)
            s3 = sum(d**3 for d in degree)
            s4 = sum(d**4 for d in degree)
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
            neighbor_square = sum(
                sum(degree[u] for u in neighbors) ** 2
                for neighbors in adj
            )
            numerators = (
                3 * s2 + t**3 - 6 * t**2 - t,
                -24 * edge_product
                + 12 * t * s2
                - 4 * s3
                + t**4
                - 14 * t**3
                + 23 * t**2
                - 2 * t,
                60 * neighbor_square
                - 120 * t * edge_product
                + 60 * edge_product_sum
                + 30 * t**2 * s2
                - 90 * t * s2
                + 55 * s2
                - 20 * t * s3
                - 30 * s3
                + 5 * s4
                + t**5
                - 25 * t**4
                + 125 * t**3
                - 115 * t**2
                - 6 * t,
            )
            denominators = (6, 24, 120)
            polynomial = independence_poly(n, adj)
            expected = tuple(
                polynomial[k] if k < len(polynomial) else 0
                for k in (3, 4, 5)
            )
            assert all(
                numerator % denominator == 0
                for numerator, denominator in zip(numerators, denominators)
            )
            got = tuple(
                numerator // denominator
                for numerator, denominator in zip(numerators, denominators)
            )
            assert got == expected, (n, got, expected, adj)


def verify_four_or_more_nonleaf_large() -> None:
    """Universal symbolic certificate for h>=4 and t=n-1>=14."""
    (
        t,
        s2,
        s3,
        s4,
        edge_product,
        edge_product_sum,
        neighbor_square,
        i3,
        i4,
        i5,
    ) = rank3_degree_formulas()
    q, r3, u4, core_product, joint = sp.symbols("Q R U C M")
    substitutions = {
        s2: q + 3 * t - 1,
        s3: r3 + 3 * q + 4 * t - 2,
        s4: u4 + 4 * r3 + 6 * q + 5 * t - 3,
        edge_product: core_product + q + 2 * t - 1,
    }
    # Combine L and P21 before substitution because they have equal
    # coefficients in i5.
    raw = sp.expand(144 * (4 * i4**2 + i3 * i4 - 5 * i3 * i5))
    raw = sp.expand(raw.subs(neighbor_square, joint - edge_product_sum))
    excess = sp.factor(raw.subs(substitutions))
    n3_factor = 3 * q + t**3 - 6 * t**2 + 8 * t - 3
    assert sp.factor(sp.diff(excess, joint) + 60 * n3_factor) == 0
    assert sp.factor(sp.diff(excess, u4) + 5 * n3_factor) == 0

    # For x_v=d_v-1 and s=t-1:
    #   L+P21 <= R+(t+7)Q+8C+10t-6,
    #   U <= (t-4)R when at least four x_v are positive.
    lower = sp.factor(
        excess.subs(
            {
                joint: r3 + (t + 7) * q + 8 * core_product + 10 * t - 6,
                u4: (t - 4) * r3,
            },
            simultaneous=True,
        )
    )

    # On the relaxed moment domain Q<=(t-1)^2, C>=t-2,
    # R>=Q^2/(t-1), this lower polynomial is increasing in C and R.
    d_c = sp.factor(sp.diff(lower, core_product))
    c_floor = (t - 1) * (3 * t**3 - 29 * t**2 + 46 * t - 8)
    assert sp.factor(
        d_c
        - 24
        * (
            c_floor
            + 48 * core_product
            + 8 * r3
            + 9 * (t - 1) * ((t - 1) ** 2 - q)
        )
    ) == 0

    d_r = sp.factor(sp.diff(lower, r3))
    r_floor = (t - 1) * (7 * t**3 - 56 * t**2 + 84 * t - 32)
    assert sp.factor(
        d_r
        - (
            r_floor
            + 192 * core_product
            + 32 * r3
            + (51 * t - 186) * ((t - 1) ** 2 - q)
        )
    ) == 0

    w = sp.symbols("w")
    assert_shift_nonnegative(c_floor, {t: w + 14}, (w,))
    assert_shift_nonnegative(r_floor, {t: w + 14}, (w,))

    reduced = sp.factor(
        lower.subs(
            {core_product: t - 2, r3: q**2 / (t - 1)},
            simultaneous=True,
        )
    )
    numerator = sp.factor(reduced * (t - 1) ** 2)
    assert not numerator.has(core_product, r3, u4, joint)

    # The quartic numerator is globally strongly convex in Q.  Its exact
    # curvature minimum and one Newton lower bound are positive for t>=14.
    curvature = sp.diff(numerator, q, 2)
    q_vertex = 3 * (t - 1) * (17 * t - 62) / 64
    curvature_floor = sp.factor(curvature.subs(q, q_vertex))
    curvature_quartic = 896 * t**4 - 2427 * t**3 + 13519 * t**2 - 78784 * t + 42220
    assert curvature_floor == (t - 1) * curvature_quartic / 64
    assert_shift_nonnegative(curvature_quartic, {t: w + 14}, (w,))

    q0 = (t - 1) * (9 * t - 50) / 14
    newton_floor = sp.factor(
        numerator.subs(q, q0)
        - sp.diff(numerator, q).subs(q, q0) ** 2 / (2 * curvature_floor)
    )
    newton_numerator, newton_denominator = sp.fraction(newton_floor)
    assert_shift_nonnegative(newton_numerator, {t: w + 14}, (w,))
    assert_shift_nonnegative(newton_denominator, {t: w + 14}, (w,))


def finite_base() -> list[dict[str, object]]:
    """Exact h>=4 base cases not covered by the t>=14 certificate."""
    expected = {
        10: (22, 5415),
        11: (138, 16318),
        12: (476, 42161),
        13: (1270, 94984),
        14: (3122, 197692),
    }
    rows: list[dict[str, object]] = []
    for n in range(8, 15):
        checked = 0
        minimum: tuple[int, str] | None = None
        for _, adj, raw in trees_geng_raw(n):
            if sum(len(neighbors) >= 2 for neighbors in adj) < 4:
                continue
            polynomial = independence_poly(n, adj)
            alpha = len(polynomial) - 1
            if alpha < 7:
                continue
            checked += 1
            deficit = (
                4 * polynomial[4] ** 2
                + polynomial[3] * polynomial[4]
                - 5 * polynomial[3] * polynomial[5]
            )
            if minimum is None or deficit < minimum[0]:
                minimum = (deficit, raw.decode())
        if checked:
            assert minimum is not None and minimum[0] > 0
            assert (checked, minimum[0]) == expected[n]
            rows.append(
                {
                    "n": n,
                    "checked": checked,
                    "minimum_deficit": minimum[0],
                    "witness_graph6": minimum[1],
                }
            )
    assert {row["n"] for row in rows} == set(expected)
    return rows


def main() -> None:
    verify_degree_formulas_against_dp()
    verify_two_nonleaf_core()
    verify_three_nonleaf_core()
    verify_four_or_more_nonleaf_large()
    rows = finite_base()
    print({"symbolic_certificates": "passed", "finite_base": rows})


if __name__ == "__main__":
    main()
