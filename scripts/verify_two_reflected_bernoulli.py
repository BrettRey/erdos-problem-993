#!/usr/bin/env python3
"""Exact audit harness for the two-reflected-Bernoulli theorem.

For a Poisson-binomial sum X with parameters in (0, 1/2] and two
independent Bernoulli variables Y_1, Y_2 with parameters in [0, 1/2], the
proof in

    notes/literature/two_reflected_bernoulli_effective_drop_2026-07-10.md

gives an effective-drop lower bound at the first strict descent of
X-Y_1-Y_2.  This script is a deterministic falsification artifact, not part
of the proof.  It uses exact rational arithmetic for the finite grid and
SymPy exact polynomials for the derivative and Bernstein certificates.
"""

from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp


CORE_PROBABILITIES = (
    Fraction(1, 10),
    Fraction(1, 6),
    Fraction(1, 4),
    Fraction(1, 3),
    Fraction(2, 5),
    Fraction(1, 2),
)
HEAVY_PROBABILITIES = (
    Fraction(1, 4),
    Fraction(1, 3),
    Fraction(2, 5),
    Fraction(1, 2),
)
DUST_PROBABILITIES = (
    Fraction(1, 100),
    Fraction(1, 50),
    Fraction(1, 20),
    Fraction(1, 10),
)
Q_VALUES = (
    Fraction(0),
    Fraction(1, 1000),
    Fraction(1, 100),
    Fraction(1, 10),
    Fraction(1, 5),
    Fraction(1, 3),
    Fraction(2, 5),
    Fraction(1, 2),
)


MATRIX_2 = (
    (Fraction(0), Fraction(1, 2), Fraction(3, 2), Fraction(3), Fraction(4)),
    (Fraction(1, 2), Fraction(15, 8), Fraction(55, 12), Fraction(9), Fraction(14)),
    (Fraction(3, 2), Fraction(55, 12), Fraction(191, 18), Fraction(62, 3), Fraction(100, 3)),
    (Fraction(3), Fraction(9), Fraction(62, 3), Fraction(40), Fraction(64)),
    (Fraction(4), Fraction(14), Fraction(100, 3), Fraction(64), Fraction(96)),
)
MATRIX_1 = (
    (Fraction(0), Fraction(2), Fraction(5), Fraction(9), Fraction(12)),
    (Fraction(2), Fraction(15, 2), Fraction(16), Fraction(29), Fraction(46)),
    (Fraction(5), Fraction(16), Fraction(581, 18), Fraction(57), Fraction(272, 3)),
    (Fraction(9), Fraction(29), Fraction(57), Fraction(98), Fraction(152)),
    (Fraction(12), Fraction(46), Fraction(272, 3), Fraction(152), Fraction(224)),
)
MATRIX_0 = (
    (Fraction(0), Fraction(2), Fraction(4), Fraction(6), Fraction(8)),
    (Fraction(2), Fraction(6), Fraction(31, 3), Fraction(15), Fraction(20)),
    (Fraction(4), Fraction(31, 3), Fraction(52, 3), Fraction(25), Fraction(100, 3)),
    (Fraction(6), Fraction(15), Fraction(25), Fraction(36), Fraction(48)),
    (Fraction(8), Fraction(20), Fraction(100, 3), Fraction(48), Fraction(64)),
)


def convolve_bernoulli(pmf: Sequence[Fraction], p: Fraction) -> list[Fraction]:
    if not 0 <= p <= 1:
        raise AssertionError(f"Bernoulli parameter outside [0, 1]: {p}")
    out = [Fraction(0) for _ in range(len(pmf) + 1)]
    for k, value in enumerate(pmf):
        out[k] += value * (1 - p)
        out[k + 1] += value * p
    return out


def poisson_binomial_pmf(parameters: Sequence[Fraction]) -> list[Fraction]:
    pmf = [Fraction(1)]
    for p in parameters:
        pmf = convolve_bernoulli(pmf, p)
    return pmf


def first_strict_descent(pmf: Sequence[Fraction], support_start: int = 0) -> int:
    for index in range(1, len(pmf)):
        if pmf[index] < pmf[index - 1]:
            return support_start + index
    raise AssertionError("finite positive pmf has no strict descent")


def at(seq: Sequence[Fraction], index: int) -> Fraction:
    return seq[index] if 0 <= index < len(seq) else Fraction(0)


def a_at(a: Sequence[Fraction], k: int) -> Fraction:
    return at(a, k)


def c_at(c: Sequence[Fraction], k: int) -> Fraction:
    """Read a signed pmf stored on support -2, -1, ..., m."""

    return at(c, k + 2)


def signed_pmf_two(
    a: Sequence[Fraction], q1: Fraction, q2: Fraction
) -> list[Fraction]:
    """Return the pmf of X-Y_1-Y_2 on support -2, ..., m."""

    if not 0 <= q1 <= Fraction(1, 2) or not 0 <= q2 <= Fraction(1, 2):
        raise AssertionError(f"q outside [0, 1/2]: {(q1, q2)}")
    out = [Fraction(0) for _ in range(len(a) + 2)]
    w0 = (1 - q1) * (1 - q2)
    w1 = q1 * (1 - q2) + (1 - q1) * q2
    w2 = q1 * q2
    for j, value in enumerate(a):
        out[j + 2] += w0 * value
        out[j + 1] += w1 * value
        out[j] += w2 * value
    return out


def tilde_at(a: Sequence[Fraction], e: Fraction, f: Fraction, k: int) -> Fraction:
    return a_at(a, k) + e * a_at(a, k + 1) + f * a_at(a, k + 2)


def effective_drop(c: Sequence[Fraction], descent: int) -> Fraction:
    center = c_at(c, descent)
    if center <= 0:
        raise AssertionError("effective drop has zero center coefficient")
    return 1 - c_at(c, descent - 1) * c_at(c, descent + 1) / center**2


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "exact": fraction_text(value),
        "decimal": f"{float(value):.17g}",
    }


def parameter_blocks(parameters: Sequence[Fraction]) -> list[dict[str, Any]]:
    counts = Counter(parameters)
    return [
        {"count": counts[p], "p": fraction_text(p)}
        for p in sorted(counts)
    ]


def sympy_rational(value: Fraction) -> sp.Rational:
    return sp.Rational(value.numerator, value.denominator)


def bernstein_polynomial(
    matrix: Sequence[Sequence[Fraction]], u: sp.Symbol, v: sp.Symbol
) -> sp.Expr:
    return sp.expand(
        sum(
            sympy_rational(matrix[i][j])
            * sp.binomial(4, i)
            * u**i
            * (1 - u) ** (4 - i)
            * sp.binomial(4, j)
            * v**j
            * (1 - v) ** (4 - j)
            for i in range(5)
            for j in range(5)
        )
    )


def symbolic_certificate() -> dict[str, Any]:
    """Reconstruct every polynomial identity displayed in the proof note."""

    k, u, v, x, y = sp.symbols("K u v x y")
    e = u + v
    f = u * v

    am, a0, a1, a2, a3 = sp.symbols("a_m a_0 a_1 a_2 a_3")
    center = a0 + e * a1 + f * a2
    previous = am + e * a0 + f * a1
    following = a1 + e * a2 + f * a3
    n0 = a0**2 - am * a1
    n1 = a1**2 - a0 * a2
    n2 = a2**2 - a1 * a3
    c01 = a0 * a1 - am * a2
    c12 = a1 * a2 - a0 * a3
    long_gap = a1**2 - am * a3
    six_term = (
        n0
        + (u**2 + v**2) * n1
        + f * long_gap
        + f**2 * n2
        + e * c01
        + e * f * c12
    )
    if sp.expand(center**2 - previous * following - six_term) != 0:
        raise AssertionError("symbolic six-term curvature identity failed")

    denominator = (1 + e * x + f * x * y) ** 2
    numerator = (
        1 / k
        + 2 * e * x / (k + 1)
        + (u**2 + v**2) * x**2 / (k + 1)
        + 2 * f * (2 * k + 1) * x**2 / ((k + 1) * (k + 2))
        + 2 * f * e * x**2 * y / (k + 2)
        + f**2 * x**2 * y**2 / (k + 2)
    )
    g = numerator / denominator

    h = (
        2 * k**2 * u * v * x**2
        - k**2 * u * v * x * y
        + k**2 * u * x
        + k**2 * v * x
        + k**2
        + k * u**2 * x**2
        - k * u * v * x * y
        + 3 * k * u * x
        + k * v**2 * x**2
        + 3 * k * v * x
        + 3 * k
        + 2
    )
    dy_claim = (
        -2
        * u
        * v
        * x
        * h
        / (k * (k + 1) * (k + 2) * (1 + e * x + f * x * y) ** 3)
    )
    if sp.cancel(sp.diff(g, y) - dy_claim) != 0:
        raise AssertionError("symbolic y-derivative identity failed")

    diagonal = g.subs(y, x)
    j = (
        k * u * v * (u - v) ** 2 * x**3
        + 3 * k * u * v * (u + v) * x**2
        + (8 * k + 4) * u * v * x
        + (k + 2) * (u + v)
    )
    diagonal_claim = (
        -2
        * j
        / (k * (k + 1) * (k + 2) * (1 + u * x) ** 3 * (1 + v * x) ** 3)
    )
    if sp.cancel(sp.diff(diagonal, x) - diagonal_claim) != 0:
        raise AssertionError("symbolic diagonal derivative identity failed")

    paid = (
        k + 4 * u / (1 + u) ** 2 + 4 * v / (1 + v) ** 2
    ) * g.subs({x: 1, y: 1}) - 1
    claimed_denominator = (
        k * (k + 1) * (k + 2) * (1 + u) ** 4 * (1 + v) ** 4
    )
    paid_numerator = sp.cancel(paid * claimed_denominator)
    paid_polynomial = sp.Poly(paid_numerator, k)
    if paid_polynomial.degree() != 2:
        raise AssertionError(
            f"variance-payment numerator has degree {paid_polynomial.degree()}"
        )

    matrices = {2: MATRIX_2, 1: MATRIX_1, 0: MATRIX_0}
    for power, matrix in matrices.items():
        coefficient = paid_polynomial.coeff_monomial(k**power)
        reconstructed = bernstein_polynomial(matrix, u, v)
        if sp.expand(coefficient - reconstructed) != 0:
            raise AssertionError(
                f"degree-(4,4) Bernstein reconstruction failed for C_{power}"
            )
        if any(entry < 0 for row in matrix for entry in row):
            raise AssertionError(f"negative Bernstein coefficient in C_{power}")

    all_entries = [
        entry
        for matrix in matrices.values()
        for row in matrix
        for entry in row
    ]
    return {
        "six_term_curvature_identity": "passed",
        "y_derivative_identity": "passed",
        "diagonal_derivative_identity": "passed",
        "variance_payment_numerator_degree": paid_polynomial.degree(),
        "bernstein_reconstructions": 3,
        "bernstein_coefficients_checked": len(all_entries),
        "bernstein_minimum": fraction_text(min(all_entries)),
        "bernstein_maximum": fraction_text(max(all_entries)),
        "polynomial_engine": f"sympy {sp.__version__}",
    }


def curvature_lower_bound(
    k_capital: int,
    u: Fraction,
    v: Fraction,
    x: Fraction,
    y: Fraction,
) -> Fraction:
    e = u + v
    f = u * v
    numerator = (
        Fraction(1, k_capital)
        + 2 * e * x / (k_capital + 1)
        + (u**2 + v**2) * x**2 / (k_capital + 1)
        + 2
        * f
        * (2 * k_capital + 1)
        * x**2
        / ((k_capital + 1) * (k_capital + 2))
        + 2 * f * e * x**2 * y / (k_capital + 2)
        + f**2 * x**2 * y**2 / (k_capital + 2)
    )
    return numerator / (1 + e * x + f * x * y) ** 2


def verify_curvature_components(
    a: Sequence[Fraction],
    c: Sequence[Fraction],
    q1: Fraction,
    q2: Fraction,
    k: int,
) -> dict[str, Fraction]:
    """Check the six-term identity and all six component lower bounds."""

    k_capital = k + 1
    if k_capital < 1:
        raise AssertionError(f"curvature index must be nonnegative, got {k}")
    u = q1 / (1 - q1)
    v = q2 / (1 - q2)
    e = u + v
    f = u * v

    am = a_at(a, k - 1)
    a0 = a_at(a, k)
    a1 = a_at(a, k + 1)
    a2 = a_at(a, k + 2)
    a3 = a_at(a, k + 3)
    n0 = a0**2 - am * a1
    n1 = a1**2 - a0 * a2
    n2 = a2**2 - a1 * a3
    c01 = a0 * a1 - am * a2
    c12 = a1 * a2 - a0 * a3
    long_gap = a1**2 - am * a3

    lower_bounds = (
        ("N0", n0, a0**2 / k_capital),
        ("N1", n1, a1**2 / (k_capital + 1)),
        ("N2", n2, a2**2 / (k_capital + 2)),
        ("C01", c01, 2 * a0 * a1 / (k_capital + 1)),
        ("C12", c12, 2 * a1 * a2 / (k_capital + 2)),
        (
            "L",
            long_gap,
            2
            * (2 * k_capital + 1)
            * a1**2
            / ((k_capital + 1) * (k_capital + 2)),
        ),
    )
    for label, actual, lower in lower_bounds:
        if actual < lower:
            raise AssertionError(
                f"{label} lower bound failed at k={k}: {actual} < {lower}"
            )

    center = tilde_at(a, e, f, k)
    previous = tilde_at(a, e, f, k - 1)
    following = tilde_at(a, e, f, k + 1)
    lhs = center**2 - previous * following
    rhs = (
        n0
        + (u**2 + v**2) * n1
        + f * long_gap
        + f**2 * n2
        + e * c01
        + e * f * c12
    )
    if lhs != rhs:
        raise AssertionError(
            f"six-term curvature identity failed at k={k}: {lhs} != {rhs}"
        )
    if center <= 0 or a0 <= 0 or a1 <= 0:
        raise AssertionError(f"nonpositive normalization term at k={k}")

    normalization = (1 - q1) * (1 - q2)
    for index in (k - 1, k, k + 1):
        if c_at(c, index) != normalization * tilde_at(a, e, f, index):
            raise AssertionError(f"normalized/actual coefficient mismatch at k={index}")

    delta = lhs / center**2
    if delta != effective_drop(c, k):
        raise AssertionError(f"tilde/actual effective-drop mismatch at k={k}")
    x = a1 / a0
    y = a2 / a1
    lower = curvature_lower_bound(k_capital, u, v, x, y)

    weighted_component_lower = (
        a0**2 / k_capital
        + (u**2 + v**2) * a1**2 / (k_capital + 1)
        + f
        * 2
        * (2 * k_capital + 1)
        * a1**2
        / ((k_capital + 1) * (k_capital + 2))
        + f**2 * a2**2 / (k_capital + 2)
        + e * 2 * a0 * a1 / (k_capital + 1)
        + e * f * 2 * a1 * a2 / (k_capital + 2)
    )
    if lower != weighted_component_lower / center**2:
        raise AssertionError(f"normalized G reconstruction failed at k={k}")
    if delta < lower:
        raise AssertionError(f"curvature lower bound failed at k={k}")

    return {
        "delta": delta,
        "lower": lower,
        "x": x,
        "y": y,
        "u": u,
        "v": v,
    }


def verify_shifted_newton(
    parameters: Sequence[Fraction],
    c: Sequence[Fraction],
    q1: Fraction,
    q2: Fraction,
    s: int,
) -> Fraction:
    shifted = poisson_binomial_pmf((*parameters, 1 - q1, 1 - q2))
    if first_strict_descent(shifted) != s:
        raise AssertionError("left-branch shifted law has the wrong descent")
    for j, value in enumerate(shifted):
        if value != c_at(c, j - 2):
            raise AssertionError(f"shift identity failed at shifted index {j}")
    delta = 1 - shifted[s - 1] * shifted[s + 1] / shifted[s] ** 2
    if delta != effective_drop(c, s - 2):
        raise AssertionError("shifted/direct effective-drop mismatch")
    lower = Fraction(1, s + 1)
    if delta < lower:
        raise AssertionError("left-branch shifted Newton bound failed")
    return lower


def verify_row(
    parameters: Sequence[Fraction],
    a: Sequence[Fraction],
    variance_x: Fraction,
    s: int,
    q1: Fraction,
    q2: Fraction,
) -> dict[str, Any]:
    c = signed_pmf_two(a, q1, q2)
    d = first_strict_descent(c, support_start=-2)
    t1 = q1 / (1 - q1)
    t2 = q2 / (1 - q2)
    e = t1 + t2
    f = t1 * t2
    aa = lambda index: a_at(a, index)
    cap_a = aa(s - 2) - aa(s - 3)
    cap_b = aa(s - 1) - aa(s - 2)
    cap_c = aa(s - 1) - aa(s)
    cap_e = aa(s) - aa(s + 1)
    delta0 = cap_a + e * cap_b - f * cap_c
    delta1 = cap_b - e * cap_c - f * cap_e
    expected = s - 2 if delta0 < 0 else s - 1 if delta1 < 0 else s
    if d != expected or d not in (s - 2, s - 1, s):
        raise AssertionError(
            f"selection failed: S={s}, D={d}, expected={expected}, "
            f"q={(q1, q2)}, delta={(delta0, delta1)}, params={parameters}"
        )

    normalization = (1 - q1) * (1 - q2)
    diff0 = c_at(c, s - 2) - c_at(c, s - 3)
    diff1 = c_at(c, s - 1) - c_at(c, s - 2)
    if diff0 != normalization * delta0 or diff1 != normalization * delta1:
        raise AssertionError("shoulder-difference identity failed")
    if delta0 == 0 and c_at(c, s - 2) != c_at(c, s - 3):
        raise AssertionError("delta0 equality did not create its plateau")
    if delta1 == 0 and c_at(c, s - 1) != c_at(c, s - 2):
        raise AssertionError("delta1 equality did not create its plateau")

    curvature = verify_curvature_components(a, c, q1, q2, d)
    variance_y = q1 * (1 - q1) + q2 * (1 - q2)
    variance_total = variance_x + variance_y

    if d == s - 2:
        proof_lower = verify_shifted_newton(parameters, c, q1, q2, s)
        payment = Fraction(s + 1)
    else:
        if not 0 <= curvature["y"] <= curvature["x"] < 1:
            raise AssertionError(
                f"descent-ratio region failed: x={curvature['x']}, "
                f"y={curvature['y']}, D={d}, S={s}"
            )
        payment = Fraction(d + 1) + 4 * variance_y
        if payment * curvature["lower"] < 1:
            raise AssertionError(
                f"variance-payment inequality failed: payment={payment}, "
                f"G={curvature['lower']}"
            )
        proof_lower = curvature["lower"]

    if 4 * variance_total < payment:
        raise AssertionError(
            f"localization did not pay denominator: 4V={4 * variance_total}, "
            f"payment={payment}"
        )
    proof_scaled = variance_total * proof_lower
    scaled = variance_total * curvature["delta"]
    if proof_scaled < Fraction(1, 4):
        raise AssertionError(
            f"proof-chain quarter bound failed: {proof_scaled}, "
            f"S={s}, D={d}, q={(q1, q2)}"
        )
    if scaled < Fraction(1, 4):
        raise AssertionError(
            f"actual quarter bound failed: {scaled}, "
            f"S={s}, D={d}, q={(q1, q2)}"
        )

    return {
        "descent": d,
        "selection": f"S{d - s:+d}" if d != s else "S",
        "delta0_equal": delta0 == 0,
        "delta1_equal": delta1 == 0,
        "delta": curvature["delta"],
        "curvature_lower": curvature["lower"],
        "proof_lower": proof_lower,
        "variance_total": variance_total,
        "proof_scaled": proof_scaled,
        "scaled": scaled,
        "payment": payment,
    }


def core_profiles(min_m: int, max_m: int) -> Iterable[tuple[Fraction, ...]]:
    for m in range(min_m, max_m + 1):
        yield from itertools.combinations_with_replacement(CORE_PROBABILITIES, m)


def grouped_profiles() -> Iterable[tuple[Fraction, ...]]:
    for heavy_p, dust_p, heavy_count, dust_count in itertools.product(
        HEAVY_PROBABILITIES,
        DUST_PROBABILITIES,
        range(2, 9),
        (1, 2, 4, 8),
    ):
        yield tuple(sorted((heavy_p,) * heavy_count + (dust_p,) * dust_count))


def anchored_checks() -> dict[str, Any]:
    """Exercise both plateau conventions with exact, human-sized examples."""

    threshold_parameters = (Fraction(11, 31),) * 5
    threshold_a = poisson_binomial_pmf(threshold_parameters)
    threshold_s = first_strict_descent(threshold_a)
    threshold_variance = sum(
        (p * (1 - p) for p in threshold_parameters), Fraction(0)
    )
    threshold_row = verify_row(
        threshold_parameters,
        threshold_a,
        threshold_variance,
        threshold_s,
        Fraction(20, 119),
        Fraction(0),
    )
    if not threshold_row["delta1_equal"] or threshold_row["descent"] != threshold_s:
        raise AssertionError("delta1 plateau anchor failed")

    fair_parameters = (Fraction(1, 2),) * 5
    fair_a = poisson_binomial_pmf(fair_parameters)
    fair_s = first_strict_descent(fair_a)
    fair_variance = Fraction(5, 4)
    fair_row = verify_row(
        fair_parameters,
        fair_a,
        fair_variance,
        fair_s,
        Fraction(1, 2),
        Fraction(1, 2),
    )
    if not fair_row["delta0_equal"] or fair_row["descent"] != fair_s - 1:
        raise AssertionError("delta0 plateau anchor failed")

    return {
        "delta1_plateau": {
            "parameters": [{"count": 5, "p": "11/31"}],
            "q": ["20/119", "0"],
            "S": threshold_s,
            "D": threshold_row["descent"],
            "scaled": fraction_record(threshold_row["scaled"]),
        },
        "delta0_plateau": {
            "parameters": [{"count": 5, "p": "1/2"}],
            "q": ["1/2", "1/2"],
            "S": fair_s,
            "D": fair_row["descent"],
            "scaled": fraction_record(fair_row["scaled"]),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-m", type=int, default=4)
    parser.add_argument("--max-m", type=int, default=8)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/two_reflected_bernoulli_verify_2026-07-10.json"
        ),
    )
    args = parser.parse_args()
    if args.min_m < 1 or args.max_m < args.min_m:
        parser.error("require 1 <= min-m <= max-m")

    symbolic = symbolic_certificate()
    anchors = anchored_checks()
    seen: set[tuple[Fraction, ...]] = set()
    generated = Counter()
    eligible = Counter()
    selections = Counter()
    rows = 0
    delta0_equal_rows = 0
    delta1_equal_rows = 0
    payment_checks = 0
    shifted_newton_checks = 0
    minimum: tuple[Fraction, dict[str, Any]] | None = None

    for family, profiles in (
        ("core_multisets", core_profiles(args.min_m, args.max_m)),
        ("heavy_dust_blocks", grouped_profiles()),
    ):
        for parameters in profiles:
            generated[family] += 1
            if parameters in seen:
                continue
            seen.add(parameters)
            a = poisson_binomial_pmf(parameters)
            variance_x = sum(
                (p * (1 - p) for p in parameters), Fraction(0)
            )
            if variance_x < 1:
                continue
            eligible[family] += 1
            s = first_strict_descent(a)
            if s < 2 or s > len(parameters) - 1:
                raise AssertionError(f"support localization failed: m={len(parameters)}, S={s}")
            if s + 1 > 4 * variance_x:
                raise AssertionError(
                    f"one-sided localization failed: S={s}, V_X={variance_x}"
                )

            for q1, q2 in itertools.combinations_with_replacement(Q_VALUES, 2):
                row = verify_row(parameters, a, variance_x, s, q1, q2)
                rows += 1
                selections[row["selection"]] += 1
                delta0_equal_rows += int(row["delta0_equal"])
                delta1_equal_rows += int(row["delta1_equal"])
                if row["selection"] == "S-2":
                    shifted_newton_checks += 1
                else:
                    payment_checks += 1

                if minimum is None or row["scaled"] < minimum[0]:
                    minimum = (
                        row["scaled"],
                        {
                            "family_first_seen": family,
                            "parameters": parameter_blocks(parameters),
                            "variance_x": fraction_record(variance_x),
                            "q": [fraction_text(q1), fraction_text(q2)],
                            "S": s,
                            "D": row["descent"],
                            "selection": row["selection"],
                            "delta": fraction_record(row["delta"]),
                            "curvature_lower": fraction_record(
                                row["curvature_lower"]
                            ),
                            "proof_lower": fraction_record(row["proof_lower"]),
                            "variance_total": fraction_record(row["variance_total"]),
                            "scaled": fraction_record(row["scaled"]),
                            "proof_scaled": fraction_record(row["proof_scaled"]),
                        },
                    )

    if minimum is None:
        raise AssertionError("grid produced no eligible profiles")
    if selections["S"] == 0 or selections["S-1"] == 0:
        raise AssertionError(f"grid missed an observed descent branch: {selections}")

    output = {
        "status": "passed",
        "arithmetic": "fractions.Fraction and exact SymPy polynomials",
        "conditions": {
            "x_parameters": "0 < p_i <= 1/2",
            "variance_x": ">= 1",
            "q_i": "0 <= q_i <= 1/2",
            "q_values": [fraction_text(q) for q in Q_VALUES],
        },
        "symbolic_certificate": symbolic,
        "grid": {
            "core_probability_values": [
                fraction_text(p) for p in CORE_PROBABILITIES
            ],
            "core_m_range": [args.min_m, args.max_m],
            "heavy_probability_values": [
                fraction_text(p) for p in HEAVY_PROBABILITIES
            ],
            "dust_probability_values": [
                fraction_text(p) for p in DUST_PROBABILITIES
            ],
            "heavy_count_range": [2, 8],
            "dust_counts": [1, 2, 4, 8],
            "q_pairs": "combinations with replacement",
        },
        "counts": {
            "profiles_generated": dict(generated),
            "profiles_eligible_and_unique": dict(eligible),
            "unique_profiles_seen": len(seen),
            "rows": rows,
            "selection_cases": dict(selections),
            "delta0_equal_rows": delta0_equal_rows,
            "delta1_equal_rows": delta1_equal_rows,
            "curvature_decomposition_checks": rows,
            "individual_newton_and_long_gap_checks": 6 * rows,
            "normalized_G_reconstructions": rows,
            "variance_payment_checks": payment_checks,
            "shifted_newton_checks": shifted_newton_checks,
            "proof_chain_quarter_checks": rows,
            "actual_quarter_checks": rows,
            "failures": 0,
        },
        "minimum_scaled_effective_drop": minimum[1],
        "anchored_checks": anchors,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")

    print(f"passed {rows:,} exact rows across {sum(eligible.values()):,} profiles")
    print(f"selection cases: {dict(selections)}")
    print(
        "minimum (Var(X)+Var(Y))*Delta = "
        f"{fraction_text(minimum[0])} = {float(minimum[0]):.17g}"
    )
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
