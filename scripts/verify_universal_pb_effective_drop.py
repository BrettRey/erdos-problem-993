#!/usr/bin/env python3
"""Replayable exact certificate for the universal effective-drop theorem.

The proof certificate is symbolic and exact.  It checks the endpoint-aware
Hillion--Johnson recurrence reduction, the forced mass window, the finite
Sturm cells, and the analytic Bernstein tail.  A separate finite
Poisson--binomial grid is only a falsification scan; no proof step depends on
that scan.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations_with_replacement
import json
from pathlib import Path
from typing import Any, Iterable

import sympy as sp


DEFAULT_OUT = Path(
    "results/universal_pb_effective_drop_certificate_2026-07-10.json"
)
SCAN_PROBABILITIES = (
    Fraction(1, 10),
    Fraction(1, 5),
    Fraction(1, 3),
    Fraction(1, 2),
    Fraction(2, 3),
    Fraction(4, 5),
    Fraction(9, 10),
)


def require_zero(expression: sp.Expr, label: str) -> None:
    if sp.factor(expression) != 0:
        raise AssertionError(f"{label} did not simplify to zero: {expression}")


def exact_string(value: Any) -> str:
    if isinstance(value, Fraction):
        return f"{value.numerator}/{value.denominator}"
    return sp.sstr(sp.factor(value))


def sturm_variations(sequence: Iterable[sp.Expr], symbol: sp.Symbol, at: int) -> int:
    signs: list[int] = []
    for expression in sequence:
        value = sp.sign(sp.cancel(expression.subs(symbol, at)))
        if value == 0:
            continue
        if value not in (-1, 1):
            raise AssertionError(f"undetermined Sturm sign at {at}: {value}")
        signs.append(int(value))
    return sum(left != right for left, right in zip(signs, signs[1:]))


def power_to_bernstein(
    polynomial: sp.Expr,
    symbol: sp.Symbol,
    degree: int,
) -> list[sp.Expr]:
    poly = sp.Poly(sp.expand(polynomial), symbol)
    if poly.degree() > degree:
        raise AssertionError("polynomial degree exceeds requested Bernstein degree")
    power = [poly.coeff_monomial(symbol**index) for index in range(degree + 1)]
    return [
        sp.factor(
            sum(
                power[index]
                * sp.binomial(bernstein_index, index)
                / sp.binomial(degree, index)
                for index in range(bernstein_index + 1)
            )
        )
        for bernstein_index in range(degree + 1)
    ]


def check_recurrence_and_mass_algebra() -> dict[str, Any]:
    x, y, d, j, r, h = sp.symbols("x y d j r H", positive=True)

    # The two cubic inequalities are precisely the two one-step bounds on
    # reciprocal curvature.  Their differences factor without assumptions.
    forward_difference = x - y * (1 - x)
    backward_difference = y - x * (1 - y)
    require_zero(
        forward_difference - x * y * (1 / y - 1 / x + 1),
        "forward reciprocal-curvature translation",
    )
    require_zero(
        backward_difference - x * y * (1 / x - 1 / y + 1),
        "backward reciprocal-curvature translation",
    )
    require_zero(
        x - (1 - x) - (2 * x - 1),
        "endpoint curvature is at least one half",
    )

    induction_step = (d / (1 - r * d)) / (1 - d / (1 - r * d))
    require_zero(
        induction_step - d / (1 - (r + 1) * d),
        "curvature-window induction",
    )
    require_zero(
        1 - d / (1 - j * d) - (1 - (j + 1) * d) / (1 - j * d),
        "curvature-factor telescope",
    )

    a_d = (1 - 2 * d) / (1 - d)
    d_of_h = 1 / (h + 1)
    require_zero(
        a_d.subs(d, d_of_h) - (h - 1) / h,
        "H parameterization of crossing factor",
    )

    right_step = a_d * (1 - r * d)
    left_step = (1 - (r + 2) * d) / (1 - d)
    ratio_step = sp.factor((right_step / left_step).subs(d, d_of_h))
    expected_ratio = (h - 1) * (h + 1 - r) / ((h + 1) * (h - r - 1))
    require_zero(ratio_step - expected_ratio, "R/L ratio recurrence")
    require_zero(
        expected_ratio - 1 - 2 * r / ((h + 1) * (h - r - 1)),
        "R/L positivity gap",
    )

    # Bonferroni is used only while the partial sum is at most one.  This
    # identity is the induction step for prod(1-x_i) >= 1-sum(x_i).
    partial_sum, next_term = sp.symbols("S z", nonnegative=True)
    require_zero(
        (1 - partial_sum) * (1 - next_term)
        - (1 - partial_sum - next_term)
        - partial_sum * next_term,
        "Bonferroni product induction",
    )

    return {
        "status": "passed",
        "endpoint_convention": "delta_0 = delta_n = 1",
        "reciprocal_curvature_law": "abs(1/delta_(j+1)-1/delta_j) <= 1",
        "window_induction": "delta_(D+r) <= d/(1-r*d)",
        "endpoint_exclusion_condition": "(r+1)*d < 1",
        "crossing_factor_a": "(1-2*d)/(1-d)",
        "right_mass_bound": "R_r=a^r*product(j=1..r-1,1-j*d)",
        "left_mass_bound": "L_r=(1-d)^(-r)*product(j=2..r+1,1-j*d)",
        "R_over_L_step": "(H-1)*(H+1-r)/((H+1)*(H-r-1))",
        "R_over_L_step_minus_one": "2*r/((H+1)*(H-r-1))",
    }


def window_weights(h: sp.Symbol, count: int) -> tuple[dict[int, sp.Expr], dict[int, sp.Expr]]:
    a = (h - 1) / h
    right = {
        index: sp.factor(
            a**index
            * sp.prod(1 - sp.Rational(j, 1) / (h + 1) for j in range(1, index))
        )
        for index in range(1, count + 1)
    }
    left = {
        index: sp.factor(
            sp.prod(1 - sp.Rational(j, 1) / h for j in range(1, index + 1))
        )
        for index in range(1, count + 1)
    }
    return right, left


def asymmetric_window_form(h: sp.Symbol, count: int) -> sp.Expr:
    right, left = window_weights(h, count)
    total_mass = 1 + sum(right.values()) + sum(left.values())
    first_moment = sum(
        index * (right[index] - left[index]) for index in range(1, count + 1)
    )
    second_moment = sum(
        index * index * (right[index] + left[index])
        for index in range(1, count + 1)
    )
    return sp.factor(total_mass * second_moment - first_moment**2)


def scalar_target(h: sp.Symbol) -> sp.Expr:
    return (h + 1) * (3 * h + 4) / 4


def check_asymmetric_cell() -> dict[str, Any]:
    h = sp.symbols("H")
    difference = sp.cancel(asymmetric_window_form(h, 3) - scalar_target(h))
    numerator, denominator = sp.fraction(difference)
    expected_numerator = (
        -3 * h**10
        - 16 * h**9
        + 750 * h**8
        - 3676 * h**7
        + 6613 * h**6
        - 5460 * h**5
        + 800 * h**4
        + 4696 * h**3
        - 7176 * h**2
        + 4272 * h
        - 912
    )
    require_zero(numerator - expected_numerator, "asymmetric numerator")
    require_zero(denominator - 4 * h**5 * (h + 1) ** 3, "asymmetric denominator")

    polynomial = sp.Poly(numerator, h)
    sequence = sp.sturm(polynomial.as_expr(), h)
    left_variations = sturm_variations(sequence, h, 3)
    right_variations = sturm_variations(sequence, h, 4)
    root_count = int(sp.polys.polytools.count_roots(polynomial, 3, 4))
    if root_count != 0 or left_variations != right_variations:
        raise AssertionError("asymmetric cell has a real root on [3,4]")

    values = {
        "H=3": sp.factor(difference.subs(h, 3)),
        "H=7/2": sp.factor(difference.subs(h, sp.Rational(7, 2))),
        "H=4": sp.factor(difference.subs(h, 4)),
    }
    if not all(value > 0 for value in values.values()):
        raise AssertionError("asymmetric cell is not positive at its checkpoints")

    return {
        "status": "passed",
        "interval": "3 < H < 4",
        "window_radius": 3,
        "numerator": exact_string(expected_numerator),
        "denominator": exact_string(denominator),
        "degree": polynomial.degree(),
        "sturm_sequence_length": len(sequence),
        "sturm_variations_left": left_variations,
        "sturm_variations_right": right_variations,
        "root_count_closed_interval": root_count,
        "checkpoint_values": {key: exact_string(value) for key, value in values.items()},
    }


def symmetric_window_form(h: sp.Symbol, count: int) -> sp.Expr:
    weight = sp.Integer(1)
    weights: list[tuple[int, sp.Expr]] = []
    for index in range(1, count + 1):
        weight = sp.cancel(weight * (h - index) / h)
        weights.append((index, weight))
    total_mass = 1 + 2 * sum(value for _, value in weights)
    second_moment = 2 * sum(index * index * value for index, value in weights)
    return sp.factor(total_mass * second_moment)


def check_compact_sturm_cells() -> list[dict[str, Any]]:
    h = sp.symbols("H")
    rows: list[dict[str, Any]] = []
    for cell in range(4, 16):
        difference = sp.cancel(symmetric_window_form(h, cell) - scalar_target(h))
        numerator, denominator = sp.fraction(difference)
        require_zero(
            denominator - 4 * h ** (2 * cell),
            f"compact-cell denominator m={cell}",
        )
        polynomial = sp.Poly(numerator, h)
        expected_degree = 2 * cell + 2
        if polynomial.degree() != expected_degree:
            raise AssertionError(
                f"unexpected degree in cell {cell}: {polynomial.degree()}"
            )
        sequence = sp.sturm(polynomial.as_expr(), h)
        left_variations = sturm_variations(sequence, h, cell)
        right_variations = sturm_variations(sequence, h, cell + 1)
        root_count = int(
            sp.polys.polytools.count_roots(polynomial, cell, cell + 1)
        )
        midpoint = sp.Rational(2 * cell + 1, 2)
        endpoint_values = (
            polynomial.eval(cell),
            polynomial.eval(midpoint),
            polynomial.eval(cell + 1),
        )
        if root_count != 0 or left_variations != right_variations:
            raise AssertionError(f"Sturm root in compact cell m={cell}")
        if not all(value > 0 for value in endpoint_values):
            raise AssertionError(f"nonpositive compact-cell checkpoint m={cell}")
        rows.append(
            {
                "m": cell,
                "interval": f"[{cell},{cell + 1}]",
                "degree": expected_degree,
                "denominator": f"4*H^{2 * cell}",
                "sturm_sequence_length": len(sequence),
                "sturm_variations_left": left_variations,
                "sturm_variations_right": right_variations,
                "root_count_closed_interval": root_count,
                "left_numerator_positive": True,
                "midpoint_numerator_positive": True,
                "right_numerator_positive": True,
            }
        )
    return rows


def check_analytic_tail() -> dict[str, Any]:
    h, radius, t, u = sp.symbols("H R t u")
    index = sp.symbols("r", integer=True, positive=True)
    c = 1 - index * (index + 1) / (2 * h)
    sum_zero = sp.factor(1 + 2 * sp.summation(c, (index, 1, radius)))
    sum_two = sp.factor(
        2 * sp.summation(index**2 * c, (index, 1, radius))
    )
    expected_zero = (
        6 * h * radius
        + 3 * h
        - radius**3
        - 3 * radius**2
        - 2 * radius
    ) / (3 * h)
    expected_two = (
        radius
        * (radius + 1)
        * (
            40 * h * radius
            + 20 * h
            - 12 * radius**3
            - 33 * radius**2
            - 17 * radius
            + 2
        )
        / (60 * h)
    )
    require_zero(sum_zero - expected_zero, "analytic-tail S0 sum")
    require_zero(sum_two - expected_two, "analytic-tail T0 sum")

    tail_difference = sp.factor(sum_zero * sum_two - scalar_target(h))
    scaled = sp.expand(180 * h**2 * tail_difference)

    radius_five = sp.expand(scaled.subs({radius: 5, h: 16 + 5 * t}))
    radius_five_bernstein = power_to_bernstein(radius_five, t, 4)
    expected_radius_five = [
        sp.Integer(424800),
        sp.Integer(1350000),
        sp.Integer(2254950),
        sp.Rational(11439225, 4),
        sp.Integer(2800350),
    ]
    for actual, expected in zip(radius_five_bernstein, expected_radius_five):
        require_zero(actual - expected, "R=5 Bernstein coefficient")

    lower = radius * (radius + 1) / 2
    triangular_cell = sp.expand(
        scaled.subs(h, lower + (radius + 1) * t)
    )
    bernstein = power_to_bernstein(triangular_cell, t, 4)
    shifted = [sp.factor(value.subs(radius, u + 6)) for value in bernstein]

    q_polynomials = [
        121 * u**4 + 2474 * u**3 + 17431 * u**2 + 46014 * u + 25400,
        121 * u**5
        + 3442 * u**4
        + 37967 * u**3
        + 199405 * u**2
        + 480475 * u
        + 384510,
        121 * u**6
        + 4410 * u**5
        + 65863 * u**4
        + 512764 * u**3
        + 2172926 * u**2
        + 4668316 * u
        + 3829440,
        121 * u**5
        + 3684 * u**4
        + 43735 * u**3
        + 249961 * u**2
        + 671859 * u
        + 644400,
        121 * u**4 + 2958 * u**3 + 25579 * u**2 + 88782 * u + 91440,
    ]
    prefactors = [
        (u + 6) ** 2 * (u + 7) ** 2 / 16,
        (u + 6) * (u + 7) ** 2 / 16,
        (u + 7) ** 2 / 16,
        (u + 7) ** 2 * (u + 8) / 16,
        (u + 7) ** 2 * (u + 8) ** 2 / 16,
    ]
    for index_value, (actual, prefactor, polynomial) in enumerate(
        zip(shifted, prefactors, q_polynomials)
    ):
        require_zero(
            actual - prefactor * polynomial,
            f"analytic Bernstein factor index={index_value}",
        )
        coefficients = sp.Poly(polynomial, u).all_coeffs()
        if not coefficients or not all(coefficient > 0 for coefficient in coefficients):
            raise AssertionError(
                f"shifted Bernstein polynomial {index_value} is not coefficient-positive"
            )

    return {
        "status": "passed",
        "bonferroni_bound": "b_r >= c_r = 1-r(r+1)/(2H)",
        "S0": exact_string(expected_zero),
        "T0": exact_string(expected_two),
        "scaled_polynomial": "G=180*H^2*(S0*T0-target)",
        "R=5_interval": "[16,21]",
        "R=5_bernstein_coefficients": [
            exact_string(value) for value in radius_five_bernstein
        ],
        "R>=6_cell": "H=R(R+1)/2+(R+1)t, 0<=t<=1",
        "R>=6_shift": "u=R-6>=0",
        "R>=6_bernstein_prefactors": [exact_string(value) for value in prefactors],
        "R>=6_positive_polynomials": [
            exact_string(value) for value in q_polynomials
        ],
    }


def check_max_atom_and_scalar_algebra() -> dict[str, Any]:
    mass, d = sp.symbols("M d", positive=True)
    u = sp.symbols("u", real=True)
    uniform_cell_variance = sp.integrate(u**2, (u, -sp.Rational(1, 2), sp.Rational(1, 2)))
    bathtub_variance = mass * sp.integrate(
        u**2,
        (u, -1 / (2 * mass), 1 / (2 * mass)),
    )
    require_zero(uniform_cell_variance - sp.Rational(1, 12), "unit-cell variance")
    require_zero(
        bathtub_variance - 1 / (12 * mass**2),
        "bounded-density bathtub variance",
    )

    p0, p1, x0, x1, x2 = sp.symbols("p0 p1 x0 x1 x2")
    p2 = 1 - p0 - p1
    mean = p0 * x0 + p1 * x1 + p2 * x2
    variance = (
        p0 * (x0 - mean) ** 2
        + p1 * (x1 - mean) ** 2
        + p2 * (x2 - mean) ** 2
    )
    pairwise = (
        p0 * p1 * (x0 - x1) ** 2
        + p0 * p2 * (x0 - x2) ** 2
        + p1 * p2 * (x1 - x2) ** 2
    )
    require_zero(variance - pairwise, "pairwise variance identity")

    scalar_threshold = (3 + d) / (4 * d**2)
    target = 1 / (4 * d)
    require_zero(
        target / scalar_threshold - d / (3 + d),
        "pairwise branch mass threshold",
    )
    require_zero(
        ((3 + d) / d - 1) / 12 - target,
        "max-atom branch reaches quarter target",
    )

    return {
        "status": "passed",
        "smoothed_variance": "Var(X+U)=Var(X)+1/12",
        "bathtub_minimum": "1/(12*M^2)",
        "max_atom_bound": "Var(X)>=(M^(-2)-1)/12",
        "pairwise_bound": "Var(X)>=M^2*A",
        "scalar_threshold": exact_string(scalar_threshold),
        "quarter_target": exact_string(target),
        "branch_mass_squared": "d/(3+d)",
    }


def convolve_bernoulli(
    probabilities: tuple[Fraction, ...],
) -> list[Fraction]:
    coefficients = [Fraction(1)]
    for probability in probabilities:
        updated = [Fraction(0)] * (len(coefficients) + 1)
        for index, coefficient in enumerate(coefficients):
            updated[index] += coefficient * (1 - probability)
            updated[index + 1] += coefficient * probability
        coefficients = updated
    return coefficients


def exact_pbd_row(probabilities: tuple[Fraction, ...]) -> dict[str, Any] | None:
    coefficients = convolve_bernoulli(probabilities)
    degree = len(probabilities)
    variance = sum(probability * (1 - probability) for probability in probabilities)

    ratios = {
        index: coefficients[index] / coefficients[index - 1]
        for index in range(1, degree + 1)
    }
    descent = next((index for index in range(1, degree + 1) if ratios[index] < 1), None)
    if descent is None:
        # The first descent is the zero coefficient outside the support, where
        # the effective quotient is undefined.  It is not part of this theorem.
        return None

    deltas = {0: Fraction(1), degree: Fraction(1)}
    for index in range(1, degree):
        deltas[index] = 1 - ratios[index + 1] / ratios[index]
    for index in range(degree):
        left = deltas[index]
        right = deltas[index + 1]
        if right * (1 - left) > left or left * (1 - right) > right:
            raise AssertionError(
                f"HJ recurrence failed on exact PBD grid at n={degree}, j={index}"
            )

    drop = deltas[descent]
    scaled = variance * drop
    if variance >= 1 and scaled < Fraction(1, 4):
        raise AssertionError(
            f"finite exact PBD counterexample: p={probabilities}, V*delta={scaled}"
        )
    return {
        "n": degree,
        "variance": variance,
        "descent": descent,
        "drop": drop,
        "scaled": scaled,
        "probabilities": probabilities,
    }


def run_exact_pbd_scan(max_n: int) -> dict[str, Any]:
    processed = 0
    eligible = 0
    terminal_outside_support = 0
    best: dict[str, Any] | None = None
    for degree in range(4, max_n + 1):
        for probabilities in combinations_with_replacement(SCAN_PROBABILITIES, degree):
            processed += 1
            row = exact_pbd_row(probabilities)
            if row is None:
                terminal_outside_support += 1
                continue
            if row["variance"] < 1:
                continue
            eligible += 1
            if best is None or row["scaled"] < best["scaled"]:
                best = row
    if best is None:
        raise AssertionError("finite exact PBD grid had no eligible rows")
    return {
        "status": "no counterexample on this finite grid",
        "role": "falsification scan only; not used by the proof certificate",
        "probability_grid": [exact_string(value) for value in SCAN_PROBABILITIES],
        "min_n": 4,
        "max_n": max_n,
        "vectors_processed": processed,
        "eligible_variance_at_least_one": eligible,
        "terminal_outside_support_skipped": terminal_outside_support,
        "minimum_scaled_drop": exact_string(best["scaled"]),
        "minimum_scaled_drop_decimal": float(best["scaled"]),
        "minimizer": {
            "n": best["n"],
            "probabilities": [exact_string(value) for value in best["probabilities"]],
            "variance": exact_string(best["variance"]),
            "descent": best["descent"],
            "drop": exact_string(best["drop"]),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--scan-max-n", type=int, default=9)
    parser.add_argument("--skip-scan", action="store_true")
    args = parser.parse_args()
    if args.scan_max_n < 4:
        parser.error("--scan-max-n must be at least 4")

    recurrence = check_recurrence_and_mass_algebra()
    max_atom = check_max_atom_and_scalar_algebra()
    asymmetric = check_asymmetric_cell()
    compact = check_compact_sturm_cells()
    tail = check_analytic_tail()
    scan = None if args.skip_scan else run_exact_pbd_scan(args.scan_max_n)

    result = {
        "certificate": "universal Poisson-binomial effective-drop theorem",
        "certificate_date": "2026-07-10",
        "arithmetic": "exact SymPy rationals/polynomials and fractions.Fraction",
        "sympy_version": sp.__version__,
        "theorem_scope": {
            "variance_assumption": "Var >= 1",
            "endpoint_convention": "delta_0=delta_n=1 with f outside support zero",
            "cubic_recurrence": "delta_(j+/-1)*(1-delta_j) <= delta_j",
            "conclusion": "Var*delta_D >= 1/4 at the first supported strict descent",
        },
        "proof_certificates": {
            "recurrence_and_mass_window": recurrence,
            "max_atom_and_scalar_reduction": max_atom,
            "asymmetric_cell": asymmetric,
            "compact_sturm_cells": compact,
            "analytic_tail": tail,
        },
        "finite_exact_pbd_scan": scan,
        "status": "passed",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "compact_cells": len(compact),
        "scan_vectors": None if scan is None else scan["vectors_processed"],
        "minimum_scan_scaled_drop": None if scan is None else scan["minimum_scaled_drop"],
        "out": str(args.out),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
