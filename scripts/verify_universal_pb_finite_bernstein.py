#!/usr/bin/env python3
"""Exact Bernstein replacement for the universal PBD finite Sturm cells.

For each finite scalar cell in the universal effective-drop proof, substitute
``H = m + t`` and convert the resulting numerator to the full-degree
Bernstein basis on ``0 <= t <= 1``.  Every coefficient is checked exactly for
strict positivity, and the Bernstein reconstruction is checked symbolically.

The optional Lean emitter writes exact coefficient data and ``ring`` theorem
stubs.  It is deliberately independent of any Aristotle submission.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import sympy as sp


DEFAULT_OUT = Path(
    "results/universal_pb_finite_bernstein_certificate_2026-07-10.json"
)


def power_to_bernstein(
    polynomial: sp.Expr,
    symbol: sp.Symbol,
    degree: int,
) -> list[sp.Rational]:
    """Return coefficients in the degree-``degree`` Bernstein basis."""
    poly = sp.Poly(sp.expand(polynomial), symbol, domain=sp.QQ)
    if poly.degree() > degree:
        raise AssertionError("polynomial degree exceeds Bernstein degree")
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


def bernstein_reconstruction(
    coefficients: list[sp.Rational],
    symbol: sp.Symbol,
) -> sp.Expr:
    degree = len(coefficients) - 1
    return sp.expand(
        sum(
            coefficient
            * sp.binomial(degree, index)
            * symbol**index
            * (1 - symbol) ** (degree - index)
            for index, coefficient in enumerate(coefficients)
        )
    )


def window_weights(
    h: sp.Symbol,
    count: int,
) -> tuple[dict[int, sp.Expr], dict[int, sp.Expr]]:
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


def symmetric_window_form(h: sp.Symbol, count: int) -> sp.Expr:
    weight = sp.Integer(1)
    weights: list[tuple[int, sp.Expr]] = []
    for index in range(1, count + 1):
        weight = sp.cancel(weight * (h - index) / h)
        weights.append((index, weight))
    total_mass = 1 + 2 * sum(value for _, value in weights)
    second_moment = 2 * sum(index * index * value for index, value in weights)
    return sp.factor(total_mass * second_moment)


def scalar_target(h: sp.Symbol) -> sp.Expr:
    return (h + 1) * (3 * h + 4) / 4


def exact_string(value: sp.Expr) -> str:
    value = sp.Rational(value)
    if value.q == 1:
        return str(value.p)
    return f"{value.p}/{value.q}"


@dataclass(frozen=True)
class Cell:
    name: str
    left_endpoint: int
    degree: int
    numerator_in_t: sp.Expr
    coefficients: tuple[sp.Rational, ...]

    @property
    def digest(self) -> str:
        payload = ",".join(exact_string(value) for value in self.coefficients)
        return sha256(payload.encode("ascii")).hexdigest()

    def summary(self) -> dict[str, Any]:
        numerators = [abs(int(value.p)) for value in self.coefficients]
        denominators = [int(value.q) for value in self.coefficients]
        return {
            "name": self.name,
            "interval": f"[{self.left_endpoint},{self.left_endpoint + 1}]",
            "degree": self.degree,
            "coefficient_count": len(self.coefficients),
            "all_strictly_positive": all(value > 0 for value in self.coefficients),
            "minimum_coefficient": exact_string(min(self.coefficients)),
            "maximum_numerator_digits": max(len(str(value)) for value in numerators),
            "maximum_denominator_digits": max(len(str(value)) for value in denominators),
            "coefficient_sha256": self.digest,
        }


def make_cell(
    name: str,
    left_endpoint: int,
    difference: sp.Expr,
    expected_denominator: sp.Expr,
    expected_degree: int,
    h: sp.Symbol,
    t: sp.Symbol,
) -> Cell:
    numerator, denominator = sp.fraction(sp.cancel(difference))
    if sp.factor(denominator - expected_denominator) != 0:
        raise AssertionError(f"{name}: unexpected denominator {denominator}")
    if sp.sign(denominator.subs(h, sp.Rational(2 * left_endpoint + 1, 2))) != 1:
        raise AssertionError(f"{name}: denominator is not positive at midpoint")
    numerator_in_t = sp.expand(numerator.subs(h, left_endpoint + t))
    degree = sp.Poly(numerator_in_t, t, domain=sp.QQ).degree()
    if degree != expected_degree:
        raise AssertionError(
            f"{name}: expected degree {expected_degree}, found {degree}"
        )
    coefficients = tuple(power_to_bernstein(numerator_in_t, t, degree))
    if not coefficients or not all(value > 0 for value in coefficients):
        raise AssertionError(f"{name}: nonpositive Bernstein coefficient")
    reconstructed = bernstein_reconstruction(list(coefficients), t)
    if sp.expand(reconstructed - numerator_in_t) != 0:
        raise AssertionError(f"{name}: Bernstein reconstruction failed")
    return Cell(name, left_endpoint, degree, numerator_in_t, coefficients)


def build_cells() -> list[Cell]:
    h, t = sp.symbols("H t")
    cells = [
        make_cell(
            "asymmetric_3_4",
            3,
            asymmetric_window_form(h, 3) - scalar_target(h),
            4 * h**5 * (h + 1) ** 3,
            10,
            h,
            t,
        )
    ]
    cells.extend(
        make_cell(
            f"symmetric_{cell}_{cell + 1}",
            cell,
            symmetric_window_form(h, cell) - scalar_target(h),
            4 * h ** (2 * cell),
            2 * cell + 2,
            h,
            t,
        )
        for cell in range(4, 16)
    )
    return cells


def lean_rational(value: sp.Rational) -> str:
    if value.q == 1:
        return f"({value.p} : ℝ)"
    return f"(({value.p} : ℝ) / {value.q})"


def lean_polynomial(polynomial: sp.Expr, symbol: sp.Symbol) -> str:
    poly = sp.Poly(sp.expand(polynomial), symbol, domain=sp.QQ)
    terms: list[str] = []
    for (power,), coefficient in poly.terms():
        scalar = lean_rational(abs(coefficient))
        monomial = scalar if power == 0 else f"{scalar} * t ^ {power}"
        if not terms:
            terms.append(monomial if coefficient > 0 else f"-{monomial}")
        else:
            terms.append((" + " if coefficient > 0 else " - ") + monomial)
    return "".join(terms) if terms else "(0 : ℝ)"


def emit_lean(cells: list[Cell], destination: Path) -> None:
    t = sp.symbols("t")
    lines = [
        "/- Generated by scripts/verify_universal_pb_finite_bernstein.py. -/",
        "import Mathlib",
        "",
        "namespace Formal.UniversalPBFiniteBernstein",
        "",
    ]
    for cell in cells:
        identifier = cell.name
        coefficient_values = ", ".join(
            exact_string(value) for value in cell.coefficients
        )
        lines.extend(
            [
                f"/-- Bernstein coefficients for {cell.name}; degree {cell.degree}. -/",
                f"def {identifier}Coefficients : List ℚ := [{coefficient_values}]",
                "",
                f"theorem {identifier}_coefficients_positive :",
                f"    ∀ c ∈ {identifier}Coefficients, 0 < c := by",
                f"  norm_num [{identifier}Coefficients]",
                "",
                f"theorem {identifier}_bernstein_identity (t : ℝ) :",
                f"    {lean_polynomial(cell.numerator_in_t, t)} =",
            ]
        )
        basis_terms = []
        for index, coefficient in enumerate(cell.coefficients):
            effective = sp.Rational(coefficient) * sp.binomial(cell.degree, index)
            basis_terms.append(
                f"{lean_rational(effective)} * t ^ {index} * "
                f"(1 - t) ^ {cell.degree - index}"
            )
        lines.append("      " + " +\n      ".join(basis_terms) + " := by")
        lines.extend(["  ring", ""])
    lines.extend(["end Formal.UniversalPBFiniteBernstein", ""])
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--emit-lean",
        type=Path,
        help="Optionally emit exact Lean coefficient data and ring theorem stubs.",
    )
    args = parser.parse_args()

    cells = build_cells()
    asymmetric = cells[0]
    compact = cells[1:]
    total_coefficients = sum(len(cell.coefficients) for cell in cells)
    combined_payload = "\n".join(cell.digest for cell in cells)
    result = {
        "certificate": "Bernstein replacement for universal PBD finite scalar cells",
        "certificate_date": "2026-07-10",
        "arithmetic": "exact SymPy rationals/polynomials",
        "sympy_version": sp.__version__,
        "basis": "binomial(n,i)*t^i*(1-t)^(n-i) on 0<=t<=1",
        "substitution": "H=left_endpoint+t",
        "asymmetric_cells": 1,
        "compact_symmetric_cells": len(compact),
        "total_cells": len(cells),
        "total_coefficients": total_coefficients,
        "all_reconstructions_exact": True,
        "all_coefficients_strictly_positive": all(
            value > 0 for cell in cells for value in cell.coefficients
        ),
        "combined_cell_digest_sha256": sha256(
            combined_payload.encode("ascii")
        ).hexdigest(),
        "cells": [cell.summary() for cell in cells],
        "lean_emitter": {
            "available": True,
            "contents": "List ℚ coefficient data, positivity proofs, and ring identities",
        },
        "status": "passed",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if args.emit_lean is not None:
        emit_lean(cells, args.emit_lean)

    print(
        json.dumps(
            {
                "status": result["status"],
                "cells": result["total_cells"],
                "coefficients": total_coefficients,
                "all_positive": result["all_coefficients_strictly_positive"],
                "out": str(args.out),
                "lean": None if args.emit_lean is None else str(args.emit_lean),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
