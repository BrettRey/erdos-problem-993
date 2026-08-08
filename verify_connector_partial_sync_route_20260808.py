#!/usr/bin/env python3
"""Exact audits for the connector-subdivided double-star theorem.

The theorem is proved algebraically in the accompanying note.  These finite
checks independently audit its formulas, coefficient expansions, and tree DP.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import comb
from pathlib import Path

from indpoly import independence_poly, is_log_concave, is_unimodal
from verify_double_star_log_concavity_20260808 import (
    binomial_poly,
    poly_add,
    poly_mul,
)
from verify_two_hub_ratio_order_obstructions_20260808 import (
    partial_synchronization_failure,
    synchronization_failures,
)


def shift(poly: list[int], amount: int = 1) -> list[int]:
    return [0] * amount + poly


def choose(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0


def weak_synchronization_failure(
    left: list[int],
    right: list[int],
) -> dict[str, int] | None:
    degree = max(len(left), len(right)) - 1
    a = left + [0] * (degree + 2 - len(left))
    b = right + [0] * (degree + 2 - len(right))
    for k in range(degree + 1):
        a_previous = 0 if k == 0 else a[k - 1]
        b_previous = 0 if k == 0 else b[k - 1]
        lhs = 2 * a[k] * b[k]
        rhs = a_previous * b[k + 1] + a[k + 1] * b_previous
        if lhs < rhs:
            return {"k": k, "lhs": lhs, "rhs": rhs}
    return None


def contraction_star_poly(p: int, q: int) -> list[int]:
    """Return C=(1+x)^(p+q)+x from the base proof."""
    return poly_add(binomial_poly(p + q), [0, 1])


def v_basis_coefficient(k: int, i: int, j: int) -> int:
    """Coefficient of binom(p,i)binom(q,j) in the hard V_k gap."""
    r = i + j
    h = 2 * k - r
    return (
        choose(r, k) * choose(k, h)
        - choose(r, k + 1) * choose(k + 1, h)
        + (choose(i, k - 1) + choose(j, k - 1)) * choose(k - 1, h - 1)
        - (choose(i, k - 2) + choose(j, k - 2)) * choose(k - 2, h - 1)
    )


def w_basis_coefficient(k: int, i: int, j: int) -> int:
    """Coefficient of binom(p,i)binom(q,j) in the hard W_k gap."""
    r = i + j
    h = 2 * k - 1 - r
    return (
        choose(r, k - 1) * choose(k - 1, h)
        - choose(r, k - 2) * choose(k - 2, h)
        + 2
        * (choose(i, k - 1) + choose(j, k - 1))
        * choose(k - 1, h - 1)
        - (choose(i, k) + choose(j, k)) * choose(k, h - 1)
        - (choose(i, k - 2) + choose(j, k - 2)) * choose(k - 2, h - 1)
    )


def base_gap(p: int, q: int, k: int, kind: str) -> int:
    """Evaluate one synchronization gap directly from its coefficients."""
    base = connector_formula(p, q, 0)
    contraction = contraction_star_poly(p, q)
    degree = max(len(base), len(contraction), k + 2) + 2
    b = base + [0] * (degree - len(base))
    c = contraction + [0] * (degree - len(contraction))

    def coefficient(poly: list[int], index: int) -> int:
        return 0 if index < 0 else poly[index]

    if kind == "v":
        return (
            coefficient(c, k) * coefficient(b, k)
            - coefficient(c, k + 1) * coefficient(b, k - 1)
        )
    if kind == "w":
        return (
            2 * coefficient(c, k - 1) * coefficient(b, k)
            - coefficient(c, k - 2) * coefficient(b, k + 1)
            - coefficient(c, k) * coefficient(b, k - 1)
        )
    raise ValueError(f"unknown gap kind: {kind}")


def finite_difference_basis_coefficient(kind: str, k: int, i: int, j: int) -> int:
    """Extract a product-binomial-basis coefficient by finite differences."""
    return sum(
        (-1) ** (i - u + j - v)
        * choose(i, u)
        * choose(j, v)
        * base_gap(u, v, k, kind)
        for u in range(i + 1)
        for v in range(j + 1)
    )


def audit_finite_difference_basis(max_k: int) -> int:
    """Independently derive the claimed basis coefficients for small k."""
    checked = 0
    for k in range(4, max_k + 1):
        for total in range(2 * k + 2):
            for i in range(total + 1):
                j = total - i
                expected_v = v_basis_coefficient(k, i, j)
                expected_w = w_basis_coefficient(k, i, j)
                assert finite_difference_basis_coefficient("v", k, i, j) == expected_v
                assert finite_difference_basis_coefficient("w", k, i, j) == expected_w
                checked += 2
    return checked


def audit_binomial_basis(max_k: int) -> int:
    """Audit the general coefficient identities and proof bounds exactly."""
    checked = 0
    for k in range(4, max_k + 1):
        for r in range(k, 2 * k + 1):
            h = 2 * k - r
            assert k * (k - 1) >= 2 * h * (k - h + 1)
            for i in range(r + 1):
                j = r - i
                direct = v_basis_coefficient(k, i, j)
                simplified = Fraction(choose(r, k) * choose(k, h), k - h + 1)
                simplified += Fraction(choose(k - 1, h - 1), k - 1) * (
                    choose(i, k - 2) * (2 - j)
                    + choose(j, k - 2) * (2 - i)
                )
                negative_bound = (
                    max(j - 2, 0) * choose(i, k - 2)
                    + max(i - 2, 0) * choose(j, k - 2)
                )
                assert negative_bound <= 2 * choose(r, k)
                assert simplified.denominator == 1
                assert direct == simplified.numerator
                assert direct >= 0
                checked += 1

        for r in range(k - 1, 2 * k):
            h = 2 * k - 1 - r
            assert k * (k - 1) >= 2 * h
            for i in range(r + 1):
                j = r - i
                direct = w_basis_coefficient(k, i, j)
                denominator = k - h + 1
                simplified = Fraction(
                    2 * choose(r, k - 1) * choose(k - 1, h),
                    denominator,
                )
                if h >= 1:
                    simplified += Fraction(
                        choose(k - 1, h - 1),
                        (k - 1) * denominator,
                    ) * (
                        choose(i, k - 2)
                        * (2 * (i - k + 2) - j * (j - 1))
                        + choose(j, k - 2)
                        * (2 * (j - k + 2) - i * (i - 1))
                    )
                negative_bound = (
                    j * (j - 1) * choose(i, k - 2)
                    + i * (i - 1) * choose(j, k - 2)
                )
                assert negative_bound <= 4 * choose(r, k)
                assert simplified.denominator == 1
                assert direct == simplified.numerator
                assert direct >= 0
                checked += 1
    return checked


def path_poly(order: int) -> list[int]:
    """Return I(P_order), with P_-2=0 and P_-1=P_0=1."""
    if order == -2:
        return [0]
    if order in (-1, 0):
        return [1]
    previous_previous = [1]
    previous = [1]
    for _ in range(1, order + 1):
        previous_previous, previous = previous, poly_add(
            previous,
            shift(previous_previous),
        )
    return previous


def connector_formula(p: int, q: int, internal_vertices: int) -> list[int]:
    """Independence polynomial for two leaf-stars joined through a path."""
    first = poly_mul(binomial_poly(p + q), path_poly(internal_vertices))
    second = shift(
        poly_mul(
            poly_add(binomial_poly(p), binomial_poly(q)),
            path_poly(internal_vertices - 1),
        )
    )
    third = shift(path_poly(internal_vertices - 2), 2)
    return poly_add(poly_add(first, second), third)


def connector_adj(p: int, q: int, internal_vertices: int) -> list[list[int]]:
    """Build the connector tree explicitly for independent tree-DP replay."""
    order = p + q + internal_vertices + 2
    adj = [[] for _ in range(order)]

    def add_edge(left: int, right: int) -> None:
        adj[left].append(right)
        adj[right].append(left)

    right_hub = internal_vertices + 1
    for vertex in range(right_hub):
        add_edge(vertex, vertex + 1)
    next_vertex = right_hub + 1
    for _ in range(p):
        add_edge(0, next_vertex)
        next_vertex += 1
    for _ in range(q):
        add_edge(right_hub, next_vertex)
        next_vertex += 1
    return adj


def verify(
    max_leaves: int,
    ratio_max_leaves: int,
    max_internal: int,
    dp_max_leaves: int,
    dp_max_internal: int,
    basis_max_k: int,
    basis_finite_difference_max_k: int,
) -> dict[str, object]:
    dp_cases = 0
    for p in range(dp_max_leaves + 1):
        for q in range(dp_max_leaves + 1):
            for internal in range(dp_max_internal + 1):
                adj = connector_adj(p, q, internal)
                assert independence_poly(len(adj), adj) == connector_formula(
                    p,
                    q,
                    internal,
                )
                dp_cases += 1

    ratio_base_cases = 0
    for p in range(ratio_max_leaves + 1):
        for q in range(ratio_max_leaves + 1):
            base = connector_formula(p, q, 0)
            contraction = contraction_star_poly(p, q)
            assert not synchronization_failures(contraction, base)
            assert weak_synchronization_failure(shift(contraction), base) is None
            ratio_base_cases += 1

    full_partial_base_cases = 0
    recurrence_cases = 0
    for p in range(max_leaves + 1):
        for q in range(max_leaves + 1):
            polynomials = [
                connector_formula(p, q, internal)
                for internal in range(max_internal + 1)
            ]
            assert partial_synchronization_failure(polynomials[1], polynomials[0]) is None
            assert partial_synchronization_failure(
                polynomials[1],
                shift(polynomials[0]),
            ) is None
            full_partial_base_cases += 1
            for internal in range(2, max_internal + 1):
                assert polynomials[internal] == poly_add(
                    polynomials[internal - 1],
                    shift(polynomials[internal - 2]),
                )
                assert is_log_concave(polynomials[internal])
                assert is_unimodal(polynomials[internal])
                recurrence_cases += 1

    basis_coefficients_checked = audit_binomial_basis(basis_max_k)
    finite_difference_coefficients_checked = audit_finite_difference_basis(
        basis_finite_difference_max_k
    )

    # Ordinary synchronization does not propagate: an exact route obstruction
    # used to prevent accidental promotion of the stronger claim.
    p = q = 2
    internal = 3
    previous = connector_formula(p, q, internal - 1)
    shifted_previous_previous = shift(connector_formula(p, q, internal - 2))
    sync_failure = synchronization_failures(previous, shifted_previous_previous)
    assert sync_failure
    assert connector_formula(p, q, internal) == poly_add(
        previous,
        shifted_previous_previous,
    )

    return {
        "theorem_scope": (
            "The algebraic proof covers all p,q,t >= 0 for double stars with "
            "unit pendant arms and t internal hub-connector vertices."
        ),
        "computation_role": (
            "Finite exact audit of formulas and general coefficient identities; "
            "the finite ranges are not the proof of the theorem."
        ),
        "max_leaves_each_side": max_leaves,
        "ratio_audit_max_leaves_each_side": ratio_max_leaves,
        "max_internal_connector_vertices": max_internal,
        "full_partial_base_parameter_cases": full_partial_base_cases,
        "ratio_base_parameter_cases": ratio_base_cases,
        "recurrence_parameter_cases": recurrence_cases,
        "tree_dp_cases": dp_cases,
        "binomial_basis_max_k": basis_max_k,
        "binomial_basis_coefficients_checked": basis_coefficients_checked,
        "basis_finite_difference_max_k": basis_finite_difference_max_k,
        "finite_difference_coefficients_checked": (
            finite_difference_coefficients_checked
        ),
        "ordinary_synchronization_obstruction": {
            "p": p,
            "q": q,
            "internal_vertices": internal,
            "summand_failure": sync_failure[0],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-leaves", type=int, default=30)
    parser.add_argument("--ratio-max-leaves", type=int, default=100)
    parser.add_argument("--max-internal", type=int, default=20)
    parser.add_argument("--dp-max-leaves", type=int, default=5)
    parser.add_argument("--dp-max-internal", type=int, default=8)
    parser.add_argument("--basis-max-k", type=int, default=60)
    parser.add_argument("--basis-finite-difference-max-k", type=int, default=10)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/connector_partial_sync_route_20260808.json"),
    )
    args = parser.parse_args()
    report = verify(
        args.max_leaves,
        args.ratio_max_leaves,
        args.max_internal,
        args.dp_max_leaves,
        args.dp_max_internal,
        args.basis_max_k,
        args.basis_finite_difference_max_k,
    )
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"event": "verified", "out": str(args.out), **report}))


if __name__ == "__main__":
    main()
