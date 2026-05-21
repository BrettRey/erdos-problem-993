"""Crude Route-2 perturbation certificate for fixed-r spider lanes.

This complements ``fixed_r_hubon_mode_certificate.py``.  It assumes the
hub-off Route-2 reserve is at least ``1/(R a)`` and checks that the full
hub-on perturbation is smaller than that reserve for all ``a >= threshold``.

The perturbation is split into two deliberately conservative pieces:

* fugacity shift: if T = max_k G_k/F_m, then
  ``|lambda-lambda_0| <= 4T`` and the hub-off mean changes by at most
  ``2 N^2 T`` using ``d mu / d lambda = Var/lambda`` and ``lambda >= 1/2``;
* hub-on mixture in B=S(2^(a-1),r): for ``lambda in [1/2,2]``,
  ``G_B(lambda)/F_B(lambda)`` is bounded by
  ``C_r (3/4)^(a-1)``.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import ceil

import sympy as sp

from fixed_r_huboff_certificate import stabilized_shifts
from fixed_r_hubon_mode_certificate import best_witness_term, shifted_positive
from route2_spider_lane_scan import path_polys


t = sp.symbols("t", integer=True, positive=True)


def eval_poly_fraction(poly: list[int], x: Fraction) -> Fraction:
    total = Fraction(0)
    power = Fraction(1)
    for coeff in poly:
        total += coeff * power
        power *= x
    return total


def tail_step_decreasing_expr(a, n, r: int, power: int):
    """Polynomial positivity equivalent to (a+r)^power*T decreasing by a+=3."""
    return sp.combsimp(
        (a + r) ** power * (a + 1) * (a + 2) * (a + 3)
        - 2 * (a + r + 3) ** power * (n + 1) * (n + 2) * (a - n + 1)
    )


def tail_weighted_step_decreasing_expr(a, n, weight):
    """Polynomial positivity equivalent to weight(a)*T(a) decreasing by a+=3."""
    return sp.combsimp(
        weight * (a + 1) * (a + 2) * (a + 3)
        - 2 * weight.subs(a, a + 3) * (n + 1) * (n + 2) * (a - n + 1)
    )


def check_r(r: int, threshold: int, reserve_denom: int) -> bool:
    paths = path_polys(max(2, r))
    pr = paths[r]
    prm1 = paths[r - 1]
    hub_on_sum = sum(prm1)
    shifts = stabilized_shifts(r, threshold, paths)
    mixture_constant = (
        2
        * eval_poly_fraction(prm1, Fraction(2))
        / eval_poly_fraction(pr, Fraction(1, 2))
    )

    print(f"r={r}, threshold a>={threshold}, reserve denominator R={reserve_denom}")
    print(f"stabilized shifts D_q=3m-2a: {shifts}")
    print(f"mixture constant C_r <= {float(mixture_constant):.6e}")

    all_ok = True
    for residue, shift in shifts.items():
        a0 = threshold
        while a0 % 3 != residue:
            a0 += 1
        threshold_t = ceil((threshold - residue) / 3)
        m0 = (2 * a0 + shift) // 3
        s, tail_bound = best_witness_term(pr, hub_on_sum, a0, m0)
        n0 = m0 - s
        degree_bound = a0 + r
        lambda_shift_bound = 2 * degree_bound * degree_bound * tail_bound
        mixture_bound = (
            degree_bound
            * mixture_constant
            * Fraction(3, 4) ** (a0 - 1)
        )
        perturb_bound = lambda_shift_bound + mixture_bound
        reserve = Fraction(1, reserve_denom * a0)
        tail_small_ok = tail_bound < Fraction(1, 2)
        threshold_ok = perturb_bound < reserve

        a = 3 * t + residue
        m = sp.simplify((2 * a + shift) / 3)
        n = sp.simplify(m - s)
        # The reserve target is proportional to 1/a, so certify decrease of
        # a times each perturbation term, not just the perturbation itself.
        lambda_decreasing_ok = shifted_positive(
            tail_weighted_step_decreasing_expr(a, n, a * (a + r) ** 2),
            threshold_t,
        )
        mixture_decreasing_ok = (a0 + 3) * (a0 + r + 3) * 27 < a0 * (a0 + r) * 64
        all_ok = (
            all_ok
            and tail_small_ok
            and threshold_ok
            and lambda_decreasing_ok
            and mixture_decreasing_ok
        )

        print(
            f"  q={residue}: witness s={s}, n0={n0}, "
            f"lambda_bound={float(lambda_shift_bound):.3e}, "
            f"mixture_bound={float(mixture_bound):.3e}, "
            f"reserve={float(reserve):.3e}, "
            f"tail_small_ok={tail_small_ok}, "
            f"threshold_ok={threshold_ok}, "
            f"lambda_decreasing_ok={lambda_decreasing_ok}, "
            f"mixture_decreasing_ok={mixture_decreasing_ok}"
        )

    return all_ok


def main() -> None:
    ap = argparse.ArgumentParser(description="Fixed-r Route-2 perturbation certificate.")
    ap.add_argument("--r", type=int, required=True)
    ap.add_argument("--threshold", type=int, default=200)
    ap.add_argument("--reserve-denom", type=int, default=1000)
    args = ap.parse_args()

    ok = check_r(args.r, args.threshold, args.reserve_denom)
    assert ok


if __name__ == "__main__":
    main()
