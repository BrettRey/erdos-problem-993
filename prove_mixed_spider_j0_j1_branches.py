#!/usr/bin/env python3
"""Attempted proof of combined-margin positivity on reduced mixed-spider branches.

Branches handled:
  1) j=0 family: T = S(2^k)
  2) j=1 family: T = S(2^k,1)

For these families, combined bound equals focused tie margin:
  combined = mu(T, lambda_m) - (m-1),
  where lambda_m = i_{m-1}/i_m and m = leftmost mode at lambda=1.

This script checks the exact finite base cases and verifies the closed-form
inequality templates used in the analytic bounds.
"""

from __future__ import annotations

import math
from fractions import Fraction


def coeff_j0(k: int, t: int) -> int:
    a = math.comb(k, t) * (1 << t) if 0 <= t <= k else 0
    b = math.comb(k, t - 1) if 1 <= t <= k + 1 else 0
    return a + b


def coeff_j1(k: int, t: int) -> int:
    a = math.comb(k, t) * (1 << t) if 0 <= t <= k else 0
    b = 0
    if 1 <= t <= k + 1:
        b = math.comb(k, t - 1) * ((1 << (t - 1)) + 1)
    return a + b


def margin_from_coeffs(coeffs: list[int], m: int) -> Fraction:
    im1 = coeffs[m - 1]
    im = coeffs[m]
    lam = Fraction(im1, im)
    z = Fraction(0)
    mu_num = Fraction(0)
    p = Fraction(1)
    for idx, c in enumerate(coeffs):
        w = Fraction(c) * p
        z += w
        mu_num += idx * w
        p *= lam
    mu = mu_num / z
    return mu - Fraction(m - 1)


def mode_j0(k: int) -> int:
    return (2 * k + 1) // 3


def mode_j1(k: int) -> int:
    return (2 * k + 3) // 3


def a_rF_minus_third_j1(k: int) -> Fraction:
    """A(r_F) - 1/3 for the j=1 branch.

    r_F = F_{m-1}/F_m where F_t are coefficients from (1+2x)^k(1+x).
    """
    m = mode_j1(k)
    rF = Fraction(m * (2 * k - m + 3), 2 * (k - m + 2) * (2 * k - m + 2))
    A = Fraction(2 * k) * rF / (1 + 2 * rF) + rF / (1 + rF) - Fraction(m - 1)
    return A - Fraction(1, 3)


def j1_closed_form_minus_third(k: int) -> Fraction:
    """Closed forms for A(r_F)-1/3 by residue class k mod 3."""
    r = k % 3
    if r == 0:
        t = k // 3
        num = 68 * t**3 + 102 * t**2 + 44 * t + 6
        den = 3 * (192 * t**4 + 424 * t**3 + 330 * t**2 + 106 * t + 12)
        return Fraction(num, den)
    if r == 1:
        t = (k - 1) // 3
        num = 196 * t**3 + 582 * t**2 + 528 * t + 152
        den = 3 * (192 * t**4 + 776 * t**3 + 1134 * t**2 + 708 * t + 160)
        return Fraction(num, den)
    t = (k - 2) // 3
    num = 132 * t**3 + 492 * t**2 + 588 * t + 228
    den = 3 * (192 * t**4 + 984 * t**3 + 1860 * t**2 + 1536 * t + 468)
    return Fraction(num, den)


def check_finite_exact() -> None:
    # j=0 base cases not covered by the k>=5 template.
    for k in [2, 3, 4]:
        m = mode_j0(k)
        coeffs = [coeff_j0(k, t) for t in range(k + 2)]
        margin = margin_from_coeffs(coeffs, m)
        if margin <= 0:
            raise AssertionError(f"j=0 finite base failed at k={k}: margin={margin}")

    # j=1 base cases not covered by the k>=7 template.
    for k in [1, 2, 3, 4, 5, 6]:
        m = mode_j1(k)
        coeffs = [coeff_j1(k, t) for t in range(k + 2)]
        margin = margin_from_coeffs(coeffs, m)
        if margin <= 0:
            raise AssertionError(f"j=1 finite base failed at k={k}: margin={margin}")


def check_j0_template() -> None:
    # Lower bound used for k>=5:
    #   margin >= 1/3 * ( (k+2)/(k+1) - (2k-5)*(2/3)^k ).
    # The second term's factor is maximized at k=5 since
    #   a_{k+1}/a_k = ((2k-3)/(2k-5))*2/3 < 1 for k>=5.
    a5 = Fraction(2 * 5 - 5) * Fraction(2, 3) ** 5  # = 5*(2/3)^5
    if a5 >= 1:
        raise AssertionError("Unexpected: a5 >= 1")

    for k in range(5, 200):
        lhs = Fraction(k + 2, k + 1) - Fraction(2 * k - 5) * Fraction(2, 3) ** k
        if lhs <= 0:
            raise AssertionError(f"j=0 template bound failed at k={k}: lhs={lhs}")


def check_j1_template() -> None:
    # Verify closed forms equal direct A(r_F)-1/3 values and are positive.
    for k in range(1, 500):
        direct = a_rF_minus_third_j1(k)
        closed = j1_closed_form_minus_third(k)
        if direct != closed:
            raise AssertionError(
                f"j=1 closed form mismatch at k={k}: direct={direct}, closed={closed}"
            )
        if closed <= 0:
            raise AssertionError(f"j=1 A(r_F)-1/3 not positive at k={k}: {closed}")

    # For k>=7, correction term is strictly below 1/3:
    #   (m-2)*(2/3)^(k-1) <= ((2k-3)/3)*(2/3)^(k-1) < 1/3.
    # b_{k+1}/b_k = ((2k-1)/(2k-3))*2/3 < 1 for k>=3, so max at k=7.
    b7 = Fraction(2 * 7 - 3, 3) * Fraction(2, 3) ** 6
    if not (b7 < Fraction(1, 3)):
        raise AssertionError(f"Unexpected: b7 >= 1/3 ({b7})")
    for k in range(7, 200):
        b = Fraction(2 * k - 3, 3) * Fraction(2, 3) ** (k - 1)
        if not (b < Fraction(1, 3)):
            raise AssertionError(f"j=1 correction bound failed at k={k}: b={b}")


def main() -> None:
    check_finite_exact()
    check_j0_template()
    check_j1_template()

    print("Finite exact checks: PASS")
    print("j=0 template inequality checks: PASS")
    print("j=1 closed-form and correction checks: PASS")
    print(
        "Result: reduced-branch proof templates for j=0 and j=1 are internally consistent and positive."
    )


if __name__ == "__main__":
    main()
