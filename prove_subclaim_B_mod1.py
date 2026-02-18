#!/usr/bin/env python3
"""Sub-claim B check/proof scaffold for mixed spiders (k ≡ 1 mod 3).

Claim (Sub-claim B):
  For k >= 4 with k ≡ 1 (mod 3),
    margin(k,0) < margin(k,1),
  where margin(k,j) = mu(S(2^k,1^j), lambda_m) - (m-1).

This script:
1) computes exact base cases t=1..5 (k=3t+1),
2) evaluates the derived analytic lower bound for t>=6.
"""

from __future__ import annotations

from fractions import Fraction

from conjecture_a_mixed_spider_exact_margin import compute_margin


def g_core(t: int) -> Fraction:
    """G(t) = A1(r_F) - A0(a) for k=3t+1."""
    num = 66 * t**4 + 227 * t**3 + 269 * t**2 + 134 * t + 24
    den = 288 * t**5 + 1356 * t**4 + 2477 * t**3 + 2196 * t**2 + 948 * t + 160
    return Fraction(num, den)


def delta0_bound(t: int) -> Fraction:
    """Upper bound on A0(lambda0)-A0(a), valid for t>=2."""
    num = (3 * t + 1) * (2 * t * t - t - 2) * (2 * t + 1)
    den = (3 * t + 2) ** 2 * (t + 2) * (2 ** (2 * t + 1))
    return Fraction(num, den)


def corr1_bound(t: int) -> Fraction:
    """Upper bound on p1*g1 term via p1<= (2/3)^(k-1), g1<k/3-1."""
    # k = 3t+1 => k-1 = 3t
    return Fraction(3 * t - 2, 3) * (Fraction(2, 3) ** (3 * t))


def analytic_lb(t: int) -> Fraction:
    """Derived lower bound for margin(k,1)-margin(k,0), t>=2."""
    return g_core(t) - delta0_bound(t) - corr1_bound(t)


def exact_gap_t(t: int) -> Fraction:
    """Exact gap margin(k,1)-margin(k,0) at k=3t+1."""
    k = 3 * t + 1
    m0 = compute_margin(k, 0).margin
    m1 = compute_margin(k, 1).margin
    return m1 - m0


def main() -> None:
    print("Sub-claim B (k ≡ 1 mod 3): margin(k,0) < margin(k,1)")
    print("k = 3t+1")

    print("\nExact base cases:")
    for t in range(1, 6):
        k = 3 * t + 1
        gap = exact_gap_t(t)
        print(f"  t={t:2d}, k={k:2d}, gap={float(gap):.12f} ({'PASS' if gap > 0 else 'FAIL'})")

    print("\nAnalytic bound check (t=6..500):")
    min_lb = None
    min_t = None
    for t in range(6, 501):
        lb = analytic_lb(t)
        if min_lb is None or lb < min_lb:
            min_lb = lb
            min_t = t
        if lb <= 0:
            raise AssertionError(f"Lower bound failed at t={t}: {lb}")
    if min_lb is None or min_t is None:
        raise AssertionError("Unexpected empty scan")
    print(f"  PASS: min lower bound = {float(min_lb):.12f} at t={min_t}")

    print("\nLarge-t helper inequalities (checked t=6..500):")
    # t*(8/27)^t <= 1/t^3
    bad = 0
    for t in range(6, 501):
        lhs = Fraction(t, 1) * (Fraction(8, 27) ** t)
        rhs = Fraction(1, t**3)
        if lhs > rhs:
            bad += 1
            break
    print(f"  t*(8/27)^t <= 1/t^3: {'PASS' if bad == 0 else 'FAIL'}")

    # t^4 <= 4^t
    bad = 0
    for t in range(4, 501):
        if t**4 > 4**t:
            bad += 1
            break
    print(f"  t^4 <= 4^t (t>=4): {'PASS' if bad == 0 else 'FAIL'}")


if __name__ == "__main__":
    main()

