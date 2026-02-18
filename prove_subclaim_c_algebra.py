#!/usr/bin/env python3
"""Algebraic proof attempt for Sub-claim C on mixed spiders.

Sub-claim C (from notes/mixed_spider_combined_2026-02-18.md):
  For k >= 3 with k == 0 or 2 (mod 3),
    margin(k,1) < margin(k,0),
  where margin(k,j) is the tie-fugacity margin of S(2^k,1^j).

We decompose
  margin = A + E,  E = r(B-A)/(1+r).

For fixed k (residue 0 or 2), define
  A_gap = A(k,0) - A(k,1).
Then
  margin(k,1) - margin(k,0)
  = -(A_gap) + (E1 - E0)
  <= -(A_gap) - E0,  provided E1 <= 0.

So a sufficient algebraic lane is:
  1) E1 <= 0,
  2) A_gap > -E0.

This script checks that lane exactly (Fraction arithmetic), alongside the direct
inequality itself and a coarse main-term bound A_gap >= 1/(4k).

The lambda formulas are closed forms specialized to residues k mod 3 in
Sub-claim C (k == 0 or 2 mod 3), using p = 4^t.
"""

from __future__ import annotations

import argparse
import math
from fractions import Fraction


def mode_j0(k: int) -> int:
    return (2 * k + 1) // 3


def mode_j1(k: int) -> int:
    return (2 * k + 3) // 3


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


def lambda_j0_closed(k: int) -> Fraction:
    if k % 3 == 0:
        # k = 3t, m0 = 2t.
        t = k // 3
        p = 4**t
        num = Fraction(p, 2) + Fraction(2 * t - 1, t + 2)
        den = Fraction(1, 1) + Fraction(p * (t + 1), 2 * t)
        return num / den
    if k % 3 == 2:
        # k = 3t+2, m0 = 2t+1.
        t = (k - 2) // 3
        p = 4**t
        num = Fraction(p, 1) + Fraction(2 * t, t + 3)
        den = Fraction(1, 1) + Fraction(2 * p * (t + 2), 2 * t + 1)
        return num / den
    raise ValueError("Sub-claim C only applies to k % 3 in {0,2}")


def lambda_j1_closed(k: int) -> Fraction:
    if k % 3 == 0:
        # k = 3t, m1 = 2t+1.
        t = k // 3
        p = 4**t
        num = (2 * t + 1) * (p * (2 * t + 1) + 2 * t)
        den = (t + 1) * (p * (4 * t + 1) + 2 * t + 1)
        return Fraction(num, den)
    if k % 3 == 2:
        # k = 3t+2, m1 = 2t+2.
        t = (k - 2) // 3
        p = 4**t
        num = p * (4 * t + 5) + (2 * t + 1)
        den = (t + 2) * (4 * p + 1)
        return Fraction(num, den)
    raise ValueError("Sub-claim C only applies to k % 3 in {0,2}")


def lambda_from_coeffs_j0(k: int) -> Fraction:
    m = mode_j0(k)
    return Fraction(coeff_j0(k, m - 1), coeff_j0(k, m))


def lambda_from_coeffs_j1(k: int) -> Fraction:
    m = mode_j1(k)
    return Fraction(coeff_j1(k, m - 1), coeff_j1(k, m))


def r_term(k: int, j: int, lam: Fraction) -> Fraction:
    # r = lambda * (1+lambda)^(k-j) / (1+2lambda)^k
    return (
        lam
        * (Fraction(1, 1) + lam) ** (k - j)
        / (Fraction(1, 1) + 2 * lam) ** k
    )


def terms_j0(k: int, lam: Fraction) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    m = mode_j0(k)
    a = Fraction(2 * k) * lam / (Fraction(1, 1) + 2 * lam) - Fraction(m - 1)
    b = Fraction(1, 1) + Fraction(k) * lam / (Fraction(1, 1) + lam) - Fraction(m - 1)
    r = r_term(k, 0, lam)
    margin = a + r * (b - a) / (Fraction(1, 1) + r)
    return a, b, r, margin


def terms_j1(k: int, lam: Fraction) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    m = mode_j1(k)
    a = (
        Fraction(2 * k) * lam / (Fraction(1, 1) + 2 * lam)
        + lam / (Fraction(1, 1) + lam)
        - Fraction(m - 1)
    )
    b = Fraction(1, 1) + Fraction(k) * lam / (Fraction(1, 1) + lam) - Fraction(m - 1)
    r = r_term(k, 1, lam)
    margin = a + r * (b - a) / (Fraction(1, 1) + r)
    return a, b, r, margin


def frac_to_float(x: Fraction) -> float:
    return x.numerator / x.denominator


def main() -> None:
    ap = argparse.ArgumentParser(description="Algebraic proof attempt for mixed-spider Sub-claim C.")
    ap.add_argument("--k-min", type=int, default=3)
    ap.add_argument("--k-max", type=int, default=2000)
    ap.add_argument(
        "--validate-formulas",
        type=int,
        default=120,
        help="Validate closed lambda formulas against coefficient ratios up to this k.",
    )
    args = ap.parse_args()

    if args.k_min < 3:
        raise ValueError("k-min must be >= 3")
    if args.k_max < args.k_min:
        raise ValueError("k-max must be >= k-min")

    formula_fails: list[tuple[int, str]] = []
    if args.validate_formulas > 0:
        for k in range(3, min(args.validate_formulas, args.k_max) + 1):
            if k % 3 == 1:
                continue
            lam0_c = lambda_j0_closed(k)
            lam1_c = lambda_j1_closed(k)
            lam0_e = lambda_from_coeffs_j0(k)
            lam1_e = lambda_from_coeffs_j1(k)
            if lam0_c != lam0_e:
                formula_fails.append((k, "j0"))
            if lam1_c != lam1_e:
                formula_fails.append((k, "j1"))

    bad_subclaim_c: list[dict[str, object]] = []
    bad_sign_b0: list[int] = []
    bad_sign_b1: list[int] = []
    bad_a_gap: list[dict[str, object]] = []
    bad_e1_sign: list[int] = []
    bad_a_beats_minus_e0: list[dict[str, object]] = []
    bad_minus_e0_quarterk: list[dict[str, object]] = []

    finite_witness: dict[int, Fraction] = {}

    total_checked = 0
    for k in range(args.k_min, args.k_max + 1):
        if k % 3 == 1:
            continue
        total_checked += 1

        lam0 = lambda_j0_closed(k)
        lam1 = lambda_j1_closed(k)
        a0, b0, _, margin0 = terms_j0(k, lam0)
        a1, b1, _, margin1 = terms_j1(k, lam1)

        diff = margin1 - margin0
        if diff >= 0:
            bad_subclaim_c.append(
                {
                    "k": k,
                    "margin_j0": margin0,
                    "margin_j1": margin1,
                    "diff": diff,
                }
            )

        if k in (3, 5, 6, 8):
            finite_witness[k] = diff

        if k >= 6 and (b0 - a0) > 0:
            bad_sign_b0.append(k)
        if k >= 5 and (b1 - a1) > 0:
            bad_sign_b1.append(k)

        a_gap = a0 - a1
        if k >= 6 and a_gap <= Fraction(1, 4 * k):
            bad_a_gap.append({"k": k, "a_gap": a_gap, "bound": Fraction(1, 4 * k)})

        e0 = margin0 - a0
        e1 = margin1 - a1
        if k >= 5 and e1 > 0:
            bad_e1_sign.append(k)
        if k >= 6 and (-e0) >= Fraction(1, 4 * k):
            bad_minus_e0_quarterk.append(
                {
                    "k": k,
                    "minus_e0": -e0,
                    "bound": Fraction(1, 4 * k),
                }
            )
        if k >= 6 and not (a_gap > -e0):
            bad_a_beats_minus_e0.append(
                {
                    "k": k,
                    "a_gap": a_gap,
                    "minus_e0": -e0,
                    "e1": e1,
                    "margin_diff": diff,
                }
            )

    print("Sub-claim C algebraic proof attempt")
    print(f"k range checked: {args.k_min}..{args.k_max} (residues 0,2 only)")
    print(f"total k checked: {total_checked}")
    print()

    print("Closed-form lambda validation")
    print(f"failures: {len(formula_fails)}")
    if formula_fails:
        print(f"first failures: {formula_fails[:8]}")
    print()

    print("Direct target inequality")
    print(f"margin(k,1) < margin(k,0) failures: {len(bad_subclaim_c)}")
    if bad_subclaim_c:
        first = bad_subclaim_c[0]
        print(f"first failure: k={first['k']}, diff={first['diff']}")
    print()

    print("Finite base witnesses (exact diff = margin(k,1)-margin(k,0))")
    for k in sorted(finite_witness):
        d = finite_witness[k]
        print(f"k={k:2d}: diff={d} ({frac_to_float(d):+.12f})")
    print()

    print("Proof-template diagnostics")
    print(f"B0-A0 <= 0 for k>=6 failures: {len(bad_sign_b0)}")
    print(f"B1-A1 <= 0 for k>=5 failures: {len(bad_sign_b1)}")
    print(f"A-gap >= 1/(4k) for k>=6 failures: {len(bad_a_gap)}")
    print(f"E1 <= 0 for k>=5 failures: {len(bad_e1_sign)}")
    print(f"-E0 < 1/(4k) for k>=6 failures: {len(bad_minus_e0_quarterk)}")
    print(f"A-gap > -E0 for k>=6 failures: {len(bad_a_beats_minus_e0)}")


if __name__ == "__main__":
    main()
