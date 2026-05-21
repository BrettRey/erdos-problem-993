"""Symbolic fixed-r certificate data for the lane S(2^a,8).

This does not claim the full Route-2 proof.  It extracts the pieces needed for
that proof:

1. boundary difference signs for the proposed full-polynomial mode shifts;
2. positivity of the hub-off Route-2 reserve Delta^0_{a,8} - 1/6.

The remaining proof obligation is to turn these boundary/asymptotic facts into
a full mode-and-perturbation argument for the complete polynomial

    P_8(x)(1+2x)^a + x P_7(x)(1+x)^a.
"""

from __future__ import annotations

import sympy as sp


t = sp.symbols("t", integer=True, positive=True)
x = sp.symbols("x")

P8 = [1, 8, 21, 20, 5]
P7 = [1, 7, 15, 10, 1]
SHIFTS = {0: 9, 1: 7, 2: 8}


def off_coeff(a, k):
    return sp.combsimp(
        sum(pj * 2 ** (k - j) * sp.binomial(a, k - j) for j, pj in enumerate(P8))
    )


def full_coeff(a, k):
    off = sum(pj * 2 ** (k - j) * sp.binomial(a, k - j) for j, pj in enumerate(P8))
    on = sum(qj * sp.binomial(a, k - 1 - j) for j, qj in enumerate(P7))
    return sp.combsimp(off + on)


def positive_coefficients(expr) -> bool:
    num, den = sp.fraction(sp.factor(expr))
    return all(c > 0 for c in sp.Poly(num, t).coeffs()) and all(
        c > 0 for c in sp.Poly(den, t).coeffs()
    )


def main() -> None:
    path_poly = sum(c * x**i for i, c in enumerate(P8))

    for residue, shift in SHIFTS.items():
        a = 3 * t + residue
        m = (2 * a + shift) / 3
        m = sp.simplify(m)
        print(f"\nresidue q={residue}, a=3t+{residue}, m={m}, shift={shift}")

        inc = sp.factor(sp.combsimp(full_coeff(a, m) - full_coeff(a, m - 1)))
        dec = sp.factor(sp.combsimp(full_coeff(a, m + 1) - full_coeff(a, m)))
        print("full boundary c_m - c_{m-1}:")
        print(inc)
        print("full boundary c_{m+1} - c_m:")
        print(dec)

        lam0 = sp.factor(sp.combsimp(off_coeff(a, m - 1) / off_coeff(a, m)))
        b = a - 1
        huboff_mean = (
            2 * b * lam0 / (1 + 2 * lam0)
            + lam0 * sp.diff(path_poly, x).subs(x, lam0) / path_poly.subs(x, lam0)
        )
        reserve = sp.factor(sp.combsimp(huboff_mean - (m - sp.Rational(3, 2)) - sp.Rational(1, 6)))
        print("hub-off lambda0:")
        print(lam0)
        print("hub-off reserve Delta0 - 1/6:")
        print(reserve)
        print(f"reserve numerator/denominator have positive coefficients: {positive_coefficients(reserve)}")


if __name__ == "__main__":
    main()
