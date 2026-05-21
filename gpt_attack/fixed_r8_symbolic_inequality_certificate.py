"""Symbolic inequality checks for the S(2^a,8) asymptotic proof.

For a >= 200, write a = 3t + q and use the stabilized hub-off mode m.  This
script verifies, by shifting t to the threshold and checking positive
coefficients, the three polynomial inequalities used in the proof skeleton:

1. left hub-off margin  (F_m - F_{m-1})/F_m >= 1/(10a);
2. right hub-off margin (F_m - F_{m+1})/F_m >= 1/(10a);
3. hub-off Route-2 reserve Delta^0_{a,8} - 1/6 >= 1/(4a).
4. hub-off fugacity lambda_0 stays in [3/4, 2].

The script is intentionally narrow: it certifies these symbolic inequalities,
not the whole Route-2 proof.
"""

from __future__ import annotations

import sympy as sp


t = sp.symbols("t", integer=True, positive=True)
u = sp.symbols("u", integer=True, nonnegative=True)
x = sp.symbols("x")

P8 = [1, 8, 21, 20, 5]
SHIFTS = {0: 9, 1: 7, 2: 8}
THRESHOLDS = {0: 67, 1: 67, 2: 66}


def off_coeff(a, k):
    return sp.combsimp(
        sum(pj * 2 ** (k - j) * sp.binomial(a, k - j) for j, pj in enumerate(P8))
    )


def shifted_num_has_positive_coefficients(expr, threshold: int) -> bool:
    num, den = sp.fraction(sp.factor(sp.combsimp(expr)))
    shifted_num = sp.Poly(sp.expand(num.subs(t, u + threshold)), u)
    shifted_den = sp.Poly(sp.expand(den.subs(t, u + threshold)), u)
    return all(c > 0 for c in shifted_num.coeffs()) and all(
        c > 0 for c in shifted_den.coeffs()
    )


def main() -> None:
    path_poly = sum(c * x**i for i, c in enumerate(P8))
    all_ok = True

    for residue, shift in SHIFTS.items():
        threshold = THRESHOLDS[residue]
        a = 3 * t + residue
        m = sp.simplify((2 * a + shift) / 3)
        fm = off_coeff(a, m)
        left_margin = sp.combsimp((fm - off_coeff(a, m - 1)) / fm)
        right_margin = sp.combsimp((fm - off_coeff(a, m + 1)) / fm)
        lam0 = sp.combsimp(off_coeff(a, m - 1) / fm)
        b = a - 1
        huboff_mean = (
            2 * b * lam0 / (1 + 2 * lam0)
            + lam0 * sp.diff(path_poly, x).subs(x, lam0) / path_poly.subs(x, lam0)
        )
        reserve = sp.combsimp(huboff_mean - (m - sp.Rational(3, 2)) - sp.Rational(1, 6))

        checks = {
            "left_margin_ge_1_over_10a": left_margin - sp.Rational(1, 10) / a,
            "right_margin_ge_1_over_10a": right_margin - sp.Rational(1, 10) / a,
            "reserve_ge_1_over_4a": reserve - sp.Rational(1, 4) / a,
            "lambda0_ge_3_over_4": lam0 - sp.Rational(3, 4),
            "lambda0_ge_1_over_2": lam0 - sp.Rational(1, 2),
            "lambda0_le_2": sp.Rational(2) - lam0,
        }

        print(f"\nresidue q={residue}, threshold t>={threshold} (a>=200 slice)")
        for name, expr in checks.items():
            ok = shifted_num_has_positive_coefficients(expr, threshold)
            all_ok = all_ok and ok
            print(f"  {name}: {ok}")

    assert all_ok


if __name__ == "__main__":
    main()
