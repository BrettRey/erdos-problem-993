#!/usr/bin/env python3
"""Exact check of the range H >= 16 in Proposition 3.1 (scalar inequality).

Companion to verify_universal_pb_finite_bernstein.py, which certifies the
compact range 3 < H <= 16. This script reproduces every symbolic step of the
manuscript's subsection "The range H >= 16":

  1. the closed forms of tilde S_J and tilde T_J from the lower bounds
     lambda_r = 1 - r(r+1)/(2H);
  2. the expansion of N_J(H) = H^2 (tilde S_J tilde T_J - Q(H)) as an
     explicit quartic in H;
  3. the degree-four Bernstein coefficients of N_5 on 16 <= H <= 21;
  4. for J = u + 6 >= 6, the identity
         beta_i(J) = mu_i(J) pi_i(u) / 2880,   i = 0, ..., 4,
     between the Bernstein coefficients in t of
     N_J(J(J+1)/2 + (J+1) t) and the multipliers mu_i and polynomials
     pi_i printed in the manuscript, as a polynomial identity in u;
  5. strict positivity of every coefficient of every pi_i.

All arithmetic is exact (SymPy rationals). Exit status 0 means every check
passed. Requires SymPy (reference version 1.14.0).
"""
from __future__ import annotations

import sys

import sympy as sp

H, t, u, J = sp.symbols("H t u J")


def bernstein(poly_t: sp.Expr, d: int) -> list[sp.Expr]:
    """Degree-d Bernstein coefficients on [0, 1] of a polynomial in t."""
    a = sp.Poly(sp.expand(poly_t), t).all_coeffs()[::-1]
    a = a + [0] * (d + 1 - len(a))
    return [
        sp.simplify(sum(a[j] * sp.binomial(i, j) / sp.binomial(d, j) for j in range(i + 1)))
        for i in range(d + 1)
    ]


def main() -> int:
    ok = True

    def check(label: str, cond: bool) -> None:
        nonlocal ok
        print(("PASS " if cond else "FAIL ") + label)
        ok = ok and bool(cond)

    Q = (3 * H + 4) * (H + 1) / 4
    delta = 1 / (H + 1)
    check("Q(H) equals (3+delta)/(4 delta^2) at delta = 1/(H+1)",
          sp.simplify(Q - (3 + delta) / (4 * delta**2)) == 0)

    # 1. Closed forms, verified for J = 1..40 against direct sums.
    sig3 = J**2 * (J + 1) ** 2 / 4
    sig4 = J * (J + 1) * (2 * J + 1) * (3 * J**2 + 3 * J - 1) / 30
    S_closed = 2 * J + 1 - J * (J + 1) * (J + 2) / (3 * H)
    T_closed = J * (J + 1) * (2 * J + 1) / 3 - (sig3 + sig4) / H
    closed_ok = True
    for j in range(1, 41):
        lam = [1 - sp.Rational(r * (r + 1), 2) / H for r in range(1, j + 1)]
        S_dir = 1 + 2 * sum(lam)
        T_dir = 2 * sum((r + 1 - 1) ** 2 * lam[r - 1] for r in range(1, j + 1))
        closed_ok &= sp.simplify(S_dir - S_closed.subs(J, j)) == 0
        closed_ok &= sp.simplify(T_dir - T_closed.subs(J, j)) == 0
    check("closed forms of tilde S_J, tilde T_J (J = 1..40)", closed_ok)

    # 2. Explicit quartic N_J(H).
    N = sp.expand(H**2 * (S_closed * T_closed - Q))
    alpha = 2 * J + 1
    beta = J * (J + 1) * (J + 2) / 3
    gamma = J * (J + 1) * (2 * J + 1) / 3
    Sig = sig3 + sig4
    N_explicit = (
        -sp.Rational(3, 4) * H**4
        - sp.Rational(7, 4) * H**3
        + (alpha * gamma - 1) * H**2
        - (alpha * Sig + beta * gamma) * H
        + beta * Sig
    )
    check("N_J(H) equals the explicit quartic", sp.simplify(N - N_explicit) == 0)

    # 3. J = 5 cell, H = 16 + 5t.
    N5 = N.subs(J, 5).subs(H, 16 + 5 * t)
    b5 = bernstein(N5, 4)
    expected5 = [sp.Integer(2360), sp.Integer(7500), sp.Rational(25055, 2),
                 sp.Rational(254205, 16), sp.Rational(31115, 2)]
    check("J = 5 Bernstein coefficients match the manuscript", b5 == expected5)
    check("J = 5 Bernstein coefficients strictly positive", all(c > 0 for c in b5))

    # 4. Uniform family J = u + 6.
    NJ = sp.expand(N.subs(H, J * (J + 1) / 2 + (J + 1) * t))
    bJ = bernstein(NJ, 4)
    mu = [J**2 * (J + 1) ** 2, J * (J + 1) ** 2, (J + 1) ** 2,
          (J + 1) ** 2 * (J + 2), (J + 1) ** 2 * (J + 2) ** 2]
    pi = [
        121 * u**4 + 2474 * u**3 + 17431 * u**2 + 46014 * u + 25400,
        121 * u**5 + 3442 * u**4 + 37967 * u**3 + 199405 * u**2 + 480475 * u + 384510,
        121 * u**6 + 4410 * u**5 + 65863 * u**4 + 512764 * u**3
        + 2172926 * u**2 + 4668316 * u + 3829440,
        121 * u**5 + 3684 * u**4 + 43735 * u**3 + 249961 * u**2 + 671859 * u + 644400,
        121 * u**4 + 2958 * u**3 + 25579 * u**2 + 88782 * u + 91440,
    ]
    for i in range(5):
        lhs = sp.expand(bJ[i].subs(J, u + 6))
        rhs = sp.expand(mu[i].subs(J, u + 6) * pi[i] / 2880)
        check(f"beta_{i}(J) = mu_{i}(J) pi_{i}(u) / 2880 identically in u",
              sp.simplify(lhs - rhs) == 0)

    # 5. Positivity of the pi_i coefficients (hence pi_i(u) > 0 for u >= 0).
    for i, p in enumerate(pi):
        coeffs = sp.Poly(p, u).all_coeffs()
        check(f"pi_{i} has strictly positive coefficients", all(c > 0 for c in coeffs))

    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
