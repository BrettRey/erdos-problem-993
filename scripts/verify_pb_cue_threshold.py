#!/usr/bin/env python3
"""Certify the CUE threshold in Example 1.6 with exact rational arithmetic.

V_N = N/4 - (2/pi^2) S_N with S_N = sum_{d odd, 1<=d<N} (N-d)/d^2 (rational).
Rational bounds pi_lo < pi < pi_hi come from Machin's formula
pi = 16 arctan(1/5) - 4 arctan(1/239) and the alternating arctangent series.
For N > 4:  V_N < 1  iff  (N/4 - 1) pi^2 < 2 S_N,  certified by pi_hi;
            V_N > 1  iff  (N/4 - 1) pi^2 > 2 S_N,  certified by pi_lo.
Also prints rigorous enclosures of V_N, 1/(4 V_N) and the ULC bound at
N = 10000. Standard library only. Exit status 0 means every check passed.
"""
from __future__ import annotations

from fractions import Fraction as F
import sys


def arctan_bounds(inv_x: int, terms: int) -> tuple[F, F]:
    """Bounds on arctan(1/inv_x) from the alternating series (terms >= 2)."""
    x = F(1, inv_x)
    s = F(0)
    for k in range(terms):
        s += F((-1) ** k, 2 * k + 1) * x ** (2 * k + 1)
    nxt = F((-1) ** terms, 2 * terms + 1) * x ** (2 * terms + 1)
    return (s, s + nxt) if nxt > 0 else (s + nxt, s)


def pi_bounds(terms: int = 30) -> tuple[F, F]:
    a_lo, a_hi = arctan_bounds(5, terms)
    b_lo, b_hi = arctan_bounds(239, terms)
    return 16 * a_lo - 4 * b_hi, 16 * a_hi - 4 * b_lo


def s_n(n: int) -> F:
    return sum((F(n - d, d * d) for d in range(1, n, 2)), F(0))


def v_bounds(n: int, pi_lo: F, pi_hi: F) -> tuple[F, F]:
    s = s_n(n)
    return F(n, 4) - 2 * s / pi_lo**2, F(n, 4) - 2 * s / pi_hi**2


def main() -> int:
    pi_lo, pi_hi = pi_bounds()
    assert pi_lo < pi_hi and float(pi_hi - pi_lo) < 1e-30
    ok = True
    for n in (1996, 1998):
        lo, hi = v_bounds(n, pi_lo, pi_hi)
        print(f"V_{n} in [{float(lo):.15f}, {float(hi):.15f}]")
    lo96, hi96 = v_bounds(1996, pi_lo, pi_hi)
    lo98, hi98 = v_bounds(1998, pi_lo, pi_hi)
    cert = hi96 < 1 < lo98
    print("certified V_1996 < 1 < V_1998:", cert)
    ok &= cert
    lo, hi = v_bounds(10000, pi_lo, pi_hi)
    print(f"V_10000 in [{float(lo):.12f}, {float(hi):.12f}]")
    print(f"1/(4 V_10000) in [{float(1 / (4 * hi)):.12f}, {float(1 / (4 * lo)):.12f}]")
    ulc = F(10001, (5000 + 2) * 5000)
    print(f"ULC bound at D, N = 10000: {float(ulc):.12f}")
    ok &= 0.2149 < float(1 / (4 * hi)) and float(1 / (4 * lo)) < 0.2150 and abs(float(ulc) - 0.0004) < 1e-6
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
