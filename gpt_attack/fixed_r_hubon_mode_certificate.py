"""Hub-on mode-domination certificate for fixed-r spider lanes.

For S(2^a,r), write

    I_{a,r}(x) = F_a(x) + G_a(x)
               = P_r(x)(1+2x)^a + x P_{r-1}(x)(1+x)^a.

Given a hub-off margin certificate

    F_m - F_k >= F_m/(D a)        for k != m,

this script certifies the simple sufficient condition

    2 max_k G_k / F_m < 1/(D a).

The bound is deliberately elementary.  Since

    G_k <= P_{r-1}(1) 2^a

and any chosen term p_s x^s in P_r gives

    F_m >= p_s 2^(m-s) binom(a,m-s),

we prove the resulting upper bound at the threshold and prove it decreases
when a advances by 3 in the same residue class.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import ceil, comb

import sympy as sp

from fixed_r_huboff_certificate import stabilized_shifts
from route2_spider_lane_scan import path_polys


t = sp.symbols("t", integer=True, positive=True)
u = sp.symbols("u", integer=True, nonnegative=True)


def shifted_positive(expr, threshold_t: int) -> bool:
    shifted = sp.Poly(sp.expand(sp.together(expr).subs(t, u + threshold_t)), u)
    return all(c > 0 for c in shifted.coeffs())


def best_witness_term(
    path_poly: list[int],
    hub_on_sum: int,
    a: int,
    m: int,
) -> tuple[int, Fraction]:
    best_s = -1
    best_bound: Fraction | None = None
    for s, coeff in enumerate(path_poly):
        n = m - s
        if coeff <= 0 or not (0 <= n <= a):
            continue
        bound = Fraction(hub_on_sum * 2 ** (a - n), coeff * comb(a, n))
        if best_bound is None or bound < best_bound:
            best_s = s
            best_bound = bound
    if best_bound is None:
        raise ValueError("no valid witness term")
    return best_s, best_bound


def check_r(r: int, threshold: int, margin_denom: int) -> bool:
    paths = path_polys(max(2, r))
    pr = paths[r]
    prm1 = paths[r - 1]
    hub_on_sum = sum(prm1)
    shifts = stabilized_shifts(r, threshold, paths)
    print(f"r={r}, threshold a>={threshold}, margin denominator D={margin_denom}")
    print(f"stabilized shifts D_q=3m-2a: {shifts}")

    all_ok = True
    for residue, shift in shifts.items():
        a0 = threshold
        while a0 % 3 != residue:
            a0 += 1
        threshold_t = ceil((threshold - residue) / 3)
        m0 = (2 * a0 + shift) // 3
        s, tail_bound = best_witness_term(pr, hub_on_sum, a0, m0)
        mode_perturb = 2 * tail_bound
        target = Fraction(1, margin_denom * a0)
        threshold_ok = mode_perturb < target

        a = 3 * t + residue
        m = sp.simplify((2 * a + shift) / 3)
        n = sp.simplify(m - s)
        # H(a) = a * P_{r-1}(1) 2^(a-n) / (p_s binom(a,n)).
        # For a -> a+3 in the same residue class, n -> n+2.
        # Certify H(a+3)/H(a) < 1 after shifting t to the threshold.
        step_margin = sp.combsimp(
            a * (a + 1) * (a + 2)
            - 2 * (n + 1) * (n + 2) * (a - n + 1)
        )
        decreasing_ok = shifted_positive(step_margin, threshold_t)
        all_ok = all_ok and threshold_ok and decreasing_ok

        print(
            f"  q={residue}: witness s={s}, n0={m0 - s}, "
            f"2G/F <= {float(mode_perturb):.3e}, "
            f"target={float(target):.3e}, "
            f"threshold_ok={threshold_ok}, decreasing_ok={decreasing_ok}"
        )

    return all_ok


def main() -> None:
    ap = argparse.ArgumentParser(description="Fixed-r hub-on mode domination certificate.")
    ap.add_argument("--r", type=int, required=True)
    ap.add_argument("--threshold", type=int, default=200)
    ap.add_argument("--margin-denom", type=int, default=1000)
    args = ap.parse_args()

    ok = check_r(args.r, args.threshold, args.margin_denom)
    assert ok


if __name__ == "__main__":
    main()
