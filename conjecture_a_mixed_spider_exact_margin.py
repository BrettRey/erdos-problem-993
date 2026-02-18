#!/usr/bin/env python3
"""Exact rational margin computation for mixed spiders S(2^k, 1^j).

Tree family: T = S(2^k, 1^j) with
  I_T(x) = (1+2x)^k (1+x)^j + x(1+x)^k

Coefficients:
  i_n(T) = sum_{i=0}^{n} C(k,i)*2^i*C(j,n-i) + C(k,n-1)

Goal: prove margin = mu(T, lambda_m) - (m-1) > 1/3 for all k >= 1, j >= 0.

This script uses exact rational arithmetic (fractions.Fraction) to:
1. Compute mode m exactly (argmax of i_n)
2. Compute lambda_m = Fraction(i_{m-1}, i_m) exactly
3. Compute margin = mu(T, lambda_m) - (m-1) exactly as a Fraction
4. Check monotonicity and recursion structure for induction

Key outputs:
- Table of exact margins for small (k,j)
- Check: margin > 1/3 for all checked cases
- Check: margin(k,0) is monotone decreasing in k
- Tip-leaf recursion: margin(k,j) = alpha*c1 + beta*c2 (exact)
- c1 = mu(S(2^{k-1},1^{j+1}), lambda_m^T) - (m-1): checks if >= 0
- c2 = mu(S(2^{k-1},1^j), lambda_m^T) - (m-2): checks if >= 0
"""

from __future__ import annotations

import argparse
import json
import math
import time
from fractions import Fraction
from functools import lru_cache
from typing import NamedTuple


# ---------------------------------------------------------------------------
# Exact coefficient computation
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def binom(n: int, k: int) -> int:
    """Exact binomial coefficient."""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)


def mixed_spider_coeffs(k: int, j: int, deg: int | None = None) -> list[int]:
    """Exact integer coefficients of I(S(2^k, 1^j); x).

    I = (1+2x)^k * (1+x)^j + x*(1+x)^k
    i_n = sum_{i=0}^{n} C(k,i)*2^i*C(j,n-i) + C(k,n-1)
    """
    max_deg = 1 + 2 * k + j if deg is None else deg  # conservative upper bound
    coeffs = [0] * (max_deg + 1)

    # Part A: (1+2x)^k * (1+x)^j
    for n in range(max_deg + 1):
        s = 0
        for i in range(n + 1):
            c_ki = binom(k, i)
            if c_ki == 0:
                continue
            c_jni = binom(j, n - i)
            if c_jni == 0:
                continue
            s += c_ki * (2 ** i) * c_jni
        coeffs[n] += s

    # Part B: x*(1+x)^k
    for n in range(1, max_deg + 1):
        coeffs[n] += binom(k, n - 1)

    # Strip trailing zeros
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()

    return coeffs


def product_graph_coeffs(k: int, j: int) -> list[int]:
    """Coefficients of (1+2x)^k * (1+x)^j (no hub-included term)."""
    max_deg = k + j
    coeffs = [0] * (max_deg + 1)
    for n in range(max_deg + 1):
        s = 0
        for i in range(n + 1):
            c_ki = binom(k, i)
            if c_ki == 0:
                continue
            c_jni = binom(j, n - i)
            if c_jni == 0:
                continue
            s += c_ki * (2 ** i) * c_jni
        coeffs[n] = s
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


def mode_and_lambda(coeffs: list[int]) -> tuple[int, Fraction]:
    """Return (mode, lambda_m) where lambda_m = i_{m-1} / i_m."""
    m = 0
    best = coeffs[0]
    for i in range(1, len(coeffs)):
        if coeffs[i] > best:
            best = coeffs[i]
            m = i
    if m == 0:
        return 0, Fraction(0)
    lam = Fraction(coeffs[m - 1], coeffs[m])
    return m, lam


def poly_eval(coeffs: list[int], x: Fraction) -> Fraction:
    """Evaluate polynomial with integer coefficients at exact rational x."""
    result = Fraction(0)
    xn = Fraction(1)
    for c in coeffs:
        result += c * xn
        xn *= x
    return result


def poly_mean(coeffs: list[int], x: Fraction) -> Fraction:
    """Compute mu = sum(n * c_n * x^n) / sum(c_n * x^n) at exact rational x."""
    num = Fraction(0)
    den = Fraction(0)
    xn = Fraction(1)
    for n, c in enumerate(coeffs):
        val = c * xn
        den += val
        num += n * val
        xn *= x
    if den == 0:
        raise ZeroDivisionError("I(T; x) = 0 at given x")
    return num / den


# ---------------------------------------------------------------------------
# Margin computation
# ---------------------------------------------------------------------------

class MarginData(NamedTuple):
    k: int
    j: int
    n: int
    mode: int
    lam: Fraction
    margin: Fraction  # mu(T, lam) - (m-1)
    margin_float: float
    above_third: bool  # margin > 1/3


def compute_margin(k: int, j: int) -> MarginData:
    """Compute exact tie-fugacity margin for S(2^k, 1^j)."""
    coeffs = mixed_spider_coeffs(k, j)
    m, lam = mode_and_lambda(coeffs)
    if m == 0:
        return MarginData(k, j, 1 + 2 * k + j, 0, Fraction(0), Fraction(0), 0.0, True)
    mu = poly_mean(coeffs, lam)
    margin = mu - (m - 1)
    return MarginData(
        k=k, j=j, n=1 + 2 * k + j, mode=m,
        lam=lam, margin=margin, margin_float=float(margin),
        above_third=(margin > Fraction(1, 3))
    )


# ---------------------------------------------------------------------------
# Tip-leaf recursion decomposition
# ---------------------------------------------------------------------------

class TipLeafData(NamedTuple):
    k: int
    j: int
    margin: Fraction
    c1: Fraction   # mu(A, lam_T) - (m-1),  A = S(2^{k-1}, 1^{j+1})
    c2: Fraction   # mu(B, lam_T) - (m-2),  B = S(2^{k-1}, 1^j)
    alpha: Fraction  # I(A, lam_T) / I(T, lam_T)
    beta: Fraction   # lam_T * I(B, lam_T) / I(T, lam_T)
    combined: Fraction  # alpha*c1 + beta*c2 (should equal margin)
    identity_err: Fraction  # combined - margin (should be 0)


def tip_leaf_decomp(k: int, j: int) -> TipLeafData | None:
    """Exact tip-leaf decomposition for S(2^k, 1^j) (requires k >= 1)."""
    if k < 1:
        return None

    coeffs_T = mixed_spider_coeffs(k, j)
    m, lam = mode_and_lambda(coeffs_T)
    if m == 0:
        return None

    # A = S(2^{k-1}, 1^{j+1}), B = S(2^{k-1}, 1^j)
    coeffs_A = mixed_spider_coeffs(k - 1, j + 1)
    coeffs_B = mixed_spider_coeffs(k - 1, j)

    iT = poly_eval(coeffs_T, lam)
    iA = poly_eval(coeffs_A, lam)
    iB = poly_eval(coeffs_B, lam)

    alpha = iA / iT
    beta = lam * iB / iT

    mu_A = poly_mean(coeffs_A, lam)
    mu_B = poly_mean(coeffs_B, lam)

    c1 = mu_A - (m - 1)
    c2 = mu_B - (m - 2)

    combined = alpha * c1 + beta * c2
    margin = poly_mean(coeffs_T, lam) - (m - 1)

    return TipLeafData(
        k=k, j=j, margin=margin, c1=c1, c2=c2,
        alpha=alpha, beta=beta,
        combined=combined, identity_err=combined - margin
    )


def unit_leaf_decomp(k: int, j: int) -> dict | None:
    """Exact unit-leaf decomposition for S(2^k, 1^j) (requires j >= 1)."""
    if j < 1:
        return None

    coeffs_T = mixed_spider_coeffs(k, j)
    m, lam = mode_and_lambda(coeffs_T)
    if m == 0:
        return None

    # A = S(2^k, 1^{j-1}), B = (1+2x)^k * (1+x)^{j-1}
    coeffs_A = mixed_spider_coeffs(k, j - 1)
    coeffs_B = product_graph_coeffs(k, j - 1)

    iT = poly_eval(coeffs_T, lam)
    iA = poly_eval(coeffs_A, lam)
    iB = poly_eval(coeffs_B, lam)

    alpha = iA / iT
    beta = lam * iB / iT

    mu_A = poly_mean(coeffs_A, lam)
    mu_B = poly_mean(coeffs_B, lam)

    c1 = mu_A - (m - 1)
    c2 = mu_B - (m - 2)

    combined = alpha * c1 + beta * c2
    margin = poly_mean(coeffs_T, lam) - (m - 1)

    return {
        "k": k, "j": j, "m": m, "lam": lam,
        "margin": margin, "c1": c1, "c2": c2,
        "alpha": alpha, "beta": beta,
        "combined": combined, "identity_err": combined - margin,
    }


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description="Exact margin analysis for S(2^k, 1^j).")
    ap.add_argument("--k-max", type=int, default=40)
    ap.add_argument("--j-max", type=int, default=20)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--check-monotone-k", action="store_true",
                    help="Check if margin(k,0) is monotone decreasing in k.")
    ap.add_argument("--check-tip-recursion", action="store_true",
                    help="Check tip-leaf recursion identity and sign patterns.")
    ap.add_argument("--table", action="store_true",
                    help="Print exact margin table for small (k,j).")
    args = ap.parse_args()

    print(f"Exact rational margin analysis for S(2^k, 1^j), k=1..{args.k_max}, j=0..{args.j_max}")
    t0 = time.time()

    # --- Basic margin scan ---
    total = 0
    below_third = 0
    min_margin_float = float("inf")
    min_margin_witness = None

    for k in range(1, args.k_max + 1):
        for j in range(0, args.j_max + 1):
            d = compute_margin(k, j)
            total += 1
            if not d.above_third:
                below_third += 1
                print(f"  FAIL margin<=1/3: k={k}, j={j}, margin={float(d.margin):.8f}")
            if d.margin_float < min_margin_float:
                min_margin_float = d.margin_float
                min_margin_witness = d

    print(f"\nTotal checked: {total}, margin<=1/3: {below_third}")
    if min_margin_witness is not None:
        d = min_margin_witness
        print(f"Min margin: {float(d.margin):.8f} = {d.margin} at (k={d.k}, j={d.j})")
    print(f"Wall: {time.time()-t0:.2f}s")

    # --- Monotonicity of margin(k,0) ---
    if args.check_monotone_k:
        print("\n--- Monotonicity of margin(k,0) in k ---")
        prev = None
        mono_fails = 0
        for k in range(1, args.k_max + 1):
            d = compute_margin(k, 0)
            if prev is not None and d.margin > prev.margin:
                mono_fails += 1
                print(f"  NON-MONOTONE: k={k}: margin={float(d.margin):.8f} > margin(k-1)={float(prev.margin):.8f}")
            if args.verbose:
                print(f"  k={k:3d}: margin={float(d.margin):.8f} = {d.margin}")
            prev = d
        print(f"Monotone decreasing fails: {mono_fails} (0 = monotone)")

    # --- Tip-leaf recursion check ---
    if args.check_tip_recursion:
        print("\n--- Tip-leaf recursion: margin = alpha*c1 + beta*c2 ---")
        identity_fails = 0
        c1_neg_count = 0
        c2_neg_count = 0
        combined_neg_count = 0

        for k in range(1, min(args.k_max + 1, 31)):
            for j in range(0, min(args.j_max + 1, 21)):
                d = tip_leaf_decomp(k, j)
                if d is None:
                    continue
                if d.identity_err != 0:
                    identity_fails += 1
                    print(f"  IDENTITY FAIL: k={k}, j={j}, err={d.identity_err}")
                if d.c1 < 0:
                    c1_neg_count += 1
                    if args.verbose:
                        print(f"  c1<0: k={k}, j={j}, c1={float(d.c1):.6f}, c2={float(d.c2):.6f}, "
                              f"combined={float(d.combined):.6f}")
                if d.c2 < 0:
                    c2_neg_count += 1
                    print(f"  c2<0: k={k}, j={j}, c2={float(d.c2):.6f}")
                if d.combined <= Fraction(0):
                    combined_neg_count += 1
                    print(f"  COMBINED<=0: k={k}, j={j}, combined={float(d.combined):.6f}")

        print(f"Identity failures: {identity_fails}")
        print(f"c1 < 0: {c1_neg_count}")
        print(f"c2 < 0: {c2_neg_count}")
        print(f"combined <= 0: {combined_neg_count}")

    # --- Margin table ---
    if args.table:
        print("\n--- Exact margin table (float) ---")
        print(f"{'j\\k':>4}", end="")
        for k in range(1, min(args.k_max + 1, 11)):
            print(f"  k={k:3d} ", end="")
        print()
        for j in range(0, min(args.j_max + 1, 11)):
            print(f"j={j:2d} ", end="")
            for k in range(1, min(args.k_max + 1, 11)):
                d = compute_margin(k, j)
                print(f"  {float(d.margin):.4f} ", end="")
            print()

    # --- Pure spider: exact margin formula check ---
    print("\n--- Pure spider S(2^k): margin and 1/3-gap ---")
    print(f"{'k':>5} {'mode':>5} {'margin':>12} {'margin-1/3':>14} {'m*margin':>12}")
    for k in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40]:
        if k > args.k_max:
            continue
        d = compute_margin(k, 0)
        gap = d.margin - Fraction(1, 3)
        print(f"{k:5d} {d.mode:5d} {float(d.margin):12.8f} {float(gap):14.8f} {float(d.mode * d.margin):12.8f}")

    # --- Mixed spider: c2/c1 ratio to understand anti-correlation ---
    print("\n--- Tip-leaf: c2 / |c1| when c1 < 0 ---")
    min_ratio = float("inf")
    min_ratio_witness = None
    for k in range(1, min(args.k_max + 1, 31)):
        for j in range(0, min(args.j_max + 1, 21)):
            d = tip_leaf_decomp(k, j)
            if d is None or d.c1 >= 0:
                continue
            ratio = float(d.c2 / (-d.c1))
            if ratio < min_ratio:
                min_ratio = ratio
                min_ratio_witness = (k, j, d)
    if min_ratio_witness is not None:
        k, j, d = min_ratio_witness
        print(f"Min c2/|c1| = {min_ratio:.4f} at (k={k}, j={j}), "
              f"c1={float(d.c1):.6f}, c2={float(d.c2):.6f}, alpha={float(d.alpha):.4f}, beta={float(d.beta):.4f}")
        print(f"  Need: alpha*|c1| <= beta*c2, i.e., (alpha/beta) <= c2/|c1|")
        print(f"  alpha/beta = {float(d.alpha/d.beta):.4f}")
    else:
        print("No c1<0 cases in tip-leaf scan.")

    print("\nDone.")


if __name__ == "__main__":
    main()
