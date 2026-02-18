#!/usr/bin/env python3
"""LC-defect approach to Sub-claim A: prove T1 >= |T2| for the F>=0 condition.

For mixed spiders S(2^k,1^j), the lambda-monotonicity condition F>=0 decomposes as:
  F = T1 + T2 >= 0
where:
  A=a_m, B=a_{m-1}, C=a_{m-2}, D=a_{m+1} (coefficients of A_j=(1+2x)^k(1+x)^j at mode m)
  T1 = (A^2 - BD) + (AC - B^2)          [always positive from data]
  T2 = G*a'_m - H*a'_{m+1}              [always negative from data]
    with G=C(k,m-1), H=C(k,m-2), a'_t = [x^t](1+2x)^k(1+x)^{j+2}

This script:
1. Verifies the LC-defect ratio Δ(m)/Δ(m-1) > 1 at the mode (j=0 case, closed form).
2. Checks the ratio T1_j / T1_0 and T2_j / T1_j for general j.
3. Tries to bound |T2| / T1 away from 1 to prove F > 0.

Key algebraic fact (j=0): Δ(t) = 4^t * C(k,t)^2 * (k+1) / [(k-t+1)(t+1)]
  So Δ(m)/Δ(m-1) = 4*(k-m+1)*(k-m+2) / [m*(m+1)].
  At m = round(2k/3): ratio = 4*(k/3)*(k/3+1) / [(2k/3)*(2k/3+1)] ~ 1 + 15/(2k) > 1.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from fractions import Fraction
from typing import Any


def binom(n: int, r: int) -> int:
    if r < 0 or r > n:
        return 0
    if r == 0 or r == n:
        return 1
    r = min(r, n - r)
    result = 1
    for i in range(r):
        result = result * (n - i) // (i + 1)
    return result


def build_aj_coeffs(k: int, j: int) -> list[int]:
    """Coefficients of A_j(x) = (1+2x)^k * (1+x)^j."""
    n = k + j + 1
    a = [0] * n
    for ell in range(k + 1):
        c_k_ell = binom(k, ell) * (2**ell)
        for s in range(j + 1):
            c_j_s = binom(j, s)
            t = ell + s
            if t < n:
                a[t] += c_k_ell * c_j_s
    # Trim trailing zeros
    while a and a[-1] == 0:
        a.pop()
    return a


def build_g_coeffs(k: int, deg: int) -> list[int]:
    """Coefficients of E(x) = x*(1+x)^k, i.e. g_t = C(k,t-1)."""
    g = [0] * (deg + 1)
    for t in range(1, min(k + 2, deg + 1)):
        g[t] = binom(k, t - 1)
    return g


def find_mode_and_lambda(a: list[int], g: list[int]) -> tuple[int, Fraction]:
    """Mode and tie-fugacity of A+G (exact Fraction)."""
    deg = max(len(a), len(g)) - 1

    def coeff(t: int) -> int:
        v = (a[t] if t < len(a) else 0) + (g[t] if t < len(g) else 0)
        return v

    m = 0
    best = coeff(0)
    for t in range(1, deg + 1):
        v = coeff(t)
        if v > best:
            best = v
            m = t

    if m == 0:
        return 0, Fraction(0)
    im1 = coeff(m - 1)
    im = coeff(m)
    return m, Fraction(im1, im)


def compute_F_components(k: int, j: int) -> dict[str, Any]:
    """Compute T1, T2, F for (k,j) using exact integer arithmetic.

    Returns dict with keys: k, j, m, T1, T2, F, ratio_T2_T1, ok.
    """
    a_j = build_aj_coeffs(k, j)
    a_j2 = build_aj_coeffs(k, j + 2)  # A_{j+2}
    g = build_g_coeffs(k, max(len(a_j), len(a_j2)) + 2)

    m, lam = find_mode_and_lambda(a_j, g)

    if m < 2 or m + 1 >= len(a_j):
        return {"k": k, "j": j, "m": m, "ok": False, "skip": "mode too small/large"}

    # Coefficients of A_j at positions m-2..m+1
    A = a_j[m]
    B = a_j[m - 1]
    C = a_j[m - 2]
    D = a_j[m + 1] if m + 1 < len(a_j) else 0

    # LC defects
    lc_m = A * A - B * D    # >= 0 if LC at m
    lc_m1 = A * C - B * B  # <= 0 if LC at m-1 (a_{m-1}^2 >= a_{m-2}*a_m)

    T1 = lc_m + lc_m1

    # A_{j+2} coefficients at m and m+1
    Ap_m = a_j2[m] if m < len(a_j2) else 0
    Ap_m1 = a_j2[m + 1] if m + 1 < len(a_j2) else 0

    G = binom(k, m - 1)   # C(k, m-1)
    H = binom(k, m - 2)   # C(k, m-2)

    T2 = G * Ap_m - H * Ap_m1

    F = T1 + T2

    ratio = Fraction(T2, T1) if T1 != 0 else None

    # Also compute Δ(m)/Δ(m-1) for the A_j sequence LC defect
    # Δ(t) for A_j is harder; for j=0 it's the closed form.
    # For general j, compute it explicitly.
    # delta_t = a[t]^2 - a[t-1]*a[t+1]
    delta_m = A * A - B * D        # = lc_m >= 0
    delta_m_minus1 = B * B - C * A  # a_{m-1}^2 - a_{m-2}*a_m

    return {
        "k": k,
        "j": j,
        "m": m,
        "lc_m": lc_m,
        "lc_m1": lc_m1,
        "T1": T1,
        "T2": T2,
        "F": F,
        "ratio_T2_T1": float(ratio) if ratio is not None else None,
        "delta_m": delta_m,
        "delta_m1": delta_m_minus1,
        "delta_ratio": float(Fraction(delta_m, delta_m_minus1)) if delta_m_minus1 != 0 else None,
        "ok": F >= 0,
    }


def delta_ratio_j0_closed(k: int, m: int) -> Fraction:
    """Closed-form Δ(m)/Δ(m-1) for A_0 = (1+2x)^k.

    Δ(t) = 4^t * C(k,t)^2 * (k+1) / [(k-t+1)*(t+1)]
    Δ(m)/Δ(m-1) = 4*(k-m+1)*(k-m+2) / [m*(m+1)].
    """
    return Fraction(4 * (k - m + 1) * (k - m + 2), m * (m + 1))


def verify_j0_closed_form(k_max: int) -> dict[str, Any]:
    """Verify that Δ(m)/Δ(m-1) > 1 at mode for j=0 case."""
    failures = []
    min_excess = None
    min_rec = None

    for k in range(6, k_max + 1):
        a0 = build_aj_coeffs(k, 0)
        g = build_g_coeffs(k, len(a0) + 2)
        m, _ = find_mode_and_lambda(a0, g)

        if m < 2 or m + 1 >= len(a0):
            continue

        ratio = delta_ratio_j0_closed(k, m)

        # Also verify the closed form matches the explicit computation
        A = a0[m]; B = a0[m - 1]; C = a0[m - 2]
        D = a0[m + 1] if m + 1 < len(a0) else 0
        delta_m_explicit = Fraction(A * A - B * D)
        delta_m1_explicit = Fraction(B * B - C * A)

        ratio_explicit = delta_m_explicit / delta_m1_explicit if delta_m1_explicit != 0 else None

        excess = ratio - 1
        rec = {"k": k, "m": m, "delta_ratio": float(ratio), "excess": float(excess)}

        if ratio <= 1:
            failures.append(rec)

        if min_excess is None or excess < min_excess:
            min_excess = excess
            min_rec = rec

    return {
        "k_max": k_max,
        "failures": failures,
        "min_excess": float(min_excess) if min_excess is not None else None,
        "min_rec": min_rec,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k-max", type=int, default=200)
    ap.add_argument("--j-max", type=int, default=40)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    t0 = time.time()
    print(f"LC-defect approach: k=6..{args.k_max}, j=0..{args.j_max}", flush=True)

    # Part 1: Verify closed-form Δ(m)/Δ(m-1) > 1 at mode for j=0
    print("\n=== Part 1: j=0 closed-form Δ(m)/Δ(m-1) > 1 ===", flush=True)
    r0 = verify_j0_closed_form(args.k_max)
    print(f"j=0 failures: {len(r0['failures'])}", flush=True)
    print(f"j=0 min excess over 1: {r0['min_excess']:.6f} at {r0['min_rec']}", flush=True)

    # Prove algebraically: Δ(m)/Δ(m-1) = 4(k-m+1)(k-m+2)/[m(m+1)] > 1
    # iff 4(k-m+1)(k-m+2) > m(m+1)
    # At mode m of I(k,0) (which also = mode of A_0 for j=0):
    # m = floor((2k+1)/3) for pure spiders.
    # Check: 4(k-m+1)(k-m+2) > m(m+1) for m = floor((2k+1)/3)?
    print("\n=== Part 1b: Algebraic check 4(k-m+1)(k-m+2) > m(m+1) ===", flush=True)
    algebraic_fails = []
    for k in range(6, args.k_max + 1):
        a0 = build_aj_coeffs(k, 0)
        g = build_g_coeffs(k, len(a0) + 2)
        m, _ = find_mode_and_lambda(a0, g)
        if m < 1:
            continue
        lhs = 4 * (k - m + 1) * (k - m + 2)
        rhs = m * (m + 1)
        if lhs <= rhs:
            algebraic_fails.append({"k": k, "m": m, "lhs": lhs, "rhs": rhs})
    print(f"Algebraic inequality fails: {len(algebraic_fails)}", flush=True)
    if algebraic_fails:
        print(f"  Examples: {algebraic_fails[:3]}", flush=True)

    # Part 2: Compute F=T1+T2 for general (k,j) and track ratio |T2|/T1
    print("\n=== Part 2: General (k,j) F-components ===", flush=True)
    f_fails = 0
    t1_fails = 0
    t2_pos = 0
    min_F = None
    min_F_rec = None
    max_ratio = None  # max |T2|/T1
    max_ratio_rec = None
    total = 0

    # Also track the delta_ratio (Δ(m)/Δ(m-1)) for general j
    min_delta_ratio = None
    min_dr_rec = None

    for k in range(6, args.k_max + 1):
        for j in range(0, args.j_max + 1):
            rec = compute_F_components(k, j)
            if rec.get("skip"):
                continue
            total += 1

            T1 = rec["T1"]
            T2 = rec["T2"]
            F = rec["F"]

            if T1 <= 0:
                t1_fails += 1
            if T2 > 0:
                t2_pos += 1
            if F < 0:
                f_fails += 1

            if min_F is None or F < min_F:
                min_F = F
                min_F_rec = rec

            ratio = abs(T2) / T1 if T1 > 0 else None
            if ratio is not None and (max_ratio is None or ratio > max_ratio):
                max_ratio = ratio
                max_ratio_rec = rec

            # Δ(m)/Δ(m-1) for A_j
            dr = rec.get("delta_ratio")
            if dr is not None and (min_delta_ratio is None or dr < min_delta_ratio):
                min_delta_ratio = dr
                min_dr_rec = rec

        if (k - 6) % 50 == 0 or k == args.k_max:
            print(
                f"k={k}: total={total}, F_fails={f_fails}, T1_fails={t1_fails}, "
                f"T2_pos={t2_pos}, max|T2|/T1={max_ratio:.6f}",
                flush=True,
            )

    print(f"\nSummary:", flush=True)
    print(f"  total pairs: {total}", flush=True)
    print(f"  F < 0 failures: {f_fails}", flush=True)
    print(f"  T1 <= 0 failures: {t1_fails}", flush=True)
    print(f"  T2 > 0 cases: {t2_pos}", flush=True)
    print(f"  max |T2|/T1 = {max_ratio:.8f}", flush=True)
    if max_ratio_rec:
        print(f"    at (k,j)=({max_ratio_rec['k']},{max_ratio_rec['j']}), T1={max_ratio_rec['T1']}, T2={max_ratio_rec['T2']}", flush=True)
    if min_F_rec:
        print(f"  min F = {min_F} at (k,j)=({min_F_rec['k']},{min_F_rec['j']})", flush=True)
    print(f"  min Δ(m)/Δ(m-1) = {min_delta_ratio:.6f}", flush=True)
    if min_dr_rec:
        print(f"    at (k,j)=({min_dr_rec['k']},{min_dr_rec['j']}), m={min_dr_rec['m']}", flush=True)

    print(f"\nwall={time.time()-t0:.1f}s", flush=True)

    # Key algebraic claim for j=0:
    print("\n=== Part 3: Algebraic proof sketch for j=0 ===", flush=True)
    print("Δ(m)/Δ(m-1) = 4(k-m+1)(k-m+2)/[m(m+1)] > 1", flush=True)
    print("iff 4(k-m+1)(k-m+2) > m(m+1)", flush=True)
    print("At m = floor((2k+1)/3), let k=3t+r (r=0,1,2):", flush=True)
    for r in range(3):
        print(f"  r={r}: checking k=6..{args.k_max} with k≡{r} mod 3:", flush=True)
        min_margin_alg = None
        for t in range(2, (args.k_max - r) // 3 + 1):
            k = 3 * t + r
            if k < 6:
                continue
            a0 = build_aj_coeffs(k, 0)
            g0 = build_g_coeffs(k, len(a0) + 2)
            m, _ = find_mode_and_lambda(a0, g0)
            margin_alg = 4 * (k - m + 1) * (k - m + 2) - m * (m + 1)
            if min_margin_alg is None or margin_alg < min_margin_alg:
                min_margin_alg = margin_alg
                min_k = k
                min_m = m
        if min_margin_alg is not None:
            print(f"    min margin = {min_margin_alg} at k={min_k}, m={min_m}", flush=True)

    if args.out:
        os.makedirs(os.path.dirname(args.out) if os.path.dirname(args.out) else ".", exist_ok=True)
        # Don't write large data; just write summary
        summary = {
            "params": vars(args),
            "j0_closed_form": r0,
            "general": {
                "total": total,
                "F_fails": f_fails,
                "T1_fails": t1_fails,
                "T2_pos": t2_pos,
                "max_ratio_T2_T1": max_ratio,
                "max_ratio_rec": {k: (str(v) if isinstance(v, Fraction) else v) for k, v in (max_ratio_rec or {}).items()},
                "min_F": min_F,
                "min_delta_ratio": min_delta_ratio,
            },
            "wall_s": time.time() - t0,
        }
        with open(args.out, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
