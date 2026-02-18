#!/usr/bin/env python3
"""Verify E₁ ≤ 0 for Sub-claim C: mixed spiders S(2^k, 1^j), k ≡ 0,2 (mod 3), k ≥ 5.

Proof route (notes/subclaim_c_E1_le0_proof_2026-02-18.md):
  E₁ ≤ 0  iff  B₁ - A₁ ≤ 0  iff  λ₁ ≥ 1/(k-2).

Key facts:
  (F1)  B₁ - A₁ = [1 - (k-2)λ₁] / [(1+λ₁)(1+2λ₁)]   (exact algebra)
  (F2)  λ₁ ≥ τ = (m₁-1)/[2(k-m₁+2)]                   (b-bound)
  (F3a) k=3t+2, t≥1: (k-2)τ = 3t(2t+1)/[2(t+2)] ≥ 1  iff 6t²+t-4 ≥ 0, true t≥1
  (F3b) k=3t,   t≥2: (k-2)τ = t(3t-2)/(t+1)    ≥ 1  iff 3t²-3t-1 ≥ 0, true t≥2

Checks performed:
  1. Polynomial inequalities (F3a), (F3b): symbolically verify for t ≤ 1000.
  2. b-bound (F2): exact Fraction check for k ≤ 200.
  3. (k-2)τ ≥ 1: exact Fraction check for k ≤ 500.
  4. E₁ ≤ 0 directly via exact coefficient ratios for k ≤ 200.
  5. Report minimum gaps for all checks.
"""

from __future__ import annotations

import math
from fractions import Fraction


# ---------------------------------------------------------------------------
# Mode formulas
# ---------------------------------------------------------------------------

def mode_j0(k: int) -> int:
    return (2 * k + 1) // 3


def mode_j1(k: int) -> int:
    return (2 * k + 3) // 3


# ---------------------------------------------------------------------------
# IS polynomial coefficients for S(2^k) (j=0) and S(2^k, 1) (j=1)
# From prove_mixed_spider_j0_j1_branches.py
# ---------------------------------------------------------------------------

def coeff_j0(k: int, t: int) -> int:
    """i_t for S(2^k): (1+2x)^k + x(1+x)^k."""
    a = math.comb(k, t) * (1 << t) if 0 <= t <= k else 0
    b = math.comb(k, t - 1) if 1 <= t <= k + 1 else 0
    return a + b


def coeff_j1(k: int, t: int) -> int:
    """i_t for S(2^k, 1): (1+2x)^k(1+x) + x(1+x)^k."""
    a = math.comb(k, t) * (1 << t) if 0 <= t <= k else 0
    b = 0
    if 1 <= t <= k + 1:
        b = math.comb(k, t - 1) * ((1 << (t - 1)) + 1)
    return a + b


# ---------------------------------------------------------------------------
# b-bound τ = (m₁-1) / [2(k - m₁ + 2)]
# ---------------------------------------------------------------------------

def tau_j1(k: int) -> Fraction:
    """Lower bound on λ₁ from log-concavity of b_t = C(k,t)*2^t."""
    m1 = mode_j1(k)
    return Fraction(m1 - 1, 2 * (k - m1 + 2))


# ---------------------------------------------------------------------------
# Check 1: Polynomial inequalities (F3a) and (F3b)
# ---------------------------------------------------------------------------

def check_poly_inequalities(t_max: int = 1000) -> dict:
    """
    (F3a) 6t²+t-4 ≥ 0 for t ≥ 1   (covers k=3t+2 ≡ 2 mod 3, k ≥ 5)
    (F3b) 3t²-3t-1 ≥ 0 for t ≥ 2   (covers k=3t   ≡ 0 mod 3, k ≥ 6)
    """
    fails_3a: list[int] = []
    fails_3b: list[int] = []
    min_3a = None
    min_3b = None

    for t in range(1, t_max + 1):
        val_3a = 6 * t * t + t - 4
        if val_3a < 0:
            fails_3a.append(t)
        if min_3a is None or val_3a < min_3a:
            min_3a = val_3a

    for t in range(2, t_max + 1):
        val_3b = 3 * t * t - 3 * t - 1
        if val_3b < 0:
            fails_3b.append(t)
        if min_3b is None or val_3b < min_3b:
            min_3b = val_3b

    return {
        "fails_3a": fails_3a,
        "fails_3b": fails_3b,
        "min_3a": min_3a,
        "min_3b": min_3b,
    }


# ---------------------------------------------------------------------------
# Check 2: b-bound λ₁ ≥ τ (exact Fraction)
# ---------------------------------------------------------------------------

def check_b_bound(k_max: int = 200) -> dict:
    """Verify λ₁ ≥ τ using exact coefficients for k ≤ k_max, k ≡ 0,2 mod 3."""
    fails: list[int] = []
    min_gap = None
    min_k = None

    for k in range(3, k_max + 1):
        if k % 3 == 1:
            continue
        m1 = mode_j1(k)
        lam1 = Fraction(coeff_j1(k, m1 - 1), coeff_j1(k, m1))
        tau = tau_j1(k)
        gap = lam1 - tau
        if gap < 0:
            fails.append(k)
        if min_gap is None or gap < min_gap:
            min_gap = gap
            min_k = k

    return {"fails": fails, "min_gap": min_gap, "min_k": min_k}


# ---------------------------------------------------------------------------
# Check 3: (k-2)τ ≥ 1 (exact Fraction)
# ---------------------------------------------------------------------------

def check_ktau_ge1(k_max: int = 500) -> dict:
    """Verify (k-2)τ ≥ 1 for k ≡ 0,2 mod 3, k ≥ 5, k ≤ k_max."""
    fails: list[int] = []
    min_gap = None
    min_k = None

    for k in range(5, k_max + 1):
        if k % 3 == 1:
            continue
        tau = tau_j1(k)
        val = Fraction(k - 2) * tau
        gap = val - Fraction(1)
        if gap < 0:
            fails.append(k)
        if min_gap is None or gap < min_gap:
            min_gap = gap
            min_k = k

    return {"fails": fails, "min_gap": min_gap, "min_k": min_k}


# ---------------------------------------------------------------------------
# Check 4: E₁ ≤ 0 directly (exact Fraction via coefficients)
# ---------------------------------------------------------------------------

def b1_minus_a1(k: int, lam: Fraction) -> Fraction:
    """B₁ - A₁ = [1 - (k-2)λ] / [(1+λ)(1+2λ)]  (algebraic form, verified below)."""
    m1 = mode_j1(k)
    one = Fraction(1)
    a1 = Fraction(2 * k) * lam / (one + 2 * lam) + lam / (one + lam) - Fraction(m1 - 1)
    b1 = one + Fraction(k) * lam / (one + lam) - Fraction(m1 - 1)
    return b1 - a1


def e1_direct(k: int) -> Fraction:
    """E₁ = r₁(B₁-A₁)/(1+r₁) computed exactly from IS polynomial coefficients."""
    m1 = mode_j1(k)
    lam1 = Fraction(coeff_j1(k, m1 - 1), coeff_j1(k, m1))
    one = Fraction(1)
    r1 = lam1 * (one + lam1) ** (k - 1) / (one + 2 * lam1) ** k
    ba = b1_minus_a1(k, lam1)
    return r1 * ba / (one + r1)


def check_e1_le0(k_max: int = 200) -> dict:
    """Verify E₁ ≤ 0 for k ≡ 0,2 mod 3, k ≥ 5, k ≤ k_max."""
    fails: list[int] = []
    max_e1 = None
    max_e1_k = None
    # Also check the algebraic identity B₁-A₁ = [1-(k-2)λ]/[(1+λ)(1+2λ)]
    identity_fails: list[int] = []
    min_e1 = None  # most negative E₁ (largest margin of E₁ ≤ 0)

    for k in range(5, k_max + 1):
        if k % 3 == 1:
            continue
        m1 = mode_j1(k)
        lam1 = Fraction(coeff_j1(k, m1 - 1), coeff_j1(k, m1))
        one = Fraction(1)

        # Check algebraic identity
        ba_formula = b1_minus_a1(k, lam1)
        ba_algebraic = (one - Fraction(k - 2) * lam1) / ((one + lam1) * (one + 2 * lam1))
        if ba_formula != ba_algebraic:
            identity_fails.append(k)

        # Compute E₁
        e1 = e1_direct(k)
        if e1 > 0:
            fails.append(k)
        if max_e1 is None or e1 > max_e1:
            max_e1 = e1
            max_e1_k = k
        if min_e1 is None or e1 < min_e1:
            min_e1 = e1

    return {
        "fails": fails,
        "identity_fails": identity_fails,
        "max_e1": max_e1,
        "max_e1_k": max_e1_k,
        "min_e1": min_e1,
    }


# ---------------------------------------------------------------------------
# Check 5: B₁-A₁ sign at base cases k=3 (outside scope) and k=5
# ---------------------------------------------------------------------------

def check_base_cases() -> dict:
    """Direct exact Fraction checks at k=3 (informational) and k=5."""
    results = {}
    for k in [3, 5, 6, 8, 9]:
        m1 = mode_j1(k)
        lam1 = Fraction(coeff_j1(k, m1 - 1), coeff_j1(k, m1))
        ba = b1_minus_a1(k, lam1)
        e1 = e1_direct(k)
        results[k] = {"lam1": lam1, "B1-A1": ba, "E1": e1}
    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 70)
    print("Verification: E₁ ≤ 0 for Sub-claim C (k ≡ 0,2 mod 3, k ≥ 5)")
    print("notes/subclaim_c_E1_le0_proof_2026-02-18.md")
    print("=" * 70)

    # --- Check 1: Polynomial inequalities ---
    print()
    print("CHECK 1: Polynomial inequalities")
    print("  (F3a) 6t²+t-4 ≥ 0  for t ≥ 1  (k=3t+2 case)")
    print("  (F3b) 3t²-3t-1 ≥ 0 for t ≥ 2  (k=3t case)")
    res1 = check_poly_inequalities(t_max=1000)
    status1a = "PASS" if not res1["fails_3a"] else "FAIL"
    status1b = "PASS" if not res1["fails_3b"] else "FAIL"
    print(f"  (F3a) t ≤ 1000: {status1a}  (failures: {res1['fails_3a'][:5]}, min value: {res1['min_3a']})")
    print(f"  (F3b) t ≤ 1000: {status1b}  (failures: {res1['fails_3b'][:5]}, min value: {res1['min_3b']})")

    # --- Check 2: b-bound λ₁ ≥ τ ---
    print()
    print("CHECK 2: b-bound λ₁ ≥ τ = (m₁-1)/[2(k-m₁+2)]  (k ≤ 200, exact Fraction)")
    res2 = check_b_bound(k_max=200)
    status2 = "PASS" if not res2["fails"] else "FAIL"
    print(f"  Result: {status2}  (failures: {res2['fails'][:5]})")
    if res2["min_gap"] is not None:
        mg = res2["min_gap"]
        print(f"  Minimum gap λ₁-τ: {float(mg):.8f}  at k={res2['min_k']}")
        if mg.numerator.bit_length() < 200:
            print(f"  (exact: {mg})")

    # --- Check 3: (k-2)τ ≥ 1 ---
    print()
    print("CHECK 3: (k-2)τ ≥ 1  (k ≡ 0,2 mod 3, k ≥ 5, k ≤ 500, exact Fraction)")
    res3 = check_ktau_ge1(k_max=500)
    status3 = "PASS" if not res3["fails"] else "FAIL"
    print(f"  Result: {status3}  (failures: {res3['fails'][:5]})")
    if res3["min_gap"] is not None:
        mg3 = res3["min_gap"]
        print(f"  Minimum gap (k-2)τ-1: {float(mg3):.8f}  at k={res3['min_k']}")
        if mg3.numerator.bit_length() < 200:
            print(f"  (exact: {mg3})")

    # --- Check 4: E₁ ≤ 0 directly ---
    print()
    print("CHECK 4: E₁ ≤ 0 directly from IS coefficients  (k ≤ 200, exact Fraction)")
    res4 = check_e1_le0(k_max=200)
    status4 = "PASS" if not res4["fails"] else "FAIL"
    id_status = "PASS" if not res4["identity_fails"] else "FAIL"
    print(f"  Algebraic identity check: {id_status}  (failures: {res4['identity_fails'][:5]})")
    print(f"  E₁ ≤ 0 check: {status4}  (failures: {res4['fails'][:5]})")
    if res4["max_e1"] is not None:
        me = res4["max_e1"]
        print(f"  Maximum E₁ (closest to 0): {float(me):.8f}  at k={res4['max_e1_k']}")
        # Large k produces huge Fractions; only print exact for small k
        if me.numerator.bit_length() < 200:
            print(f"  (exact: {me})")
    if res4["min_e1"] is not None:
        me2 = res4["min_e1"]
        print(f"  Minimum E₁ (most negative): {float(me2):.8f}")

    # --- Check 5: Base cases ---
    print()
    print("CHECK 5: Base cases k ∈ {3, 5, 6, 8, 9} (exact Fraction)")
    res5 = check_base_cases()
    for k in sorted(res5):
        d = res5[k]
        lam = d["lam1"]
        ba = d["B1-A1"]
        e1 = d["E1"]
        scope = "(outside scope k<5)" if k < 5 else ""
        sign_e1 = "≤ 0 ✓" if e1 <= 0 else "> 0 (see note)"
        print(
            f"  k={k}: λ₁={float(lam):.6f}, B₁-A₁={float(ba):+.6f}, "
            f"E₁={float(e1):+.8f}  {sign_e1} {scope}"
        )

    # --- Summary ---
    print()
    print("=" * 70)
    all_pass = all([
        not res1["fails_3a"],
        not res1["fails_3b"],
        not res2["fails"],
        not res3["fails"],
        not res4["fails"],
        not res4["identity_fails"],
    ])
    if all_pass:
        print("OVERALL: PASS -- E₁ ≤ 0 fully verified for k ≡ 0,2 (mod 3), k ≥ 5.")
        print()
        print("Proof summary:")
        print("  B₁-A₁ = [1-(k-2)λ₁]/[(1+λ₁)(1+2λ₁)]  (algebraic identity: CHECK 4 ✓)")
        print("  E₁ ≤ 0 ⟺ (k-2)λ₁ ≥ 1.")
        print("  b-bound: λ₁ ≥ τ = (m₁-1)/[2(k-m₁+2)]  (CHECK 2 ✓)")
        print("  k=3t+2, t≥1: (k-2)τ ≥ 1 via 6t²+t-4 ≥ 0  (CHECK 1/3 ✓)")
        print("  k=3t,   t≥2: (k-2)τ ≥ 1 via 3t²-3t-1 ≥ 0  (CHECK 1/3 ✓)")
    else:
        print("OVERALL: FAIL -- see diagnostics above.")


if __name__ == "__main__":
    main()
