#!/usr/bin/env python3
"""Rigorous proof that the Cluster Slack bound holds.

CLAIM: f(x) = x + 2 - 2·sqrt(x/(1-x)) < 1  for all x ∈ (1/3, 1/2].

PROOF (no numerical scanning, pure algebra + one sign check):

Step 1.  Substitute t = sqrt(x/(1-x)).
         When x ∈ (1/3, 1/2], we have t ∈ (1/√2, 1].
         Then x = t²/(1+t²), and

           f(x) - 1 = g(x) = (1 + 2t² - 2t - 2t³) / (1 + t²).

         The denominator 1+t² > 0, so sign(g) = sign(N(t)) where

           N(t) = -2t³ + 2t² - 2t + 1.

Step 2.  Show N(t) < 0 for all t ∈ (1/√2, 1].

         N'(t) = -6t² + 4t - 2.

         Discriminant of N'(t) = 16 - 48 = -32 < 0.
         Leading coefficient of N' is -6 < 0.
         Therefore N'(t) < 0 for ALL real t.

         ==> N is strictly decreasing on all of ℝ.

Step 3.  Evaluate N at the left endpoint:

           N(1/√2) = -2/(2√2) + 2/2 - 2/√2 + 1
                    = -1/√2 + 1 - 2/√2 + 1
                    = 2 - 3/√2
                    = 2 - 3√2/2
                    ≈ -0.1213.

         Since N(1/√2) < 0 and N is strictly decreasing,
         N(t) < N(1/√2) < 0 for all t > 1/√2.

Step 4.  Therefore g(x) = N(t)/(1+t²) < 0,
         hence f(x) = g(x) + 1 < 1.  □

The maximum of f on the interval is:
  sup f = lim_{x→1/3+} f(x) = 1/3 + 2 - √2 ≈ 0.9191,
giving slack ≈ 0.0809.
"""
import math
from fractions import Fraction

def N(t):
    """N(t) = -2t³ + 2t² - 2t + 1"""
    return -2*t**3 + 2*t**2 - 2*t + 1

def N_prime(t):
    """N'(t) = -6t² + 4t - 2"""
    return -6*t**2 + 4*t - 2

def f(x):
    """f(x) = x + 2 - 2·sqrt(x/(1-x))"""
    return x + 2 - 2*math.sqrt(x/(1-x))

def main():
    print("=" * 65)
    print("  RIGOROUS PROOF: Cluster Slack Bound")
    print("=" * 65)
    print()
    print("Claim: f(x) = x + 2 - 2√(x/(1-x)) < 1  for x ∈ (1/3, 1/2].")
    print()

    # ── Step 1: Substitution ──
    print("STEP 1: Substitution t = √(x/(1-x))")
    print("  x ∈ (1/3, 1/2]  ⟹  t ∈ (1/√2, 1]")
    print("  x = t²/(1+t²)")
    print()
    print("  f(x) - 1 = N(t)/(1+t²)")
    print("  where N(t) = -2t³ + 2t² - 2t + 1")
    print()

    # Verify substitution at t=1 (x=1/2)
    t_test = 1.0
    x_test = t_test**2 / (1 + t_test**2)
    lhs = f(x_test) - 1
    rhs = N(t_test) / (1 + t_test**2)
    print(f"  Verify at t=1 (x=1/2):  f(x)-1 = {lhs:.6f},  N(t)/(1+t²) = {rhs:.6f}")
    assert abs(lhs - rhs) < 1e-12, "Substitution check failed!"
    print("  ✓ Substitution verified.")
    print()

    # ── Step 2: N'(t) < 0 everywhere ──
    print("STEP 2: N(t) is strictly decreasing")
    print("  N'(t) = -6t² + 4t - 2")
    print()
    disc = 4**2 - 4*(-6)*(-2)  # = 16 - 48 = -32
    print(f"  Discriminant of N'(t) = 4² - 4·(-6)·(-2) = 16 - 48 = {disc}")
    assert disc < 0, "Discriminant should be negative!"
    print(f"  Leading coefficient = -6 < 0")
    print()
    print("  Negative leading coeff + negative discriminant")
    print("  ⟹  N'(t) < 0 for ALL real t")
    print("  ⟹  N is strictly decreasing on ℝ.")
    print()

    # Spot-check
    for t in [0.0, 0.5, 0.707, 1.0, 2.0]:
        print(f"    N'({t:.3f}) = {N_prime(t):.4f}", end="")
        assert N_prime(t) < 0, f"N'({t}) should be negative!"
        print("  < 0  ✓")
    print()

    # ── Step 3: N(1/√2) < 0 ──
    print("STEP 3: Evaluate N at the left endpoint")
    t0 = 1 / math.sqrt(2)
    N_t0 = N(t0)

    # Exact arithmetic: N(1/√2) = 2 - 3/√2 = 2 - 3√2/2
    exact_val = 2 - 3*math.sqrt(2)/2
    print(f"  N(1/√2) = -2·(1/√2)³ + 2·(1/√2)² - 2·(1/√2) + 1")
    print(f"          = -1/√2 + 1 - 2/√2 + 1")
    print(f"          = 2 - 3/√2")
    print(f"          = 2 - 3√2/2")
    print(f"          ≈ {exact_val:.10f}")
    print()
    assert abs(N_t0 - exact_val) < 1e-12, "Exact value mismatch!"
    assert N_t0 < 0, "N(1/√2) should be negative!"
    print(f"  2 - 3√2/2 < 0  ⟺  4 < 3√2  ⟺  16 < 18.  TRUE.  ✓")
    print()

    # Also check N(1)
    N_1 = N(1.0)
    print(f"  N(1) = -2 + 2 - 2 + 1 = {N_1:.0f}")
    assert N_1 == -1, "N(1) should be -1!"
    print(f"  ✓ N(1) = -1 < 0")
    print()

    # ── Step 4: Conclusion ──
    print("STEP 4: Conclusion")
    print("  N is strictly decreasing (Step 2)")
    print("  N(1/√2) < 0 (Step 3)")
    print("  ⟹  N(t) ≤ N(1/√2) < 0  for all t ≥ 1/√2")
    print("  ⟹  g(x) = N(t)/(1+t²) < 0  for x ∈ (1/3, 1/2]")
    print("  ⟹  f(x) = g(x) + 1 < 1.  □")
    print()

    sup_f = Fraction(1, 3) + 2 - float(math.sqrt(2))  # not exact but close
    sup_f_exact = 1/3 + 2 - math.sqrt(2)
    slack = 1 - sup_f_exact
    print(f"  sup f = lim_{{x→1/3+}} f(x) = 1/3 + 2 - √2 ≈ {sup_f_exact:.6f}")
    print(f"  Slack = 1 - sup f ≈ {slack:.6f}")
    print()
    print("=" * 65)
    print("  PROOF IS COMPLETE.  No sampling. No numerics.")
    print("  Only: substitution, discriminant sign, one evaluation.")
    print("=" * 65)

if __name__ == "__main__":
    main()
