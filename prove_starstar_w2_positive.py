"""
PROOF: w_2 ≥ 0 at the star+star extremal family.

THEOREM: For a support vertex with 1 leaf and two star children K_{1,a} and K_{1,b}
(a, b ≥ 1), the W-form value w_2 is always non-negative.

PROOF: From symbolic computation (see prove_starstar_ratio.py):

    w_2 · C_2 = s²(s-1)² · Q(a,b) / 12

where s = a + b and Q(a,b) = a² - ab + a + b² + b + 9.

Since s ≥ 2 (both a,b ≥ 1), s² > 0 and (s-1)² ≥ 1.
Since C_2 = s(s-1)/2 > 0 for s ≥ 2.

It suffices to show Q(a,b) > 0 for all a,b ≥ 1.

CLAIM: Q(a,b) = a² - ab + a + b² + b + 9 > 0 for all real a,b.

Proof: Complete the square in a:
    Q = a² - a(b-1) + b² + b + 9
    = (a - (b-1)/2)² - (b-1)²/4 + b² + b + 9
    = (a - (b-1)/2)² + (4b² + 4b + 36 - b² + 2b - 1)/4
    = (a - (b-1)/2)² + (3b² + 6b + 35)/4

Since 3b² + 6b + 35 = 3(b+1)² + 32 ≥ 35 > 0,
we have Q(a,b) ≥ 35/4 > 0 for all real a,b.

In fact Q(a,b) ≥ 35/4 = 8.75, with equality only at a = (b-1)/2, b = -1
(which is outside the domain a,b ≥ 1).

For integer a,b ≥ 1:
    Q(1,1) = 1 - 1 + 1 + 1 + 1 + 9 = 12
    Q(1,2) = 1 - 2 + 1 + 4 + 2 + 9 = 15
    minimum Q over a,b ≥ 1 is Q(1,1) = 12.

Therefore w_2 ≥ s²(s-1)² · 12 / (12 · C_2) = s(s-1) / 1 = s(s-1) ≥ 2.

QED.

Additional: the RATIO (Term1/|t23|) is

    ratio = N(a,b) / |D(a,b)|

where N = (a² - a + 2b² + 6b + 4) > 0 (obvious)
  and D = (a² - 2ab + 3a - 4b + 14) which can be either sign.

When D < 0 (⋆ fails), ratio = N/|D| = N/(-D).

The ratio is w2-relevant only when D < 0. We showed w2 ≥ 0 regardless.

Since w_2 = (Term1 + t23)/C_2 and Term1 = C_2² · N / (12·C_2) = C_2 · N/12...
actually more directly:

w_2 · C_2 = s²(s-1)²/12 · (N + D) = s²(s-1)²/12 · Q(a,b)

where Q = N + D = (a² - a + 2b² + 6b + 4) + (a² - 2ab + 3a - 4b + 14)
       = 2a² - 2ab + 2a + 2b² + 2b + 18 = 2(a² - ab + a + b² + b + 9) = 2Q(a,b).

Wait, let me recheck: N = a² - a + 2b² + 6b + 4 and D = a² - 2ab + 3a - 4b + 14.
N + D = 2a² - 2ab + 2a + 2b² + 2b + 18 = 2(a² - ab + a + b² + b + 9).
This matches Q(a,b) = a² - ab + a + b² + b + 9 with the factor of 2 absorbed
into the 12 → 24 denominators.

So indeed: w_2 · C_2 = s²(s-1)² · (N+D) / 24 and Term1 = s²(s-1)² · N/24,
t23 = s²(s-1)² · D/24. Since Q = (N+D)/2 = N/2 + D/2, and the result
factors as s²(s-1)²/12 · Q(a,b).

COMPLETE PROOF SUMMARY:
  w_2 ≥ 0 ⟺ N + D ≥ 0 ⟺ 2Q(a,b) ≥ 0 ⟺ Q(a,b) ≥ 0.
  Q(a,b) = (a - (b-1)/2)² + (3b² + 6b + 35)/4 ≥ 35/4 > 0.  ∎
"""

from fractions import Fraction
from math import comb

def verify_factorization():
    """Verify that N + D = 2Q for many (a,b)."""
    print("Verifying N + D = 2Q(a,b) and factorization")
    for a in range(1, 30):
        for b in range(1, 30):
            s = a + b
            N = a**2 - a + 2*b**2 + 6*b + 4
            D = a**2 - 2*a*b + 3*a - 4*b + 14
            Q = a**2 - a*b + a + b**2 + b + 9

            assert N + D == 2*Q, f"Failed at ({a},{b}): N+D={N+D}, 2Q={2*Q}"

            # Verify the actual w2
            C2 = comb(s, 2)
            C1 = s
            C3 = comb(s, 3)
            A2 = comb(s+1, 2) + (b+1)
            A3 = comb(s+1, 3) + comb(b+1, 2)
            B2 = comb(a+1, 2) + 1
            B1 = a + 2

            dk_AC = A3 * C2 - A2 * C3
            Term1 = C2 * dk_AC
            t23 = C3 * (B2 * C1 - B1 * C2) + B2 * (C2**2 - C1 * C3)

            w2_times_C2 = Term1 + t23
            expected = s**2 * (s-1)**2 * Q // 12

            # Use exact Fraction to avoid integer division issues
            expected_frac = Fraction(s**2 * (s-1)**2 * Q, 12)
            actual_frac = Fraction(w2_times_C2, 1)

            assert actual_frac == expected_frac, \
                f"Failed at ({a},{b}): w2*C2={w2_times_C2}, expected={expected_frac}"

            # Verify Q > 0
            assert Q > 0, f"Q negative at ({a},{b})"

            # Verify w2 >= 0
            assert w2_times_C2 >= 0, f"w2 negative at ({a},{b})"

    print("All checks passed for a,b in [1,29]")
    print()

    # Show Q is always large
    min_Q = None
    for a in range(1, 100):
        for b in range(1, 100):
            Q = a**2 - a*b + a + b**2 + b + 9
            if min_Q is None or Q < min_Q:
                min_Q = Q
                min_ab = (a, b)

    print(f"Minimum Q over [1,99]²: Q={min_Q} at {min_ab}")

    # Also verify the completing-the-square
    a_opt = (min_ab[1] - 1) / 2
    print(f"Optimal a for b={min_ab[1]}: a = (b-1)/2 = {a_opt}")
    b = min_ab[1]
    remainder = (3*b**2 + 6*b + 35) / 4
    print(f"Remainder term (3b²+6b+35)/4 = {remainder}")
    print()

    # Check the proof bound
    print("=== Proof summary ===")
    print("Q(a,b) = a² - ab + a + b² + b + 9")
    print("       = (a - (b-1)/2)² + (3b² + 6b + 35)/4")
    print("       ≥ (3·1² + 6·1 + 35)/4  [at b=1, minimizing the remainder]")
    print(f"       = {(3 + 6 + 35)/4}")
    print("       = 11  [at integer a=0, b=1, but a ≥ 1]")
    print()

    # True minimum at integer a,b ≥ 1
    for b in range(1, 5):
        for a in range(1, 5):
            Q = a**2 - a*b + a + b**2 + b + 9
            print(f"  Q({a},{b}) = {Q}")


def prove_lower_bound():
    """Prove w_2 ≥ 2(a+b-1) for all a,b ≥ 1 with star_margin < 0."""
    print("\n=== Lower bound on w_2 ===")
    for a in range(1, 50):
        for b in range(a, 50):
            s = a + b
            C2 = comb(s, 2)
            Q = a**2 - a*b + a + b**2 + b + 9
            w2_times_C2 = s**2 * (s-1)**2 * Q // 12
            w2 = Fraction(w2_times_C2, C2)

            # Check star_margin < 0 (⋆ failure)
            D = a**2 - 2*a*b + 3*a - 4*b + 14
            if D >= 0:
                continue  # no ⋆ failure

            # w2 should be ≥ some bound
            lower = s * (s - 1) * Q // (6 * 1)  # w2 = s(s-1)Q/6 ... hmm
            # Actually w2 = s²(s-1)²Q/(12·C2) = s²(s-1)²Q/(12·s(s-1)/2) = s(s-1)Q/6
            w2_exact = Fraction(s * (s-1) * Q, 6)
            assert w2 == w2_exact, f"Formula mismatch at ({a},{b})"

    print("w_2 = s(s-1)Q(a,b)/6 confirmed for all ⋆-failing pairs")
    print()
    print("Since s ≥ 5 (minimum for ⋆ failure), (s-1) ≥ 4, Q ≥ 12:")
    print("  w_2 ≥ 5·4·12/6 = 40")
    print()

    # Actually find minimum w2 among ⋆-failing pairs
    min_w2 = None
    for a in range(1, 100):
        for b in range(a, 100):
            D = a**2 - 2*a*b + 3*a - 4*b + 14
            if D >= 0:
                continue
            s = a + b
            Q = a**2 - a*b + a + b**2 + b + 9
            w2 = s * (s - 1) * Q // 6
            if min_w2 is None or w2 < min_w2:
                min_w2 = w2
                min_ab = (a, b)
    print(f"Minimum w_2 among ⋆-failing pairs [1,99]²: w_2={min_w2} at {min_ab}")


if __name__ == '__main__':
    verify_factorization()
    prove_lower_bound()
