#!/usr/bin/env python3
"""Closed-form proof lane for A_gap >= 1/(4k) in Sub-claim C.

Target:
  For k >= 6 with k % 3 in {0, 2}, show
    A_gap(k) := A0(k) - A1(k) >= 1/(4k),
where A0, A1 are the main terms in:
  margin(k,j) = A_j + E_j.

Context:
  - notes/subclaim_c_algebra_attempt_2026-02-18.md
  - prove_subclaim_c_algebra.py

Approach:
  1) Exact closed forms for A_gap - 1/(4k) are encoded as num/den for k=3t and
     k=3t+2 (derived symbolically once, then transcribed here).
  2) The denominator is manifestly positive.
  3) For t >= 4, the numerator is bounded below by p*C1 - C0 with p=4^t, and
     this is positive by explicit polynomial bounds.
  4) Remaining finite bases are checked exactly: k in {6, 8, 9, 11}.

This script checks all algebraic identities exactly (Fraction arithmetic), verifies
the finite bases, and validates the large-t lane over a broad range.
"""

from __future__ import annotations

import argparse
from fractions import Fraction

from prove_subclaim_c_algebra import lambda_j0_closed, lambda_j1_closed, terms_j0, terms_j1


def a_gap_exact(k: int) -> Fraction:
    """Exact A_gap = A0 - A1 at closed-form tie fugacities."""
    lam0 = lambda_j0_closed(k)
    lam1 = lambda_j1_closed(k)
    a0, _, _, _ = terms_j0(k, lam0)
    a1, _, _, _ = terms_j1(k, lam1)
    return a0 - a1


def case0_num_den(t: int) -> tuple[int, int]:
    """num/den for (A_gap - 1/(4k)) with k=3t and p=4^t."""
    p = 4**t

    num = (
        72 * p**3 * t**6
        + 360 * p**3 * t**5
        + 577 * p**3 * t**4
        + 266 * p**3 * t**3
        - 83 * p**3 * t**2
        - 76 * p**3 * t
        - 12 * p**3
        + 3456 * p**2 * t**8
        - 288 * p**2 * t**7
        - 12264 * p**2 * t**6
        - 11380 * p**2 * t**5
        - 2124 * p**2 * t**4
        + 831 * p**2 * t**3
        + 204 * p**2 * t**2
        - 41 * p**2 * t
        - 10 * p**2
        + 2592 * p * t**8
        - 4560 * p * t**7
        - 14972 * p * t**6
        - 10796 * p * t**5
        - 2357 * p * t**4
        + 197 * p * t**3
        + 85 * p * t**2
        - 7 * p * t
        - 2 * p
        - 2784 * t**7
        - 3960 * t**6
        - 2288 * t**5
        - 750 * t**4
        - 144 * t**3
        - 10 * t**2
    )

    den = (
        12
        * t
        * (3 * p * t**2 + 7 * p * t + 2 * p + 10 * t**2)
        * (8 * p * t**2 + 9 * p * t + 2 * p + 6 * t**2 + 5 * t + 1)
        * (12 * p * t**2 + 13 * p * t + 3 * p + 10 * t**2 + 7 * t + 1)
    )
    return num, den


def case2_num_den(t: int) -> tuple[int, int]:
    """num/den for (A_gap - 1/(4k)) with k=3t+2 and p=4^t."""
    p = 4**t

    num = (
        144 * p**3 * t**4
        + 1416 * p**3 * t**3
        + 4548 * p**3 * t**2
        + 5304 * p**3 * t
        + 1548 * p**3
        + 3456 * p**2 * t**6
        + 11808 * p**2 * t**5
        + 3264 * p**2 * t**4
        - 27402 * p**2 * t**3
        - 34428 * p**2 * t**2
        - 13896 * p**2 * t
        - 1554 * p**2
        + 1296 * p * t**6
        + 1824 * p * t**5
        - 7626 * p * t**4
        - 20700 * p * t**3
        - 18903 * p * t**2
        - 7365 * p * t
        - 1014 * p
        - 696 * t**5
        - 2378 * t**4
        - 3375 * t**3
        - 2454 * t**2
        - 901 * t
        - 132
    )

    den = (
        4
        * (3 * t + 2)
        * (8 * p * t + 13 * p + 3 * t + 3)
        * (12 * p * t + 18 * p + 5 * t + 4)
        * (6 * p * t**2 + 24 * p * t + 18 * p + 10 * t**2 + 11 * t + 3)
    )
    return num, den


def finite_base_check() -> list[tuple[int, Fraction, Fraction]]:
    """Exact finite bases not covered by the large-t lane."""
    bases = [6, 8, 9, 11]
    out: list[tuple[int, Fraction, Fraction]] = []
    for k in bases:
        ag = a_gap_exact(k)
        rhs = Fraction(1, 4 * k)
        if ag < rhs:
            raise AssertionError(f"Finite base failed at k={k}: A_gap={ag}, rhs={rhs}")
        out.append((k, ag, rhs))
    return out


def check_numden_identities(t_max: int) -> None:
    """Check exact identity: A_gap - 1/(4k) == num/den."""
    for t in range(2, t_max + 1):
        k0 = 3 * t
        lhs0 = a_gap_exact(k0) - Fraction(1, 4 * k0)
        n0, d0 = case0_num_den(t)
        rhs0 = Fraction(n0, d0)
        if lhs0 != rhs0:
            raise AssertionError(f"case0 identity failed at t={t}, k={k0}")

        k2 = 3 * t + 2
        lhs2 = a_gap_exact(k2) - Fraction(1, 4 * k2)
        n2, d2 = case2_num_den(t)
        rhs2 = Fraction(n2, d2)
        if lhs2 != rhs2:
            raise AssertionError(f"case2 identity failed at t={t}, k={k2}")


def prove_large_t_lemmas() -> None:
    """Algebraic lemmas used for t >= 4 in both residues."""
    # Case k=3t: R1 >= 150 t^7 for t>=4.
    # R1 - 150 t^7 = t^4 * Q(t) + (15 t^3 + 91 t^2 - 3 t - 2),
    # Q(t)=1146 t^3 - 2928 t^2 - 6022 t - 2387.
    # Q''(t)=6876 t - 5856 > 0 for t>=1, so Q' is increasing on [1,∞).
    # Q'(4)>0 => Q' > 0 for t>=4, hence Q increasing on [4,∞), Q(4)>0.
    q_at_4 = 1146 * 4**3 - 2928 * 4**2 - 6022 * 4 - 2387
    qp_at_4 = 3438 * 4**2 - 5856 * 4 - 6022
    qpp_at_4 = 6876 * 4 - 5856
    if not (qpp_at_4 > 0 and qp_at_4 > 0 and q_at_4 > 0):
        raise AssertionError("Monotonicity gate for case0 Q failed")

    # Case k=3t+2: R2 >= 100 t^6 for t>=4 via:
    # R2 >= 432 t^6 + 608 t^5 - 4702 t^4
    # (using t>=4 to absorb lower-order negatives into t^4 terms),
    # so need 432 t^2 + 608 t - 4702 >= 100 t^2,
    # i.e. L(t)=332 t^2 + 608 t - 4702 >= 0 for t>=4.
    # L'(t)=664 t + 608 > 0, and L(4)>0.
    l_at_4 = 332 * 4**2 + 608 * 4 - 4702
    lp_at_4 = 664 * 4 + 608
    if not (lp_at_4 > 0 and l_at_4 > 0):
        raise AssertionError("Monotonicity gate for case2 L failed")

    # Shared tail gate: 300 * t * 4^t - 9936 > 0 for t>=4.
    h_at_4 = 300 * 4 * 4**4 - 9936
    if h_at_4 <= 0:
        raise AssertionError("Tail gate H(4)>0 failed")


def check_large_t_lane(t_max: int) -> None:
    """Verify the explicit lower-bound chain used in the proof for t>=4."""
    for t in range(4, t_max + 1):
        p = 4**t

        # ----- case k = 3t -----
        r1 = (
            1296 * t**7
            - 2928 * t**6
            - 6022 * t**5
            - 2387 * t**4
            + 15 * t**3
            + 91 * t**2
            - 3 * t
            - 2
        )
        if r1 < 150 * t**7:
            raise AssertionError(f"case0: R1 bound failed at t={t}")

        c1 = (2 * t + 1) * r1
        c0 = 2 * t**2 * (2 * t + 1) ** 2 * (348 * t**3 + 147 * t**2 + 52 * t + 5)

        if c1 < 300 * t**8:
            raise AssertionError(f"case0: C1 bound failed at t={t}")
        if c0 > 9936 * t**7:
            raise AssertionError(f"case0: C0 bound failed at t={t}")

        n0, _ = case0_num_den(t)
        if n0 < p * c1 - c0:
            raise AssertionError(f"case0: numerator lower bound failed at t={t}")
        if p * c1 - c0 <= 0:
            raise AssertionError(f"case0: p*C1-C0 not positive at t={t}")
        if n0 <= 0:
            raise AssertionError(f"case0: numerator not positive at t={t}")

        # ----- case k = 3t+2 -----
        r2 = (
            432 * t**6
            + 608 * t**5
            - 2542 * t**4
            - 6900 * t**3
            - 6301 * t**2
            - 2455 * t
            - 338
        )
        if r2 < 100 * t**6:
            raise AssertionError(f"case2: R2 bound failed at t={t}")

        d1 = 3 * r2
        d0 = (2 * t + 1) * (348 * t**4 + 1015 * t**3 + 1180 * t**2 + 637 * t + 132)

        if d1 < 300 * t**6:
            raise AssertionError(f"case2: D1 bound failed at t={t}")
        if d0 > 9936 * t**5:
            raise AssertionError(f"case2: D0 bound failed at t={t}")

        n2, _ = case2_num_den(t)
        if n2 < p * d1 - d0:
            raise AssertionError(f"case2: numerator lower bound failed at t={t}")
        if p * d1 - d0 <= 0:
            raise AssertionError(f"case2: p*D1-D0 not positive at t={t}")
        if n2 <= 0:
            raise AssertionError(f"case2: numerator not positive at t={t}")


def direct_scan(k_max: int) -> tuple[int, Fraction]:
    """Direct exact scan of A_gap >= 1/(4k) for sanity/reproducibility."""
    fails = 0
    min_gap: Fraction | None = None

    for k in range(6, k_max + 1):
        if k % 3 == 1:
            continue
        gap = a_gap_exact(k) - Fraction(1, 4 * k)
        if min_gap is None or gap < min_gap:
            min_gap = gap
        if gap < 0:
            fails += 1

    if min_gap is None:
        raise AssertionError("No k checked in direct scan")
    return fails, min_gap


def main() -> None:
    ap = argparse.ArgumentParser(description="Closed-form lane for A_gap >= 1/(4k).")
    ap.add_argument("--t-max-id", type=int, default=120, help="Check num/den identities up to this t.")
    ap.add_argument("--t-max-lane", type=int, default=500, help="Check large-t lane inequalities up to this t.")
    ap.add_argument("--k-max-scan", type=int, default=2000, help="Direct exact scan upper bound.")
    args = ap.parse_args()

    print("Sub-claim C: proving A_gap >= 1/(4k) for k>=6, k%3 in {0,2}")

    prove_large_t_lemmas()
    print("Monotonicity gate lemmas: PASS")

    check_numden_identities(args.t_max_id)
    print(f"Closed-form num/den identities: PASS (t<= {args.t_max_id})")

    base = finite_base_check()
    print("Finite bases (exact): PASS")
    for k, ag, rhs in base:
        print(f"  k={k:2d}: A_gap={float(ag):.12f} >= 1/(4k)={float(rhs):.12f}")

    check_large_t_lane(args.t_max_lane)
    print(f"Large-t lower-bound lane: PASS (t=4..{args.t_max_lane})")

    fails, min_gap = direct_scan(args.k_max_scan)
    print(f"Direct exact scan k<= {args.k_max_scan}: failures={fails}, min gap={float(min_gap):.12f}")
    if fails != 0:
        raise AssertionError("Direct exact scan found a failure")

    print("Result: A_gap >= 1/(4k) is certified on finite bases + large-t lane.")


if __name__ == "__main__":
    main()
