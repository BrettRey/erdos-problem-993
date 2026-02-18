#!/usr/bin/env python3
"""Attempt to close Sub-claim C analytically.

Sub-claim C:
  For k >= 3 with k % 3 in {0,2},
    margin(k,1) < margin(k,0),
for mixed spiders S(2^k,1^j), where margin(k,j) is the tie-fugacity margin.

This script combines:
  1) the proved E1 <= 0 lane (see notes/subclaim_c_E1_le0_proof_2026-02-18.md),
  2) a new analytic lower bound A_gap > 1/(4k) for k >= 14,
  3) a new analytic upper bound -E0 < 1/(4k) for k >= 14,
  4) exact finite checks on the remaining small k in {3,5,6,8,9,11,12}.

Notation:
  A_gap = A0 - A1,
  diff  = margin(k,1) - margin(k,0) = -A_gap + (E1 - E0).
"""

from __future__ import annotations

from fractions import Fraction

from prove_subclaim_c_algebra import lambda_j0_closed, lambda_j1_closed, terms_j0, terms_j1


def agap_minus_quarterk_res0_num_den(t: int) -> tuple[int, int]:
    """Exact numerator/denominator for D0 = A_gap - 1/(4k), k=3t."""
    p = 4**t
    c3 = (t + 2) * (3 * t + 1) * (24 * t**4 + 64 * t**3 + 27 * t**2 - 17 * t - 6)
    c2 = (3 * t + 1) * (
        1152 * t**7
        - 480 * t**6
        - 3928 * t**5
        - 2484 * t**4
        + 120 * t**3
        + 237 * t**2
        - 11 * t
        - 10
    )
    c1 = (2 * t + 1) * (
        1296 * t**7
        - 2928 * t**6
        - 6022 * t**5
        - 2387 * t**4
        + 15 * t**3
        + 91 * t**2
        - 3 * t
        - 2
    )
    c0 = -2 * t**2 * (2 * t + 1) ** 2 * (348 * t**3 + 147 * t**2 + 52 * t + 5)
    num = p**3 * c3 + p**2 * c2 + p * c1 + c0
    den = (
        12
        * t
        * (3 * p * t**2 + 7 * p * t + 2 * p + 10 * t**2)
        * (8 * p * t**2 + 9 * p * t + 2 * p + 6 * t**2 + 5 * t + 1)
        * (12 * p * t**2 + 13 * p * t + 3 * p + 10 * t**2 + 7 * t + 1)
    )
    return num, den


def agap_minus_quarterk_res2_num_den(t: int) -> tuple[int, int]:
    """Exact numerator/denominator for D2 = A_gap - 1/(4k), k=3t+2."""
    p = 4**t
    d3 = 12 * (t + 3) * (12 * t**3 + 82 * t**2 + 133 * t + 43)
    d2 = 6 * (
        576 * t**6
        + 1968 * t**5
        + 544 * t**4
        - 4567 * t**3
        - 5738 * t**2
        - 2316 * t
        - 259
    )
    d1 = 3 * (
        432 * t**6
        + 608 * t**5
        - 2542 * t**4
        - 6900 * t**3
        - 6301 * t**2
        - 2455 * t
        - 338
    )
    d0 = -(2 * t + 1) * (348 * t**4 + 1015 * t**3 + 1180 * t**2 + 637 * t + 132)
    num = p**3 * d3 + p**2 * d2 + p * d1 + d0
    den = (
        4
        * (3 * t + 2)
        * (8 * p * t + 13 * p + 3 * t + 3)
        * (12 * p * t + 18 * p + 5 * t + 4)
        * (6 * p * t**2 + 24 * p * t + 18 * p + 10 * t**2 + 11 * t + 3)
    )
    return num, den


def direct_agap_minus_quarterk(k: int) -> Fraction:
    lam0 = lambda_j0_closed(k)
    lam1 = lambda_j1_closed(k)
    a0 = terms_j0(k, lam0)[0]
    a1 = terms_j1(k, lam1)[0]
    return (a0 - a1) - Fraction(1, 4 * k)


def direct_diff(k: int) -> Fraction:
    """diff = margin(k,1)-margin(k,0)."""
    m0 = terms_j0(k, lambda_j0_closed(k))[3]
    m1 = terms_j1(k, lambda_j1_closed(k))[3]
    return m1 - m0


def validate_transcribed_fraction_identities() -> None:
    # Residue 0: k=3t, t>=2.
    for t in range(2, 61):
        k = 3 * t
        num, den = agap_minus_quarterk_res0_num_den(t)
        lhs = Fraction(num, den)
        rhs = direct_agap_minus_quarterk(k)
        if lhs != rhs:
            raise AssertionError(f"res0 identity mismatch at t={t}: lhs={lhs}, rhs={rhs}")

    # Residue 2: k=3t+2, t>=2 (k>=8).
    for t in range(2, 61):
        k = 3 * t + 2
        num, den = agap_minus_quarterk_res2_num_den(t)
        lhs = Fraction(num, den)
        rhs = direct_agap_minus_quarterk(k)
        if lhs != rhs:
            raise AssertionError(f"res2 identity mismatch at t={t}: lhs={lhs}, rhs={rhs}")


def prove_agap_gt_quarterk_for_k_ge_14() -> None:
    # ------------------------------------------------------------------
    # Residue 0: k=3t, t>=4
    # D0 numerator: p^3*c3 + p^2*c2 + p*c1 + c0, p=4^t
    #
    # For t>=4:
    #   c2 > 0, c1 > 0,
    # so N0 > p^3*c3 + c0.
    #
    # Lower bounds:
    #   c3 >= 72 t^6
    #   |c0| <= 9936 t^7
    # hence N0 > 72 t^6 (64^t - 138 t).
    # ------------------------------------------------------------------

    # c2 positivity via cubic lower envelope:
    # 1152 t^7 - 480 t^6 - 3928 t^5 - 2484 t^4 >= t^4 * g2(t),
    # g2(t)=1152 t^3 - 480 t^2 - 3928 t - 2484.
    def g2(t: int) -> int:
        return 1152 * t**3 - 480 * t**2 - 3928 * t - 2484

    def g2p(t: int) -> int:
        return 3456 * t**2 - 960 * t - 3928

    if not (g2p(4) > 0 and g2(4) > 0):
        raise AssertionError("Failed g2 positivity base at t=4")

    # g2' is increasing for t>=4 (g2''=6912t-960 > 0), so g2'>0 and g2 increasing.
    if not ((6912 * 4 - 960) > 0):
        raise AssertionError("Failed g2'' positivity base at t=4")

    # c1 positivity via cubic lower envelope:
    # 1296 t^7 - 2928 t^6 - 6022 t^5 - 2387 t^4 >= t^4 * g1(t),
    # g1(t)=1296 t^3 - 2928 t^2 - 6022 t - 2387.
    def g1(t: int) -> int:
        return 1296 * t**3 - 2928 * t**2 - 6022 * t - 2387

    def g1p(t: int) -> int:
        return 3888 * t**2 - 5856 * t - 6022

    if not (g1p(4) > 0 and g1(4) > 0):
        raise AssertionError("Failed g1 positivity base at t=4")

    # g1' increasing for t>=4 (g1''=7776t-5856 > 0), so g1'>0 and g1 increasing.
    if not ((7776 * 4 - 5856) > 0):
        raise AssertionError("Failed g1'' positivity base at t=4")

    # 64^t / t is strictly increasing for t>=1:
    # ratio at step t->t+1 is 64t/(t+1) > 1.
    if not ((64 * 4) > (4 + 1)):
        raise AssertionError("Failed monotone-ratio base for 64^t/t")
    if not ((64**4) > 138 * 4):
        raise AssertionError("Failed base inequality 64^t > 138t at t=4")

    # ------------------------------------------------------------------
    # Residue 2: k=3t+2, t>=4
    # D2 numerator: p^3*d3 + p^2*d2 + p*d1 + d0.
    #
    # For t>=4:
    #   d2 > 0, d1 > 0,
    # so N2 > p^3*d3 + d0.
    #
    # Lower bounds:
    #   d3 >= 144 t^4
    #   |d0| <= 9936 t^5
    # hence N2 > 144 t^4 (64^t - 69 t).
    # ------------------------------------------------------------------

    # d2 increasing for t>=4 (derivative lower bound >0).
    def d2(t: int) -> int:
        return (
            576 * t**6
            + 1968 * t**5
            + 544 * t**4
            - 4567 * t**3
            - 5738 * t**2
            - 2316 * t
            - 259
        )

    def d2p_lb(t: int) -> Fraction:
        return (
            Fraction(3456) * t**5
            - Fraction(2880, 4) * t**5
            - Fraction(19640, 16) * t**5
            - Fraction(9936, 64) * t**5
            - Fraction(13701, 256) * t**5
            - Fraction(11476, 1024) * t**5
            - Fraction(2316, 4096) * t**5
        )

    if not (d2p_lb(4) > 0 and d2(4) > 0):
        raise AssertionError("Failed d2 positivity base at t=4")

    # d1 increasing for t>=4 (derivative lower bound >0).
    def d1(t: int) -> int:
        return (
            432 * t**6
            + 608 * t**5
            - 2542 * t**4
            - 6900 * t**3
            - 6301 * t**2
            - 2455 * t
            - 338
        )

    def d1p_lb(t: int) -> Fraction:
        return (
            Fraction(2592) * t**5
            - Fraction(10168, 16) * t**5
            - Fraction(20700, 64) * t**5
            - Fraction(12602, 256) * t**5
            - Fraction(2455, 1024) * t**5
        )

    if not (d1p_lb(4) > 0 and d1(4) > 0):
        raise AssertionError("Failed d1 positivity base at t=4")

    # 64^t > 69t at t=4, and 64^t / t increasing => for all t>=4.
    if not ((64**4) > 69 * 4):
        raise AssertionError("Failed base inequality 64^t > 69t at t=4")

    # ------------------------------------------------------------------
    # Finite exact checks that bridge to t>=4 arguments.
    # ------------------------------------------------------------------
    # Residue 0 small t not covered above: t=2,3  (k=6,9)
    for t in (2, 3):
        num, den = agap_minus_quarterk_res0_num_den(t)
        if Fraction(num, den) <= 0:
            raise AssertionError(f"A_gap-1/(4k) failed at res0 t={t}")

    # Residue 2 small t not covered above: t=2,3  (k=8,11)
    for t in (2, 3):
        num, den = agap_minus_quarterk_res2_num_den(t)
        if Fraction(num, den) <= 0:
            raise AssertionError(f"A_gap-1/(4k) failed at res2 t={t}")


def prove_minus_e0_lt_quarterk_for_k_ge_14() -> None:
    # For j=0 branch:
    #   -E0 = r0(A0-B0)/(1+r0) < r0 * (k/6),
    # using lambda<=1 => A0-B0 = k*lambda/((1+lambda)(1+2lambda)) - 1 < k/6.
    #
    # Also lambda0 >= a_k:
    #   k=3t   : a_k = k/(k+3)
    #   k=3t+2 : a_k = (k-1)/(k+2)
    #
    # so r0 <= s_k^k where s_k = (1+a_k)/(1+2a_k):
    #   k=3t   : s_k = (2k+3)/(3k+3)
    #   k=3t+2 : s_k = (2k+1)/(3k)
    #
    # For k>=14 (k%3 in {0,2}), both satisfy s_k <= 29/42.
    # Hence -E0 < (k/6)*(29/42)^k.
    #
    # Let f(k)=k^2*(29/42)^k. Then
    #   f(k+1)/f(k) = (29/42)*((k+1)/k)^2 <= (29/42)*(15/14)^2 < 1  for k>=14,
    # so f decreases on k>=14.
    # Check f(14) < 3/2. Therefore f(k) < 3/2 for all k>=14.
    # Conclude -E0 < (1/(6k))*f(k) < (1/(6k))*(3/2)=1/(4k).

    s = Fraction(29, 42)
    ratio_bound = s * Fraction(15 * 15, 14 * 14)
    if not (ratio_bound < 1):
        raise AssertionError("Failed monotone ratio bound for f(k)")

    f14 = Fraction(14 * 14) * (s**14)
    if not (f14 < Fraction(3, 2)):
        raise AssertionError("Failed base bound f(14) < 3/2")

    # Additional finite confirmation on a moderate window.
    for k in range(14, 401):
        if k % 3 == 1:
            continue
        if k % 3 == 0:
            sk = Fraction(2 * k + 3, 3 * k + 3)
        else:
            sk = Fraction(2 * k + 1, 3 * k)
        if sk > s:
            raise AssertionError(f"s_k <= 29/42 failed at k={k}")
        if not (Fraction(k * k) * (s**k) < Fraction(3, 2)):
            raise AssertionError(f"k^2 s^k < 3/2 failed at k={k}")


def finite_small_k_diff_checks() -> None:
    # Remaining cases not covered by k>=14 analytic lane.
    small_k = [3, 5, 6, 8, 9, 11, 12]
    for k in small_k:
        d = direct_diff(k)
        if d >= 0:
            raise AssertionError(f"Sub-claim C finite base failed at k={k}: diff={d}")


def sanity_scan(k_max: int = 500) -> None:
    # Computational sanity: direct target inequality.
    for k in range(3, k_max + 1):
        if k % 3 == 1:
            continue
        d = direct_diff(k)
        if d >= 0:
            raise AssertionError(f"Sanity scan failure at k={k}: diff={d}")


def main() -> None:
    print("Sub-claim C complete lane (analytic + finite exact checks)")
    print("Target: margin(k,1) < margin(k,0) for k>=3, k%3 in {0,2}")
    print()

    print("1) Validating transcribed A_gap-1/(4k) fraction identities...")
    validate_transcribed_fraction_identities()
    print("   PASS")
    print()

    print("2) Proving A_gap > 1/(4k) for k>=14...")
    prove_agap_gt_quarterk_for_k_ge_14()
    print("   PASS")
    print()

    print("3) Proving -E0 < 1/(4k) for k>=14...")
    prove_minus_e0_lt_quarterk_for_k_ge_14()
    print("   PASS")
    print()

    print("4) Finite exact base checks for remaining small k...")
    finite_small_k_diff_checks()
    print("   PASS")
    print()

    print("5) Optional sanity scan (direct diff) up to k=500...")
    sanity_scan(500)
    print("   PASS")
    print()

    print("Conclusion:")
    print("  Combined with the proved E1<=0 theorem, this closes Sub-claim C.")


if __name__ == "__main__":
    main()
