"""
Check if the correction/c_k ratio stays bounded below 1 as a,b grow.

The key question: does max_k |correction_k|/c_k converge to a limit < 1?
If yes, then w_k >= (1 - limit) * c_k > 0 for all k, proving w_k >= 0.
"""

from math import comb
from fractions import Fraction


def max_correction_ratio(a, b):
    """Find max |correction_k|/c_k over all k for star(a)+star(b)."""
    s = a + b
    worst_ratio = Fraction(0)
    worst_k = -1

    for k in range(0, s):
        Ck = comb(s, k)
        Ckm = comb(s, k - 1) if k >= 1 else 0
        Ckp = comb(s, k + 1)

        ck = Ck * Ck - Ckm * Ckp

        if ck <= 0:
            continue

        # A = (1+x)^{s+1} + x(1+x)^{b+1}
        Ak = comb(s + 1, k) + (comb(b + 1, k - 1) if k >= 1 else 0)
        Akp = comb(s + 1, k + 1) + comb(b + 1, k)

        # B = (1+x)^{a+1} + x + x^2
        Bk = comb(a + 1, k) + (1 if k == 1 else 0) + (1 if k == 2 else 0)
        Bkm = (comb(a + 1, k - 1) if k >= 1 else 0) + (1 if k == 2 else 0) + (1 if k == 3 else 0)

        dk_AC = Akp * Ck - Ak * Ckp
        star_margin = Bk * Ck - Bkm * Ckp

        wk = dk_AC + star_margin
        correction = dk_AC - ck  # dk_AC = ck + delta_k, so correction = delta_k + star_margin... no

        # Actually: dk_AC = ck + delta_k (where delta_k = dk(D,C))
        # So correction = delta_k + star_margin
        # But dk_AC + star_margin = ck + delta_k + star_margin = ck + correction
        # correction = wk - ck
        correction = wk - ck

        if correction < 0:
            ratio = Fraction(-correction, ck)
            if ratio > worst_ratio:
                worst_ratio = ratio
                worst_k = k

    return worst_ratio, worst_k


def main():
    print("=== Correction ratio convergence as a,b grow ===")
    print(f"{'a':>6s} {'b':>6s} {'s':>6s} {'max|corr|/ck':>14s} {'k':>4s} {'1-ratio':>12s}")
    print("-" * 55)

    # Scan along the extremal line a/b ≈ √3 - 1 ≈ 0.732
    for s in [10, 20, 30, 50, 70, 100, 150, 200, 300, 500]:
        # Optimal a ≈ s * (1 - 1/√3)
        from math import sqrt
        a_opt = int(round(s * (1 - 1/sqrt(3))))
        b_opt = s - a_opt
        if a_opt < 1:
            a_opt = 1
            b_opt = s - 1

        ratio, k = max_correction_ratio(a_opt, b_opt)
        print(f"{a_opt:6d} {b_opt:6d} {s:6d} {float(ratio):14.8f} {k:4d} {1-float(ratio):12.8f}")

    # Also scan balanced a = b
    print("\n=== Balanced a = b ===")
    for a in [5, 10, 20, 30, 50, 70, 100, 150, 200]:
        ratio, k = max_correction_ratio(a, a)
        print(f"  a=b={a:4d}: max|corr|/ck = {float(ratio):.8f} at k={k}, gap = {1-float(ratio):.8f}")

    # Scan along a=1 (most asymmetric)
    print("\n=== Most asymmetric: a=1 ===")
    for b in [5, 10, 20, 50, 100, 200, 500]:
        ratio, k = max_correction_ratio(1, b)
        print(f"  a=1, b={b:4d}: max|corr|/ck = {float(ratio):.8f} at k={k}, gap = {1-float(ratio):.8f}")

    # Scan along a=2
    print("\n=== a=2 ===")
    for b in [5, 10, 20, 50, 100, 200]:
        ratio, k = max_correction_ratio(2, b)
        print(f"  a=2, b={b:4d}: max|corr|/ck = {float(ratio):.8f} at k={k}, gap = {1-float(ratio):.8f}")

    print("\n=== Does the gap vanish? ===")
    # Check: does 1 - ratio → 0 as s → ∞?
    # If gap → 0, the ratio approaches 1 and we can't use this method.
    # If gap → const > 0, we have an algebraic proof.


if __name__ == '__main__':
    main()
