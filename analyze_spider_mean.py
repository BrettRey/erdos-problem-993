#!/usr/bin/env python3
"""Analyze the mean μ for the extremal spider family S(2^k).

For S(2^k), n = 2k+1:
  I(x) = (1+2x)^k + x(1+x)^k
  I(1) = 3^k + 2^k
  I'(1) = 2k·3^{k-1} + 2^k + k·2^{k-1}
  μ = I'(1)/I(1)

Asymptotically: μ → 2k/3 = (n-1)/3 < n/3.
Gap: n/3 - μ → 1/3.

Also check S(2^k, 1^1) (one extra unit arm).
"""

from fractions import Fraction


def spider_2k_exact(k):
    """Exact mean for S(2^k) using Fraction arithmetic."""
    I1 = 3**k + 2**k
    Ip1 = 2 * k * 3**(k - 1) + 2**k + k * 2**(k - 1)
    mu = Fraction(Ip1, I1)
    n = 2 * k + 1
    return n, mu


def spider_2k_1j_exact(k, j):
    """Exact mean for S(2^k, 1^j).

    I(x) = (1+2x)^k (1+x)^j + x(1+x)^k (1)^j
          = (1+2x)^k (1+x)^j + x(1+x)^k
    Actually for spider with hub:
    I = prod I(P_{a_i}) + x prod I(P_{a_i - 1})
    P_2: I = 1+2x, P_1: I = 1+x, P_0: I = 1
    So I(S(2^k,1^j)) = (1+2x)^k (1+x)^j + x·(1+x)^k · 1^j
                      = (1+2x)^k (1+x)^j + x(1+x)^k
    """
    # I(1) = 3^k * 2^j + 2^k
    I1 = 3**k * 2**j + 2**k

    # I'(x) = 2k(1+2x)^{k-1}(1+x)^j + j(1+2x)^k(1+x)^{j-1}
    #        + (1+x)^k + kx(1+x)^{k-1}
    # I'(1) = 2k·3^{k-1}·2^j + j·3^k·2^{j-1} + 2^k + k·2^{k-1}
    if j > 0:
        Ip1 = (2 * k * 3**(k-1) * 2**j
               + j * 3**k * 2**(j-1)
               + 2**k + k * 2**(k-1))
    else:
        Ip1 = 2 * k * 3**(k-1) + 2**k + k * 2**(k-1)

    mu = Fraction(Ip1, I1)
    n = 1 + 2 * k + j
    return n, mu


def main():
    print("SPIDER S(2^k) MEAN ANALYSIS", flush=True)
    print("=" * 70, flush=True)
    print(flush=True)

    print(f"{'k':>4s} {'n':>5s} {'μ':>12s} {'n/3':>10s} "
          f"{'n/3-μ':>10s} {'μ/(n/3)':>10s} {'mode_bound':>10s} "
          f"{'thr':>5s}", flush=True)
    print("-" * 70, flush=True)

    for k in list(range(2, 31)) + [50, 100, 200, 500]:
        n, mu = spider_2k_exact(k)
        n_third = Fraction(n, 3)
        gap = n_third - mu
        ratio = float(mu / n_third)
        mode_bound = int(mu) + 1  # ceil(μ) since mode ≤ ceil(μ)
        thr = n // 3 + 1

        print(f"{k:4d} {n:5d} {float(mu):12.6f} {float(n_third):10.4f} "
              f"{float(gap):10.6f} {ratio:10.6f} {mode_bound:10d} "
              f"{thr:5d}", flush=True)

    print(flush=True)
    print("PROOF SKETCH:", flush=True)
    print("-" * 70, flush=True)
    print(flush=True)
    print("For S(2^k), n = 2k+1:", flush=True)
    print("  I(1) = 3^k + 2^k", flush=True)
    print("  I'(1) = 2k·3^{k-1} + 2^k + k·2^{k-1}", flush=True)
    print("  μ = I'(1)/I(1)", flush=True)
    print(flush=True)
    print("  n/3 - μ = [3^k + 2^{k-1}(k-4)] / [3(3^k + 2^k)]", flush=True)
    print(flush=True)
    print("  For k ≥ 5: numerator = 3^k + 2^{k-1}(k-4) > 0", flush=True)
    print("  So μ < n/3, hence mode ≤ ceil(μ) ≤ floor(n/3)+1 = thr", flush=True)
    print(flush=True)

    # Verify the gap formula
    print("VERIFY GAP FORMULA: n/3 - μ = [3^k + 2^{k-1}(k-4)] / [3(3^k+2^k)]",
          flush=True)
    for k in [2, 3, 4, 5, 10, 20]:
        n, mu = spider_2k_exact(k)
        n_third = Fraction(n, 3)
        gap_actual = n_third - mu
        gap_formula = Fraction(3**k + 2**(k-1) * (k - 4),
                               3 * (3**k + 2**k))
        match = "✓" if gap_actual == gap_formula else "✗"
        print(f"  k={k:3d}: gap={float(gap_actual):.8f}, "
              f"formula={float(gap_formula):.8f} {match}", flush=True)

    print(flush=True)

    # Also check S(2^k, 1^1)
    print("SPIDER S(2^k, 1^1) (one extra unit arm):", flush=True)
    print(f"{'k':>4s} {'n':>5s} {'μ':>12s} {'n/3':>10s} "
          f"{'n/3-μ':>10s} {'ratio':>10s}", flush=True)
    for k in list(range(2, 21)) + [50, 100]:
        n, mu = spider_2k_1j_exact(k, 1)
        n_third = Fraction(n, 3)
        gap = n_third - mu
        ratio = float(mu / n_third)
        print(f"{k:4d} {n:5d} {float(mu):12.6f} {float(n_third):10.4f} "
              f"{float(gap):10.6f} {ratio:10.6f}", flush=True)

    print(flush=True)

    # Key question: is S(2^k) the global maximizer of μ/(n/3) among
    # d_leaf ≤ 1 trees? The exhaustive data says yes.
    # But what about S(3^k)? (arms of length 3)
    print("SPIDER S(3^k) (arms of length 3):", flush=True)
    print(f"{'k':>4s} {'n':>5s} {'μ':>12s} {'n/3':>10s} "
          f"{'ratio':>10s}", flush=True)

    for k in list(range(2, 21)) + [50, 100]:
        # I(P_3) = 1 + 3x + x^2, I(P_2) = 1 + 2x
        # I(S(3^k)) = (1+3x+x^2)^k + x(1+2x)^k
        # I(1) = 5^k + 3^k
        # I'(x) = k(3+2x)(1+3x+x^2)^{k-1} + (1+2x)^k + 2kx(1+2x)^{k-1}
        # I'(1) = 5k·5^{k-1} + 3^k + 2k·3^{k-1}
        #       = k·5^k + 3^k + 2k·3^{k-1}
        #       Wait, 3+2x at x=1 = 5. k·5·5^{k-1} = k·5^k. Yes.
        I1 = 5**k + 3**k
        Ip1 = k * 5**k + 3**k + 2 * k * 3**(k - 1)
        mu = Fraction(Ip1, I1)
        n = 3 * k + 1
        n_third = Fraction(n, 3)
        ratio = float(mu / n_third)
        print(f"{k:4d} {n:5d} {float(mu):12.6f} {float(n_third):10.4f} "
              f"{ratio:10.6f}", flush=True)

    print(flush=True)
    print("ASYMPTOTIC RATIOS:", flush=True)
    print("  S(2^k): μ/(n/3) → 2k/(2k+1) → 1 from below", flush=True)
    print("  S(3^k): μ/(n/3) → 3k·5^k / ((3k+1)·(5^k+3^k)/3)", flush=True)
    print("         ≈ 3k/(3k+1) · 3/(1+(3/5)^k) → 3·(1) = 3?", flush=True)
    print("  Wait: μ → k·5^k/5^k = k. n/3 = (3k+1)/3 ≈ k.", flush=True)
    print("  So ratio → k/k = 1 also. Gap?", flush=True)
    print(flush=True)

    # Compare gaps for S(2^k) vs S(3^k) at same n
    print("GAP COMPARISON at similar n:", flush=True)
    for n_target in [31, 61, 101, 201]:
        # S(2^k): n = 2k+1
        k2 = (n_target - 1) // 2
        n2, mu2 = spider_2k_exact(k2)

        # S(3^k): n = 3k+1
        k3 = (n_target - 1) // 3
        n3 = 3 * k3 + 1
        I1_3 = 5**k3 + 3**k3
        Ip1_3 = k3 * 5**k3 + 3**k3 + 2 * k3 * 3**(k3 - 1)
        mu3 = Fraction(Ip1_3, I1_3)

        print(f"  n≈{n_target}: S(2^{k2}) n={n2} μ/(n/3)="
              f"{float(mu2/Fraction(n2,3)):.6f}; "
              f"S(3^{k3}) n={n3} μ/(n/3)="
              f"{float(mu3/Fraction(n3,3)):.6f}", flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
