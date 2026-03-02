"""
Scan w_k for ALL k at the star+star extremal family.

For the support vertex with 1 leaf + star(a) + star(b):
  C = (1+x)^s  (s = a+b)
  A = (1+x)^{s+1} + x(1+x)^{b+1}
  B = (1+x)^{a+1} + x + x^2

Check w_k for all k from 0 to deg(C)-1.
"""

from fractions import Fraction
from math import comb


def _coeff(p, k):
    return p[k] if 0 <= k < len(p) else 0


def binomial_coeffs(n):
    """Coefficients of (1+x)^n."""
    return [comb(n, k) for k in range(n + 1)]


def compute_all_wk(a, b):
    """Compute w_k for all k at the star+star support vertex."""
    s = a + b

    # C = (1+x)^s
    C = binomial_coeffs(s)

    # A = (1+x)^{s+1} + x(1+x)^{b+1}
    A1 = binomial_coeffs(s + 1)
    A2_shifted = [0] + binomial_coeffs(b + 1)  # x(1+x)^{b+1}
    deg_A = max(len(A1), len(A2_shifted))
    A = [0] * deg_A
    for i in range(len(A1)):
        A[i] += A1[i]
    for i in range(len(A2_shifted)):
        A[i] += A2_shifted[i]

    # B = (1+x)^{a+1} + x + x^2
    B1 = binomial_coeffs(a + 1)
    deg_B = max(len(B1), 3)
    B = [0] * deg_B
    for i in range(len(B1)):
        B[i] += B1[i]
    B[1] += 1  # +x
    B[2] += 1  # +x^2

    results = []
    for k in range(len(C)):
        Ck = _coeff(C, k)
        Ckm = _coeff(C, k - 1)
        Ckp = _coeff(C, k + 1)
        Ak = _coeff(A, k)
        Akp = _coeff(A, k + 1)
        Bk = _coeff(B, k)
        Bkm = _coeff(B, k - 1)

        if Ck == 0:
            continue

        # Term1 (Karlin) = C_k * (A_{k+1}*C_k - A_k*C_{k+1})
        dk_AC = Akp * Ck - Ak * Ckp
        Term1 = Ck * dk_AC

        # Term2 = C_{k+1} * (B_k*C_{k-1} - B_{k-1}*C_k)
        Term2 = Ckp * (Bk * Ckm - Bkm * Ck)

        # Term3 = B_k * (C_k^2 - C_{k-1}*C_{k+1})
        ck_C = Ck * Ck - Ckm * Ckp  # LC defect of C
        Term3 = Bk * ck_C

        # Star margin = B_k*C_k - B_{k-1}*C_{k+1}
        star_margin = Bk * Ck - Bkm * Ckp

        # w_k * C_k = Term1 + Term2 + Term3
        wk_times_Ck = Term1 + Term2 + Term3
        wk = Fraction(wk_times_Ck, Ck)

        # Ratio (only when t23 < 0)
        t23 = Term2 + Term3
        ratio = None
        if t23 < 0 and Term1 > 0:
            ratio = Fraction(Term1, -t23)

        results.append({
            'k': k,
            'wk': wk,
            'Term1': Term1,
            'Term2': Term2,
            'Term3': Term3,
            't23': t23,
            'ratio': ratio,
            'star_margin': star_margin,
            'ck_C': ck_C,
        })

    return results


def main():
    print("=== All-k scan for star+star family ===\n")

    # First, check which k values have ⋆ failures
    print("--- Which k values have ⋆ failures? ---")
    for a in range(1, 20):
        for b in range(a, 20):
            results = compute_all_wk(a, b)
            failures = [r for r in results if r['star_margin'] < 0]
            if failures:
                fail_ks = [r['k'] for r in failures]
                ratios = [float(r['ratio']) if r['ratio'] else 'N/A' for r in failures]
                if fail_ks != [2]:  # Only report if not just k=2
                    print(f"  ({a},{b}) s={a+b}: ⋆ fails at k={fail_ks}")

    print("\n--- Detailed analysis for specific (a,b) ---")
    for a, b in [(1, 4), (5, 9), (10, 20), (20, 40), (50, 100)]:
        results = compute_all_wk(a, b)
        print(f"\n(a,b) = ({a},{b}), s={a+b}:")
        print(f"  {'k':>3s} {'w_k':>12s} {'T1':>12s} {'t23':>12s} {'ratio':>8s} {'⋆':>3s}")
        for r in results:
            if r['k'] > (a + b) // 2 + 2:
                break
            ratio_str = f"{float(r['ratio']):.3f}" if r['ratio'] else "---"
            star_str = "F" if r['star_margin'] < 0 else "+"
            print(f"  {r['k']:3d} {float(r['wk']):12.1f} {r['Term1']:12d} {r['t23']:12d} {ratio_str:>8s} {star_str:>3s}")

    # Comprehensive: check ALL (a,b) and ALL k, find minimum w_k
    print("\n=== Minimum w_k across all k and (a,b) ===")
    min_wk = None
    min_info = None
    neg_count = 0
    total = 0

    for a in range(1, 80):
        for b in range(a, 80):
            results = compute_all_wk(a, b)
            for r in results:
                total += 1
                if r['wk'] < 0:
                    neg_count += 1
                    print(f"  NEGATIVE: ({a},{b}) k={r['k']}: w_k={r['wk']}")
                if min_wk is None or r['wk'] < min_wk:
                    min_wk = r['wk']
                    min_info = (a, b, r['k'], r['ratio'])

    print(f"\nTotal checks: {total}")
    print(f"Negative w_k: {neg_count}")
    print(f"Minimum w_k: {float(min_wk):.4f} ({min_wk}) at (a,b)=({min_info[0]},{min_info[1]}), k={min_info[2]}")
    if min_info[3]:
        print(f"  ratio at that point: {float(min_info[3]):.6f}")

    # Find minimum w_k for each k separately
    print("\n=== Minimum w_k by k index ===")
    for target_k in range(0, 10):
        min_wk_k = None
        min_info_k = None
        for a in range(1, 80):
            for b in range(a, 80):
                results = compute_all_wk(a, b)
                for r in results:
                    if r['k'] == target_k:
                        if min_wk_k is None or r['wk'] < min_wk_k:
                            min_wk_k = r['wk']
                            min_info_k = (a, b, r['ratio'])
        if min_wk_k is not None:
            ratio_str = f"{float(min_info_k[2]):.3f}" if min_info_k[2] else "---"
            print(f"  k={target_k}: min w_k = {float(min_wk_k):.4f} at ({min_info_k[0]},{min_info_k[1]}), ratio={ratio_str}")


if __name__ == '__main__':
    main()
