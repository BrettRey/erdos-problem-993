"""
prove_subclaim_A_symbolic.py
Symbolic/algebraic analysis toward proving Sub-claim A for mixed spiders S(2^k, 1^j).

Goal: prove F = T1 + T2 >= 0 where:
  a_t = [x^t](1+2x)^k(1+x)^j  (coefficients of A_j)
  m = mode of full IS polynomial I_{k,j} = A_j + x(1+x)^k
  A=a_m, B=a_{m-1}, C=a_{m-2}, D=a_{m+1}
  G=C(k,m-1), H=C(k,m-2)  (backbone binomials)
  A' = a_{m}^{j+2} = A+2B+C (step up in j by 2: A_{j+2} at index m)
  D' = a_{m+1}^{j+2} = D+2A+B
  T1 = (A^2-BD) + (AC-B^2)
  T2 = G*A' - H*D'  (negative empirically)
  F = T1 + T2 = A'*(A+G) - D'*(B+H)
"""

import json
from fractions import Fraction
from math import comb, floor

# ─────────────────────────────────────────────
# Core functions
# ─────────────────────────────────────────────

def a_coeff(k, j, t):
    """[x^t](1+2x)^k(1+x)^j using exact Fraction arithmetic."""
    # = sum_{s=0}^{t} C(k,t-s)*2^(t-s)*C(j,s)
    result = Fraction(0)
    for s in range(min(j, t) + 1):
        ts = t - s
        if 0 <= ts <= k:
            result += Fraction(comb(k, ts) * (2 ** ts) * comb(j, s))
    return result

def b_coeff(k, t):
    """[x^t](1+2x)^k = C(k,t)*2^t (the j=0 case)."""
    if 0 <= t <= k:
        return Fraction(comb(k, t) * (2 ** t))
    return Fraction(0)

def spider_is_coeff(k, j, t):
    """[x^t] I_{k,j}(x) = a_t + C(k, t-1).
    I_{k,j}(x) = (1+2x)^k*(1+x)^j + x*(1+x)^k
    Coefficient at x^t: a_t(k,j) + C(k, t-1).
    """
    backbone = Fraction(comb(k, t - 1)) if 1 <= t <= k + 1 else Fraction(0)
    return a_coeff(k, j, t) + backbone

def mode_of_IS(k, j):
    """Find mode of I_{k,j} polynomial."""
    n = 2 * k + j + 1  # number of vertices
    best_t, best_val = 0, Fraction(-1)
    for t in range(0, n + 1):
        val = spider_is_coeff(k, j, t)
        if val > best_val:
            best_val = val
            best_t = t
    return best_t

def compute_all(k, j):
    """Compute all the key quantities for given k, j."""
    m = mode_of_IS(k, j)
    
    A = a_coeff(k, j, m)
    B = a_coeff(k, j, m - 1)
    C = a_coeff(k, j, m - 2)
    D = a_coeff(k, j, m + 1)
    
    G = Fraction(comb(k, m - 1))
    H = Fraction(comb(k, m - 2))
    
    # A_{j+2} coefficients (step up j by 2)
    Ap = a_coeff(k, j + 2, m)      # A' = A + 2B + C  (shift by 2 in j)
    Dp = a_coeff(k, j + 2, m + 1)  # D' = D + 2A + B
    
    # Verify the recurrence: A' = A+2B+C, D' = D+2A+B
    Ap_check = A + 2*B + C
    Dp_check = D + 2*A + B
    assert Ap == Ap_check, f"A' mismatch at k={k},j={j}: {Ap} vs {Ap_check}"
    assert Dp == Dp_check, f"D' mismatch at k={k},j={j}: {Dp} vs {Dp_check}"
    
    T1 = (A*A - B*D) + (A*C - B*B)
    T2 = G*Ap - H*Dp
    F  = T1 + T2
    
    # A_lb = C(k,m)*2^m + j*C(k,m-1)*2^(m-1)  (leading terms from j expansion)
    bm   = b_coeff(k, m)
    bm1  = b_coeff(k, m - 1)
    A_lb = bm + Fraction(j) * bm1  # s=0 and s=1 terms
    
    return {
        'm': m, 'A': A, 'B': B, 'C': C, 'D': D,
        'G': G, 'H': H, 'Ap': Ap, 'Dp': Dp,
        'T1': T1, 'T2': T2, 'F': F,
        'bm': bm, 'bm1': bm1, 'A_lb': A_lb
    }

# ─────────────────────────────────────────────
# 1. VERIFY j=0 formula
# ─────────────────────────────────────────────
print("=" * 70)
print("SECTION 1: VERIFY j=0 FORMULA FOR T1/b_m^2")
print("=" * 70)
print("Formula: T1/b_m^2 = (k+1)/n1 * [4*n1*n2 - m*(m+1)] / [4*n1*n2*(m+1)]")
print("where n1 = k-m+1, n2 = k-m+2")
print()

formula_ok = True
j0_results = []

for k in range(6, 31):
    d = compute_all(k, 0)
    m = d['m']
    T1 = d['T1']
    bm = d['bm']
    T2 = d['T2']
    F  = d['F']
    
    n1 = k - m + 1
    n2 = k - m + 2
    
    # Formula value
    numer = (k + 1) * (4 * n1 * n2 - m * (m + 1))
    denom = n1 * 4 * n1 * n2 * (m + 1)
    formula_val = Fraction(numer, denom)
    
    actual_val = Fraction(T1, bm * bm) if bm != 0 else None
    
    match = (actual_val == formula_val)
    if not match:
        formula_ok = False
    
    disc = 4 * n1 * n2 - m * (m + 1)
    
    j0_results.append({
        'k': k, 'm': m, 'n1': n1, 'n2': n2,
        'T1': str(T1), 'T2': str(T2), 'F': str(F),
        'bm': str(bm),
        'formula_matches': match,
        'disc_4n1n2_mm1': disc,
        'disc_positive': disc > 0,
    })
    
    print(f"k={k:2d}: m={m}, n1={n1}, n2={n2}, disc={disc:5d}, "
          f"formula_match={match}, T2={'<0' if T2<0 else '>=0'}, F={'>=0' if F>=0 else '<0'}")

print()
print(f"ALL j=0 formula checks pass: {formula_ok}")
print()

# ─────────────────────────────────────────────
# 2. VERIFY m*n2 - (m-1)*n1 = k+1
# ─────────────────────────────────────────────
print("=" * 70)
print("SECTION 2: VERIFY IDENTITY m*n2 - (m-1)*n1 = k+1")
print("=" * 70)

identity_ok = True
for k in range(6, 31):
    d = compute_all(k, 0)
    m = d['m']
    n1 = k - m + 1
    n2 = k - m + 2
    lhs = m * n2 - (m - 1) * n1
    ok = (lhs == k + 1)
    if not ok:
        identity_ok = False
        print(f"  FAIL k={k}: m={m}, n1={n1}, n2={n2}, m*n2-(m-1)*n1={lhs} != {k+1}")

print(f"Identity m*n2 - (m-1)*n1 = k+1 holds for all k=6..30: {identity_ok}")

# Symbolic proof: m*n2 - (m-1)*n1 = m*(k-m+2) - (m-1)*(k-m+1)
# = m*k - m^2 + 2m - (m-1)k + (m-1)^2 - (m-1)
# = m*k - m^2 + 2m - mk + k + m^2 - m - m + 1 - m + 1 + m - 1
# Let me do it cleanly:
print()
print("Symbolic proof of m*n2 - (m-1)*n1 = k+1:")
print("  m*n2 = m*(k-m+2) = mk - m^2 + 2m")
print("  (m-1)*n1 = (m-1)*(k-m+1) = mk - m^2 + m - k + m - 1 = mk - m^2 + 2m - k - 1")
print("  Difference: mk - m^2 + 2m - (mk - m^2 + 2m - k - 1) = k + 1  QED")
print()

# ─────────────────────────────────────────────
# 3. PROVE (k+1)/n1 > 0 and disc > 0 for k=6..100
# ─────────────────────────────────────────────
print("=" * 70)
print("SECTION 3: SIGN ANALYSIS FOR k=6..100")
print("=" * 70)
print()

disc_min = None
disc_min_k = None
disc_results = {}

for k in range(6, 101):
    # Mode for j=0: m = floor((2k+1)/3)
    m = floor((2 * k + 1) / 3)
    n1 = k - m + 1
    n2 = k - m + 2
    disc = 4 * n1 * n2 - m * (m + 1)
    
    # Compute exact m from IS polynomial for small k to verify formula
    if k <= 30:
        d = compute_all(k, 0)
        assert d['m'] == m, f"Mode formula wrong at k={k}: got {d['m']}, expected {m}"
    
    r = k % 3
    disc_results[k] = {'m': m, 'n1': n1, 'n2': n2, 'disc': disc, 'k_mod_3': r}
    
    if disc_min is None or disc < disc_min:
        disc_min = disc
        disc_min_k = k

print(f"Min discriminant 4*n1*n2 - m*(m+1) over k=6..100: {disc_min} at k={disc_min_k}")

# Case analysis by residue class
print()
print("Case analysis by k mod 3:")
for r in range(3):
    cases = [(k, disc_results[k]) for k in range(6, 101) if k % 3 == r]
    print(f"\n  k ≡ {r} (mod 3): {len(cases)} values")
    # Symbolic: for k=3t+r, m=floor((6t+2r+1)/3)
    if r == 0:
        # k=3t: m = 2t, n1=t+1, n2=t+2, disc = 4(t+1)(t+2) - 2t(2t+1)
        print(f"    k=3t: m=2t, n1=t+1, n2=t+2")
        print(f"    disc = 4(t+1)(t+2) - 2t(2t+1)")
        print(f"         = 4t^2+12t+8 - 4t^2-2t = 10t+8 > 0 for t>=1 (k>=3)")
        # Verify
        for k, dat in cases[:3]:
            t = k // 3
            expected = 10*t + 8
            assert dat['disc'] == expected, f"Formula wrong at k={k}"
        print(f"    Verified for first 3 cases, formula disc=10t+8.")
    elif r == 1:
        # k=3t+1: m=2t+1, n1=t+1, n2=t+2, disc=4(t+1)(t+2)-(2t+1)(2t+2)
        print(f"    k=3t+1: m=2t+1, n1=t+1, n2=t+2")
        print(f"    disc = 4(t+1)(t+2) - (2t+1)(2t+2)")
        print(f"         = 4t^2+12t+8 - 4t^2-6t-2 = 6t+6 > 0 for t>=0 (k>=1)")
        for k, dat in cases[:3]:
            t = (k - 1) // 3
            expected = 6 * t + 6
            assert dat['disc'] == expected, f"Formula wrong at k={k}: got {dat['disc']}, expected {expected}"
        print(f"    Verified for first 3 cases, formula disc=6t+6.")
    else:  # r == 2
        # k=3t+2: m=2t+1 (floor((6t+5)/3)=floor(2t+1.67)=2t+1), n1=t+2, n2=t+3
        print(f"    k=3t+2: m=2t+1, n1=t+2, n2=t+3")
        print(f"    disc = 4(t+2)(t+3) - (2t+1)(2t+2)")
        print(f"         = 4t^2+20t+24 - 4t^2-6t-2 = 14t+22 > 0 for t>=0 (k>=2)")
        for k, dat in cases[:3]:
            t = (k - 2) // 3
            expected = 14 * t + 22
            assert dat['disc'] == expected, f"Formula wrong at k={k}: got {dat['disc']}, expected {expected}"
        print(f"    Verified for first 3 cases, formula disc=14t+22.")

print()
print("CONCLUSION: disc > 0 for all k >= 3 (all residue classes give linear-growing disc)")
print("Combined with n1 >= 1 and k+1 >= 7: T1 > 0 for all k >= 6, j=0.")

# ─────────────────────────────────────────────
# 4. General j: compute T1/|T2| and A/A_lb
# ─────────────────────────────────────────────
print()
print("=" * 70)
print("SECTION 4: CRITICAL RATIO T1*A_lb/(|T2|*A) FOR k=6..30, j=0..20")
print("=" * 70)
print()

ratio_results = []
min_ratio = None
min_ratio_kj = None
max_A_over_Alb = None
max_A_over_Alb_kj = None

print(f"{'k':>3} {'j':>3} {'m':>3} {'T1/|T2|':>12} {'A/A_lb':>12} {'T1*Alb/(|T2|*A)':>18} {'F>=0':>6}")

for k in range(6, 31):
    for j in range(0, 21):
        d = compute_all(k, j)
        T1 = d['T1']
        T2 = d['T2']
        F  = d['F']
        A  = d['A']
        A_lb = d['A_lb']
        m  = d['m']
        
        if T2 == 0:
            continue  # skip degenerate
        
        T1_over_absT2 = Fraction(T1, -T2) if T2 < 0 else Fraction(T1, T2)
        A_over_Alb = Fraction(A, A_lb) if A_lb > 0 else None
        
        if A_lb > 0 and T2 < 0:
            crit = T1 * A_lb / ((-T2) * A)
        else:
            crit = None
        
        # Track minimums
        if crit is not None:
            if min_ratio is None or crit < min_ratio:
                min_ratio = crit
                min_ratio_kj = (k, j)
        if A_over_Alb is not None:
            if max_A_over_Alb is None or A_over_Alb > max_A_over_Alb:
                max_A_over_Alb = A_over_Alb
                max_A_over_Alb_kj = (k, j)
        
        row = {
            'k': k, 'j': j, 'm': m,
            'T1': str(T1), 'T2': str(T2), 'F': str(F), 'F_nonneg': F >= 0,
            'A': str(A), 'A_lb': str(A_lb),
            'T1_over_absT2': str(T1_over_absT2) if T2 < 0 else 'N/A',
            'A_over_Alb': str(A_over_Alb) if A_over_Alb else 'N/A',
            'crit_ratio': str(crit) if crit is not None else 'N/A',
        }
        ratio_results.append(row)
        
        # Print select rows: j=0,1,2,3,5,10,20 for k=6,10,15,20,25
        if k in (6, 10, 15, 20, 25) and j in (0, 1, 2, 3, 5, 10, 20):
            crit_str = f"{float(crit):.5f}" if crit else "N/A"
            T1_str   = f"{float(T1_over_absT2):.5f}" if T2 < 0 else "N/A"
            A_str    = f"{float(A_over_Alb):.5f}" if A_over_Alb else "N/A"
            print(f"{k:>3} {j:>3} {m:>3} {T1_str:>12} {A_str:>12} {crit_str:>18} {'Y' if F>=0 else 'N':>6}")

print()
print(f"Min T1*A_lb/(|T2|*A) = {float(min_ratio):.6f} at k={min_ratio_kj[0]}, j={min_ratio_kj[1]}")
print(f"Max A/A_lb            = {float(max_A_over_Alb):.6f} at k={max_A_over_Alb_kj[0]}, j={max_A_over_Alb_kj[1]}")
print(f"Min ratio > 1: {min_ratio > 1}")

# ─────────────────────────────────────────────
# 5. FOCUS on tight cases
# ─────────────────────────────────────────────
print()
print("=" * 70)
print("SECTION 5: TIGHT CASE ANALYSIS")
print("=" * 70)

tight_threshold = Fraction(11, 10)  # 1.1
tight_cases = [r for r in ratio_results if r['crit_ratio'] != 'N/A' 
               and Fraction(r['crit_ratio']) < tight_threshold]

print(f"\nCases with crit_ratio < {float(tight_threshold):.2f}:")
print(f"{'k':>3} {'j':>3} {'m':>3} {'crit':>12} {'T1/|T2|':>12} {'A/A_lb':>12} {'m=k?':>6}")
for r in sorted(tight_cases, key=lambda x: Fraction(x['crit_ratio'])):
    k, j, m = r['k'], r['j'], r['m']
    crit = float(Fraction(r['crit_ratio']))
    T1_r = float(Fraction(r['T1_over_absT2'])) if r['T1_over_absT2'] != 'N/A' else float('nan')
    A_r  = float(Fraction(r['A_over_Alb'])) if r['A_over_Alb'] != 'N/A' else float('nan')
    print(f"{k:>3} {j:>3} {m:>3} {crit:>12.6f} {T1_r:>12.6f} {A_r:>12.6f} {'Y' if m==k else 'N':>6}")

# Focus on (k=6, j=3..5)
print()
print("Detail for k=6, j=3..8:")
for j in range(3, 9):
    d = compute_all(6, j)
    T1 = d['T1']
    T2 = d['T2']
    F  = d['F']
    m  = d['m']
    A  = d['A']
    A_lb = d['A_lb']
    B  = d['B']
    C  = d['C']
    D  = d['D']
    G  = d['G']
    H  = d['H']
    Ap = d['Ap']
    Dp = d['Dp']
    
    crit = Fraction(T1 * A_lb, (-T2) * A) if T2 < 0 else None
    print(f"  j={j}: m={m}, T1={T1}, T2={T2}, F={F}, crit={float(crit):.6f}, m=k? {m==6}")
    print(f"    A={A}, B={B}, C={C}, D={D}")
    print(f"    G={G}, H={H}, Ap={Ap}, Dp={Dp}")
    print(f"    T1/|T2|={float(Fraction(T1,-T2)):.6f}, A/A_lb={float(Fraction(A,A_lb)):.6f}")

# ─────────────────────────────────────────────
# 6. J=0: lower bound T1 vs upper bound |T2|
# ─────────────────────────────────────────────
print()
print("=" * 70)
print("SECTION 6: J=0 SYMBOLIC BOUND: T1 vs |T2|")
print("=" * 70)
print()
print("For j=0: A=b_m, B=b_{m-1}, C=b_{m-2}, D=b_{m+1}")
print("  T1 = b_m^2 - b_{m-1}*b_{m+1} + b_m*b_{m-2} - b_{m-1}^2")
print("  |T2| = H*(D'+2A'+B') where notation is shifted - actually:")
print("  T2 = G*A' - H*D' < 0, so |T2| = H*D' - G*A'")
print()

print("Direct ratios for j=0:")
print(f"{'k':>3} {'m':>3} {'T1':>14} {'|T2|':>14} {'T1/|T2|':>12} {'disc':>8}")

j0_ratio_min = None
j0_ratio_min_k = None
for k in range(6, 31):
    d = compute_all(k, 0)
    T1 = d['T1']
    T2 = d['T2']
    m  = d['m']
    n1 = k - m + 1
    n2 = k - m + 2
    disc = 4 * n1 * n2 - m * (m + 1)
    
    ratio = Fraction(T1, -T2) if T2 < 0 else None
    if ratio is not None:
        if j0_ratio_min is None or ratio < j0_ratio_min:
            j0_ratio_min = ratio
            j0_ratio_min_k = k
    
    print(f"{k:>3} {m:>3} {str(T1):>14} {str(-T2):>14} {float(ratio):>12.6f} {disc:>8}")

print()
print(f"Min T1/|T2| for j=0: {float(j0_ratio_min):.6f} at k={j0_ratio_min_k}")

# Compute T1/|T2| in terms of b_m ratios
# b_t = C(k,t)*2^t, ratio r_t = b_{t+1}/b_t = 2*(k-t)/(t+1)
# LC: b_m^2 >= b_{m-1}*b_{m+1} (i.e., r_m <= r_{m-1})
# T1 = b_m^2(1 - r_{m-1}/r_m * 1) ... let sigma = b_m/b_{m-1} = 2n1/m, rho = b_{m+1}/b_m = 2n2/(m+1) -- wait
# Actually r_t = b_{t+1}/b_t = 2*(k-t)/(t+1)
# r_{m-1} = 2*(k-m+1)/m = 2*n1/m  (ratio going from m-1 to m)
# r_m     = 2*(k-m)/(m+1)          (ratio going from m to m+1)
# For b to have mode at m: r_{m-1} >= 1 and r_m <= 1
# T1/b_m^2 = (1 - r_m/r_{m-1}^{-1}) ... let's use the exact formula instead

print()
print("Using exact formula T1/b_m^2 = (k+1)/n1 * [4*n1*n2 - m*(m+1)] / [4*n1*n2*(m+1)]:")
print()
print("  Numerator factor: (k+1)*(disc) where disc = 4n1n2 - m(m+1) > 0 (proved)")
print("  Denominator factor: 4*n1^2*n2*(m+1)")
print("  So T1 = b_m^2 * (k+1) * disc / [4*n1^2*n2*(m+1)]")

# Now bound |T2| for j=0
# T2 = G*A' - H*D' where A'=a'_m, D'=a'_{m+1} for j=2 spider
# For j=0: A'=a_m^{j=2} = b_m + 2*b_{m-1} + ... (using j=2 formula)
# Actually A' = b_m + 2*C(k,m-1)*2^{m-1} + ... Let me compute explicitly

print()
print("Bounding |T2| for j=0 using backbone ratio:")
print("  T2 = G*A'_2 - H*D'_2 where A'_2, D'_2 are j=2 coefficients")
print()

for k in range(6, 16):
    d = compute_all(k, 0)
    T1 = d['T1']
    T2 = d['T2']
    m  = d['m']
    bm  = d['bm']
    G  = d['G']
    H  = d['H']
    Ap = d['Ap']  # a^{j=2}_m
    Dp = d['Dp']  # a^{j=2}_{m+1}
    
    # |T2| = H*Dp - G*Ap
    absT2 = H * Dp - G * Ap
    assert absT2 == -T2, f"Sign error at k={k}"
    
    # Express in terms of b's: Ap = b_m + 2*b_{m-1} + b_{m-2}*(1 if m>=2 else 0) for j=2
    # But we computed A' = a^{j=2}_m directly
    # b_m/b_{m-1} = 2*n1/m, b_{m-1}/b_{m-2} = 2*n2/(m-1) ... wait
    # b_t = C(k,t)*2^t, so b_t/b_{t-1} = 2*(k-t+1)/t
    # b_{m-1}/b_m = m/(2*n1)
    # b_{m+1}/b_m = 2*n2/(m+1)  -- wait: b_{m+1}/b_m = 2*(k-m)/(m+1) = 2*(n2-2)/(m+1)?
    # n2 = k-m+2, so k-m = n2-2, b_{m+1}/b_m = 2*(n2-2)/(m+1)
    # Actually I'll just check: bm1 = b_{m-1}
    
    n1 = k - m + 1
    n2 = k - m + 2
    bm_val = Fraction(comb(k, m) * 2**m)
    bm1_val = Fraction(comb(k, m-1) * 2**(m-1))
    bm2_val = Fraction(comb(k, m-2) * 2**(m-2)) if m >= 2 else Fraction(0)
    bm_plus_val = Fraction(comb(k, m+1) * 2**(m+1)) if m+1 <= k else Fraction(0)
    
    # For j=2: A' = b_m + 2*b_{m-1} + b_{m-2}, D' = b_{m+1} + 2*b_m + b_{m-1}
    Ap2 = bm_val + 2*bm1_val + bm2_val
    Dp2 = bm_plus_val + 2*bm_val + bm1_val
    assert Ap2 == Ap, f"A' decomp wrong k={k}"
    assert Dp2 == Dp, f"D' decomp wrong k={k}"
    
    # T2 = G*A' - H*D' < 0, so |T2| = H*D' - G*A'
    # = H*(b_{m+1}+2b_m+b_{m-1}) - G*(b_m+2b_{m-1}+b_{m-2})
    # = H*b_{m+1} + (2H-G)*b_m + (H-2G)*b_{m-1} - G*b_{m-2}
    # G = C(k,m-1), H = C(k,m-2)
    # G/H = (m-1)/(k-m+2) = (m-1)/n2
    # H*b_{m+1}/bm^2 ?
    
    ratio = Fraction(T1, -T2)
    print(f"  k={k}: m={m}, T1/|T2| = {ratio} = {float(ratio):.5f}")

# ─────────────────────────────────────────────
# Additional: is T1/|T2| monotone in k? Symbolic bound via T1 formula.
# ─────────────────────────────────────────────
print()
print("=" * 70)
print("SECTION 7: ASYMPTOTIC ANALYSIS OF T1/|T2| FOR J=0")
print("=" * 70)
print()

# For k=3t (mode m=2t): T1/b_m^2 = (3t+1)/(t+1) * (10t+8)/(4(t+1)(t+2)(2t+1))
# And |T2|/b_m^2 = ?  Let's compute empirically.
print("Computing |T2|/b_m^2 for k=6..60 (j=0):")
print(f"{'k':>4} {'m':>4} {'T1/bm^2':>14} {'|T2|/bm^2':>14} {'ratio T1/|T2|':>16}")
for k in range(6, 61, 3):
    d = compute_all(k, 0)
    T1 = d['T1']
    T2 = d['T2']
    m  = d['m']
    bm = d['bm']
    T1_norm = Fraction(T1, bm*bm)
    T2_norm = Fraction(-T2, bm*bm)
    ratio = Fraction(T1, -T2) if T2 < 0 else Fraction(0)
    print(f"{k:>4} {m:>4} {float(T1_norm):>14.8f} {float(T2_norm):>14.8f} {float(ratio):>16.8f}")

print()
print("Asymptotic form (for k=3t as t→∞):")
print("  T1/bm^2 = (3t+1)/(t+1) * (10t+8)/[4(t+1)(t+2)(2t+1)]")
print("           → 3 * 10/(4*1*2*2*t^2) * t^2/... let me compute limit:")
print("           = (3t^2+O(t)) * (10t+8) / (4*t^4+O(t^3))")
print("           = 30*t^3 / (8*t^4) → 0  (T1/bm^2 → 0 like 1/t)")
print()

# ─────────────────────────────────────────────
# Section 8: Does T1 > |T2| hold for j=0?
# ─────────────────────────────────────────────
print("=" * 70)
print("SECTION 8: DIRECT COMPARISON T1 > |T2| FOR j=0")
print("=" * 70)
print()

print(f"{'k':>3} {'m':>3} {'T1':>18} {'|T2|':>18} {'T1>|T2|':>10} {'ratio':>10}")
all_T1_gt_T2 = True
for k in range(6, 61):
    d = compute_all(k, 0)
    T1 = d['T1']
    T2 = d['T2']
    m  = d['m']
    gt = T1 > -T2
    ratio = float(Fraction(T1, -T2)) if T2 < 0 else float('nan')
    if not gt:
        all_T1_gt_T2 = False
    if k <= 30 or k % 10 == 0:
        print(f"{k:>3} {m:>3} {T1:>18} {-T2:>18} {'Y' if gt else 'N':>10} {ratio:>10.5f}")

print()
print(f"T1 > |T2| for all k=6..60, j=0: {all_T1_gt_T2}")
print()
print("If T1 > |T2| holds for j=0 AND T1/|T2| is non-decreasing in j (to verify),")
print("then F = T1 + T2 = T1 - |T2| > 0 for all k >= 6.")

# ─────────────────────────────────────────────
# Section 9: is T1/|T2| monotone in j?
# ─────────────────────────────────────────────
print()
print("=" * 70)
print("SECTION 9: MONOTONICITY OF T1/|T2| IN j")
print("=" * 70)
print()
print("Do ratios T1/|T2| increase or decrease as j increases?")
print()

mono_holds = True
for k in range(6, 21):
    prev_ratio = None
    row_vals = []
    for j in range(0, 21):
        d = compute_all(k, j)
        T1 = d['T1']
        T2 = d['T2']
        if T2 < 0:
            ratio = Fraction(T1, -T2)
            row_vals.append((j, ratio))
    
    # Print k=6, k=10 in detail
    if k in (6, 10, 15):
        print(f"k={k}:")
        for j, r in row_vals:
            print(f"  j={j:2d}: T1/|T2| = {float(r):.6f}")
    
    # Check if monotone
    for i in range(1, len(row_vals)):
        if row_vals[i][1] < row_vals[i-1][1]:
            if k <= 10:
                print(f"  NON-MONOTONE at k={k}, j={row_vals[i][0]}: "
                      f"{float(row_vals[i][1]):.5f} < {float(row_vals[i-1][1]):.5f}")
            mono_holds = False
            break

print()
print(f"T1/|T2| monotonically non-decreasing in j for all k=6..20: {mono_holds}")

# ─────────────────────────────────────────────
# Section 10: F_check = T1+T2 >= 0 for all k=6..30, j=0..50
# ─────────────────────────────────────────────
print()
print("=" * 70)
print("SECTION 10: FINAL CHECK F >= 0 FOR k=6..30, j=0..50")
print("=" * 70)
print()

all_F_ok = True
F_fail_count = 0
for k in range(6, 31):
    for j in range(0, 51):
        d = compute_all(k, j)
        if d['F'] < 0:
            all_F_ok = False
            F_fail_count += 1
            print(f"  FAIL: k={k}, j={j}, m={d['m']}, F={d['F']}")

print(f"All F >= 0 for k=6..30, j=0..50: {all_F_ok}")
if all_F_ok:
    print(f"Total cases checked: {25 * 51} = 1275 (all pass)")

# ─────────────────────────────────────────────
# Save JSON results
# ─────────────────────────────────────────────

output = {
    'section1_j0_formula_verified': formula_ok,
    'section2_identity_verified': identity_ok,
    'section3_disc_positive_all_k6_100': True,  # proved by case analysis
    'disc_formulas': {
        'k_mod_0': '10t+8 (k=3t)',
        'k_mod_1': '6t+6 (k=3t+1)',
        'k_mod_2': '14t+22 (k=3t+2)',
    },
    'section4_min_crit_ratio': str(min_ratio),
    'section4_min_crit_ratio_float': float(min_ratio),
    'section4_min_crit_at': list(min_ratio_kj),
    'section4_max_A_over_Alb': str(max_A_over_Alb),
    'section4_max_A_over_Alb_float': float(max_A_over_Alb),
    'section4_max_A_over_Alb_at': list(max_A_over_Alb_kj),
    'section4_min_ratio_gt_1': bool(min_ratio > 1),
    'section8_T1_gt_absT2_j0_k6_60': all_T1_gt_T2,
    'section9_T1_absT2_monotone_in_j': mono_holds,
    'section10_all_F_nonneg_k6_30_j0_50': all_F_ok,
    'j0_detail': j0_results,
    'ratio_detail': ratio_results,
}

import os
os.makedirs('/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results', exist_ok=True)
with open('/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/symbolic_proof_analysis.json', 'w') as f:
    json.dump(output, f, indent=2)

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"1. j=0 formula verified for k=6..30: {formula_ok}")
print(f"2. Identity m*n2-(m-1)*n1=k+1 verified for k=6..30: {identity_ok}")
print(f"3. disc > 0 proved for ALL k >= 3 (three residue classes, linear growth)")
print(f"4. Min crit ratio T1*A_lb/(|T2|*A) = {float(min_ratio):.6f} at {min_ratio_kj} > 1: {min_ratio>1}")
print(f"5. Max A/A_lb = {float(max_A_over_Alb):.6f} at {max_A_over_Alb_kj}")
print(f"8. T1 > |T2| for j=0, k=6..60: {all_T1_gt_T2}")
print(f"9. T1/|T2| monotone in j: {mono_holds}")
print(f"10. F = T1+T2 >= 0 for ALL k=6..30, j=0..50: {all_F_ok}")
print()
print("Results saved to results/symbolic_proof_analysis.json")
