#!/usr/bin/env python3
"""Closed-form excess for the geometric LC distribution at the tie point.

Setup: p_{M-1} = p_M = p* (tied modes)
Left decay:  p_{M-1-j} = (s')^j · p*  for j = 1,2,...
Right decay: p_{M+j}   = r^j · p*     for j = 1,2,...

(This is the EXACT geometric case, the worst case for LC.)

Normalization (infinite tails for simplicity):
  2p* + p* s'/(1-s') + p* r/(1-r) = 1
  p* [2 + s'/(1-s') + r/(1-r)] = 1
  p* [(2(1-s')(1-r) + s'(1-r) + r(1-s')) / ((1-s')(1-r))] = 1
  numerator = 2(1-s')(1-r) + s'(1-r) + r(1-s')
            = 2 - 2s' - 2r + 2s'r + s' - s'r + r - rs'
            = 2 - s' - r  (the cross terms cancel!)
  
  So p* = (1-s')(1-r) / (2-s'-r)

Mean:
  mean = Σk k·p_k = (M-1)p* + M·p* + Σ_{j≥1} (M-1-j)(s')^j p* + Σ_{j≥1} (M+j) r^j p*
       = p* [(2M-1) + (M-1)s'/(1-s') - s'/(1-s')² + M r/(1-r) + r/(1-r)²]

Excess = M - mean:
  = M - p* [(2M-1) + (M-1)s'/(1-s') - s'/(1-s')² + M r/(1-r) + r/(1-r)²]
  = M(1 - p*[2 + s'/(1-s') + r/(1-r)]) + p*[1 - s'/(1-s') + s'/(1-s')² - r/(1-r)²]
  
  Since p*[2 + s'/(1-s') + r/(1-r)] = 1:
  first part = M(1-1) = 0

  Wait! That means M cancels out! Let me recompute.

  excess = M - p* [(2M-1) + (M-1)s'/(1-s') - s'/(1-s')² + M r/(1-r) + r/(1-r)²]
  
  Group M terms:
  = M [1 - p*(2 + s'/(1-s') + r/(1-r))]  +  p* [1 - s'/(1-s') + s'/(1-s')² + r/(1-r)²]
  
  Hmm wait, let me expand carefully:
  = M - p*(2M-1) - p*(M-1)*s'/(1-s') + p*s'/(1-s')² - p*M r/(1-r) - p*r/(1-r)²
  = M - 2Mp* + p* - p*(M-1)*s'/(1-s') + p*s'/(1-s')² - p*M r/(1-r) - p*r/(1-r)²
  = M[1 - 2p* - p*s'/(1-s') - p*r/(1-r)] + p*[1 + s'/(1-s') + s'/(1-s')² - r/(1-r)²]
  
  Wait, let me be more careful:
  = M - 2Mp* + p* - Mp*s'/(1-s') + p*s'/(1-s') + p*s'/(1-s')² - Mp*r/(1-r) - p*r/(1-r)²
  = M[1 - 2p* - p*s'/(1-s') - p*r/(1-r)] + p*[1 + s'/(1-s') + s'/(1-s')² - r/(1-r)²]
  
  But 2p* + p*s'/(1-s') + p*r/(1-r) = 1, so 1 - [2p* + p*s'/(1-s') + p*r/(1-r)] = 0.
  
  THE M TERMS CANCEL! 
  
  excess = p* [1 + s'/(1-s') + s'/(1-s')² - r/(1-r)²]
         = p* [1/(1-s') + s'/(1-s')² - r/(1-r)²]
         = p* [(1-s'+s')/(1-s')² - r/(1-r)²]
         = p* [1/(1-s')² - r/(1-r)²]

GORGEOUS. The excess is M-INDEPENDENT!

  excess = p* · [1/(1-s')² - r/(1-r)²]

  where p* = (1-s')(1-r)/(2-s'-r).

So: excess = (1-s')(1-r)/(2-s'-r) · [1/(1-s')² - r/(1-r)²]
           = (1-r)/(2-s'-r) · [1/(1-s') - r(1-s')/(1-r)²]

Hmm, let me simplify differently:
  excess = p* [1/(1-s')² - r/(1-r)²]
  
Let u = 1-s', v = 1-r. Then s' = 1-u, r = 1-v, u,v ∈ (0,1].
  
  p* = u·v/(u+v)  [since 2-s'-r = 2-(1-u)-(1-v) = u+v]
  
  excess = u·v/(u+v) · [1/u² - (1-v)/v²]
         = u·v/(u+v) · [v² - u²(1-v)] / (u²v²)
         = [v² - u²(1-v)] / [(u+v) u v]
         = [v² - u² + u²v] / [(u+v) u v]
         = [(v-u)(v+u) + u²v] / [(u+v) u v]
         = [(v-u) + u²v/(u+v)] / [u v / 1]    ... this is getting messy.

Let me just evaluate numerically.
"""

def compute_excess(s_prime, r):
    """Compute excess for geometric LC at tie point."""
    u = 1 - s_prime
    v = 1 - r
    
    if u <= 0 or v <= 0:
        return float('inf')
    
    pstar = u * v / (u + v)
    
    # excess = pstar * [1/u^2 - (1-v)/v^2]
    excess = pstar * (1/u**2 - (1-v)/v**2)
    # = pstar * (1/u^2 - r/v^2) where r = 1-v
    
    return excess

def main():
    print("=" * 70)
    print("  EXACT EXCESS FORMULA FOR GEOMETRIC LC AT TIE POINT")
    print("  excess = p* · [1/(1-s')² - r/(1-r)²]")
    print("  where p* = (1-s')(1-r)/(2-s'-r)")
    print()
    print("  KEY: THE EXCESS IS INDEPENDENT OF M!")
    print("=" * 70)
    print()
    
    # Table of excess for various (s', r) combinations
    print(f"{'s_prime':>7} {'r':>7} {'p*':>7} {'excess':>8} {'< 1?':>5}")
    print("-" * 40)
    
    for s_prime in [0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.49, 0.499, 0.4999]:
        for r in [0.0, 0.01, 0.1, 0.3, 0.5, 0.9]:
            u, v = 1 - s_prime, 1 - r
            pstar = u * v / (u + v)
            ex = compute_excess(s_prime, r)
            ok = "  ✓" if ex < 1 else " ✗"
            if s_prime in [0.0, 0.4, 0.49] and r in [0.0, 0.01, 0.5]:
                print(f"{s_prime:7.4f} {r:7.4f} {pstar:7.4f} {ex:8.5f} {ok}")
    
    print()
    
    # Find the WORST CASE (maximum excess over all s',r ∈ [0,1))
    print("  SEARCHING FOR MAXIMUM EXCESS:")
    max_ex = 0
    best_sr = (0, 0)
    
    for s100 in range(0, 100):
        s_prime = s100 / 100.0
        for r100 in range(0, 100):
            r = r100 / 100.0
            ex = compute_excess(s_prime, r)
            if ex > max_ex:
                max_ex = ex
                best_sr = (s_prime, r)
    
    print(f"  Max excess = {max_ex:.6f} at s'={best_sr[0]:.2f}, r={best_sr[1]:.2f}")
    
    # Refine
    s0, r0 = best_sr
    for ds in range(-10, 11):
        for dr in range(-10, 11):
            s_prime = s0 + ds * 0.001
            r = r0 + dr * 0.001
            if 0 <= s_prime < 1 and 0 <= r < 1:
                ex = compute_excess(s_prime, r)
                if ex > max_ex:
                    max_ex = ex
                    best_sr = (s_prime, r)
    
    print(f"  Refined: max excess = {max_ex:.8f} at s'={best_sr[0]:.4f}, r={best_sr[1]:.4f}")
    
    # Further refinement
    s0, r0 = best_sr
    for ds in range(-100, 101):
        for dr in range(-100, 101):
            s_prime = s0 + ds * 0.0001
            r = r0 + dr * 0.0001
            if 0 <= s_prime < 1 and 0 <= r < 1:
                ex = compute_excess(s_prime, r)
                if ex > max_ex:
                    max_ex = ex
                    best_sr = (s_prime, r)
    
    print(f"  Fine: max excess = {max_ex:.10f} at s'={best_sr[0]:.6f}, r={best_sr[1]:.6f}")
    
    # The limit as r → 0
    print()
    print("  LIMIT r → 0:")
    for s_prime in [0.3, 0.4, 0.45, 0.48, 0.49, 0.495, 0.499, 0.4999]:
        ex = compute_excess(s_prime, 0.0001)  # r ≈ 0
        print(f"    s'={s_prime:.4f}: excess = {ex:.8f}")
    
    # At r=0: excess = p* · 1/(1-s')² where p* = (1-s')·1/(2-s') = (1-s')/(2-s')
    # excess = (1-s')/(2-s') · 1/(1-s')² = 1/[(2-s')(1-s')]
    print()
    print("  At r=0 exactly: excess = 1 / [(2-s')(1-s')]")
    print("  This is maximized as s' → 1:")
    print("  lim_{s'→1} 1/[(2-s')(1-s')] = 1/[1·0] → ∞ ???")
    print()
    print("  WAIT — but s' < 1 is forced by the LC condition!")
    print("  Need to check: what constrains s' from approaching 1?")
    print()
    
    # Actually: for r=0, the normalization gives:
    # 2p* + p*s'/(1-s') = 1
    # p*(2 + s'/(1-s')) = 1
    # p*(2-2s'+s')/(1-s') = 1
    # p*(2-s')/(1-s') = 1
    # p* = (1-s')/(2-s')
    #
    # For this to be a valid distribution, we need p* > 0 and S_L < 1.
    # S_L = p*·s'/(1-s') = s'/(2-s')
    # 2p* + S_L = (2-2s'+s')/(2-s') = (2-s')/(2-s') = 1 ✓
    #
    # No constraint on s' other than s' ∈ [0, 1)!
    
    # So at r=0: excess = 1/[(2-s')(1-s')]
    # At s'=0: excess = 1/2
    # At s'=1/2: excess = 1/(3/2 · 1/2) = 1/(3/4) = 4/3 > 1 !!!
    
    print("  CRITICAL: at r=0, s'=0.5: excess = 1/(1.5·0.5) = 4/3 = 1.333...")
    print("  BUT this is the GEOMETRIC case, not actual LC sequences!")
    print()
    print("  The geometric sequence with s'=0.5 and r=0 is:")
    print("  ..., (1/2)^3, (1/2)^2, 1/2, 1, 1, 0, 0, ...")
    print("  = ..., 1/8, 1/4, 1/2, 1, 1, 0, 0, ...")
    print()
    
    # Check: is this sequence LC?
    p = [1/8, 1/4, 1/2, 1, 1, 0]
    print(f"  Sequence: {p}")
    for k in range(1, len(p)-1):
        lc = p[k]**2 >= p[k-1]*p[k+1]
        print(f"    k={k}: {p[k]}² = {p[k]**2:.4f} ≥ {p[k-1]*p[k+1]:.4f}? {lc}")
    
    # The last check: p[4]² = 1 ≥ p[3]·p[5] = 1·0 = 0 ✓
    # But p[3]² = 1 ≥ p[2]·p[4] = 0.5·1 = 0.5 ✓
    # So it IS LC!
    
    print()
    print("  THE GEOMETRIC SEQUENCE WITH s'=0.5, r=0 IS LC!")
    print("  But its excess = 4/3 > 1.")
    print()
    print("  HOWEVER: this is an INFINITE sequence (infinite left tail).")
    print("  For FINITE sequences supported on {0,...,d} with d finite,")
    print("  what happens?")
    print()
    
    # Finite version
    for M in [2, 5, 10, 20, 50, 100]:
        d = M + 1  # M is at positions M-1, M; right tail is just M
        # p = [s'^(M-1-k) for k=0..M-2, 1, 1, 0]
        # where s' = 0.5
        s_prime = 0.5
        poly = []
        for k in range(M-1):
            poly.append(s_prime ** (M-1-k))
        poly.append(1.0)  # position M-1
        poly.append(1.0)  # position M
        
        total = sum(poly)
        mean_val = sum(k * poly[k] for k in range(len(poly))) / total
        mode = max(range(len(poly)), key=lambda k: poly[k])
        excess = mode - mean_val
        
        print(f"    M={M}: support 0..{d}, mode={mode}, mean={mean_val:.4f}, excess={excess:.4f}")
    
    print()
    print("  So for FINITE geometric LC sequences, the excess CONVERGES to 4/3 > 1!")
    print("  This means the property mean ∈ [M-1,M] does NOT hold for all LC sequences!")
    print()
    print("  BUT WAIT — we TESTED 10K LC polys and found 0 violations!")
    print("  Those were products of (1+ti·x), which are a SUBSET of LC sequences.")
    print("  The geometric seq with s'=0.5, r=0 is NOT a product of linear factors!")
    print()
    
    # Verify: is [1/8, 1/4, 1/2, 1, 1] the IS poly of some tree?
    # IS polys have a_0 = 1 always (empty set).
    # This sequence has a_0 = 1/8, which is not 1. So NOT an IS poly as-is.
    # But we can normalize: multiply by 8 → [1, 2, 4, 8, 8]
    print("  Normalized: [1, 2, 4, 8, 8]")
    print("  Is this an IS poly of some graph?")
    print("  a_0 = 1 ✓ (empty set)")
    
    # Check: a_1 = 2 means 2 vertices. But a_4 = 8 means 8 IS of size 4.
    # This would require at most C(n,4) IS. With n ≥ 4 vertices.
    # But a_1 = 2 means n = 2!!! Can't have IS of size 4 with 2 vertices.
    # So this is NOT an IS poly.
    
    print("  a_1 = 2 → graph has 2 vertices. But a_4 = 8 → impossible!")
    print("  NOT an IS poly. Safe for our application.")
    print()
    
    print("  CONCLUSION:")
    print("  The mean ∈ [M-1,M] property at tie points does NOT hold for")
    print("  all LC sequences. It fails for geometric sequences with s' ≥ 1/2.")
    print("  BUT it holds for IS polynomials of all graphs tested (n ≤ 9).")
    print("  The IS poly structure provides additional constraints beyond LC.")
    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
