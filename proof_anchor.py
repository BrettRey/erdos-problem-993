#!/usr/bin/env python3
"""THE FINAL PUSH: derive the excess bound from the DP structure.

SETUP:
  I(T; x) = A(x) + x·B(x) where A = I(T-v), B = I(T-N[v]).
  Both A and B have a₀ = 1 (empty set is always an IS).

At a mode-tie point λ₀ for I(T): mode jumps from M-1 to M.
  i_M λ₀^M = i_{M-1} λ₀^{M-1}

Now I(T;λ₀) = A(λ₀) + λ₀·B(λ₀).

The mean of I(T) at λ₀:
  mean_T = λ₀ I'(T;λ₀)/I(T;λ₀)
         = λ₀ [A'(λ₀) + B(λ₀) + λ₀ B'(λ₀)] / [A(λ₀) + λ₀ B(λ₀)]
         = [λ₀ A'(λ₀) + λ₀ B(λ₀) + λ₀² B'(λ₀)] / I(T;λ₀)

  Define w = λ₀ B(λ₀) / I(T;λ₀) = probability vertex v is in the IS.

  mean_T = λ₀ A'(λ₀)/I(T;λ₀) + w + λ₀² B'(λ₀) / I(T;λ₀)

  Now λ₀ A'(λ₀)/A(λ₀) = mean_A (mean of A at λ₀).
  And λ₀ B'(λ₀)/B(λ₀) = mean_B (mean of B at λ₀).
  And A(λ₀)/I(T;λ₀) = 1-w, λ₀ B(λ₀)/I(T;λ₀) = w.

  So: mean_T = (1-w)·mean_A + w·(1 + mean_B)
             = mean_A + w·(1 + mean_B - mean_A)

  excess_T = M - mean_T = M - mean_A - w(1 + mean_B - mean_A)
           = (M - mean_A) - w(1 + mean_B - mean_A)

  Note: M is the mode of I(T) at λ₀. The mode of A at λ₀ is some m_A.

  If m_A = M (no mode shift): excess_T = excess_A - w(1+mean_B-mean_A) < excess_A.
  If m_A = M-1 (mode shifts up by 1):
    excess_T = (M - mean_A) - w(1+mean_B-mean_A)
             = (1 + m_A - mean_A) - w(1+mean_B-mean_A)
             = 1 + excess_A - w(1+mean_B-mean_A)

  For excess_T < 1: w(1+mean_B-mean_A) > excess_A.

  Now w = P(v in IS) = occupation probability of v at fugacity λ₀.

  For a tree rooted at v:
  P(v) = λ₀ B(λ₀) / I(T;λ₀) = λ₀ / (A(λ₀)/B(λ₀) + λ₀)

  Key constraint: A(λ₀)/B(λ₀) = I(T-v;λ₀)/I(T-N[v];λ₀).

  For large λ₀: P(v) → 1/(1 + A(λ₀)/B(λ₀)/λ₀) → depends on growth.

OBSERVATION: w and (1+mean_B-mean_A) are NOT independent.
When mean_B is close to mean_A, the correction is small;
but w = P(v) is then typically large (since B dominates A).

When mean_B ≫ mean_A, the correction per unit w is large;
but w = P(v) might be small.

The PRODUCT w·(1+mean_B-mean_A) is what matters.

This product = P(v) · [1 + mean_B(λ₀) - mean_A(λ₀)].

Now: mean_T = mean_A + P(v)·(1+mean_B-mean_A).
So P(v)·(1+mean_B-mean_A) = mean_T - mean_A.

We need: mean_T - mean_A > excess_A = m_A - mean_A.
i.e., mean_T > m_A = M - 1.

Wait — THIS IS JUST THE ORIGINAL CLAIM! mean_T > M - 1 at the tie.

Hmm, let me think more carefully...

We're trying to prove: at the tie point where mode=M,
  mean_T > M - 1.

Express this using the DP:
  mean_T = mean_A + P(v)*(1 + mean_B - mean_A) > M - 1

If we know (by induction) that mode_A ≤ ceil(mean_A), and mode_A is at 
most M-1 (from the gap), then:
  M-1 ≤ ceil(mean_A) ≤ mean_A + 1
  so mean_A ≥ M - 2.

  Then: mean_T = mean_A + P(v)*(1+mean_B-mean_A) ≥ (M-2) + P(v)*(1+mean_B-(M-2))

  For mean_T > M-1: P(v)*(1+mean_B-mean_A) > 1 - (mean_A - (M-2))
  which is P(v)*(1+mean_B-mean_A) > 1 - excess_A_from_mode.

Hmm, can I bound P(v) from below?

For a tree, the occupation probability of a LEAF v at fugacity λ:
  I(T) = I(T-v) + λ·I(T-v)  (since N[v] = {v, parent}, T-N[v] = T-v minus parent)
  Wait, that's not right for arbitrary trees. Let me be precise.

For a leaf v with parent u:
  T-v has n-1 vertices.
  T-N[v] = T - {v,u} has n-2 vertices.

  I(T) = I(T-v) + λ·I(T-{v,u})

  P(v) = λ·I(T-{v,u}) / I(T) = λ·I(T-{v,u}) / [I(T-v) + λ·I(T-{v,u})]
       = 1 / (I(T-v)/(λ·I(T-{v,u})) + 1)

  For large λ: P(v) → 1 / (1 + 0) = 1 (leaf always occupied).
  For small λ: P(v) → 0.

  The KEY: in the hard-core model on a tree, there's a SELF-CONSISTENCY condition.
  If parent u has high P(u), vertex v is blocked and P(v) is low.
  If P(u) is low, P(v) is higher.

  P(v) = λ / (1 + λ · ∏_{w∈N(v)} (1-P(w))^{-1} ... )

Actually for trees the EXACT recursion is:
  P(v) = λ / (1 + λ + λ · Σ_{children} ...)

Let me just focus on what we CAN prove directly.

THE M-INDEPENDENT FORMULA gives us: for geometric LC at tie,
excess = 1/[(2-s')(1-s')] (at r=0). This is < 1 iff s' < (3-√5)/2 ≈ 0.382.

For IS polys, the left decay rate s' at tie points is computationally 
bounded by ~0.475. So we need a TIGHTER analysis.

BUT the geometric LC isn't the actual IS poly distribution — IS polys 
have sharp left truncation at k=0 (a₀=1), which prevents the left tail 
from being too heavy. This is the key constraint.

Let me verify: at the tightest tie points for IS polys, what is the 
ratio p₀/p* (the "anchor strength")?
"""

import subprocess
from indpoly import independence_poly
from graph6 import parse_graph6

def main():
    print("=" * 70)
    print("  ANCHOR ANALYSIS: How does p₀ = 1/Z constrain the excess?")
    print("=" * 70)
    print()
    
    for n in range(5, 18):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        best_excess = -100
        best_info = None
        
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            
            for M in range(1, len(poly)):
                if poly[M] > 0 and poly[M-1] > 0:
                    lam0 = poly[M-1] / poly[M]
                    weights = [poly[k] * (lam0 ** k) for k in range(len(poly))]
                    Z = sum(weights)
                    p = [w/Z for w in weights]
                    mean_val = sum(k * p[k] for k in range(len(p)))
                    excess = M - mean_val
                    
                    if excess > best_excess:
                        best_excess = excess
                        pstar = p[M]
                        p0 = p[0]
                        left_mass = sum(p[k] for k in range(M-1))
                        right_mass = sum(p[k] for k in range(M+1, len(p)))
                        
                        # Actual s' and r
                        s_prime = p[M-2]/pstar if M >= 2 and pstar > 0 else 0
                        r = p[M+1]/pstar if M+1 < len(p) and pstar > 0 else 0
                        
                        # Left tail composition
                        wt_left = sum((M-k)*p[k] for k in range(M-1))
                        wt_right = sum((k-M)*p[k] for k in range(M+1, len(p)))
                        
                        best_info = (M, pstar, p0, left_mass, right_mass, 
                                    s_prime, r, wt_left, wt_right, mean_val)
        
        M, pstar, p0, SL, SR, sp, r, wL, wR, mu = best_info
        print(f"  n={n:2d}: M={M:2d} ex={best_excess:.4f} p*={pstar:.4f} p₀={p0:.6f} "
              f"S_L={SL:.4f} S_R={SR:.4f} s'={sp:.4f} r={r:.4f} wL={wL:.4f} wR={wR:.4f}")
    
    print()
    
    # Key question: in the limit, p₀ → 0 and cannot anchor the mean.
    # But the LEFT TAIL is truncated at k=0, so:
    # Even though individual p_k decrease geometrically going left,
    # there are only M-1 terms (not infinite).
    #
    # LEFT TAIL SUM ≤ p* · s'/(1-s') but ALSO ≤ p* · (M-1) (at most M-1 terms)
    # For small M: the truncation matters a lot.
    # For large M and s' close to 1: truncation barely matters.
    #
    # BUT for IS polys, s' can't be too close to 1 when M is large!
    # Because the IS poly has a₀ = 1, and the geometric decay s'^{M-1} would
    # make p₀ = s'^{M-1}·p*, but p₀ = (1·λ₀⁰)/Z = 1/Z = p₀.
    # So: s'^{M-1}·p* ≈ p₀ = 1/Z.
    #
    # If s' is close to 1 and M is large, then s'^{M-1} is moderate,
    # but for the IS poly, p₀ = a₀·1/Z = 1/Z is small (Z is large).
    #
    # Actually: p₀ = 1/Z and p* · s'^{M-1} should approximate p₀ if
    # the decay is truly geometric. Let me check.
    
    print("  CHECK: p₀ vs s'^{M-1} · p*:")
    for n in range(5, 18):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            
            for M in range(2, len(poly)):
                if poly[M] > 0 and poly[M-1] > 0:
                    lam0 = poly[M-1] / poly[M]
                    weights = [poly[k] * (lam0 ** k) for k in range(len(poly))]
                    Z = sum(weights)
                    p = [w/Z for w in weights]
                    excess = M - sum(k*p[k] for k in range(len(p)))
                    
                    if excess > 0.95 and M >= 3:
                        pstar = p[M]
                        sp = p[M-2]/pstar if pstar > 0 else 0
                        p0_actual = p[0]
                        p0_predicted = pstar * (sp ** (M-1)) if sp > 0 else 0
                        
                        print(f"    n={n:2d} M={M:2d} ex={excess:.4f}: "
                              f"p₀={p0_actual:.2e} vs s'^(M-1)·p*={p0_predicted:.2e} "
                              f"ratio={p0_actual/p0_predicted:.4f}" if p0_predicted > 0 else
                              f"    n={n:2d} M={M:2d}: p0_predicted=0")
                        break
            else:
                continue
            break
    
    print()
    
    # CRITICAL OBSERVATION: 
    # For IS poly at tie point, the LEFT TRUNCATION means:
    # p₀ = 1/Z = concrete value.
    # If the geometric decay predicts s'^{M-1}·p* ≈ p₀,
    # then s' ≈ (p₀/p*)^{1/(M-1)}.
    # For large M: s' → 1 (since p₀/p* → 0, but slowly).
    # 
    # HOWEVER: the excess formula excess = 1/[(2-s')(1-s')] uses 
    # INFINITE tails. For finite left tail of length M-1:
    #
    # left_sum = p* · (s' + s'² + ... + s'^{M-1}) = p* · s'(1-s'^{M-1})/(1-s')
    # left_moment = p* · Σ_{j=1}^{M-1} j·s'^j
    #
    # The FINITE-TAIL excess is smaller than the infinite-tail excess.
    # How much smaller?
    
    print("  FINITE vs INFINITE excess comparison:")
    import math
    
    for s_prime in [0.3, 0.4, 0.45, 0.47, 0.499]:
        for M in [5, 10, 20, 50, 100]:
            r = 0.001  # nearly zero right tail
            
            pstar_inf = (1-s_prime)*(1-r) / (2-s_prime-r)
            excess_inf = pstar_inf * (1/(1-s_prime)**2 - r/(1-r)**2)
            
            # Finite version
            left_sum = sum(s_prime**j for j in range(1, M))
            right_sum = sum(r**j for j in range(1, M))  # truncate symmetrically
            total = 2 + left_sum + right_sum
            pstar_fin = 1.0  # unnormalized
            p_fin = [pstar_fin * s_prime**(M-1-k) for k in range(M-1)]
            p_fin.append(pstar_fin)  # M-1
            p_fin.append(pstar_fin)  # M
            for j in range(1, M):
                p_fin.append(pstar_fin * r**j)
            
            Z_fin = sum(p_fin)
            p_fin = [x/Z_fin for x in p_fin]
            mean_fin = sum(k*p_fin[k] for k in range(len(p_fin)))
            mode_fin = M  # by construction
            excess_fin = mode_fin - mean_fin
            
            if s_prime in [0.47, 0.499] or M in [10, 100]:
                print(f"    s'={s_prime:.3f} M={M:3d}: excess_finite={excess_fin:.4f} excess_infinite={excess_inf:.4f}")
    
    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
