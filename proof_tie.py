#!/usr/bin/env python3
"""Exploit the tie condition at mode-jump points.

At λ₀ = i_{M-1}/i_M, the mode jumps from M-1 to M.
At this point: i_{M-1} λ₀^{M-1} = i_M λ₀^M (exact tie).

Let p_k = i_k λ₀^k / Z where Z = I(T; λ₀).

The tie says p_{M-1} = p_M, and both are maximum.

For excess < 1, we need:
  M - mean < 1
  where mean = Σ k p_k

Equivalently: mean > M - 1.

Now, since p_M = p_{M-1} (tied maximum):
  mean ≥ (M-1) p_{M-1} + M p_M = (M-1) p_M + M p_M = (2M-1) p_M

So: mean ≥ (2M-1) p_M

For mean > M-1, it suffices to show:
  (2M-1) p_M > M - 1
  p_M > (M-1)/(2M-1)

For M ≥ 2: (M-1)/(2M-1) < 1/2.

So: if p_M > 1/2, we're done. But p_M < 1/2 for large trees (many coefficients).

Hmm, this bound is too weak for individual p_M. Let me be smarter.

APPROACH: Use the tie + unimodality of the weighted coefficients.

At the jump point, the distribution {p_k} has two adjacent maxima at M-1 and M.
If the distribution is unimodal (or bimodal with adjacent modes), then:
  p_{M-2} ≤ p_{M-1} = p_M ≥ p_{M+1}

This means the distribution is peaked around M-1 and M.

The mean is pulled by the tails away from M.
The left tail pulls mean down, the right tail pulls mean up.

mean - (M - 1/2) = Σ_k (k - M + 1/2) p_k
  = Σ_{k<M-1} (k - M + 1/2) p_k + (-1/2) p_{M-1} + (1/2) p_M 
    + Σ_{k>M} (k - M + 1/2) p_k
  = Σ_{k<M-1} (k - M + 1/2) p_k + Σ_{k>M} (k - M + 1/2) p_k

Since p_{M-1} = p_M, the two middle terms cancel to 0.

mean - (M - 1/2) = Σ_{k<M-1} (k - M + 1/2) p_k + Σ_{k>M} (k - M + 1/2) p_k

For k < M-1: k - M + 1/2 ≤ -3/2 (since k ≤ M-2)
For k > M: k - M + 1/2 ≥ 3/2 (since k ≥ M+1)

Hmm, this still doesn't simplify. But the POINT of the tie is that the 
distribution is symmetric-ish around M-1/2. If the right tail exceeds the 
left tail (weighted), mean > M-1/2 > M-1, and we're done.

For trees, IS polys tend to be right-skewed (heavy right tail). 
Let me check: at the mode-jump point, is the right tail heavier than the left?

ALTERNATIVELY: Use the tree recurrence I(T;λ) = I(T-v;λ) + λ·I(T-N[v];λ).
At the jump point, this gives a relationship between the two sub-trees
and the tie condition.

Let me first check the skewness at jump points.
"""

import subprocess
from indpoly import independence_poly
from graph6 import parse_graph6

def main():
    print("=" * 80)
    print("  TIE-CONDITION ANALYSIS AT MODE JUMPS")
    print("=" * 80)
    print()
    
    # For each tree, find the TIGHTEST mode jump (largest excess at jump)
    # and analyze the distribution at that point
    
    print(f"{'n':>3} {'M':>3} {'excess':>8} {'p_M':>8} {'S_L':>8} {'S_R':>8} "
          f"{'wt_L':>8} {'wt_R':>8} {'asym':>8}")
    print("-" * 78)
    
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
                    total_w = sum(weights)
                    if total_w > 0:
                        p = [w/total_w for w in weights]
                        mean_val = sum(k * p[k] for k in range(len(p)))
                        excess = M - mean_val
                        
                        if excess > best_excess:
                            best_excess = excess
                            # Left/right tail analysis
                            S_L = sum(p[k] for k in range(M-1))  # mass below M-1
                            S_R = sum(p[k] for k in range(M+1, len(p)))  # mass above M
                            # Weighted distance from center (M-1/2)
                            wt_L = sum((M - 0.5 - k) * p[k] for k in range(M-1))
                            wt_R = sum((k - M + 0.5) * p[k] for k in range(M+1, len(p)))
                            asym = wt_R - wt_L  # positive means right-skewed
                            
                            best_info = (M, p[M], S_L, S_R, wt_L, wt_R, asym, line)
        
        M, pM, SL, SR, wt_L, wt_R, asym, _ = best_info
        print(f"{n:3d} {M:3d} {best_excess:8.4f} {pM:8.4f} {SL:8.4f} {SR:8.4f} "
              f"{wt_L:8.4f} {wt_R:8.4f} {asym:+8.4f}")
    
    print()
    print("  asym = wt_R - wt_L: positive means right tail heavier")
    print("  For excess < 1: need mean > M-1, i.e., excess < 1")
    print("  mean = (M-1/2) + (wt_R - wt_L) = (M-1/2) + asym")
    print("  So excess = M - mean = M - (M-1/2) - asym = 1/2 - asym")
    print()
    
    # WAIT: excess = 1/2 - asym? Let me verify.
    # mean = Σk p_k = (M-1)p_{M-1} + M p_M + Σ_{k<M-1} k p_k + Σ_{k>M} k p_k
    # Since p_{M-1} = p_M:
    # mean = (2M-1)p_M + Σ_{k<M-1} k p_k + Σ_{k>M} k p_k
    # 
    # Actually, let me compute it more carefully:
    # mean = Σ k p_k
    # M - mean = Σ (M-k) p_k = Σ_{k<M} (M-k) p_k - Σ_{k>M} (k-M) p_k
    #          = p_{M-1} + Σ_{k<M-1} (M-k) p_k - Σ_{k>M} (k-M) p_k
    #          = p_M + Σ_{k<M-1} (M-k) p_k - Σ_{k>M} (k-M) p_k  (since p_{M-1}=p_M)
    
    # So: excess = p_M + (left-tail weighted sum) - (right-tail weighted sum)
    
    # For excess < 1:
    # p_M + left_sum - right_sum < 1
    # Since p_M + p_{M-1} + S_L + S_R = 1 and p_{M-1}=p_M:
    # 2p_M + S_L + S_R = 1
    # p_M = (1 - S_L - S_R) / 2
    
    # excess = (1-S_L-S_R)/2 + left_sum - right_sum
    # left_sum = Σ_{k<M-1} (M-k) p_k ≤ M * S_L
    # right_sum = Σ_{k>M} (k-M) p_k ≥ S_R
    
    # excess ≤ (1-S_L-S_R)/2 + M*S_L - S_R
    #        = 1/2 - S_L/2 - S_R/2 + M*S_L - S_R
    #        = 1/2 + (M-1/2)*S_L - 3S_R/2
    
    # For this to be < 1: (M-1/2)*S_L < 1/2 + 3S_R/2
    
    # Hmm, still depends on M*S_L, which can be large.
    
    # But the MONOTONICITY constraint from the tie is crucial:
    # At the tie point, p_k is maximized at M-1 and M.
    # If IS polys have the property that i_k λ^k is unimodal (for trees),
    # then p_{M-2} ≤ p_{M-1} = p_M ≥ p_{M+1}.
    # This constrains the tail masses!
    
    # In fact, for p_k to be unimodal with peak at {M-1, M}:
    # p_0 ≤ p_1 ≤ ... ≤ p_{M-1} = p_M ≥ p_{M+1} ≥ ... ≥ p_α
    
    # Under this unimodality:
    # S_L = Σ_{k<M-1} p_k ≤ (M-1) * p_{M-2} ≤ (M-1) * p_M
    # S_R = Σ_{k>M} p_k ≤ (α-M) * p_{M+1} ≤ (α-M) * p_M
    
    print("  If p_k is unimodal (which IS polys of trees should be for all λ):")
    print("  then the tail constraints give a tighter analysis.")
    print()
    
    # Let me check: is i_k λ^k unimodal for all λ for tree IS polys?
    # This is equivalent to: for all λ, the IS poly evaluated at λ has
    # the same unimodality as the original.
    # 
    # If the IS poly is log-concave, then i_k λ^k is also log-concave
    # (since multiplying by λ^k preserves LC). So yes, for LC polys,
    # unimodality is preserved for all λ.
    #
    # But we haven't proved LC! However, we can check it numerically.
    
    print("  CHECK: Is i_k λ^k unimodal at the tightest jump points?")
    
    for n in [8, 12, 16]:
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            
            for M in range(1, len(poly)):
                if poly[M] > 0 and poly[M-1] > 0:
                    lam0 = poly[M-1] / poly[M]
                    weights = [poly[k] * (lam0 ** k) for k in range(len(poly))]
                    total_w = sum(weights)
                    p = [w/total_w for w in weights]
                    
                    mean_val = sum(k * p[k] for k in range(len(p)))
                    excess = M - mean_val
                    
                    if excess > 0.9:
                        # Check unimodality
                        is_unimodal = True
                        found_peak = False
                        for j in range(1, len(p)):
                            if not found_peak:
                                if p[j] < p[j-1]:
                                    found_peak = True
                            else:
                                if p[j] > p[j-1]:
                                    is_unimodal = False
                                    break
                        
                        print(f"    n={n}, M={M}, excess={excess:.4f}: unimodal={is_unimodal}")
    
    print()
    print("=" * 80)
    
    # KEY INSIGHT from the tie condition:
    # At the tie, p_{M-1} = p_M. These are the two largest values.
    # excess = p_M + Σ_{k<M-1}(M-k)p_k - Σ_{k>M}(k-M)p_k
    # 
    # A LOWER BOUND on mean comes from the tie constraint:
    # mean ≥ (M-1)·p_{M-1} + M·p_M = (2M-1)·p_M
    # 
    # So excess ≤ M - (2M-1)p_M = M(1-2p_M) + p_M = M - 2Mp_M + p_M
    # 
    # For excess < 1: M(1-2p_M) + p_M < 1
    # i.e., p_M > (M-1)/(2M-1) = 1/2 - 1/(2(2M-1))
    # 
    # So we need p_M > 1/2 - ε. But for large trees, p_M << 1/2.
    # This bound isn't tight enough on its own.
    
    # BETTER: Use ALL coefficients, not just the tied pair.
    # The key structural constraint from trees is the DP recurrence.

if __name__ == "__main__":
    main()
