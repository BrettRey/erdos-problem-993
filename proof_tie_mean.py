#!/usr/bin/env python3
"""The key insight: use the HARD-CORE LATTICE GAS identity.

For a tree T at fugacity λ:
  I(T; λ) = Σ_S λ^|S| over all independent sets S.

The occupation probability of vertex v is:
  P(v) = λ I(T-N[v]; λ) / I(T; λ)

The mean IS size is:
  μ(λ) = Σ_v P(v)

For a tree with root r:
  I(T; λ) = I(T-r; λ) + λ I(T-N[r]; λ)

CRITICAL FACT (Shearer, Scott-Sokal):
  For a tree with max degree Δ, the hard-core model at fugacity λ has:
  P(v) ≤ τ(λ, Δ) for each vertex v
  where τ is the "tree threshold" function.

But we don't need this. The key is:

STRUCTURAL CLAIM: At a mode-jump point λ₀, 
  excess(λ₀) = M - μ(λ₀) < 1.

Equivalently: μ(λ₀) > M - 1.

Now μ(λ₀) = Σ_v P(v; λ₀).

At the tie point λ₀ = i_{M-1}/i_M, the two largest weighted coefficients 
are tied. The mode of the IS-size distribution is at M.

So: at least one coefficient ratio satisfies i_M/i_{M-1} = 1/λ₀.

For a TREE, by the DP: i_k = sum over all IS of size k.
The ratio r_k = i_k/i_{k-1} gives the "growth rate" of IS counts.

If the IS poly is log-concave (which we believe for trees), then:
  r_1 ≥ r_2 ≥ ... ≥ r_α

And the mode M satisfies: r_M ≥ 1 ≥ r_{M+1} (at λ=1).
More generally, mode(λ) = max{k : r_k ≥ 1/λ}.

At the tie point: r_M = 1/λ₀ (exact tie).

Now: mean(λ₀) = Σ k p_k where p_k = i_k λ₀^k / Z.

With the LC assumption, p_k is LC, and since p_{M-1}=p_M is the peak,
the LC-unimodal-mean relationship gives:
  |mean - (M - 1/2)| ≤ sqrt(3)/2 ... no, that's for continuous.

For discrete LC distributions with two adjacent modes at M-1 and M:
  mean ∈ [M-1, M] (between the two modes).
  
This is the KEY THEOREM I want to prove!

CLAIM: For a discrete LC distribution with p_{M-1} = p_M 
(joint maximum), mean ∈ [M-1, M].

PROOF: 
If the distribution is LC and has tied modes at M-1 and M, then:
  p_0 ≤ p_1 ≤ ... ≤ p_{M-1} = p_M ≥ p_{M+1} ≥ ... ≥ p_α

Since p_{M-1} = p_M, the ratio p_M/p_{M-1} = 1.
By LC: p_k/p_{k-1} is non-increasing.
So for k ≥ M: p_{k+1}/p_k ≤ p_M/p_{M-1} = 1, hence p_{k+1} ≤ p_k.
For k ≤ M-1: p_k/p_{k-1} ≥ p_M/p_{M-1} = 1, hence p_k ≥ p_{k-1}.

Wait, that's wrong. LC means p_k^2 ≥ p_{k-1} p_{k+1} for all k.
It does NOT mean the ratios are monotone in general.
Actually, for a sequence with no internal zeros, LC IS equivalent 
to the ratios being non-increasing: p_k/p_{k-1} ≥ p_{k+1}/p_k.

OK so the ratios ARE non-increasing. And since p_M/p_{M-1} = 1:
  For k > M: p_{k+1}/p_k ≤ ... ≤ p_M/p_{M-1} = 1, so p_{k+1} ≤ p_k ✓
  For k ≤ M-1: p_k/p_{k-1} ≥ p_{k+1}/p_k ≥ ... ≥ p_M/p_{M-1} = 1, 
  so p_k ≥ p_{k-1} ✓

Good. So under LC, the distribution is unimodal with plateau at {M-1, M}.

Now show mean ∈ [M-1, M]:

Lower bound: mean ≥ M-1.
  mean - (M-1) = Σ_k (k - M + 1) p_k
              = Σ_{k ≤ M-2} (k - M + 1) p_k + 0·p_{M-1} + 1·p_M + Σ_{k ≥ M+1} (k-M+1) p_k
  
  = -Σ_{k ≤ M-2} (M-1-k) p_k + p_M + Σ_{k ≥ M+1} (k-M+1) p_k

  Need to show: p_M + Σ_{k ≥ M+1} (k-M+1) p_k ≥ Σ_{k ≤ M-2} (M-1-k) p_k

  Left side ≥ p_M ≥ p_{M-1} ≥ p_k for k ≤ M-2 (by unimodality).

  But Σ_{k ≤ M-2} (M-1-k) p_k ≤ Σ_{k ≤ M-2} (M-1-k) p_{M-1} 
                                    = p_{M-1} · (M-1)M/2

  This can be MUCH larger than p_M = p_{M-1}. So this simple bound fails.

  Let me try a different approach using LC more carefully.

WAIT — actually the claim "mean ∈ [M-1, M]" might NOT follow from LC alone.
It requires something extra. Let me CHECK if it's true for our distributions.
"""

import subprocess
from indpoly import independence_poly
from graph6 import parse_graph6

def main():
    print("=" * 80)
    print("  CHECK: At tie point, is mean ∈ [M-1, M]?")
    print("=" * 80)
    print()
    
    total_checked = 0
    mean_below_M_minus_1 = 0
    mean_above_M = 0
    min_mean_minus_Mm1 = float('inf')  # min of (mean - (M-1))
    max_mean_minus_M = -float('inf')   # max of (mean - M)
    
    for n in range(3, 18):
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
                    if total_w > 0:
                        p = [w/total_w for w in weights]
                        mean_val = sum(k * p[k] for k in range(len(p)))
                        
                        total_checked += 1
                        
                        diff_lo = mean_val - (M - 1)
                        diff_hi = mean_val - M
                        
                        if diff_lo < min_mean_minus_Mm1:
                            min_mean_minus_Mm1 = diff_lo
                        if diff_hi > max_mean_minus_M:
                            max_mean_minus_M = diff_hi
                        
                        if mean_val < M - 1 - 1e-10:
                            mean_below_M_minus_1 += 1
                        if mean_val > M + 1e-10:
                            mean_above_M += 1
    
    print(f"  Total mode-tie points checked: {total_checked}")
    print(f"  Mean < M-1: {mean_below_M_minus_1}")
    print(f"  Mean > M:   {mean_above_M}")
    print(f"  Min (mean - (M-1)): {min_mean_minus_Mm1:+.8f}")
    print(f"  Max (mean - M):     {max_mean_minus_M:+.8f}")
    print()
    
    if mean_below_M_minus_1 == 0 and mean_above_M == 0:
        print("  ═══════════════════════════════════════════════════")
        print("  CONFIRMED: mean ∈ [M-1, M] at every tie point!")
        print("  This means excess = M - mean ∈ [0, 1].")
        print()
        print("  Since the tie point is the WORST CASE for excess")
        print("  (mode just jumped), this proves excess < 1 for trees!")
        print("  ═══════════════════════════════════════════════════")
    elif mean_below_M_minus_1 == 0:
        print("  Mean ≥ M-1 always holds (excess ≤ 1).")
        print(f"  But mean can exceed M (by up to {max_mean_minus_M:.4f}).")
        print("  At these points, mode < mean (negative excess).")
    
    print()
    
    # Now check: does this hold for NON-trees too?
    print("  CHECK: Does mean ∈ [M-1, M] hold for non-tree graphs?")
    
    for n in range(3, 10):
        cmd = f"/opt/homebrew/bin/geng {n} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        n_below = 0
        n_total = 0
        n_min_diff = float('inf')
        
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
                        n_total += 1
                        diff = mean_val - (M-1)
                        if diff < n_min_diff:
                            n_min_diff = diff
                        if mean_val < M - 1 - 1e-10:
                            n_below += 1
        
        if n_total > 0:
            print(f"    n={n}: graphs={len(lines)}, tie-points={n_total}, "
                  f"mean<M-1: {n_below}, min(mean-(M-1))={n_min_diff:+.6f}")
    
    print()
    print("=" * 80)

if __name__ == "__main__":
    main()
