#!/usr/bin/env python3
"""DIRECT PROOF via the tree DP and the left-tail truncation.

KEY INSIGHT from the geometric analysis:
- For INFINITE geometric LC, excess can exceed 1.
- For FINITE sequences starting at a₀ = 1, the LEFT TRUNCATION at k=0
  creates a massive anchor that keeps the mean high enough.

For IS polynomials: a₀ = 1 (empty set), and a₀ is LARGE relative to 
sequences with heavy left tails.

Let me formalize: 

APPROACH: "a₀ = 1 anchor" argument.

For P(x) = Σ_{k=0}^{α} aₖ xᵏ with a₀ = 1:
At tie point λ₀: p₀ = a₀/Z = 1/Z.
Mean = Σ k pₖ.
Mode = M.
Excess = M - mean.

Now: Σ pₖ = 1 means Z = Σ aₖ λ₀ᵏ.
p₀ = 1/Z.
Mean ≥ 0 (trivially, since all kpₖ ≥ 0).

But we need mean ≥ M - 1.

The anchor p₀ = 1/Z contributes 0 to the mean but contributes to normalization.
Higher k terms contribute k·pₖ to the mean.

Hmm, p₀ actually HURTS the mean (adds mass at 0, reducing the average).

Actually, let me think about this differently. The reason why IS polys 
behave well is not the anchor alone — it's the GROWTH STRUCTURE.

For an IS poly of a tree on n vertices:
- a₀ = 1
- a₁ = n  (each vertex is an IS of size 1)
- a₂ = number of non-adjacent pairs
- ...
- aₖ = number of IS of size k

These counts grow rapidly from a₀ to a₁ to ... to mode M, then decrease.
The GROWTH RATE a_{k+1}/a_k is large for small k and decreases.

At the tie point λ₀ = a_{M-1}/a_M:
p₀ = 1·λ₀⁰/Z = 1/Z
p₁ = a₁·λ₀/Z = n·λ₀/Z

The left tail starts at p₀ = 1/Z, which is TINY for large n.
p₁ = n·λ₀/Z is also small initially since λ₀ is small (when mode is large).

Actually λ₀ = a_{M-1}/a_M can be quite variable. Let me just test the 
DP approach directly.

FOR THE DP: total = dp₀ + dp₁, where dp₁ = x·P.

At the mode-jump of total from M-1 to M:
total[M] ≥ total[k] for all k.

total[M] = dp₀[M] + dp₁[M] = dp₀[M] + P[M-1]
total[M-1] = dp₀[M-1] + dp₁[M-1] = dp₀[M-1] + P[M-2]
These are equal at the tie, and both maximal.

I(T;λ) = dp₀(λ) + λ·P(λ) where dp₀ and P are IS polys of sub-trees.

At the tie point for I(T): λ_T = tie fugacity.

mean_T(λ_T) = Σ_v P_T(v; λ_T) where P_T(v) = occupation probability.

For a tree: P_T(v; λ) = λ I(T-N[v]; λ) / I(T; λ).

mean_T = Σ_v λ·I(T-N[v]) / I(T)

Now I(T) = I(T-v) + λ·I(T-N[v]) for any vertex v.
So λ·I(T-N[v]) = I(T) - I(T-v).
P(v) = 1 - I(T-v)/I(T).

mean_T = n - Σ_v I(T-v;λ)/I(T;λ)

This is nice but doesn't directly help.

Let me try the λ-deformation + induction approach:

INDUCTION: excess(T; λ) < 1 for ALL trees T and ALL λ > 0.

Base: Single vertex: I = 1+x. At λ: p₀ = 1/(1+λ), p₁ = λ/(1+λ).
  Mean = λ/(1+λ), mode = 0 if λ<1, 1 if λ>1, tie at λ=1.
  At tie: mean = 1/2, excess = 1 - 1/2 = 1/2 < 1. ✓

Inductive step: For tree T with root v and children c₁,...,cₖ:
  I(T;λ) = I(T-v;λ) + λ·I(T-N[v];λ)
  = ∏ᵢ I(Tᵢ;λ) + λ·∏ᵢ I(Tᵢ-cᵢ;λ)
  
  where Tᵢ is the subtree rooted at cᵢ.

Hmm, this is the same DP as before just expressed at general λ.

Let me instead try to prove it for the HARD-CORE MODEL directly.

For the hard-core model on a tree at fugacity λ:
The IS-size distribution has mean μ = n·P_avg where P_avg is the 
average occupation probability.

CLAIM: For trees, the variance of the IS-size distribution satisfies
  Var ≥ μ - μ²/n.

Actually, I should just compute the exact relationship between mean and
mode at the tie point and see if the tree DP gives a clean bound.

Let me try the SIMPLEST USEFUL CASE: star graph.
"""

import subprocess
from indpoly import independence_poly
from graph6 import parse_graph6

def poly_stats(poly, lam=1.0):
    """Compute mode, mean, excess at fugacity λ."""
    weights = [poly[k] * (lam ** k) for k in range(len(poly))]
    Z = sum(weights)
    p = [w/Z for w in weights]
    mode = max(range(len(p)), key=lambda k: p[k])
    mean = sum(k * p[k] for k in range(len(p)))
    return mode, mean, mode - mean

def main():
    print("=" * 70)
    print("  STAR GRAPH ANALYSIS")
    print("=" * 70)
    print()
    
    # Star graph S_n: center + n leaves.
    # IS poly: 1 + (n+1)x + C(n,2)x² + ... no wait.
    # IS of S_n: either include center (and exclude all leaves → 1 IS of size 1)
    # Wait, S_n has n+1 vertices: center + n leaves.
    # IS = ∅, {center}, or any subset of leaves (not including center).
    # a₀ = 1, a₁ = n+1 (center or any leaf), 
    # a₂ = C(n,2) (any 2 leaves),
    # a_k = C(n,k) for k ≤ n (k leaves, no center)
    # PLUS a₁ includes center.
    # Wait: a₁ = 1 (center) + n (each leaf) = n+1.
    # a₂ = C(n,2) (pairs of leaves only).
    # ...
    # a_k = C(n,k) for 2 ≤ k ≤ n.
    # And the 1-IS includes center: a₁ = n+1.
    
    # So IS poly = 1 + (n+1)x + Σ_{k=2}^{n} C(n,k) x^k
    #            = 1 + (n+1)x + (1+x)^n - 1 - nx
    #            = (1+x)^n + x
    
    # Check: for n=3 (star with 4 vertices):
    # (1+x)^3 + x = 1 + 3x + 3x² + x³ + x = 1 + 4x + 3x² + x³
    # IS of size 0: {} → 1 ✓
    # IS of size 1: center, leaf1, leaf2, leaf3 → 4 ✓
    # IS of size 2: {l1,l2}, {l1,l3}, {l2,l3} → 3 ✓
    # IS of size 3: {l1,l2,l3} → 1 ✓
    
    print("  Star S_n: IS poly = (1+x)^n + x")
    print()
    
    from math import comb
    
    for n in [3, 5, 10, 20, 50, 100]:
        # Poly: a_0=1, a_1=n+1, a_k=C(n,k) for k≥2
        poly = [comb(n, k) for k in range(n+1)]
        poly[1] += 1  # add x term
        
        # Find tie points
        for M in range(1, n+1):
            if poly[M] > 0 and poly[M-1] > 0:
                lam0 = poly[M-1] / poly[M]
                mode, mean, excess = poly_stats(poly, lam0)
                
                if excess > 0.8:
                    print(f"    S_{n}: M={M}, λ₀={lam0:.4f}, mean={mean:.4f}, excess={excess:.4f}")
    
    print()
    
    # For the star, the tie-point excess max should be findable analytically.
    # IS poly of S_n = (1+x)^n + x.
    # The interesting tie is at M where a_{M-1}/a_M is close to 1.
    
    # For (1+x)^n: a_k = C(n,k). The mode is at k = n/2.
    # The tie points are where C(n,M-1)/C(n,M) = λ₀.
    # C(n,M-1)/C(n,M) = M/(n-M+1).
    # So λ₀ = M/(n-M+1).
    
    # (Ignoring the extra +x term for large n, since a_1 = n+1 ≈ C(n,1)=n)
    
    # At the tie for pure binomial (1+x)^n:
    # p_k ∝ C(n,k)·(M/(n-M+1))^k
    # = C(n,k)·M^k/(n-M+1)^k
    
    # mean of the binomial at λ₀ = nλ₀/(1+λ₀) = n·M/(n+1)
    # mode = M (the jump point)
    # excess = M - nM/(n+1) = M/(n+1)
    
    # For the binomial mode jump at M = ⌊(n+1)/2⌋:
    # excess ≈ n/2 / (n+1) ≈ 1/2
    
    print("  For pure binomial (1+x)^n:")
    print("  At tie point for jump to M: excess = M/(n+1)")
    print("  Max excess at M ≈ n/2: excess ≈ 1/2")
    print()
    
    # Now the KEY: for the product (1+t₁x)(1+t₂x)...(1+tₙx):
    # This is the IS poly of a FOREST of edges (matching).
    # Each factor (1+tᵢx) corresponds to an isolated vertex or edge.
    
    # For a tree on n vertices:
    # I(T;x) is NOT a product of linear factors in general.
    # But by the DP: I(T;x) = I(T-v;x) + x·I(T-N[v];x)
    
    # The key structural constraint is that both I(T-v) and I(T-N[v])
    # are IS polys of smaller trees/forests.
    
    # For the INDUCTION to work, we need:
    # If I(T₁) and I(T₂) both have excess < 1 at all λ,
    # then I(T₁)·I(T₂) has excess < 1 at all λ (for products/forests).
    # AND: I(T₁) + λ·I(T₂) has excess < 1 at all λ (for the merge step).
    
    print("  CHECK: Does the product preserve excess < 1?")
    
    for n in range(5, 14):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        max_product_excess = -100
        for line in lines:
            tn, adj = parse_graph6(line)
            poly1 = independence_poly(tn, adj)
            
            # Multiply by I(P_3) = [1,3,1] (a small example)
            poly2 = [1, 3, 1]
            
            product = [0] * (len(poly1) + len(poly2) - 1)
            for i, a in enumerate(poly1):
                for j, b in enumerate(poly2):
                    product[i+j] += a * b
            
            # Check excess at all tie points
            for M in range(1, len(product)):
                if product[M] > 0 and product[M-1] > 0:
                    lam0 = product[M-1] / product[M]
                    mode, mean, excess = poly_stats(product, lam0)
                    if excess > max_product_excess:
                        max_product_excess = excess
        
        print(f"    n={n}: max product excess = {max_product_excess:.6f}")
    
    print()
    
    # Now the REAL test: does the DP merge preserve excess < 1?
    print("  CHECK: Does the DP merge preserve excess < 1 at all λ?")
    
    for n in range(3, 14):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        max_merge_excess = -100
        
        for line in lines[:200]:  # sample
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            
            for M in range(1, len(poly)):
                if poly[M] > 0 and poly[M-1] > 0:
                    lam0 = poly[M-1] / poly[M]
                    mode, mean, excess = poly_stats(poly, lam0)
                    if excess > max_merge_excess:
                        max_merge_excess = excess
        
        print(f"    n={n}: max tree excess at any tie = {max_merge_excess:.6f}")
    
    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
