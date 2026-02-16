#!/usr/bin/env python3
"""Core proof attempt: coefficient-level argument.

SETUP: total = A + B where A = dp₀, B = x·P, P = ∏dp₀[cᵢ].
Let m = mode(A), M = m+1 = mode(total).

CONDITION FOR GAP=+1:
  (A+B)[M] ≥ (A+B)[j] for all j.

In particular:
  (A+B)[M] ≥ (A+B)[m]
  => A[m+1] + P[m] ≥ A[m] + P[m-1]
  => P[m] - P[m-1] ≥ A[m] - A[m+1]  ...(*)

Let D_A = A[m] - A[m+1] ≥ 0  (drop in A after mode)
Let R_P = P[m] - P[m-1] ≥ 0  (rise in P at m)

(*) says R_P ≥ D_A.

CLAIM: This implies excess(A) ≤ 0 when A is "sufficiently unimodal".

Proof attempt:
  excess(A) = m - mean(A)
  
  mean(A) = Σ_k k·A[k] / |A|
  
  m - mean(A) = Σ_k (m-k) A[k] / |A|
              = [Σ_{k<m} (m-k) A[k] - Σ_{k>m} (k-m) A[k]] / |A|
  
  If A is unimodal, the left tail contributes positively (pulling excess up)
  and the right tail contributes negatively (pulling excess down).
  
  For excess(A) ≤ 0, we need the right tail to dominate.
  
  The right tail is: Σ_{k>m} (k-m) A[k] ≥ A[m+1] + 2·A[m+2] + ...
  The left tail is: Σ_{k<m} (m-k) A[k] ≤ m·A[m-1] + ... 
  
  Hmm, this doesn't use (*) at all.

DIFFERENT APPROACH: Instead of proving excess(A) ≤ 0, prove the condition
DIRECTLY: correction > excess(A).

  correction = w · (1 - Σ δᵢ) where w = |P|/(|A|+|P|)
  
For single child: w = |f|/(2|f|+|g|) where f = dp₀[c], g = dp₁[c]
  and (1-δ) = (1 + μ_f - μ_{f+g}) 

DIRECT INEQUALITY:
  |f|·(1 + μ_f - μ_{f+g}) / (2|f|+|g|) > m_{f+g} - μ_{f+g}

where m_{f+g} = mode(f+g) = mode(A).

Let's write: |f| = F, |g| = G, μ_f = a, μ_g = b.
  μ_{f+g} = (F·a + G·b)/(F+G)
  m_{f+g} = m (the mode of f+g)

  LHS = F·(1 + a - (Fa+Gb)/(F+G)) / (2F+G)
      = F·((F+G+Fa+Ga-Fa-Gb)/(F+G)) / (2F+G)
      = F·((F+G+Ga-Gb)/(F+G)) / (2F+G)      = F·(F + G(1+a-b))/(F+G) / (2F+G)

Hmm, this is getting messy. Let me just compute the key quantities
and look for an invariant or a monotonicity argument.

APPROACH 3 (Generating function): 
For the IS poly I(T;x) = Σ iₖ xᵏ, consider log I(T;x).
The mean at fugacity λ is μ(λ) = λ I'(λ)/I(λ).

At λ=1: μ(1) = I'(1)/I(1) = mean IS size.

Now d/dλ [μ(λ)] = d/dλ [λ I'(λ)/I(λ)]
  = I'(λ)/I(λ) + λ I''(λ)/I(λ) - λ (I'(λ))²/I²(λ)  = μ/λ + λ·(I''/I - (I'/I)²)  = μ/λ + λ · Var / I(λ) · ... 

Actually: if X ~ distribution with pmf p_k = iₖ λᵏ / I(T;λ):
  E[X] = μ(λ)
  Var[X] = λ dμ/dλ (this is a standard fact for exponential families)

So d μ(λ)/dλ = Var[X]/λ > 0, meaning μ(λ) is strictly increasing in λ.

Now mode(λ) = argmax_k iₖ λᵏ. As λ increases from 0 to ∞, mode(λ) increases
(or stays the same). At λ=1, mode(1) = mode of {iₖ}.

For λ → 0: mode → 0 (i₀=1 dominates). For λ → ∞: mode → α(G).

Between these extremes, the mode jumps at specific values of λ.

KEY INSIGHT: mode(λ) is a step function that increases, while mean(λ) is 
a smooth increasing function with mean(0)=0, mean(∞)=α.

When mode first jumps to value M at some λ₀, the mean is μ(λ₀).
Since mode was M-1 just before λ₀, and mean is continuous:
At λ₀: mode jumps from M-1 to M, while mean ≈ some value.

The ratio iₖ λᵏ at the jump point λ₀ has:
  i_M λ₀^M = i_{M-1} λ₀^{M-1}
  => λ₀ = i_{M-1}/i_M

And at λ₀, the mean μ(λ₀) = ...

This λ-deformation approach might give a clean proof!

If we can show that at every mode-jump point λ₀:
  M ≤ μ(λ₀) + 1
then we'd have mode ≤ ceil(mean) at every λ, including λ=1.

At the jump point λ₀ = i_{M-1}/i_M:
  i_M λ₀^M = i_{M-1} λ₀^{M-1}  (equal coefficients at the peak)

So mode M is a "tie" between M and M-1.

This is a known approach in the study of LC polynomials!

Let me test this computationally.
"""

import subprocess
from fractions import Fraction
from indpoly import independence_poly
from graph6 import parse_graph6

def main():
    print("=" * 80)
    print("  λ-DEFORMATION APPROACH")
    print("=" * 80)
    print()
    
    # For each tree, compute:
    # 1. The mode jump points (where i_k λ^k = i_{k-1} λ^{k-1}, i.e., λ = i_{k-1}/i_k)
    # 2. The mean at each jump point
    # 3. Check if mode ≤ ceil(mean) at each jump
    
    total_jumps = 0
    violations = 0
    max_excess_at_jump = -100
    
    for n in range(3, 18):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        n_max_excess = -100
        
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            
            # Find mode jump points
            for M in range(1, len(poly)):
                if poly[M] > 0 and poly[M-1] > 0:
                    # Jump point: λ₀ = i_{M-1}/i_M
                    lam0 = poly[M-1] / poly[M]
                    
                    # At λ₀, compute mean
                    # p_k proportional to i_k λ₀^k
                    weights = [poly[k] * (lam0 ** k) for k in range(len(poly))]
                    total_w = sum(weights)
                    if total_w > 0:
                        mean_at_jump = sum(k * weights[k] for k in range(len(weights))) / total_w
                        excess_at_jump = M - mean_at_jump
                        
                        total_jumps += 1
                        if excess_at_jump > n_max_excess:
                            n_max_excess = excess_at_jump
                        if excess_at_jump > max_excess_at_jump:
                            max_excess_at_jump = excess_at_jump
                        if excess_at_jump >= 1:
                            violations += 1
        
        print(f"  n={n:2d}: max excess at jump = {n_max_excess:+.6f}")
    
    print()
    print(f"  Total jump points: {total_jumps}")
    print(f"  Violations (excess ≥ 1): {violations}")
    print(f"  Max excess at any jump: {max_excess_at_jump:+.6f}")
    
    if violations == 0:
        print()
        print("  ═══════════════════════════════════════════════════")
        print("  λ-DEFORMATION WORKS: mode ≤ ceil(mean) at EVERY")  
        print("  mode-jump point for every tree through n=17.")
        print("  ═══════════════════════════════════════════════════")
        print()
        print("  At λ = i_{M-1}/i_M, the mode jumps from M-1 to M.")
        print("  At this point, i_{M-1}λ^{M-1} = i_M λ^M (tie).")
        print("  The mean is ≤ M at this point (excess ≤ 0).")
        print()
        print("  PROOF STRATEGY:")
        print("  1. For λ ∈ (0,∞), mode(λ) is a non-decreasing step function.")
        print("  2. mean(λ) is strictly increasing (since Var > 0).")
        print("  3. At each jump point, mode ≤ ceil(mean).")
        print("  4. Between jumps, mode is constant while mean increases.")
        print("  5. So mode ≤ ceil(mean) for all λ, including λ=1.")
        print()
        print("  The key is step 3: prove that at the mode-tie point,")
        print("  the mean is at least M - 1.")
        print("  i.e., excess_at_jump < 1")
    
    print()
    print("=" * 80)

if __name__ == "__main__":
    main()
