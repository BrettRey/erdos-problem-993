#!/usr/bin/env python3
"""PRODUCT EXCESS DEEP DIVE.

Key question: what is the tightest bound on excess(∏ IS-polys)?

We know excess(∏) ≤ 0.612 empirically for trees ≤ 9.
But we need a PROOF.

Approach: study what happens to excess under convolution.
For two distributions P, Q with modes M_P, M_Q:
  mode(P*Q) ≤ M_P + M_Q (subadditivity)
  mean(P*Q) = mean(P) + mean(Q)
  excess(P*Q) ≤ excess(P) + excess(Q) (worst case of subadditivity)

BUT: the exact mode of P*Q is usually MUCH closer to M_P + M_Q 
than the sum excess(P) + excess(Q) suggests, because the 
convolution concentrates around the sum of means.

For UNIMODAL distributions:
  Wintner (1938): convolution of unimodal distributions is unimodal.
  Medgyessy: the mode of the convolution of two symmetric unimodal 
  distributions is at the sum of the modes.

For NON-SYMMETRIC unimodal:
  The mode can shift away from the sum of modes, but by how much?

Actually, I think a simpler approach works.

THEOREM (folklore/Darroch extension): For the convolution of 
distributions from "product-type" polynomials (where P(x) = ∏(1+rᵢx)), 
excess ≤ ? 

Hmm, but tree IS polys are NOT products of linear factors.

ALTERNATIVE: Just verify the product bound for larger trees and 
check if it holds.

Actually, let me re-examine the question. The REAL question is:

We have: excess(T) = M - mean(T) with the DP formula:
  excess(T) = 1 + excess(dp₀) - correction  [at gap=+1]
  excess(T) < excess(dp₀)                    [at gap≤0]

And: correction = w·(μ₁ - μ₀) where w = P(root).

At gap=+1: excess(T) = 1 + excess(dp₀) - correction = 1 - (correction - excess(dp₀))

Define: slack = correction - excess(dp₀) = w·(1 + μ_B - μ_A) - (M₀ - μ_A)
       where M₀ = mode(dp₀), μ_A = mean(dp₀), μ_B = mean(dp₁)-1.

Wait: μ₁ = mean(dp₁) = 1 + mean(∏ dp₀[cᵢ]) (shifted by x).
And μ₀ = mean(dp₀) = mean(∏ total[cᵢ]).

correction = w·(μ₁ - μ₀) = w·(1 + mean(∏ dp₀[cᵢ]) - mean(∏ total[cᵢ]))

Since mean(∏ dp₀[cᵢ]) = Σ mean(dp₀[cᵢ]) and similarly for total:
correction = w·(1 + Σ mean(dp₀[cᵢ]) - Σ mean(total[cᵢ]))
            = w·(1 - Σ [mean(total[cᵢ]) - mean(dp₀[cᵢ])])

And mean(total[cᵢ]) - mean(dp₀[cᵢ]) = δᵢ > 0 (dp₁ always shifts mean up).

So correction = w·(1 - Σ δᵢ).

And: excess(dp₀) = mode(∏ total[cᵢ]) - mean(∏ total[cᵢ])
                  ≤ Σ mode(total[cᵢ]) - Σ mean(total[cᵢ])
                  = Σ excess(total[cᵢ])

BUT the inequality mode(∏) ≤ Σ mode is NOT tight in general.
The mode of a convolution is typically between Σ⌊mean⌋ and Σ⌈mean⌉.

Hmm, let me try a completely different tactic.

THE λ-DEFORMATION PROOF (attempt at rigor):

At any fugacity λ > 0, define:
  p_k(λ) = i_k λ^k / I(T;λ)  (IS-size distribution)
  μ(λ) = Σ k·p_k(λ)           (mean)
  mode(λ) = argmax_k p_k(λ)   (mode)

FACTS:
  1. μ(λ) is continuous and strictly increasing (derivative = Var/λ > 0).
  2. mode(λ) is non-decreasing, piecewise constant (step function).
  3. At λ=0: μ=0, mode=0.
  4. As λ→∞: μ→α(T) (IS number), mode→α(T).

At every mode-jump point λ₀ (where mode jumps from M-1 to M):
  p_{M-1}(λ₀) = p_M(λ₀) (tie).
  
We need to prove: μ(λ₀) > M - 1.

Now at λ₀, using the tree DP for trees rooted at any vertex v:
  I(T;λ) = I(T-v;λ) + λ·I(T-N[v];λ)

  μ_T(λ) = (1-P(v))·μ_{T-v}(λ) + P(v)·(1 + μ_{T-N[v]}(λ))

where P(v;λ) = λ·I(T-N[v];λ)/I(T;λ).

By induction: μ_{T-v}(λ) > mode_{T-v}(λ) - 1 (if the claim holds for T-v).
Similarly: μ_{T-N[v]}(λ) > mode_{T-N[v]}(λ) - 1 (for T-N[v]).

But: mode_{T-v}(λ) and mode_{T-N[v]}(λ) can be anything relative to M.

THE TRICK: Choose v to be a LEAF of T. Then:
  T-v is a tree on n-1 vertices.
  T-N[v] = T minus v and its parent = forest on n-2 vertices (or tree if v 
  is the only child of its parent).

For a LEAF v with parent u:
  I(T) = I(T-v) + λ·I(T-{v,u}-E_u)

Wait, T-N[v] = T minus the closed neighborhood of v = T minus {v, parent(v)}.
This is a forest (removing a leaf and its parent from a tree).

Actually T-N[v] is T minus v and all neighbors of v. Since v is a leaf, 
N(v) = {parent(v)}. So N[v] = {v, parent(v)}, and T-N[v] = T - {v, parent}.

T - {v, parent} is a forest (possibly disconnected).

Hmm, this complicates the induction because T-N[v] might not be a tree.
But it IS a FOREST, and the IS poly of a forest is the PRODUCT of 
IS polys of its components (trees). So the induction applies to each 
component separately.

FOR A FOREST F = T₁ ∪ ... ∪ Tₖ:
  μ_F = Σ μ_{Tᵢ}
  mode_F ≤ Σ mode_{Tᵢ}
  By induction: μ_{Tᵢ} > mode_{Tᵢ} - 1.
  So: μ_F > Σ mode_{Tᵢ} - k ≥ mode_F - k.

That gives μ_F > mode_F - k, which is WEAKER than μ_F > mode_F - 1.
The issue is that the product excess can be as large as the SUM of 
individual excesses.

BUT WAIT: we know empirically that product excess ≤ 0.61, not k times anything.
And theoretically, for the sum of independent rv's, 
Darroch's result gives |mode - mean| ≤ max_i p_i for INTEGRAL distributions.
Not quite right for non-integral.

ACTUALLY: There's a theorem by Darroch (1964) that says for the 
convolution of ANY unimodal distributions, |mode - mean| < 1.
Wait, is that right? Let me check.

Darroch's original result is specifically for sums of independent 
Bernoullis. But there might be extensions.

Let me check: is IS-poly an EXPONENTIAL FAMILY? Yes!
P(|S|=k; λ) = i_k λ^k / I(T;λ) is an exponential family in λ.
For exponential families, the mode and mean are always close 
(within 1 for integer-valued distributions).

Wait — is mode within 1 of mean for ALL exponential families?

For Poisson(μ): mode = ⌊μ⌋ or ⌈μ⌉ - 1. So |mode - mean| < 1 ✓
For Binomial(n,p): mode = ⌊(n+1)p⌋ or ⌊(n+1)p⌋-1. |mode-mean|<1 ✓
For NegBin: similar.

Is there a GENERAL theorem that for any exponential family 
{p(k;θ) = h(k) exp(θk - A(θ))} with k ∈ ℤ≥0, |mode - mean| < 1?

The IS-size distribution at fugacity λ is:
  p(k) = i_k λ^k / Σ i_j λ^j  = i_k exp(k log λ) / Z
  
This is a one-parameter exponential family with natural parameter θ = log λ
and sufficient statistic T = k.

For ANY one-parameter exponential family with integer-valued outcome 
k ∈ {0,1,...,d}, is |mode(θ) - mean(θ)| < 1?

I believe the answer is YES if the distribution is LOG-CONCAVE 
(which it IS for IS polys by the Alavi-Malde-Schwenk or 
Hamidoune theorem for trees and claw-free).

Wait, IS polys of trees ARE log-concave? Yes, this was proved by 
Alavi, Malde, Schwenk, and Erdős (1987): IS polys of trees are 
log-concave (in the coefficients).

And for log-concave distributions on integers:
|mode - mean| < 1.

Is this a known theorem? I think YES.

BAGNOLI & BERGSTROM (2005): For a log-concave probability mass function
on the integers, |mode - mean| < 1.

Wait, I need to check this carefully. The statement I need is:
"For a log-concave distribution on {0,1,...,d}, mode ∈ {⌊μ⌋, ⌈μ⌉}."

This would give mode ≤ ⌈μ⌉, which is exactly our claim!

Let me verify this computationally with a decisive test.
"""

import numpy as np

def is_log_concave(p):
    """Check if p is log-concave (p_k^2 >= p_{k-1}*p_{k+1})."""
    for k in range(1, len(p)-1):
        if p[k] > 0 and p[k-1] >= 0 and p[k+1] >= 0:
            if p[k]**2 < p[k-1]*p[k+1] - 1e-12:
                return False
    return True

def main():
    print("=" * 70)
    print("  TEST: For LC distributions, mode ∈ {⌊μ⌋, ⌈μ⌉}?")
    print("=" * 70)
    print()
    
    # Generate random LC distributions
    import random
    random.seed(42)
    
    violations = 0
    total = 0
    
    for trial in range(100000):
        # Generate LC distribution: start with a log-concave sequence
        # Method: generate log p_k as concave sequence
        d = random.randint(3, 20)
        # Concave log-probabilities
        peak = random.randint(1, d-1)
        slope_left = random.uniform(0.1, 2.0)
        slope_right = random.uniform(0.1, 2.0)
        
        log_p = []
        for k in range(d+1):
            if k <= peak:
                log_p.append(-slope_left * (peak - k))
            else:
                log_p.append(-slope_right * (k - peak))
        
        # Add small noise that preserves concavity
        p = [np.exp(lp) for lp in log_p]
        Z = sum(p)
        p = [x/Z for x in p]
        
        if not is_log_concave(p):
            continue
        
        total += 1
        mode = max(range(len(p)), key=lambda k: p[k])
        mean = sum(k*p[k] for k in range(len(p)))
        
        import math
        if mode != math.floor(mean) and mode != math.ceil(mean):
            violations += 1
            if violations <= 5:
                print(f"  VIOLATION: d={d}, mode={mode}, mean={mean:.4f}, ⌊μ⌋={math.floor(mean)}, ⌈μ⌉={math.ceil(mean)}")
    
    print(f"\n  Total LC distributions tested: {total}")
    print(f"  Violations: {violations}")
    if violations == 0:
        print("  ✓ mode ∈ {⌊μ⌋, ⌈μ⌉} for all tested LC distributions!")
    
    print()
    
    # Now test with specific counterexample from before: geometric LC
    print("  TEST with geometric LC (the known counterexample family):")
    import math
    for s_prime in [0.3, 0.4, 0.45, 0.49, 0.5, 0.6, 0.8, 0.95]:
        for M in [5, 10, 20]:
            # Geometric LC: p_k = C * s'^|M-k| for k near M, with tie at M-1,M
            # Build finite distribution
            p = []
            for k in range(2*M + 5):
                if k <= M-1:
                    p.append(s_prime ** (M-1-k))
                elif k == M:
                    p.append(1.0)  # tie: p_M = p_{M-1} = 1.0 unnormalized
                else:
                    # Right tail: let it decay faster
                    r = s_prime * 0.5  # right decay
                    p.append(r ** (k - M))
            
            Z = sum(p)
            p = [x/Z for x in p]
            
            lc = is_log_concave(p)
            mode = max(range(len(p)), key=lambda k: p[k])
            mean = sum(k*pk for k, pk in enumerate(p))
            
            in_range = mode == math.floor(mean) or mode == math.ceil(mean)
            
            if not in_range:
                print(f"  s'={s_prime:.2f} M={M}: LC={lc}, mode={mode}, mean={mean:.4f}, "
                      f"⌊μ⌋={math.floor(mean)}, ⌈μ⌉={math.ceil(mean)} {'✗' if not in_range else '✓'}")
    
    print()
    
    # Test with the EXACT geometric counterexample from proof_symbolic.py
    # s'=0.5, r=0 (no right tail)
    print("  EXACT geometric counterexample (s'=0.5, r=0):")
    for M in [5, 10, 20, 50]:
        p = []
        for k in range(M+2):
            if k <= M-1:
                p.append(0.5 ** (M-1-k))
            elif k == M:
                p.append(1.0)
            else:
                p.append(0.0)
        
        Z = sum(p)
        p = [x/Z for x in p]
        lc = is_log_concave(p)
        mode = max(range(len(p)), key=lambda k: p[k])
        mean = sum(k*pk for k, pk in enumerate(p))
        excess = mode - mean
        in_range = mode == math.floor(mean) or mode == math.ceil(mean)
        print(f"    M={M:3d}: LC={lc}, mode={mode}, mean={mean:.4f}, excess={excess:.4f}, "
              f"in_range={in_range}")
    
    print()
    
    # KEY THEOREM CHECK:
    # For a LC distribution p on {0,...,d}: mode ∈ {⌊μ⌋, ⌈μ⌉}?
    # 
    # Equivalent to: |mode - mean| < 1.
    #
    # This is KNOWN to be TRUE for:
    # - Poisson distributions
    # - Binomial distributions  
    # - Negative binomial distributions
    # - Products of 1-parameter exponential families with LC base measure
    #
    # For GENERAL LC distributions... I believe this is TRUE but need reference.
    #
    # PROOF SKETCH: If p is LC on integers and mode = M:
    # p_M ≥ p_k for all k.
    # mean = Σ k·p_k = M - excess.
    # 
    # Claim: excess = M - mean ∈ [0, 1).
    #
    # Lower bound (excess ≥ 0 when mode is the smallest peak):
    # Not always true. Mode can be below mean (excess < 0).
    # We just need: |excess| < 1.
    
    print("  ATTEMPTING PROOF: For LC dist, |mode - mean| < 1")
    print()
    print("  KEY STEP: Bounding excess = mode - mean for LC distributions.")
    print()
    print("  For a unimodal distribution p on {0,...,d}:")
    print("  Let p_M = max probability. At the MODE M:")
    print("  p_M ≥ p_k for all k.")
    print("  Excess = M - mean = Σ_{k≠M} (M-k) p_k")
    print("         = Σ_{k<M} (M-k) p_k - Σ_{k>M} (k-M) p_k")
    print()
    print("  By LC: p_{M-j}/p_M ≤ s^j for some s ≤ 1 (non-increasing ratios)")
    print("  Similarly right tail.")
    print()
    print("  Left contribution: Σ_{j=1}^{M} j · p_{M-j} ≤ p_M · Σ j·s^j = p_M·s/(1-s)²")
    print("  Right contribution: Σ_{j=1}^{d-M} j · p_{M+j} ≥ 0")
    print()
    print("  If the right tail exactly balances the left, excess can approach p_M/(1-s)²")
    print("  which can be > 1 for s > 0.382.")
    print()
    print("  BUT: p_M ≤ 1/(1 + s/(1-s) + r/(1-r))) = (1-s)(1-r)/(2-s-r)")
    print("  So excess ≤ (1-r)/[(2-s-r)(1-s)] = 1/[(2-s)(1-s)] when r=0")
    print("  This is > 1 for s > 0.382. ✗")
    print()
    print("  Conclusion: LC ALONE does NOT guarantee |mode - mean| < 1.")
    print("  Need additional structure (e.g., IS polynomial structure).")
    
    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
