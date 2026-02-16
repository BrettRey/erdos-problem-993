#!/usr/bin/env python3
"""THEOREM: For a log-concave polynomial P(x) = Σ aₖ xᵏ with nonneg coefficients
and no internal zeros, the mode M satisfies M ≤ ⌈mean⌉.

Equivalently: mode - mean < 1.

PROOF via λ-deformation:

Let fₖ(λ) = aₖ λᵏ. The weighted distribution p_k(λ) = f_k(λ) / P(λ).

Key facts:
  1. mean(λ) = Σ k p_k(λ) = λ P'(λ)/P(λ) is strictly increasing in λ.
     (Its derivative is Var(λ)/λ > 0.)
  
  2. mode(λ) = argmax_k aₖ λᵏ is non-decreasing in λ.
     (As λ increases, higher k terms grow faster.)
  
  3. If {aₖ} is LC, then {aₖ λᵏ} is also LC for any λ > 0.
     (LC is preserved under multiplication by geometric sequences.)

  4. At λ=0: mode(0) = 0, mean(0) = 0.
     At λ→∞: mode(∞) = deg(P), mean(∞) = deg(P).

  5. mode(λ) = max{k : aₖ λᵏ ≥ aₖ₋₁ λᵏ⁻¹} = max{k : aₖ/aₖ₋₁ ≥ 1/λ}
     Since ratios rₖ = aₖ/aₖ₋₁ are non-increasing (by LC):
     mode(λ) = max{k : r_k ≥ 1/λ}.

  6. The mode jumps from M-1 to M at λ₀ = 1/r_M = a_{M-1}/a_M.

AT THE JUMP POINT λ₀:
  {p_k(λ₀)} is LC (by fact 3) with two adjacent maxima at M-1 and M.
  
  CLAIM: For an LC distribution with adjacent modes at M-1, M:
    mean ∈ [M-1, M].

  This is the key lemma. Let me prove it.

PROOF OF CLAIM:
  By LC, the ratios p_k/p_{k-1} are non-increasing.
  Since p_M = p_{M-1}: ratio at M is exactly 1.
  For k > M: p_k/p_{k-1} ≤ 1, so the sequence is non-increasing.
  For k < M: p_k/p_{k-1} ≤ p_{k+1}/p_k ≤ ... ≤ p_M/p_{M-1} = 1.

  Wait, that's wrong. LC says p_k² ≥ p_{k-1} p_{k+1}, which gives
  p_{k+1}/p_k ≤ p_k/p_{k-1}. I.e., ratios are NON-INCREASING.

  So:
  r_M = p_M/p_{M-1} = 1
  r_{M+1} = p_{M+1}/p_M ≤ r_M = 1    → p_{M+1} ≤ p_M ✓
  r_{M-1} = p_{M-1}/p_{M-2} ≥ r_M = 1 → p_{M-1} ≥ p_{M-2} ✓

  So the distribution is unimodal with plateau {M-1, M}.

  Now bound the mean:
  
  UPPER BOUND (mean ≤ M):
    Define S_k = Σ_{j≥k} p_j (tail sum). Then mean = Σ_{k≥1} S_k.
    Equivalently: mean = Σ_k k·p_k.
    
    M - mean = Σ_k (M-k) p_k = Σ_{k<M} (M-k) p_k - Σ_{k>M} (k-M) p_k
    
    Since the distribution is peaked at {M-1, M} and decreases both ways:
    The right tail Σ_{k>M} (k-M) p_k pulls mean UP (toward M or above).
    
    Actually, it's cleaner to use the generating function approach.
    
    Since p_M = p_{M-1}, the symmetry center is at M - 1/2.
    
    For LC distributions with plateau at {M-1, M}:
    Σ_{k<M-1} p_k ≤ Σ_{k>M} p_k  ← IS THIS TRUE?
    
    No, not necessarily. Left tail can be heavier.

  Let me try a DIRECT proof.

  LOWER BOUND (mean ≥ M - 1):
    mean - (M-1) = Σ_k (k - M + 1) p_k
    
    Split: 
    = Σ_{k≤M-2} (k-M+1) p_k + 0·p_{M-1} + 1·p_M + Σ_{k≥M+1} (k-M+1) p_k
    
    = -Σ_{j=1}^{M-1} j·p_{M-1-j} + p_M + Σ_{j=1}^{α-M} j·p_{M+j}
    
    Need: p_M + Σ_{j≥1} j·p_{M+j} ≥ Σ_{j=1}^{M-1} j·p_{M-1-j}
    
    Since p_{M-1} = p_M, pair up terms:
    For j ≥ 1: j·p_{M+j} vs j·p_{M-1-j}
    
    By LC: p_{M+j}/p_M ≤ (p_{M+1}/p_M)^j   (LC gives geometric decay)
    Similarly: p_{M-1-j}/p_{M-1} ≤ (p_{M-2}/p_{M-1})^j
    
    Since p_{M+1}/p_M ≤ 1 and p_{M-2}/p_{M-1} ≤ 1:
    Both tails decay at least geometrically.
    
    The KEY is comparing the decay rates. Let:
    r = p_{M+1}/p_M (right decay rate, ≤ 1)
    s = p_{M-2}/p_{M-1} (left decay rate, ≤ 1)
    
    By LC, r·s ≤ 1 (from p_{M-1}² ≥ p_{M-2}·p_M and p_M² ≥ p_{M-1}·p_{M+1}).
    Actually: p_{M-1} = p_M. From p_M² ≥ p_{M-1}·p_{M+1} → p_M ≥ p_{M+1}, i.e., r ≤ 1. ✓
    From p_{M-1}² ≥ p_{M-2}·p_M → p_{M-1} ≥ p_{M-2}, i.e., s ≤ 1. ✓
    
    From the LC condition at M-1: p_{M-1}² ≥ p_{M-2} p_M
    → 1 ≥ s (since p_M = p_{M-1}).
    
    From LC at M: p_M² ≥ p_{M-1} p_{M+1}
    → p_M ≥ p_{M+1}, i.e., r ≤ 1. ✓
    
    But we DON'T get r ≤ s or r ≥ s from these.
    
    HOWEVER: from LC applied to the whole range:
    The ratios p_k/p_{k-1} are non-increasing.
    
    r_{M+1} = p_{M+1}/p_M ≤ r_M = 1
    r_{M-1} = p_{M-1}/p_{M-2} ≥ r_M = 1
    
    So: p_{M-2} ≤ p_{M-1} = p_M ≥ p_{M+1}
    
    But we need to compare p_{M+j} with p_{M-1-j} for j ≥ 1.
    
    Since ratios are non-increasing:
    p_{M+j}/p_{M+j-1} ≤ p_{M+1}/p_M = r
    So p_{M+j} ≤ r^j · p_M
    
    Similarly: 
    p_{M-1-j}/p_{M-j} ≤ p_{M-1}/p_M = ... wait, 
    p_k/p_{k-1} is non-increasing, so for k < M:
    p_k/p_{k-1} ≥ p_{k+1}/p_k
    
    At k = M-1: p_{M-1}/p_{M-2} ≥ p_M/p_{M-1} = 1
    At k = M-2: p_{M-2}/p_{M-3} ≥ p_{M-1}/p_{M-2} ≥ 1
    ...
    So for k < M-1: p_k ≤ p_{k+1}, i.e., the left side is increasing.
    
    More precisely: p_{M-1-j} ≤ p_{M-j} for j ≥ 1.
    And the decay rate going LEFT is:
    p_{M-1-j}/p_{M-j} ≤ p_{M-j}/p_{M-j+1} ≤ ... ≤ p_{M-1}/p_M = 1
    
    Hmm, this says the left tail is also bounded by p_M. Not helpful.
    
    Let me try a DIFFERENT APPROACH entirely.

    USE THE VARIANCE-COVARIANCE IDENTITY:
    
    For the LC distribution at λ₀:
    Var = Σ k² p_k - μ²
    
    Since d(mean)/dλ = Var/λ > 0, mean is increasing.
    
    The mean at λ₀ is between mean(0)=0 and mean(∞)=d.
    Mode jumps to M at λ₀.
    
    Just before λ₀ (at λ₀⁻): mode was M-1, and mean was close to μ₀.
    Just after λ₀ (at λ₀⁺): mode is M, and mean ≈ μ₀ (continuous).
    
    Before the jump: excess = (M-1) - μ < 1 by induction hypothesis.
    So μ > M - 2.
    
    At the jump: mode goes to M, mean stays ≈ μ.
    So excess ≈ M - μ < M - (M-2) = 2.
    
    That only gives excess < 2. We need < 1.
    
    AH BUT WAIT: the induction hypothesis for the PREVIOUS jump is:
    At the previous jump (from M-2 to M-1, at some λ₁ < λ₀):
      excess at λ₁ < 1, so mean(λ₁) > M - 2.
    
    Between λ₁ and λ₀: mean increased (since Var/λ > 0).
    So mean(λ₀) > mean(λ₁) > M - 2.
    
    NEW INFORMATION: Between jumps, how much does mean increase?
    
    mean(λ₀) - mean(λ₁) = ∫_{λ₁}^{λ₀} Var(λ)/λ dλ
    
    This is POSITIVE. But can we bound it from below?
    
    YES! Between jumps, the mode is constant at M-1.
    The distribution is peaked at M-1.
    Var(λ) ≥ some minimum value.
    
    Actually, at λ₁: mean > M-2, mode = M-1, so excess < 1.
    As λ increases, mean increases through [M-2, ...].
    The mode stays at M-1 until λ₀.
    
    At λ₀: mean has increased to some value > M-2.
    
    We need mean(λ₀) > M - 1 for the new excess to be < 1.
    
    So we need: mean(λ₀) - mean(λ₁) > 1 - excess(λ₁)?
    Not necessarily, since mean(λ₁) could be close to M-2.
    
    THIS APPROACH NEEDS A LOWER BOUND ON ∫ Var/λ dλ.

Let me try yet another angle.

DIRECT LC PROOF (algebraic):

For an LC sequence with p_{M-1} = p_M, show Σk·p_k ≥ M - 1.

Σk·p_k = Σ_{k=0}^{α} k·p_k ≥ (M-1)·p_{M-1} + M·p_M + Σ_{k>M} k·p_k
       ≥ (2M-1)·p_M + (M+1)·S_R

where S_R = Σ_{k>M} p_k.

So mean ≥ (2M-1) p_M + (M+1) S_R.

For mean ≥ M-1:
  (2M-1) p_M + (M+1) S_R ≥ M - 1
  
Since 2p_M + S_L + S_R = 1 (where S_L = Σ_{k<M-1} p_k):
  p_M = (1 - S_L - S_R) / 2

  (2M-1)(1-S_L-S_R)/2 + (M+1) S_R ≥ M - 1
  (2M-1)(1-S_L-S_R) + 2(M+1) S_R ≥ 2(M-1)
  (2M-1) - (2M-1)S_L - (2M-1)S_R + 2(M+1)S_R ≥ 2M-2
  1 - (2M-1)S_L + (2M+2-2M+1)S_R ≥ 0
  1 - (2M-1)S_L + 3S_R ≥ 0

So: mean ≥ M-1 iff (2M-1)S_L ≤ 1 + 3S_R.

For large M with S_L ≈ 0.5, this needs M ≤ (1+3S_R)/(2S_L).
That's not always satisfied for generic unimodal sequences.

BUT for LC sequences, S_L is constrained!

Under LC with plateau {M-1, M}:
  p_k ≤ p_M for all k (unimodality).
  S_L = Σ_{k<M-1} p_k ≤ (M-1) p_M (at most M-1 terms, each ≤ p_M).
  
  So (2M-1) S_L ≤ (2M-1)(M-1) p_M.
  
  This can be ≫ 1 + 3S_R for large M. So this bound alone isn't tight enough.

  We need the geometric decay from LC. Under LC:
  p_{M-1-j} ≤ s^j · p_M where s = p_{M-2}/p_{M-1} ≤ 1.
  
  S_L = Σ_{j=1}^{M-1} p_{M-1-j} ≤ p_M · Σ s^j ≤ p_M · s/(1-s) (for s < 1)
  
  Wait, but we also need S_R.
  S_R = Σ_{j=1}^{α-M} p_{M+j} ≤ p_M · r/(1-r) where r = p_{M+1}/p_M ≤ 1.
  
  Under LC: p_M² ≥ p_{M-1}·p_{M+1} → p_M ≥ p_{M+1} (since p_{M-1}=p_M), r ≤ 1.
  Also: p_{M-1}² ≥ p_{M-2}·p_M → p_{M-1} ≥ p_{M-2}, s ≤ 1.
  
  AND the KEY LC constraint is:
  Ratios are non-increasing: r_{k+1} ≤ r_k.
  At k=M: r_{M+1} = p_{M+1}/p_M ≤ r_M = p_M/p_{M-1} = 1.
  At k=M-1: r_M = 1 ≤ r_{M-1} = p_{M-1}/p_{M-2} = 1/s.
  So s ≤ 1. ✓
  
  MORE: From the decreasing ratios, for j ≥ 1:
  p_{M+j}/p_{M+j-1} ≤ p_{M+1}/p_M = r
  So p_{M+j} ≤ r^j · p_M.
  
  Similarly: p_{M-j}/p_{M-j+1} ≤ p_{M-1}/p_M = 1 (going left, ratios decrease).
  Wait: ratios r_k = p_k/p_{k-1} are non-increasing.
  For k ≤ M-1: r_k ≥ r_{k+1} ≥ ... ≥ r_M = 1.
  So p_k ≥ p_{k-1} for k ≤ M-1 (increasing on left). ✓
  
  For the LEFT decay rate: going from M-1 backward:
  p_{M-2}/p_{M-1} = 1/r_{M-1} = ... hmm.
  
  r_{M-1} = p_{M-1}/p_{M-2} ≥ 1.
  Define s' = 1/r_{M-1} = p_{M-2}/p_{M-1} ≤ 1.
  
  Then p_{M-2} = s' · p_{M-1} = s' · p_M.
  
  Since ratios are non-increasing: r_{M-2} ≥ r_{M-1}
  So p_{M-2}/p_{M-3} ≥ p_{M-1}/p_{M-2}
  i.e., p_{M-3} ≤ p_{M-2}²/p_{M-1} = (s')² · p_{M-1} = (s')² · p_M
  
  Similarly: p_{M-1-j} ≤ (s')^j · p_M.
  
  BEAUTIFUL! So S_L ≤ p_M · Σ_{j=1}^{M-1} (s')^j ≤ p_M · s'/(1-s').
  
  AND: from the LC condition at position M-1:
  p_{M-1}² ≥ p_{M-2} · p_M
  p_M² ≥ s' · p_M · p_M
  p_M ≥ s' · p_M ✓ (just s' ≤ 1)
  
  From the LC condition at position M:
  p_M² ≥ p_{M-1} · p_{M+1} = p_M · r · p_M
  ∴ r ≤ 1 ✓
  
  Now the CRUCIAL LC inequality at M-1 AND M together:
  The ratio r_{M-1} = p_{M-1}/p_{M-2} = 1/s'.
  The ratio r_M = p_M/p_{M-1} = 1.  
  The ratio r_{M+1} = p_{M+1}/p_M = r.
  
  Non-increasing ratios: r ≤ 1 ≤ 1/s'.
  So r·s' ≤ 1. 
  
  THIS CONSTRAINS THE RELATIONSHIP BETWEEN LEFT AND RIGHT DECAY!
  If right tail decays fast (r small), left tail can decay slow (s' close to 1).
  If left tail decays slow (s' close to 1), right tail must decay fast (r ≤ 1/s' ≈ 1).
  
  Actually r·s' ≤ 1 just says r ≤ 1/s', which since s' ≤ 1, says r ≤ 1/s' ≥ 1.
  Not useful beyond r ≤ 1 which we already know.
  
  ACTUALLY: More precisely, from non-increasing ratios:
  r = r_{M+1} ≤ r_M = 1  AND  1 = r_M ≤ r_{M-1} = 1/s'.
  
  But is there a stronger constraint? YES:
  From r_{M+1} ≤ r_M AND r_{M-1} ≥ r_M, we get:
  r ≤ 1 ≤ 1/s', i.e., s' ≤ 1 ≤ 1/r.
  i.e., r · s' ≤ 1 · 1 = 1? No: r ≤ 1 and 1/s' ≥ 1, but r and 1/s' are not 
  directly related.
  
  WAIT: From non-increasing ratios MORE CAREFULLY:
  r_1 ≥ r_2 ≥ ... ≥ r_M = 1 ≥ r_{M+1} ≥ ...
  
  This means ALL left ratios are ≥ 1 and ALL right ratios are ≤ 1.
  
  So the left tail increases (more slowly as we get further from center)
  and the right tail decreases. The sequence is still unimodal.
  
  But crucially: p_{M-1-j} ≤ (s')^j · p_M where s' ≤ 1
  AND: p_{M+j} ≤ r^j · p_M where r ≤ 1.
  
  Now: is there a constraint between s' and r that helps?
  Actually yes! From the LC condition at M-1:
  p_{M-1}² ≥ p_{M-2}·p_M → p_M² ≥ s'·p_M·p_M → s' ≤ 1 ✓
  
  From LC at M: p_M² ≥ p_{M-1}·p_{M+1} → p_M ≥ r·p_M → r ≤ 1 ✓
  
  From LC at the "hinge" between left and right, there's a THIRD constraint.
  Consider the LC condition at position M-1/2 (between M-1 and M):
  p_{M-1}² ≥ p_{M-2}·p_M → s' ≤ 1 (already used)
  p_M² ≥ p_{M-1}·p_{M+1} → r ≤ 1 (already used)
  
  These DON'T link s' and r directly.
  
  BUT: from LC-like conditions further out:
  p_{M-2}² ≥ p_{M-3}·p_{M-1}
  (s'·p_M)² ≥ p_{M-3}·p_M
  (s')²·p_M ≥ p_{M-3}
  But we already know p_{M-3} ≤ (s')²·p_M from the geometric bound. ✓
  
  So there's NO direct link between s' and r from LC alone.
  
  CONCLUSION: We need a different proof technique.

Let me try something completely different: THE MEAN AT THE TIE IS THE 
MIDPOINT MINUS THE "NET TORQUE".

At the tie point, p_{M-1} = p_M = p* (say).
The rest sums to 1 - 2p*.

mean = (M-1)p* + Mp* + Σ_{k≠M-1,M} k·p_k
     = (2M-1)p* + Σ_{k<M-1} k·p_k + Σ_{k>M} k·p_k

mean - (M - 1/2) = Σ_{k<M-1} (k-M+1/2)·p_k + Σ_{k>M} (k-M+1/2)·p_k

I need mean ≥ M-1, i.e., mean - (M-1/2) ≥ -1/2.

Let L = Σ_{k<M-1} (M-1/2-k)·p_k (left torque, positive)
Let R = Σ_{k>M} (k-M+1/2)·p_k   (right torque, positive)

mean - (M-1/2) = R - L.

Need: R - L ≥ -1/2, i.e., L - R ≤ 1/2.

Now: L = Σ_{j≥1} (j+1/2)·p_{M-1-j} ≤ Σ (j+1/2)·(s')^j·p*
                                       = p* · [3/2 s' + 5/2 (s')² + 7/2 (s')³ + ...]
                                       = p* · [s'(3+5s'+7(s')²+...)/2]
                                       = p* · s' · [Σ_{j≥0} (2j+3)(s')^j] / 2
                                       
And R = Σ_{j≥1} (j+1/2)·p_{M+j} ≥ 3/2·p_{M+1} = 3/2·r·p*

Hmm, the bound L ≤ p*·s'/(1-s')² · something and R ≥ some function of r...

This is getting very messy. Let me just try to verify computationally
whether a simple inequality like r + s' ≥ 1 holds (which would constrain 
the coupling between left and right tails).
"""

import random
import math

def make_lc_poly(factors):
    poly = [1]
    for t in factors:
        new_poly = [0] * (len(poly) + 1)
        for i, c in enumerate(poly):
            new_poly[i] += c
            new_poly[i+1] += c * t
        poly = new_poly
    return poly

def main():
    print("=" * 80)
    print("  LC TIE-POINT: Check relationship between r and s'")
    print("=" * 80)
    print()
    
    random.seed(42)
    
    # For LC polynomials, at the tie point λ₀ = a_{M-1}/a_M:
    # Weighted seq p_k = a_k λ₀^k / Z
    # p_{M-1} = p_M = p*
    # s' = p_{M-2}/p* (left decay starting factor)
    # r  = p_{M+1}/p* (right decay starting factor)
    
    # Check what constraints exist between r and s'
    
    min_r_plus_s = 2.0
    
    data = []
    
    for trial in range(50000):
        k = random.randint(3, 25)
        factors = [random.uniform(0.01, 50) for _ in range(k)]
        poly = make_lc_poly(factors)
        d = len(poly) - 1
        
        for M in range(2, d):
            if poly[M] > 0 and poly[M-1] > 0:
                lam0 = poly[M-1] / poly[M]
                weights = [poly[j] * (lam0 ** j) for j in range(d+1)]
                Z = sum(weights)
                p = [w/Z for w in weights]
                
                # Check this is a valid tie
                if abs(p[M] - p[M-1]) > 1e-10 * max(p[M], p[M-1]):
                    continue
                
                pstar = p[M]
                
                if pstar < 1e-12:
                    continue
                
                if M >= 2 and M < d:
                    s_prime = p[M-2] / pstar
                    r = p[M+1] / pstar
                    
                    if r + s_prime < min_r_plus_s:
                        min_r_plus_s = r + s_prime
                    
                    # Compute mean
                    mean_val = sum(j * p[j] for j in range(d+1))
                    excess = M - mean_val
                    
                    if excess > 0.9:
                        data.append((excess, r, s_prime, M, d, pstar))
    
    print(f"  Min (r + s'): {min_r_plus_s:.6f}")
    print()
    
    if data:
        data.sort(reverse=True)
        print(f"  Tightest cases (excess > 0.9):")
        print(f"  {'excess':>8} {'r':>8} {'s_prime':>8} {'M':>4} {'d':>4} {'p*':>8} {'r+s':>8} {'r*s':>8}")
        for excess, r, s_prime, M, d, pstar in data[:20]:
            print(f"  {excess:8.5f} {r:8.5f} {s_prime:8.5f} {M:4d} {d:4d} {pstar:8.5f} {r+s_prime:8.5f} {r*s_prime:8.5f}")
        
        print()
        print(f"  Pattern in tight cases:")
        print(f"    r·s' product: min={min(r*s for _,r,s,_,_,_ in data[:20]):.4f}, "
              f"max={max(r*s for _,r,s,_,_,_ in data[:20]):.4f}")
        print(f"    r values:     min={min(r for _,r,_,_,_,_ in data[:20]):.4f}, "
              f"max={max(r for _,r,_,_,_,_ in data[:20]):.4f}")
        print(f"    s' values:    min={min(s for _,_,s,_,_,_ in data[:20]):.4f}, "
              f"max={max(s for _,_,s,_,_,_ in data[:20]):.4f}")

    print()
    print("=" * 80)

if __name__ == "__main__":
    main()
