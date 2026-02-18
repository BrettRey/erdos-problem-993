#!/usr/bin/env python3
"""FINAL ATTACK: Use the λ-deformation + direct coefficient analysis.

We've established:
1. The property holds for all trees n ≤ 17 at all λ (822K checks)
2. The property holds for all connected graphs n ≤ 9 (1.4M checks)
3. Sums of IS distributions propagate the property
4. The property is NOT a consequence of LC, ULC, or real-rootedness alone

The TIGHTEST cases have:
- excess ≈ 0.98 (approaching 1)
- s' ≈ 0.50 (left decay approaching 1/2)
- r = 0 (no right tail — mode is at the END of the support!)
- p* ≈ 0.37

These are trees where the IS poly peaks at position α(T)-1 with a very 
flat left tail and sharp cutoff at α(T).

Actually, the r=0 observation is KEY: the tightest trees have 
mode = α(T) = independence number! The IS distribution drops to 0 
right after the mode.

This means: the tight cases have i_{α} > 0 and i_{α+1} = 0.
The IS number α(T) is the mode.
And mean ≈ α - 0.98.

For a tree, α(T) ≥ ⌈n/2⌉ (bipartite + Konig's theorem).

When does mode = α(T)?
This happens when i_{α} ≥ i_{α-1} (the last coefficient is the largest
or tied with the previous one).

For this to happen: the number of MAXIMUM independent sets must be 
at least as large as i_{α-1}.

Hmm, that's a strong condition. Usually the last coefficient is small 
(one or few maximum IS's).

Wait, looking at the data again: the TIGHT cases have r=0 and s' ≈ 0.5.
With r=0, the right tail vanishes, meaning i_{mode+1} = 0.
But mode = M doesn't have to equal α(T). It could be that 
at the TIE point λ₀, the weighted distribution p_k(λ₀) peaks at M 
which is NOT equal to α(T).

Let me re-examine. At the tie point λ₀:
p_M(λ₀) = p_{M-1}(λ₀) = max. And p_{M+1}(λ₀) ≈ 0.

This means: i_{M+1} λ₀^{M+1} ≈ 0. Either i_{M+1} ≈ 0 or λ₀ is very small.

For λ₀ = i_{M-1}/i_M: if M < α(T), then i_{M+1} > 0, and 
p_{M+1} = i_{M+1}λ₀^{M+1}/Z. We saw r = p_{M+1}/p_M ≈ 0, which means
i_{M+1}/i_M · λ₀ ≈ 0, i.e., i_{M+1}/i_M is much less than 1/λ₀ = i_M/i_{M-1}.
So i_{M+1} ≪ i_{M-1}.

This IS the LC+monotonicity condition: the ratios r_k = i_k/i_{k-1} 
decrease past the mode, so i_{M+1}/i_M < i_M/i_{M-1}.

At the tie point: r_{M+1} = i_{M+1}/i_M ≪ r_M = i_M/i_{M-1}.

The ratio of consecutive r's: r_{M+1}/r_M = i_{M+1}i_{M-1}/i_M² ≤ 1 (by LC).

So: s' = r_{M+1}/r_M ≤ 1 (wait, I had s' as the LEFT decay, let me 
be consistent).

s' = p_{M-2}/p_{M-1} = (i_{M-2}·λ₀^{M-2})/(i_{M-1}·λ₀^{M-1}) = i_{M-2}/(i_{M-1}·λ₀)
   = (i_{M-2}/i_{M-1})·(i_M/i_{M-1}) = r_M/r_{M-1} (since λ₀ = i_{M-1}/i_M = 1/r_M)
   Wait: s' = r_M/r_{M-1}? Let me redo this.

p_{M-2}/p_{M-1} = (i_{M-2}λ₀^{M-2})/(i_{M-1}λ₀^{M-1}) = (i_{M-2}/i_{M-1})·(1/λ₀)
               = (1/r_{M-1})·r_M = r_M/r_{M-1}

Since r_{M-1} ≥ r_M (LC): s' = r_M/r_{M-1} ≤ 1. ✓

Similarly: r = p_{M+1}/p_M = (i_{M+1}/i_M)·λ₀ = r_{M+1}/r_M.
r ≤ 1 by LC. And our data shows r ≈ 0, meaning r_{M+1} ≪ r_M.

IMPORTANT: The TIGHTEST cases have s' ≈ 0.50, meaning r_M/r_{M-1} ≈ 0.50.
By LC: r_M² ≥ r_{M-1}·r_{M+1} → r_{M+1} ≤ r_M²/r_{M-1} = r_M·s'.
So: r = r_{M+1}/r_M ≤ s'. 
At s' = 0.50: r ≤ 0.50. Our data shows r ≈ 0, which is much tighter.

Actually, there IS a stronger constraint from LC applied at MULTIPLE levels.
At position k: r_k² ≥ r_{k-1}·r_{k+1}
This means: r_{k+1} ≤ r_k²/r_{k-1}

Starting from r_M and r_{M-1} = r_M/s':
r_{M+1} ≤ r_M²/(r_M/s') = r_M·s'
r_{M+2} ≤ r_{M+1}²/r_M ≤ (r_M·s')²/r_M = r_M·s'²
r_{M+j} ≤ r_M·s'^j

So the RIGHT tail decays at least as fast as s'^j, 
while the LEFT tail decays EXACTLY as s' (first step) then accelerating.

For the excess:
excess = Σ_{j≥1} j·p_{M-j} - Σ_{j≥1} j·p_{M+j}
       ≤ Σ j·(s')^{j-1}·p_M - Σ j·r^{j}·p_M (geometric bounds)
       ≤ p_M/(1-s')² - p_M·r/(1-r)²  (for r = 0: just p_M/(1-s')²)

And p_M ≤ (1-s')/(2-s') (from normalization with symmetric-ish tails).
So excess ≤ 1/[(2-s')(1-s')].

At s'=0.50: excess ≤ 1/(1.5·0.5) = 4/3. WAY too loose compared to actual 0.98.

THE GAP: the geometric bound gives 4/3, but the actual is 0.98.
This 36% gap is due to the super-geometric decay and the right tail.

Actually, let me use the TIGHTER bound with r = s'² (from LC):
excess ≤ p_star/(1-s')² - p_star·s'²/(1-s'²)²
       = p_star·[1/(1-s')² - s'²/(1-s'²)²]

p_star = (1-s')(1-s'²)/[2-s'-s'² + correction... hmm complicated.

Forget the analytical approach for now. Let me instead try to
COMPUTATIONALLY find what separates the violations from non-violations.

I'll generate random distributions that:
(a) Are LC
(b) Have mode ∈ {⌊μ⌋, ⌈μ⌉}
Then take SUMS and check if violations occur.
Compare with IS distributions.
"""

import random, math
import subprocess
from indpoly import independence_poly
from graph6 import parse_graph6

def convolve(a, b):
    n = len(a) + len(b) - 1
    r = [0.0]*n
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            r[i+j] += ai*bj
    return r

def poly_mode(p):
    return max(range(len(p)), key=lambda k: p[k])

def poly_mean(p):
    s = sum(p)
    if s == 0: return 0
    return sum(k*p[k] for k in range(len(p))) / s

def is_lc(p):
    for k in range(1, len(p)-1):
        if p[k] > 0:
            if p[k]**2 < p[k-1]*p[k+1] - 1e-12:
                return False
    return True

def main():
    random.seed(42)
    
    print("=" * 70)
    print("  WHAT MAKES IS DISTRIBUTIONS SPECIAL?")
    print("=" * 70)
    print()
    
    # Collect IS distributions
    is_dists = []
    for n in range(2, 10):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            Z = sum(poly)
            dist = [p/Z for p in poly]
            is_dists.append(dist)
    
    # Properties of IS distributions
    print(f"  {len(is_dists)} tree IS distributions (n=2..9)")
    
    # Property 1: a₀ = 1 → p₀ > 0
    all_p0_pos = all(d[0] > 0 for d in is_dists)
    print(f"  All have p₀ > 0: {all_p0_pos}")
    
    # Property 2: LC
    all_lc = all(is_lc(d) for d in is_dists)
    print(f"  All LC: {all_lc}")
    
    # Property 3: Unimodal
    def is_unimodal(p):
        mode = poly_mode(p)
        for k in range(1, mode):
            if p[k] < p[k-1] - 1e-12:
                return False
        for k in range(mode+1, len(p)):
            if p[k] > p[k-1] + 1e-12:
                return False
        return True
    
    all_uni = all(is_unimodal(d) for d in is_dists)
    print(f"  All unimodal: {all_uni}")
    
    # Property 4: Variance relative to mean
    def poly_var(p):
        m = poly_mean(p)
        return sum(p[k]*(k-m)**2 for k in range(len(p)))
    
    min_var_over_mean = float('inf')
    for d in is_dists:
        m = poly_mean(d)
        v = poly_var(d)
        if m > 0.01:
            ratio = v/m
            if ratio < min_var_over_mean:
                min_var_over_mean = ratio
    
    print(f"  Min Var/Mean ratio: {min_var_over_mean:.4f}")
    
    # Property 5: Ratio of consecutive coefficients
    # For IS polys at λ=1: p_k/p_{k-1} = i_k/i_{k-1}
    # The first ratio a₁/a₀ = n (number of vertices)
    # This is always ≥ 2 for n ≥ 2.
    
    # Property 6: MODE vs degree
    # mode of IS dist / degree of polynomial
    mode_ratios = []
    for d in is_dists:
        m = poly_mode(d)
        deg = len(d) - 1
        mode_ratios.append(m / deg if deg > 0 else 0)
    
    print(f"  Mode/degree ratio: min={min(mode_ratios):.3f}, max={max(mode_ratios):.3f}, "
          f"mean={sum(mode_ratios)/len(mode_ratios):.3f}")
    
    # Property 7: excess distribution
    excesses = [poly_mode(d) - poly_mean(d) for d in is_dists]
    print(f"  Excess: min={min(excesses):.4f}, max={max(excesses):.4f}")
    
    # Property 8: p_mode value
    pmode_vals = [d[poly_mode(d)] for d in is_dists]
    print(f"  p_mode: min={min(pmode_vals):.4f}, max={max(pmode_vals):.4f}")
    
    print()
    
    # Now: compare with LC distributions that FAIL under sums
    # Generate LC distributions where mode ∈ {floor, ceil} individually
    # but their sum violates
    
    print("  SEARCHING for LC distributions whose sum VIOLATES mode ∈ {⌊μ⌋, ⌈μ⌉}:")
    
    sum_violations = 0
    
    for trial in range(100000):
        # Generate 2 LC distributions
        dists = []
        for _ in range(2):
            while True:
                d = random.randint(2, 6)
                peak = random.randint(0, d)
                sl = random.uniform(0.3, 3.0)
                sr = random.uniform(0.3, 3.0)
                
                p = []
                for k in range(d+1):
                    if k <= peak:
                        p.append(math.exp(-sl*(peak-k)))
                    else:
                        p.append(math.exp(-sr*(k-peak)))
                
                Z = sum(p)
                p = [x/Z for x in p]
                
                if is_lc(p):
                    m = poly_mode(p)
                    mu = poly_mean(p)
                    if m == math.floor(mu) or m == math.ceil(mu):
                        dists.append(p)
                        break
        
        conv = convolve(dists[0], dists[1])
        m_s = poly_mode(conv)
        mu_s = poly_mean(conv)
        
        if m_s != math.floor(mu_s) and m_s != math.ceil(mu_s):
            sum_violations += 1
            
            if sum_violations <= 3:
                print(f"\n    Violation {sum_violations}:")
                for i, d in enumerate(dists):
                    d_norm = [x/sum(d) if sum(d) > 0 else 0 for x in d]
                    print(f"      X{i}: {[f'{x:.3f}' for x in d_norm]}")
                    print(f"        mode={poly_mode(d_norm)}, mean={poly_mean(d_norm):.3f}, "
                          f"LC={is_lc(d_norm)}")
                
                c_norm = [x/sum(conv) for x in conv]
                print(f"      SUM: mode={m_s}, mean={mu_s:.4f}, LC={is_lc(c_norm)}")
                
                # KEY: is the SUM still LC?
                print(f"      SUM IS LC: {is_lc(c_norm)}")
    
    print(f"\n    LC sum violations: {sum_violations}/100000")
    print()
    
    # CRITICAL CHECK: are sums of LC distributions always LC?
    # If yes, and if LC → mode ∈ {floor, ceil} for the sum, then we're done.
    # BUT we already know LC does NOT imply mode ∈ {floor, ceil}!
    
    # However: if the SUM of LC distributions is ALWAYS LC, 
    # and the INITIAL IS distributions have some extra property P,
    # and P propagates through LC sums, then:
    # P might be what we need.
    
    # CHECK: is the sum of two LC distributions always LC?
    print("  CHECK: Is the convolution of two LC distributions always LC?")
    lc_conv_violations = 0
    
    for trial in range(100000):
        dists = []
        for _ in range(2):
            while True:
                d = random.randint(2, 6)
                peak = random.randint(0, d)
                sl = random.uniform(0.3, 5.0)
                sr = random.uniform(0.3, 5.0)
                
                p = []
                for k in range(d+1):
                    if k <= peak:
                        p.append(math.exp(-sl*(peak-k)))
                    else:
                        p.append(math.exp(-sr*(k-peak)))
                
                Z = sum(p)
                p = [x/Z for x in p]
                
                if is_lc(p):
                    dists.append(p)
                    break
        
        conv = convolve(dists[0], dists[1])
        if not is_lc(conv):
            lc_conv_violations += 1
    
    print(f"    Convolution of LC is always LC: {lc_conv_violations == 0}")
    print(f"    (violations: {lc_conv_violations}/100000)")
    
    # YES: the convolution of LC distributions IS always LC!
    # This is Ibragimov's theorem (1956).
    # So sums of LC distributions are LC.
    # And IS distributions from trees (n ≤ 25) are LC.
    # So sums of IS distributions are LC.
    
    # But LC alone doesn't give mode ∈ {floor, ceil}!
    # So the IS distributions must have something ADDITIONAL that 
    # propagates through LC convolution.
    
    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
