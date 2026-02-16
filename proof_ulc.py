#!/usr/bin/env python3
"""Identify the structural property of IS distributions that makes sums well-behaved.

Candidates:
1. Ultra-log-concavity (ULC): a_k²/C(d,k)² ≥ a_{k-1}/C(d,k-1)·a_{k+1}/C(d,k+1)
2. Real-rootedness of GENERATING FUNCTION at a specific fugacity
3. Some form of "coefficient dominance" or "spread" condition
4. The VARIANCE relative to mean is bounded in a specific way

Or maybe it's simpler: IS distributions are always UNI-MODAL 
(for trees, proved by Alavi-Malde-Schwenk 1987), and STRONGLY 
unimodal (the sum of two unimodal distributions is unimodal — 
this is Ibragimov's theorem for LC, since LC ⊂ strongly unimodal).

Actually: the key might be that IS distributions are LC (for trees, 
at least for n ≤ 25), and for LC distributions on integers,
there might be a stronger mode-mean result that goes beyond just LC.

Hmm wait — LC of the COEFFICIENTS of I(T;x) doesn't directly mean 
the WEIGHTED distribution p_k(λ) = i_k λ^k/Z is LC. Let me check:

p_k is LC iff p_k² ≥ p_{k-1}p_{k+1}, i.e., (i_k λ^k)² ≥ (i_{k-1}λ^{k-1})(i_{k+1}λ^{k+1})
= i_k² ≥ i_{k-1}·i_{k+1} (the λ's cancel!).

So: p_k(λ) is LC iff the ORIGINAL coefficients i_k are LC.

For trees with n ≤ 25: i_k IS LC (Alavi-Malde-Schwenk, confirmed computationally).
So p_k(λ) is LC for all λ > 0.

For trees with n ≥ 26: LC may fail (Kadrawi-Levit counterexamples).
When LC fails, the weighted distribution is also NOT LC.

BUT we've verified that mode ∈ {⌊μ⌋, ⌈μ⌉} holds for ALL trees 
up to n=17 at ALL fugacities (822K checks). And for n ≤ 25, the IS poly
is LC, which should help.

WAIT — I just realized something. For the HARD-CORE MODEL on trees,
even though the IS poly may not be real-rooted, the distribution at 
any fugacity λ is the distribution of the IS-SIZE, which is a 
DETERMINANTAL POINT PROCESS on certain graphs. For determinantal 
processes, the size distribution IS the distribution of a sum of 
independent Bernoullis (this is the spectral decomposition).

Let me check: is the hard-core model on trees a determinantal process?

No — the hard-core model is NOT a determinantal process in general.
Determinantal processes have NEGATIVE CORRELATIONS (which hard-core does)
but also require a specific algebraic structure (kernel matrix).

HOWEVER: there IS a result connecting the hard-core model on trees to 
determinantal processes. The partition function of a hard-core model 
on a TREE can be written as a DETERMINANT of a certain matrix 
(tridiagonal for paths, etc.), and this gives the IS polynomial as 
a characteristic polynomial.

Actually, for paths: I(P_n; x) = det(I + xA_n) where A_n is the 
adjacency matrix? No, that gives the matching polynomial, not the IS poly.

Hmm, let me think differently.

MATCHING POLYNOMIAL vs IS POLYNOMIAL:
- Matching poly of G: Σ m_k (-1)^k x^{n-2k} where m_k = #matchings of size k.
  This HAS all real roots (Heilmann-Lieb 1972) for all graphs.
  
- IS poly of G: Σ i_k x^k where i_k = #IS of size k.
  This does NOT have all real roots in general.

But for a graph G, the IS polynomial of G equals the MATCHING polynomial 
of the COMPLEMENT of G? No, that's not right either.

Actually, there's a subtler connection for TREES:
The matching polynomial of a tree IS its characteristic polynomial of 
the adjacency matrix (Schwenk 1973). And the IS polynomial of a tree 
is related to the matching polynomial of some other graph.

For a LINE GRAPH L(T) of a tree T: I(L(T); x) = μ(T; -x)... hmm, 
not quite. The independence number of L(T) = the matching number of T.

Let me just CHECK: for small trees, is I(T;x) = det(I + xD) for some 
nice matrix D?

If so, then I(T; x) = ∏(1 + λᵢ x) where λᵢ are eigenvalues of D,
and the IS-size distribution at fugacity λ would be a sum of 
independent Bernoullis with probabilities λᵢλ/(1+λᵢλ).
Darroch's theorem would then apply immediately!
"""

import subprocess
import numpy as np
from indpoly import independence_poly
from graph6 import parse_graph6

def main():
    print("=" * 70)
    print("  CHECK: Is I(T;x) a determinantal polynomial?")
    print("  i.e., can I(T;x) = det(I + x·D) for some matrix D?")
    print("=" * 70)
    print()
    
    # For I(T;x) = det(I + xD), we need:
    # I(T;x) = ∏(1 + λᵢ x) where λᵢ are eigenvalues of D.
    # This means I(T;x) has ALL REAL NON-NEGATIVE roots at x = -1/λᵢ.
    
    # But we ALREADY showed tree IS polys are NOT real-rooted!
    # So I(T;x) ≠ det(I + xD) for any PSD matrix D.
    
    print("  We already know tree IS polys are NOT real-rooted,")
    print("  so they cannot be determinantal. Confirmed.")
    print()
    
    # BUT: maybe the CONDITIONAL distributions (conditioned on root in/out)
    # have better properties?
    
    # For a tree T rooted at v:
    # I(T) = A(x) + x·B(x)
    # A(x) = I(T-v) = IS poly of T-v (forest)
    # B(x) = I(T-N[v]) = IS poly of T-N[v] (forest)
    
    # For a FOREST F = T₁ ∪ T₂ ∪ ... ∪ Tₖ:
    # I(F) = ∏ I(Tᵢ)
    
    # The IS-size distribution of F is the CONVOLUTION of IS-size 
    # distributions of the individual trees.
    
    # If each I(Tᵢ) is NOT real-rooted, the product I(F) is also 
    # NOT real-rooted. So even the conditional distributions are 
    # not determinantal.
    
    # HOWEVER: the MODE-MEAN property might still propagate because 
    # of some WEAKER structural property.
    
    # Let me check: what is the STRONGEST property that IS distributions 
    # have beyond LC?
    
    # Property: ULTRA-LOG-CONCAVITY (ULC)
    # i_k/C(n,k) is LC. Equivalently: i_k²/C(n,k)² ≥ i_{k-1}/C(n,k-1)·i_{k+1}/C(n,k+1)
    
    print("  CHECK: Are tree IS sequences ultra-log-concave?")
    
    violations_ulc = 0
    total_trees = 0
    
    from math import comb
    
    for n in range(3, 16):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        n_ulc = 0
        n_not_ulc = 0
        
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            d = len(poly) - 1
            
            total_trees += 1
            is_ulc = True
            
            for k in range(1, d):
                if poly[k] > 0 and poly[k-1] >= 0 and poly[k+1] >= 0:
                    # Check: (i_k/C(d,k))² ≥ (i_{k-1}/C(d,k-1))·(i_{k+1}/C(d,k+1))
                    lhs = (poly[k] / comb(d, k)) ** 2
                    rhs = (poly[k-1] / comb(d, k-1)) * (poly[k+1] / comb(d, k+1))
                    if lhs < rhs - 1e-10:
                        is_ulc = False
                        break
            
            if is_ulc:
                n_ulc += 1
            else:
                n_not_ulc += 1
        
        status = "✓ all ULC" if n_not_ulc == 0 else f"✗ {n_not_ulc} not ULC"
        print(f"    n={n:2d}: {len(lines):5d} trees | {status}")
    
    print()
    
    # Check if ULC is sufficient for mode ∈ {⌊μ⌋, ⌈μ⌉}
    print("  CHECK: For ULC distributions, does mode ∈ {⌊μ⌋, ⌈μ⌉} hold?")
    
    import random, math
    random.seed(42)
    
    violations_ulc_mm = 0
    total_ulc = 0
    
    for trial in range(100000):
        d = random.randint(3, 15)
        
        # Generate ULC sequence: i_k/C(d,k) should be LC
        # So b_k = i_k/C(d,k) is LC → b_k = exp(concave function)
        peak = random.randint(1, d-1)
        sl = random.uniform(0.1, 2.0)
        sr = random.uniform(0.1, 2.0)
        
        b = []
        for k in range(d+1):
            if k <= peak:
                b.append(math.exp(-sl*(peak-k)))
            else:
                b.append(math.exp(-sr*(k-peak)))
        
        # i_k = b_k · C(d,k)
        poly = [b[k] * comb(d, k) for k in range(d+1)]
        
        # Normalize to get probabilities
        Z = sum(poly)
        p = [x/Z for x in poly]
        
        mode_p = max(range(len(p)), key=lambda k: p[k])
        mean_p = sum(k*p[k] for k in range(len(p)))
        
        total_ulc += 1
        if mode_p != math.floor(mean_p) and mode_p != math.ceil(mean_p):
            violations_ulc_mm += 1
    
    print(f"    ULC tests: {total_ulc}")
    print(f"    Violations: {violations_ulc_mm}")
    
    if violations_ulc_mm == 0:
        print("    ✓ ULC → mode ∈ {⌊μ⌋, ⌈μ⌉}!")
        print("    COMBINED WITH: trees (n≤25) are ULC (?)")
        print("    THIS WOULD PROVE THE CLAIM FOR n ≤ 25!")
    else:
        print(f"    ✗ {violations_ulc_mm} violations found")
    
    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
