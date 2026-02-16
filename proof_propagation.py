#!/usr/bin/env python3
"""HARD-CORE COUPLING ARGUMENT.

For the hard-core model on a TREE at fugacity λ:
  S = Σ Xᵢ where Xᵢ ∈ {0,1} indicates vertex i in the IS.
  The Xᵢ are NOT independent (hard-core constraint).
  But on a tree, the distribution factorizes nicely via BP.

KEY IDEA: NEGATIVE CORRELATION.

For the hard-core model on ANY graph, the occupation variables 
satisfy NEGATIVE ASSOCIATION (Pemantle 2000 for lattice gas models):
  Cov(Xᵢ, Xⱼ) ≤ 0 for all i ≠ j.

This means: Var(S) = Σᵢ Var(Xᵢ) + Σᵢ≠ⱼ Cov(Xᵢ,Xⱼ) ≤ Σᵢ Var(Xᵢ).

Since each Xᵢ ~ Bernoulli(pᵢ): Var(Xᵢ) = pᵢ(1-pᵢ) ≤ 1/4.
So Var(S) ≤ n/4.

Now: S has mean μ = Σ pᵢ and variance σ² = Var(S) ≤ n/4.

If S were EXACTLY the sum of independent Bernoullis with the SAME 
marginal probabilities pᵢ, then Darroch's theorem would give 
|mode - mean| < 1.

But S is not independent — the hard-core constraint introduces 
negative correlations. Can we still get the mode-mean bound?

IDEA: Use STOCHASTIC DOMINATION.

Let Y = Σᵢ Yᵢ where Yᵢ ~ Bernoulli(pᵢ) independently.
By negative association: S is stochastically "more concentrated" than Y.
Specifically: P(S ≥ t) ≤ P(Y ≥ t) for all t ≥ μ (Chernoff-type upper tail).
And: P(S ≤ t) ≤ P(Y ≤ t) for all t ≤ μ (lower tail).

These are one-sided bounds. For the MODE, we need more.

Actually, negative association gives us: the generating function
G_S(z) = E[z^S] is coefficient-wise dominated by G_Y(z) = ∏(1-pᵢ+pᵢz) 
in some sense... but not exactly.

Hmm, this is delicate. Let me try a different angle.

DIRECT APPROACH VIA THE BP RECURSION ON TREES:

For a tree T, the IS-size distribution can be computed EXACTLY via the 
belief propagation. The key recursion is:

For tree T rooted at v with children c₁,...,cₖ:
  Z₀(v) = ∏ᵢ Z(cᵢ)     (v not in IS)
  Z₁(v) = λ · ∏ᵢ Z₀(cᵢ) (v in IS, children must be out)
  Z(v) = Z₀(v) + Z₁(v)

The IS-size generating function is:
  G_T(z) = Σ_S z^|S| · λ^|S| / I(T;λ)
         = I(T; zλ) / I(T; λ)

Hmm wait, G_T(z) = Σₖ pₖ z^k where pₖ = iₖ λ^k / I(T;λ).
So G_T(z) = I(T; zλ) / I(T; λ)? 
Check: G_T(z) = Σₖ (iₖ λ^k / I(T;λ)) z^k = I(T; zλ) / I(T; λ). ✓

The mean is G'(1) = λ I'(T;λ) / I(T;λ) = Σᵢ pᵢ. ✓

Now: G_T(z) = I(T; zλ) / I(T; λ).

For a tree: I(T; x) = I(T-v; x) + x·I(T-N[v]; x).
So: G_T(z) = [I(T-v; zλ) + zλ·I(T-N[v]; zλ)] / I(T; λ)
           = [I(T; λ)/I(T-v; zλ)]... hmm, this doesn't simplify easily.

Let me try using the COMPLEMENTARY APPROACH.

At the tie point λ₀ where mode jumps from M-1 to M:
  G_T(z) is the GF of a distribution with mode M (or mode M-1 depending on z).
  
  We need: mean = G'(1) > M - 1.

Actually, let me try a result on NEGATIVE LATTICE DISTRIBUTIONS.

For a negatively associated distribution on {0,...,n}:
  If the distribution is unimodal (mode at M), then |mode - mean| < 1?

By analogy with Darroch's theorem for independent Bernoullis...
Darroch's proof uses the fact that the GF has all real roots.
For negatively associated variables, the GF does NOT necessarily have all real roots.

But maybe the result still holds?

Let me test: compute the GF for hard-core on small trees and check 
mode vs mean directly for the IS-size distribution. We already 
verified this holds for all trees up to n=17. The question is: WHY?

OK let me try yet another approach: INDUCTION ON THE TREE STRUCTURE
using a CONVEX COMBINATION argument.

THEOREM (attempt): For any tree T, at any fugacity λ > 0:
  mode(T; λ) ∈ {⌊μ(T;λ)⌋, ⌈μ(T;λ)⌉}

PROOF (attempt via strong induction on n = |T|):

Base case: n=1 (single vertex). I = 1+λx. Distribution: p₀ = 1/(1+λ), p₁ = λ/(1+λ).
  mean = λ/(1+λ) ∈ (0,1).
  mode = 0 if λ<1, 1 if λ>1, both if λ=1.
  Always: mode ∈ {⌊mean⌋, ⌈mean⌉}. ✓

Inductive step: Assume the result for all trees with < n vertices.
Let T be a tree on n vertices, rooted at v with children c₁,...,cₖ.

I(T) = A(x) + x·B(x) where A = ∏ I(Tᵢ), B = ∏ I(Tᵢ - cᵢ).

At fugacity λ:
  p_k(T, λ) = (aₖ + bₖ₋₁)λᵏ / Z   (where bₖ₋₁ is from B shifted)

  Actually: I(T;x) = A(x) + x·B(x).
  i_k(T) = a_k + b_{k-1} where a_k = coeff of A, b_{k-1} = coeff of B at k-1.
  
  The distribution of IS-size is a MIXTURE:
  p_k = (a_k λ^k + b_{k-1} λ^k) / Z = (a_k λ^k + b_{k-1} λ^k) / (A(λ) + λB(λ))

  Alternatively:
  p_k = [(1-P(v)) · p_k^A + P(v) · p_{k-1}^B] 

  where P(v) = λB(λ)/(A(λ) + λB(λ)) and:
  p_k^A = a_k λ^k / A(λ) (distribution of IS-size NOT containing v)
  p_k^B = b_k λ^k / B(λ) (distribution of IS-size of T-N[v])

  So: p_k(T) = (1-P(v)) · p_k^A + P(v) · p_{k-1}^B

  The distribution of T is a mixture of:
    - Distribution A (IS-size of T not containing v) with weight 1-P(v)
    - Distribution B+1 (IS-size of T-N[v] plus 1) with weight P(v)

  Distribution A: S₀ = IS-size conditioned on v ∉ S.
    S₀ has the same distribution as the IS-size of T-v.
    Wait: not exactly. If v ∉ S, then S is an IS of T-v (since removing v 
    doesn't affect which edges exist among other vertices). But the children 
    of v are now unconstrained. So S₀ ~ IS-size of T-v at fugacity λ.
    
    Hmm: a_k = coefficient of x^k in ∏ I(Tᵢ;x). 
    But T-v is UNrooted. T-v = ∪ Tᵢ (disjoint union of subtrees).
    So IS-size of T-v at λ = Σᵢ IS-size of Tᵢ at λ (independent sum!).
    
    So: S₀ = Σᵢ Sᵢ where Sᵢ = IS-size of Tᵢ at λ. Independent!

  Distribution B+1: S₁ = 1 + IS-size of T-N[v] at λ.
    T-N[v] = T minus v and its neighbors (parents/children).
    For the root: T-N[v] = ∪ᵢ (Tᵢ - cᵢ) (disjoint union).
    So IS-size of T-N[v] = Σᵢ IS-size of (Tᵢ-cᵢ) at λ. Also independent!
    S₁ = 1 + Σᵢ S'ᵢ where S'ᵢ = IS-size of (Tᵢ-cᵢ) at λ.

  Now: S = IS-size of T at λ.
  S has the distribution: (1-P(v))·S₀ + P(v)·S₁ (mixture).
  
  mean(S) = (1-P(v))·mean(S₀) + P(v)·mean(S₁)
           = (1-P(v))·Σ mean(Sᵢ) + P(v)·(1 + Σ mean(S'ᵢ))

  By induction: for each Tᵢ (tree with fewer vertices):
    mode(Sᵢ) ∈ {⌊mean(Sᵢ)⌋, ⌈mean(Sᵢ)⌉}

  By Darroch's theorem (since Σᵢ Sᵢ is a sum of independent integer rv's,
  each from a smaller tree), if each Sᵢ satisfies mode(Sᵢ) ∈ {⌊μᵢ⌋, ⌈μᵢ⌉}...
  
  Hmm wait, that's not enough. Darroch's theorem requires the Sᵢ to be 
  Bernoulli. They're not — they're IS-size distributions of smaller trees.
  
  BUT: Σ Sᵢ IS a sum of INDEPENDENT random variables.
  For independent rv's, there are results like:
  - If each Sᵢ has mode ∈ {⌊μᵢ⌋, ⌈μᵢ⌉}, does Σ Sᵢ also?
  - NOT in general! The sum of two distributions each with mode ∈ {⌊μ⌋, ⌈μ⌉}
    does NOT necessarily have mode ∈ {⌊μ_total⌋, ⌈μ_total⌉}.
  
  Example: X ~ [0.49, 0.51] (mode=1, mean=0.51, ⌊μ⌋=0, ⌈μ⌉=1 ✓)
           Y ~ [0.49, 0.51] similarly.
           X+Y ~ [0.49*0.49, 2*0.49*0.51, 0.51*0.51] = [0.2401, 0.4998, 0.2601]
           mode = 1, mean = 1.02. ⌊μ⌋=1, ⌈μ⌉=2. mode=1. ✓
  
  Actually this still works. Let me try a worse case.
  X ~ [0.01, 0.99] (mode=1, mean=0.99)
  Y ~ [0.99, 0.01] (mode=0, mean=0.01)
  X+Y ~ [0.0099, 0.99*0.01+0.01*0.99, 0.0099] = [0.0099, 0.0198, 0.0099]
  Wait that doesn't sum to 1. 
  X+Y ~ [0.01*0.99, 0.99*0.99+0.01*0.01, 0.99*0.01] = [0.0099, 0.9802, 0.0099]
  mode=1, mean=1.0. ✓
  
  Hmm, it's hard to find counterexamples. Let me try computationally.
"""

def poly_stats(p):
    Z = sum(p)
    if Z == 0: return 0, 0, 0
    mode = max(range(len(p)), key=lambda k: p[k])
    mean = sum(k*p[k] for k in range(len(p))) / Z
    return mode, mean, mode - mean

def convolve(a, b):
    n = len(a) + len(b) - 1
    r = [0.0]*n
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            r[i+j] += ai*bj
    return r

def main():
    import random, math
    random.seed(42)
    
    print("=" * 70)
    print("  SUM OF INDEPENDENT RVs: Does mode ∈ {⌊μ⌋, ⌈μ⌉} PROPAGATE?")
    print("=" * 70)
    print()
    
    # Test: if X₁,...,Xₖ each have mode ∈ {⌊μᵢ⌋, ⌈μᵢ⌉},
    # does the SUM S = Σ Xᵢ also have mode ∈ {⌊μ_S⌋, ⌈μ_S⌉}?
    
    violations = 0
    total = 0
    
    for trial in range(100000):
        k = random.randint(2, 5)
        dists = []
        
        for _ in range(k):
            # Random distribution with mode ∈ {floor(mean), ceil(mean)}
            d = random.randint(1, 4)
            while True:
                p = [random.random() for _ in range(d+1)]
                Z = sum(p)
                p = [x/Z for x in p]
                mode, mean, _ = poly_stats(p)
                if mode == math.floor(mean) or mode == math.ceil(mean):
                    break
            dists.append(p)
        
        # Compute convolution
        result = dists[0]
        for i in range(1, k):
            result = convolve(result, dists[i])
        
        mode_s, mean_s, _ = poly_stats(result)
        
        total += 1
        if mode_s != math.floor(mean_s) and mode_s != math.ceil(mean_s):
            violations += 1
            if violations <= 5:
                print(f"  VIOLATION #{violations}: k={k}")
                for i, d in enumerate(dists):
                    m, mu, _ = poly_stats(d)
                    print(f"    X_{i}: {[f'{x:.3f}' for x in d]} mode={m} mean={mu:.3f}")
                print(f"    SUM: mode={mode_s} mean={mean_s:.4f} floor={math.floor(mean_s)} ceil={math.ceil(mean_s)}")
    
    print(f"\n  Total tests: {total}")
    print(f"  Violations: {violations}")
    
    if violations > 0:
        print("  ✗ mode ∈ {⌊μ⌋, ⌈μ⌉} does NOT propagate through independent sums!")
        print("  So we can't use induction through tree decomposition this way.")
    else:
        print("  ✓ Zero violations — propagation might hold!")
    
    print()
    
    # Test specifically with IS distributions from small trees
    print("  TEST WITH IS DISTRIBUTIONS FROM TREES:")
    import subprocess
    from indpoly import independence_poly
    from graph6 import parse_graph6
    
    tree_dists = {}
    for n in range(2, 10):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        tree_dists[n] = []
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            tree_dists[n].append(poly)
    
    violations2 = 0
    total2 = 0
    
    for _ in range(50000):
        k = random.randint(2, 4)
        choices = []
        for _ in range(k):
            n_val = random.randint(2, 7)
            p = random.choice(tree_dists[n_val])
            choices.append(p)
        
        result = choices[0]
        for i in range(1, k):
            result = convolve(result, choices[i])
        
        mode_s, mean_s, _ = poly_stats(result)
        total2 += 1
        
        if mode_s != math.floor(mean_s) and mode_s != math.ceil(mean_s):
            violations2 += 1
    
    print(f"    Tests: {total2}")
    print(f"    Violations: {violations2}")
    
    if violations2 == 0:
        print("    ✓ Propagation holds for IS distributions!")
    
    print()
    
    # NOW: the MIXTURE. The IS distribution of T is a MIXTURE of 
    # S₀ (sum of IS-sizes of subtrees) and S₁ = 1 + S₀' (sum of IS-sizes of sub-subtrees).
    # Does mode ∈ {⌊μ⌋, ⌈μ⌉} propagate through mixtures?
    
    print("  TEST: Does mode ∈ {⌊μ⌋, ⌈μ⌉} propagate through MIXTURES?")
    
    violations3 = 0
    total3 = 0
    
    for trial in range(100000):
        # Two distributions, both with mode ∈ {floor, ceil}
        while True:
            d = random.randint(2, 8)
            p1 = [random.random() for _ in range(d+1)]
            Z = sum(p1); p1 = [x/Z for x in p1]
            m1, mu1, _ = poly_stats(p1)
            if m1 == math.floor(mu1) or m1 == math.ceil(mu1):
                break
        
        while True:
            d = random.randint(2, 8)
            p2 = [random.random() for _ in range(d+1)]
            Z = sum(p2); p2 = [x/Z for x in p2]
            m2, mu2, _ = poly_stats(p2)
            if m2 == math.floor(mu2) or m2 == math.ceil(mu2):
                break
        
        # Extend to same length
        n_max = max(len(p1), len(p2))
        while len(p1) < n_max: p1.append(0)
        while len(p2) < n_max: p2.append(0)
        
        # Mixture with random weight
        w = random.random()
        mixed = [(1-w)*p1[k] + w*p2[k] for k in range(n_max)]
        
        mode_m, mean_m, _ = poly_stats(mixed)
        total3 += 1
        
        if mode_m != math.floor(mean_m) and mode_m != math.ceil(mean_m):
            violations3 += 1
            if violations3 <= 3:
                print(f"  VIOLATION #{violations3}:")
                print(f"    P1: mode={m1}, mean={mu1:.3f}")
                print(f"    P2: mode={m2}, mean={mu2:.3f}")
                print(f"    w={w:.3f}: mode={mode_m}, mean={mean_m:.4f}")
    
    print(f"\n    Total tests: {total3}")
    print(f"    Violations: {violations3}")
    
    if violations3 > 0:
        print("    ✗ Mixtures break the mode ∈ {⌊μ⌋, ⌈μ⌉} property!")
    else:
        print("    ✓ Mixtures preserve it!")
    
    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
