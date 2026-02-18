#!/usr/bin/env python3
"""STRENGTHENED INDUCTION: excess < 1 - P(root)/2.

Idea: the excess of a tree with n ≥ 2 vertices satisfies
  excess(T) < 1 - δ(T)
where δ(T) > 0 is some quantifiable slack.

If we can show δ(T) is bounded below, the induction closes.

What is the TIGHTEST relationship? Let me compute:
  slack = 1 - excess(T) for all trees
and see what correlates with it.

Candidates for δ(T):
  - 1/n (gets small)
  - P(root) (occupation prob of some vertex)
  - min_v P(v) (smallest occupation prob)
  - 1/α(T) (reciprocal of IS number)
  - Some function of the DP quantities

Let me also try the DUAL approach: instead of strengthening the 
inductive hypothesis, use a TWO-STEP induction.

TWO-STEP INDUCTION:
  Let e_k = max excess over all trees on k vertices.
  At a gap=+1 merge from child of size k to parent of size k+1:
    excess(parent) = 1 + excess(dp₀) - correction ≤ 1 + e_k - C_min
  where C_min is the minimum correction at gap=+1 merges.
  
  If C_min > e_k for all k, then excess(parent) < 1.
  
  From our data: minimum correction at gap=+1 is ~0.167 (n=2).
  And max e_k grows toward ~0.6.
  
  So 1 + 0.6 - 0.167 = 1.433 > 1. Too loose.
  
  But the correction is NOT independent of excess(dp₀)!
  When excess(dp₀) is large, the correction is ALSO large.
  
  The KEY is the SLACK = correction - excess(dp₀), 
  which we showed is always ≥ 0.38.

So: excess(parent) = 1 + excess(dp₀) - correction = 1 - slack ≤ 1 - 0.38 = 0.62.

Wait, that's ALREADY a bound! If slack ≥ 0.38, then excess ≤ 0.62 at gap=+1.
And at gap ≤ 0, excess ≤ excess(dp₀) ≤ max over subtree products.

But does excess(dp₀) stay bounded? For gap=0 merges, 
excess(total) < excess(dp₀). This means:
  excess(total) ≤ max(excess(dp₀) at gap≤0, 1-slack at gap=+1)

And excess(dp₀) = excess(∏ total[cᵢ]).

For a single-child gap=0: excess(dp₀) = excess(total[c]) which is 
the excess of a smaller tree. By induction, < 1.

For multi-child: the product concentrates, so excess(∏) ≤ Σ excess(single).
Hmm, that grows with k.

WAIT - but we MEASURED excess(dp₀) ≤ 0.61 for ALL tested trees.

Let me try to prove: excess(dp₀) ≤ 0.62 for all trees, 
using the fact that slack ≥ 0.38 at gap=+1.

This would give a SELF-CONSISTENT bound:
  At gap=+1: excess = 1 - slack ≤ 1 - 0.38 = 0.62
  At gap=0: excess ≤ excess(dp₀) < 0.62 (inherits from subtrees)
  Products of distributions with excess < 0.62: excess ≤ 0.62?

For products: mode(∏) ≤ Σ mode(Xᵢ) and mean(∏) = Σ mean(Xᵢ).
excess(∏) ≤ Σ excess(Xᵢ) < k · 0.62 → grows with k!

But empirically, excess(∏) ≤ 0.61 always. Why?

Because the MODE of a convolution doesn't reach Σ mode(Xᵢ).
The mode of a sum of independent rv's is CLOSE to the sum of means.
By the CLT, for k large, mode ≈ mean, so excess → 0.

For k=2: excess(X₁+X₂) ≤ excess(X₁) + excess(X₂)? 
Not tight. The actual excess of a convolution is less.

There's a theorem: for the sum of independent unimodal rv's,
the mode is within 1 of the mean? Not exactly.

Actually, there IS: for the sum of TWO independent unimodal distributions 
each with excess ≤ ε, the convolution has excess ≤ ε + 1? No too loose.

Ruzsa's theorem: the mode of the sum of independent rv's is 
at most max(mode(X₁), n) for some n... I don't remember.

Let me just verify the PRODUCT EXCESS numerically with exhaustive search.
"""

import subprocess
from indpoly import independence_poly
from graph6 import parse_graph6

def polymul(a, b):
    n = len(a) + len(b) - 1
    r = [0]*n
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

def poly_excess(p):
    return poly_mode(p) - poly_mean(p)

def main():
    print("=" * 70)
    print("  PRODUCT EXCESS: maximum excess of products of tree IS polys")  
    print("=" * 70)
    print()
    
    # Collect ALL tree IS polys for n ≤ 12
    tree_polys = {}
    for n in range(2, 13):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        tree_polys[n] = []
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            tree_polys[n].append(poly)
    
    # Now compute max excess of products of 2 tree IS polys
    print("  Products of 2 tree IS polys:")
    max_product_excess = -100
    
    for n1 in range(2, 10):
        for n2 in range(n1, 10):
            for p1 in tree_polys[n1]:
                for p2 in tree_polys[n2]:
                    prod = polymul(p1, p2)
                    ex = poly_excess(prod)
                    if ex > max_product_excess:
                        max_product_excess = ex
                        best_product = (n1, n2, p1[:], p2[:], ex)
    
    print(f"    Max excess: {max_product_excess:.6f}")
    if max_product_excess > 0:
        n1, n2, p1, p2, ex = best_product
        print(f"    Achieved by: n1={n1} poly={p1}, n2={n2} poly={p2}")
    print()
    
    # Products of 3
    print("  Products of 3 tree IS polys (sampled):")
    max3 = -100
    import itertools
    
    small_polys = []
    for n in range(2, 8):
        small_polys.extend(tree_polys[n])
    
    import random
    random.seed(42)
    
    for _ in range(50000):
        k = 3
        choices = random.choices(small_polys, k=k)
        prod = choices[0]
        for i in range(1, k):
            prod = polymul(prod, choices[i])
        ex = poly_excess(prod)
        if ex > max3:
            max3 = ex
    
    print(f"    Max excess (50K samples): {max3:.6f}")
    print()
    
    # Now the KEY: compute the minimum SLACK at gap=+1 merges
    # Slack = correction - excess(dp₀) = w*(1+μ_B-μ_A) - (M-μ_A)
    print("  MINIMUM SLACK at gap=+1 merges:")
    
    min_slack = float('inf')
    best_slack_info = None
    
    for n in range(3, 16):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        n_min_slack = float('inf')
        
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            
            mode = poly_mode(poly)
            mean = poly_mean(poly)
            excess = mode - mean
            
            if excess > 0.5:
                # Find the DP merge info at the root
                # Use the formula: excess = 1 + excess(dp₀) - correction
                # So slack = correction - excess(dp₀) = 1 - excess
                slack = 1 - excess
                
                if slack < n_min_slack:
                    n_min_slack = slack
                if slack < min_slack:
                    min_slack = slack
        
        if n_min_slack < float('inf'):
            print(f"    n={n:2d}: min slack = {n_min_slack:.6f}")
    
    print()
    print(f"  OVERALL minimum slack: {min_slack:.6f}")
    print(f"  This means: at gap=+1, excess ≤ 1 - {min_slack:.6f} = {1-min_slack:.6f}")
    print()
    
    # THE SELF-CONSISTENT ARGUMENT:
    # Let e* = sup excess over all trees.
    # At gap=+1: excess = 1 - slack.
    # At gap≤0: excess ≤ excess(dp₀) ≤ e* (for single child)
    #                                    ≤ e* (for multi-child, by product bound)
    # So: e* ≤ max(1 - min_slack, e*_product)
    # If 1 - min_slack < 1 and e*_product < 1, then e* < 1.
    
    # BUT: e*_product ≤ e* for single child, and products REDUCE excess.
    # So: e*_product ≤ e* (gap=0 propagates directly).
    
    # At gap=+1: excess = 1 - slack where slack ≥ 0.38 empirically.
    # At gap=0 with single child: excess ≤ excess(total[c]) = previous e*.
    # At gap=0 with multi child: excess ≤ excess(∏ total[cᵢ]) ≤ e* (or less).
    
    # So: e* = max over all merges = max(max at gap=+1, max at gap=0).
    # max at gap=+1 = max(1 - slack) ≤ 1 - min_slack.
    # max at gap=0 ≤ e* (circular!).
    
    # But at gap=0, excess STRICTLY DECREASES (excess(total) < excess(dp₀)).
    # So the maximum excess must be ACHIEVED at a gap=+1 merge!
    
    # Therefore: e* = max excess at gap=+1 = max(1 - slack).
    # And min slack ≥ 0.38 gives e* ≤ 0.62.
    
    # BUT: is this argument rigorous? The issue is:
    # At gap=0, excess(total) < excess(dp₀), but dp₀ = ∏ total[cᵢ],
    # and the excess of the product can be LARGER than any individual excess
    # if mode alignment aligns favorably.
    
    # Hmm, actually NO: for gap=0, mode(total) ≤ mode(dp₀).
    # So excess(total) = mode(total) - mean(total) ≤ mode(dp₀) - mean(total).
    # And mean(total) > mean(dp₀) (since dp₁ adds positive mass).
    # So excess(total) < mode(dp₀) - mean(dp₀) = excess(dp₀).
    
    # And excess(dp₀) = excess(∏ total[cᵢ]).
    # For a SINGLE child: this is just the child's excess. ≤ e* by induction.
    # For MULTIPLE children: the product excess ≤ Σ excess(cᵢ)... grows with k.
    
    # AH BUT: the product mode ≤ Σ mode(cᵢ), and the product mean = Σ mean(cᵢ).
    # So product excess ≤ Σ excess(cᵢ) ≤ k * e*.
    # This breaks the bound for k ≥ 2.
    
    # HOWEVER: excess(total) < excess(dp₀) at gap≤0.
    # And then dp₁ adds mass that shifts the mean UP.
    # At the PREVIOUS level, this excess(dp₀) is used.
    
    # The chain: excess at level L ≤ excess(∏ at level L-1).
    # For a single-child chain (path), this works: each step inherits or reduces.
    # For branching, the product can increase excess.
    
    # KEY FACT: excess(dp₀ at branching) = excess(∏ total[cᵢ]).
    # Each total[cᵢ] has excess < 1 by induction.
    # But ∏ can have excess > individual excesses.
    
    # UNLESS: excess(∏ total[cᵢ]) ≤ max excess(total[cᵢ]) < 1.
    # Is THIS true?
    
    print()
    print("  CRITICAL TEST: excess(A*B) ≤ max(excess(A), excess(B))?")
    
    violations = 0
    total_tests = 0
    max_ratio = 0
    
    for n1 in range(2, 10):
        for n2 in range(n1, 10):
            for p1 in tree_polys[n1]:
                for p2 in tree_polys[n2]:
                    prod = polymul(p1, p2)
                    ex_prod = poly_excess(prod)
                    ex1 = poly_excess(p1)
                    ex2 = poly_excess(p2)
                    max_ex = max(ex1, ex2)
                    total_tests += 1
                    
                    if ex_prod > max_ex + 1e-10:
                        violations += 1
                        ratio = ex_prod / max_ex if max_ex > 0 else float('inf')
                        if ratio > max_ratio:
                            max_ratio = ratio
    
    print(f"    Tests: {total_tests}")
    print(f"    Violations (excess(A*B) > max(excess(A), excess(B))): {violations}")
    if violations > 0:
        print(f"    Max ratio: {max_ratio:.4f}")
    else:
        print("    NO VIOLATIONS! Product excess ≤ max of individual excesses!")
    
    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
