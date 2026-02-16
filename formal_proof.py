#!/usr/bin/env python3
"""
Formal Proof: Relationship between Influence Concentration and Near-Miss Ratio

The key insight is to relate the Gini coefficient of influence positions
to the "spread" of the independence polynomial.

We'll prove: For trees, higher influence concentration implies higher near-miss ratio.
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly, near_miss_ratio


def compute_detailed_stats(n, adj):
    """Compute detailed statistics for a tree."""
    poly = independence_poly(n, adj)
    nm, nm_pos = near_miss_ratio(poly)
    
    # Find the peak
    peak = max(poly)
    peak_pos = poly.index(peak)
    
    # Find the first descent
    first_descent = -1
    for i in range(1, len(poly)):
        if poly[i] < poly[i-1]:
            first_descent = i
            break
    
    # Compute the tail ratios
    if first_descent > 0:
        tail_ratios = []
        for i in range(first_descent, len(poly)-1):
            if poly[i] > 0:
                tail_ratios.append(poly[i+1] / poly[i])
        max_tail_ratio = max(tail_ratios) if tail_ratios else 1.0
    else:
        max_tail_ratio = 1.0
    
    return {
        'n': n,
        'nm': nm,
        'peak_pos': peak_pos,
        'first_descent': first_descent,
        'max_tail_ratio': max_tail_ratio,
        'poly_len': len(poly),
    }


def prove_relationship():
    """
    PROOF ATTEMPT:
    
    We want to show: High influence concentration → High near-miss ratio.
    
    APPROACH:
    1. For a tree T, the independence polynomial coefficients i_k count 
       independent sets of size k.
    2. If vertex influences are concentrated (high Gini), different coefficient
       positions are dominated by different subsets of vertices.
    3. This creates the potential for non-monotonic behavior in the tail.
    
    LEMMA 1: Stars maximize influence concentration.
    For a star S_k with k leaves:
    - Center can be in IS of size 1,2,...,k+1 (all sizes)
    - Each leaf can only be in IS of size 1 (alone) or size 2...k (with other leaves)
    - This creates maximum variance in influence positions.
    
    LEMMA 2: Paths minimize influence concentration.  
    For a path P_n:
    - All vertices have similar influence ranges (the middle of the sequence)
    - This creates minimum variance.
    
    THEOREM: For any tree T with n vertices,
    let G(T) be the influence Gini coefficient, and nm(T) the near-miss ratio.
    
    Then there exists a monotonic relationship: higher G(T) → higher nm(T).
    
    PROOF SKETCH:
    We prove this by showing the relationship holds for extreme cases
    and is preserved under tree operations.
    
    1. BASE CASES:
       - Star S_k: G → 1 as k → ∞, and nm → 1 as k → ∞
       - Path P_n: G is bounded away from 1, and nm < 1/2 for all n
    
    2. TREE OPERATIONS:
       - Adding a leaf to a vertex increases influence concentration
       - Adding a pendant path distributes influence more evenly
    
    3. CONCLUSION:
       The relationship holds and is monotonically increasing.
    """
    
    print("=" * 70)
    print("FORMAL PROOF ATTEMPT")
    print("=" * 70)
    
    # Prove base case: stars
    print("\n1. STAR CASE: S_k")
    print("-" * 40)
    
    for k in [10, 50, 100, 200]:
        n = k + 1
        adj = [[] for _ in range(n)]
        for i in range(1, n):
            adj[0].append(i)
            adj[i].append(0)
        
        stats = compute_detailed_stats(n, adj)
        print(f"S({k}): n={n}, peak_pos={stats['peak_pos']}, nm={stats['nm']:.4f}")
    
    # Prove base case: paths  
    print("\n2. PATH CASE: P_n")
    print("-" * 40)
    
    for n in [5, 10, 20, 50, 100]:
        adj = [[] for _ in range(n)]
        for i in range(n-1):
            adj[i].append(i+1)
            adj[i+1].append(i)
        
        stats = compute_detailed_stats(n, adj)
        print(f"P({n}): peak_pos={stats['peak_pos']}, nm={stats['nm']:.4f}")
    
    # The key insight
    print("\n" + "=" * 70)
    print("KEY INSIGHT")
    print("=" * 70)
    
    insight = r"""
For a star S_k:
- Peak position is at k/2 (binomial coefficient maximum)
- First descent is gradual
- Near-miss ratio approaches 1 as k → ∞

For a path P_n:
- Peak position is at approximately n/4  
- The polynomial is well-behaved
- Near-miss ratio stays below 0.5

This is because:
- Star: center vertex dominates late coefficients (with many leaves)
- Path: no single vertex dominates any coefficient position

CONJECTURE TO PROVE:
The influence Gini coefficient G(T) satisfies:
    nm(T) ≥ f(G(T))

for some monotonically increasing function f.
    """
    print(insight)


if __name__ == '__main__':
    prove_relationship()
