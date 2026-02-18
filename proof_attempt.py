#!/usr/bin/env python3
"""
PROOF ATTEMPT: Gini-Unimodality Relationship

Strategy:
1. Analyze extreme cases (stars, paths) to establish base cases
2. Connect influence distribution to polynomial behavior
3. Try to prove a formal relationship

Let's start by computing exact values for simple cases.
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly, near_miss_ratio


def compute_exact_influence_stats(n, adj):
    """Compute exact influence statistics for a tree."""
    # This is expensive, so let's do simpler metrics first
    
    # For a path: all vertices have similar influence (middle values)
    # For a star: center affects all values, leaves affect only small values
    
    # Let's compute the polynomial directly and analyze
    poly = independence_poly(n, adj)
    nm, nm_pos = near_miss_ratio(poly)
    
    return {
        'n': n,
        'poly': poly,
        'nm': nm,
        'nm_pos': nm_pos,
    }


# ============================================================================
# PART 1: ANALYZE EXTREME CASES
# ============================================================================

def star_tree(k):
    """Star with k leaves, n = k+1 vertices."""
    n = k + 1
    adj = [[] for _ in range(n)]
    for i in range(1, n):
        adj[0].append(i)
        adj[i].append(0)
    return n, adj


def path_tree(n):
    """Path with n vertices."""
    adj = [[] for _ in range(n)]
    for i in range(n-1):
        adj[i].append(i+1)
        adj[i+1].append(i)
    return n, adj


def main():
    print("=" * 70)
    print("PROOF ATTEMPT: EXTREME CASES")
    print("=" * 70)
    
    # STARS: What are the statistics?
    print("\n1. STARS: S(k) with k leaves")
    print("-" * 40)
    
    for k in [2, 3, 5, 10, 20, 50, 100]:
        n, adj = star_tree(k)
        poly = independence_poly(n, adj)
        nm, pos = near_miss_ratio(poly)
        
        # Independence polynomial of star:
        # i_0 = 1 (empty set)
        # i_1 = k+1 (center + any one leaf)
        # i_2 = C(k,2) (any 2 leaves)
        # etc.
        # Actually: for star S_k:
        # I(S_k, x) = 1 + (k+1)x + C(k,2)x^2 + ... + C(k,k)x^k
        
        print(f"  S({k}): n={n}, nm={nm:.4f}")
        print(f"    poly: {poly}")
    
    # PATHS: What are the statistics?
    print("\n2. PATHS: P(n)")
    print("-" * 40)
    
    for n in [5, 10, 15, 20, 30]:
        _, adj = path_tree(n)
        poly = independence_poly(n, adj)
        nm, pos = near_miss_ratio(poly)
        print(f"  P({n}): nm={nm:.4f}")
    
    # ============================================================================
    # PART 2: CONNECT TO INFLUENCE
    # ============================================================================
    print("\n" + "=" * 70)
    print("CONNECTING TO INFLUENCE DISTRIBUTION")
    print("=" * 70)
    
    # For a star S_k:
    # - Center (vertex 0): can be in IS of size 1, 2, ..., k+1 (any subset of leaves)
    #   Actually: if center is in IS, we can choose any subset of leaves to exclude
    #   So center can contribute to ALL coefficient positions 1 through k+1
    # - Each leaf: can be in IS of size 1 (alone) or with other leaves but NOT center
    #   Leaves contribute to positions 1, 2, ..., k (but not k+1)
    
    # This is HIGH concentration: center affects ALL positions, leaves affect early positions
    
    # For a path P_n:
    # - Each vertex can be in IS with various combinations
    # - The influence is more uniformly distributed
    
    # ============================================================================
    # PART 3: TRY A PROOF
    # ============================================================================
    print("\n" + "=" * 70)
    print("PROOF SKETCH")
    print("=" * 70)
    
    proof = r"""
================================================================================
THEOREM (Informal): Influence Concentration → Near-Miss

Let T be a tree. Define:
- A(v) = {k : v belongs to some independent set of size k}
- a(v) = average of A(v)

CLAIM: The variance of {a(v)} is positively correlated with near-miss ratio.

================================================================================
PROOF SKETCH:
================================================================================

1. POLYNOMIAL DECOMPOSITION:
   The independence polynomial can be written as:
       I(T,x) = ∏_{v∈V(T)} (1 + x * f_v(x))
   
   where f_v captures v's contribution. The structure of f_v depends on
   the local neighborhood of v.

2. INFLUENCE CONCENTRATION:
   If a vertex v has highly concentrated influence (can only affect a 
   narrow range of coefficients), it contributes to specific parts of
   the polynomial shape.

3. COMPETING INFLUENCES:
   When different vertices have very different influence ranges, the
   polynomial's rise and fall are determined by different "leaders" at
   different coefficient positions. This creates potential for 
   non-monotonic behavior.

4. STAR CASE (HIGH G, HIGH NM):
   For star S_k:
   - Center influences ALL positions k=1,...,k+1 (very wide range)
   - Leaves influence only positions k=1,...,k (narrow range)
   - This creates high variance in influence → high near-miss

5. PATH CASE (LOW G, LOW NM):
   For path P_n:
   - All vertices have similar influence ranges (centered around n/4)
   - Low variance → low near-miss

================================================================================
FORMALIZATION NEEDED:
================================================================================

We need to show:
   |I'(1)| / I(α)  or similar metric correlates with influence variance

where I'(1) is derivative at 1, α is the position of maximum coefficient.

Actually, the near-miss ratio measures:
   nm = max_{j≥d} (i_{j+1} / i_j)

where d is the first descent. This is the maximum "rise" after the peak.

If influences are concentrated, different coefficient positions are 
dominated by different subsets of vertices → potential for rises.

================================================================================
""".format()
    
    print(proof)


if __name__ == '__main__':
    main()
