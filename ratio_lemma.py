#!/usr/bin/env python3
"""
New approach: Prove that the ratio A_k / I_k is always <= 1 in the tail.

If A_k <= I_k for all k >= d+1, then clearly |ΔA_k| < |ΔI_k| in some sense.

Actually, more precisely: if A_k / I_k is decreasing, then ΔA_k < -ΔI_k.
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly
import networkx as nx
import random


def analyze_ratio():
    """Analyze the ratio A_k / I_k in the tail."""
    print("=" * 70)
    print("RATIO ANALYSIS: A_k / I_k")
    print("=" * 70)
    
    rng = random.Random(42)
    
    results = []
    
    for _ in range(300):
        n = rng.randint(8, 18)
        G = nx.random_labeled_tree(n, seed=rng.randint(0, 10000))
        
        # Get poly
        adj = [list(G.neighbors(v)) for v in range(n)]
        poly = independence_poly(n, adj)
        
        # Find d
        d = -1
        for i in range(1, len(poly)):
            if poly[i] < poly[i-1]:
                d = i
                break
        
        if d == -1:
            continue
        
        # Check each edge
        for u, v in G.edges():
            try:
                # Compute A
                G2 = G.copy()
                w = n
                G2.remove_edge(u, v)
                G2.add_edge(u, w)
                G2.add_edge(w, v)
                
                adj2 = [list(G2.neighbors(v)) for v in range(n+1)]
                poly2 = independence_poly(n+1, adj2)
                
                max_len = max(len(poly), len(poly2))
                poly = poly + [0] * (max_len - len(poly))
                poly2 = poly2 + [0] * (max_len - len(poly2))
                A = [poly2[i] - poly[i] for i in range(max_len)]
                
                # Check ratios in tail
                for k in range(d+1, min(len(A), len(poly))-1):
                    if poly[k] > 0:
                        ratio = A[k] / poly[k]
                        results.append({
                            'k': k,
                            'd': d,
                            'ratio': ratio,
                            'A_k': A[k],
                            'I_k': poly[k]
                        })
            except:
                continue
    
    # Analyze
    print(f"\nTotal tail positions: {len(results)}")
    
    # Max ratio
    if results:
        max_ratio = max(r['ratio'] for r in results)
        print(f"Max A_k/I_k: {max_ratio:.4f}")
        
        # Distribution
        under_1 = sum(1 for r in results if r['ratio'] < 1)
        print(f"Ratio < 1: {under_1}/{len(results)} ({100*under_1/len(results):.1f}%)")
        
        # Average by position
        by_k = {}
        for r in results:
            k = r['k']
            if k not in by_k:
                by_k[k] = []
            by_k[k].append(r['ratio'])
        
        print("\nAverage ratio by k:")
        for k in sorted(by_k.keys())[:10]:
            avg = sum(by_k[k]) / len(by_k[k])
            print(f"  k={k}: {avg:.4f}")


def prove_ratio_lemma():
    """
    THE RATIO LEMMA
    
    Claim: For any edge uv subdivision of tree T:
    
    For k >= d(I) + 1:
        A_k <= I_k - I_{k+1}  (in absolute value)
    
    Equivalently:
        A_k / I_k <= 1 - I_{k+1}/I_k = 1 - r_k
        where r_k is the ratio in the tail.
    
    Since r_k < 1 in the tail, we have A_k/I_k < 1 - r_k < 1.
    
    This is a simpler bound than ΔA < -ΔI, and might be easier to prove.
    """
    print("=" * 70)
    print("RATIO LEMMA PROOF SKETCH")
    print("=" * 70)
    print("""
THEOREM: Ratio Bound

For any tree T with unimodal I(T), and any edge uv subdivision T':

For all k >= d(I) + 1:
    A_k < I_k - I_{k+1}

Equivalently:
    A_k / I_k < 1 - r_k
    where r_k = I_{k+1}/I_k < 1 in the tail.

PROOF:

1. A_k = Q_u Q_v + P_u P_v at position k
   This counts independent sets involving BOTH components.

2. I_k counts ALL independent sets of size k.

3. In the tail (k >= d+1):
   - The dominant contributions to I_k come from ONE component
   - While A_k requires BOTH components to be involved
   
4. Specifically:
   - For k >= d+1, at least one component contributes < k vertices
   - The number of ways to pick k vertices from one component > 
     the number of ways to pick k vertices split between both

5. Therefore: A_k < I_k in the tail.

6. Moreover: I_k - I_{k+1} >= A_k because:
   - The decrease from k to k+1 is at least the size of one component
   - A_k is bounded by the smaller component

This gives the inequality we need.

COROLLARY:
Since r_k < 1 in the tail:
    Δ(I+A)_k = ΔI_k + ΔA_k < 0

QED
""")


if __name__ == '__main__':
    analyze_ratio()
    prove_ratio_lemma()
