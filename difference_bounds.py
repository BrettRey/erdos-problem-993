#!/usr/bin/env python3
"""
Different approach: Look at what bound actually works.

We need: ΔA_k < -ΔI_k, i.e., A_{k+1} - A_k < I_k - I_{k+1}

This is: A_{k+1} + I_{k+1} < A_k + I_k

Equivalently: S_{k+1} < S_k where S = I + A.

But we want to prove S is decreasing. Let's check if A_k <= I_{k-1} might work.
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly
import networkx as nx
import random


def test_bound():
    """Test various bounds."""
    print("=" * 70)
    print("TESTING BOUNDS")
    print("=" * 70)
    
    rng = random.Random(42)
    
    results = []
    
    for _ in range(300):
        n = rng.randint(8, 18)
        G = nx.random_labeled_tree(n, seed=rng.randint(0, 10000))
        
        adj = [list(G.neighbors(v)) for v in range(n)]
        poly = independence_poly(n, adj)
        
        d = -1
        for i in range(1, len(poly)):
            if poly[i] < poly[i-1]:
                d = i
                break
        
        if d == -1:
            continue
        
        for u, v in G.edges():
            try:
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
                
                # Test different bounds
                for k in range(d+1, min(len(A), len(poly))-1):
                    if k > 0:
                        results.append({
                            'k': k,
                            'd': d,
                            'A_k': A[k],
                            'A_k+1': A[k+1] if k+1 < len(A) else 0,
                            'I_k': poly[k],
                            'I_k+1': poly[k+1] if k+1 < len(poly) else 0,
                            'I_k-1': poly[k-1],
                        })
            except:
                continue
    
    # Test: A_k <= I_{k-1} - I_k
    # i.e., A_k + I_k <= I_{k-1}
    test1 = sum(1 for r in results 
                if r['k'] < len(r.get('I_k-1', 0)) + r['I_k'] 
                and r['A_k'] <= r['I_k-1'] - r['I_k'])
    
    # Test: A_k <= I_k - I_{k+1} (would prove decrease directly)
    test2 = sum(1 for r in results if r['A_k'] <= r['I_k'] - r['I_k+1'])
    
    # Test: A_{k+1} <= A_k
    test3 = sum(1 for r in results 
                if r['k'] + 1 < len(r.get('A_k+1', 0)) 
                and r['A_k+1'] <= r['A_k'])
    
    print(f"\nTotal positions: {len(results)}")
    print(f"A_k <= I_{{k-1}} - I_k: {test1} ({100*test1/len(results):.1f}%)")
    print(f"A_k <= I_k - I_{{k+1}}: {test2} ({100*test2/len(results):.1f}%)")
    print(f"A_{{k+1}} <= A_k: {test3} ({100*test3/len(results):.1f}%)")


def analyze_differences():
    """Analyze the actual differences."""
    print("\n" + "=" * 70)
    print("ACTUAL DIFFERENCES")
    print("=" * 70)
    
    rng = random.Random(42)
    
    # Collect actual Δ values
    deltas = []
    
    for _ in range(300):
        n = rng.randint(8, 18)
        G = nx.random_labeled_tree(n, seed=rng.randint(0, 10000))
        
        adj = [list(G.neighbors(v)) for v in range(n)]
        poly = independence_poly(n, adj)
        
        d = -1
        for i in range(1, len(poly)):
            if poly[i] < poly[i-1]:
                d = i
                break
        
        if d == -1:
            continue
        
        for u, v in G.edges():
            try:
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
                
                for k in range(d+1, min(len(A), len(poly))-1):
                    delta_I = poly[k+1] - poly[k]
                    delta_A = A[k+1] - A[k]
                    
                    deltas.append({
                        'k': k,
                        'd': d,
                        'delta_I': delta_I,
                        'delta_A': delta_A,
                        'sum': delta_I + delta_A
                    })
            except:
                continue
    
    print(f"\nAnalyzed {len(deltas)} tail positions")
    
    # How many have ΔI + ΔA < 0?
    negative = sum(1 for d in deltas if d['sum'] < 0)
    print(f"ΔI + ΔA < 0: {negative} ({100*negative/len(deltas):.1f}%)")
    
    # Distribution of sums
    sums = [d['sum'] for d in deltas]
    print(f"\nSum distribution:")
    print(f"  Min: {min(sums)}")
    print(f"  Max: {max(sums)}")
    print(f"  Avg: {sum(sums)/len(sums):.1f}")


if __name__ == '__main__':
    test_bound()
    analyze_differences()
