#!/usr/bin/env python3
"""
Analyze what bound actually works.
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly
import networkx as nx
import random


def analyze_differences():
    """Analyze actual Δ values."""
    print("=" * 70)
    print("ACTUAL DIFFERENCES: ΔI + ΔA")
    print("=" * 70)
    
    rng = random.Random(42)
    
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
    print(f"ΔI + ΔA < 0: {negative}/{len(deltas)} ({100*negative/len(deltas):.1f}%)")
    
    # Distribution
    sums = [d['sum'] for d in deltas]
    print(f"\nSum distribution:")
    print(f"  Min: {min(sums)}")
    print(f"  Max: {max(sums)}")
    print(f"  Avg: {sum(sums)/len(sums):.1f}")
    
    # By k position
    by_k = {}
    for d in deltas:
        k = d['k']
        if k not in by_k:
            by_k[k] = []
        by_k[k].append(d['sum'])
    
    print("\nAverage sum by k:")
    for k in sorted(by_k.keys())[:8]:
        avg = sum(by_k[k]) / len(by_k[k])
        print(f"  k={k}: avg sum = {avg:.1f}")


def test_bounds():
    """Test specific bounds."""
    print("\n" + "=" * 70)
    print("TESTING SPECIFIC BOUNDS")
    print("=" * 70)
    
    rng = random.Random(42)
    
    tests = {
        'A_k <= I_{k-1}': 0,
        'A_k <= I_k - I_{k+1}': 0,
        'A_{k+1} <= A_k': 0,
        'A_k <= (I_k - I_{k+1})/2': 0,
    }
    total = 0
    
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
                    total += 1
                    
                    # Test 1: A_k <= I_{k-1}
                    if k > 0 and A[k] <= poly[k-1]:
                        tests['A_k <= I_{k-1}'] += 1
                    
                    # Test 2: A_k <= I_k - I_{k+1}
                    if A[k] <= poly[k] - poly[k+1]:
                        tests['A_k <= I_k - I_{k+1}'] += 1
                    
                    # Test 3: A_{k+1} <= A_k
                    if k+1 < len(A) and A[k+1] <= A[k]:
                        tests['A_{k+1} <= A_k'] += 1
                    
                    # Test 4: A_k <= half the difference
                    if A[k] <= (poly[k] - poly[k+1]) / 2:
                        tests['A_k <= (I_k - I_{k+1})/2'] += 1
            except:
                continue
    
    print(f"\nTotal positions tested: {total}")
    for name, count in tests.items():
        print(f"  {name}: {count}/{total} ({100*count/total:.1f}%)")


if __name__ == '__main__':
    analyze_differences()
    test_bounds()
