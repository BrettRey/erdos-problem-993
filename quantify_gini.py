#!/usr/bin/env python3
"""
Quantifying the Gini-Unimodality relationship.

The data shows: higher Gini of influence position → higher near-miss ratio.

Let's fit: nm ≈ 1 - C * (1 - G)^α
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly, near_miss_ratio
from targeted import make_broom
from collections import defaultdict
import math


def compute_vertex_contributions(n, adj):
    """Compute set of k values each vertex can influence."""
    parent = [-1] * n
    children = [[] for _ in range(n)]
    
    queue = [0]
    parent[0] = 0
    head = 0
    while head < len(queue):
        v = queue[head]
        head += 1
        for u in adj[v]:
            if parent[u] == -1:
                parent[u] = v
                children[v].append(u)
                queue.append(u)
    
    inc_sets = [set() for _ in range(n)]
    exc_sets = [set() for _ in range(n)]
    
    post_order = queue[::-1]
    
    for v in post_order:
        if not children[v]:
            inc_sets[v] = {1}
            exc_sets[v] = {0}
        else:
            exc_sets[v] = {0}
            for c in children[v]:
                new_set = set()
                for s1 in list(exc_sets[v]):
                    for s2 in list(exc_sets[c]):
                        new_set.add(s1 + s2)
                    for s2 in list(inc_sets[c]):
                        new_set.add(s1 + s2)
                exc_sets[v] = new_set
            
            inc_sets[v] = {1}
            for c in children[v]:
                new_set = set()
                for s1 in list(inc_sets[v]):
                    for s2 in list(exc_sets[c]):
                        new_set.add(s1 + s2)
                inc_sets[v] = new_set
    
    return [inc_sets[v] | exc_sets[v] for v in range(n)]


def compute_gini_of_influence(n, adj):
    """Compute Gini of average k-position each vertex can influence."""
    influence = compute_vertex_contributions(n, adj)
    
    avg_k_per_vertex = []
    for v in range(n):
        if influence[v]:
            avg_k_per_vertex.append(sum(influence[v]) / len(influence[v]))
        else:
            avg_k_per_vertex.append(0)
    
    sorted_avg = sorted(avg_k_per_vertex)
    N = len(sorted_avg)
    gini = sum((2*(i+1) - N - 1) * sorted_avg[i] for i in range(N))
    gini /= (N * sum(sorted_avg)) if sum(sorted_avg) > 0 else 1
    
    return gini


def main():
    print("=" * 70)
    print("QUANTIFYING THE GINI-NM RELATIONSHIP")
    print("=" * 70)
    
    # Gather data
    data = []
    for p in [3, 5, 10, 15, 20, 30, 50]:
        for s in [10, 30, 50, 100, 200, 500, 1000, 2000]:
            n, adj = make_broom(p, s)
            poly = independence_poly(n, adj)
            nm, _ = near_miss_ratio(poly)
            gini = compute_gini_of_influence(n, adj)
            
            data.append({
                'p': p, 's': s, 'n': n,
                'nm': nm,
                'gini': gini,
            })
    
    print(f"\nGathered {len(data)} data points")
    
    # Try different fits
    print("\nFITTING: nm = 1 - C * (1 - g)^α")
    
    # Linear fit in (1-G) space
    # nm = 1 - C*(1-G)  →  1-nm = C*(1-G)
    linear_fits = []
    for d in data:
        if d['gini'] < 1:
            C = (1 - d['nm']) / (1 - d['gini'])
            linear_fits.append(C)
    
    avg_C = sum(linear_fits) / len(linear_fits)
    print(f"Average C (linear fit): {avg_C:.2f}")
    
    # Check residuals
    residuals = []
    for d in data:
        predicted = 1 - avg_C * (1 - d['gini'])
        residual = d['nm'] - predicted
        residuals.append(residual)
    
    avg_residual = sum(residuals) / len(residuals)
    print(f"Average residual: {avg_residual:.4f}")
    
    # Quadratic fit: nm = 1 - C1*(1-G) - C2*(1-G)^2
    print("\nTrying quadratic: nm = 1 - C1*(1-G) - C2*(1-G)^2")
    
    # Simple approach: fit by eye
    # For G near 1: nm ≈ 1 - C1*(1-G) - C2*(1-G)^2
    # At G=0.99, 1-G=0.01, nm ≈ 0.99 → 1-nm = 0.01
    # At G=0.9, 1-G=0.1, nm varies...
    
    # Just print data for manual analysis
    print("\nData for analysis:")
    print("G       1-G      nm      1-nm    (1-nm)/(1-G)")
    for d in sorted(data, key=lambda x: x['gini']):
        g = d['gini']
        nm = d['nm']
        ratio = (1-nm)/(1-g) if g < 1 else 0
        print(f"{g:.3f}   {1-g:.3f}    {nm:.4f}   {1-nm:.4f}   {ratio:.2f}")
    
    # The key insight
    print("\n" + "=" * 70)
    print("KEY INSIGHT:")
    print("=" * 70)
    print("""
The ratio (1-nm)/(1-G) appears to converge to approximately C ≈ 4-5
as G → 1. This is consistent with the broom asymptotic:

    nm(B(p,s)) = 1 - C_broom/s + O(1/s²)
    
If we can show that 1 - G(B(p,s)) ≈ c/s for some c, then:
    nm ≈ 1 - (C_broom/c) * (1 - G)
    
This would give a THEORETICAL basis for the Gini coefficient!
""")


if __name__ == '__main__':
    main()
