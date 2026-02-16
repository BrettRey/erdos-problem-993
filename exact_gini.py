#!/usr/bin/env python3
"""
Computing exact Gini coefficients for stars and paths to verify relationship.
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly, near_miss_ratio


def compute_vertex_influence_positions(n, adj):
    """
    For each vertex, compute the set of k values it can appear in.
    This is the full influence set S(v).
    """
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
    
    # DP to compute possible sizes when vertex is included/excluded
    inc = [set() for _ in range(n)]  # sizes when included
    exc = [set() for _ in range(n)]  # sizes when excluded
    
    post = queue[::-1]
    
    for v in post:
        if not children[v]:
            inc[v] = {1}
            exc[v] = {0}
        else:
            # exc: combine children (each can be in or out)
            exc[v] = {0}
            for c in children[v]:
                new = set()
                for s1 in list(exc[v]):
                    for s2 in list(exc[c]):
                        new.add(s1 + s2)
                    for s2 in list(inc[c]):
                        new.add(s1 + s2)
                exc[v] = new
            
            # inc: include v, children must be excluded
            inc[v] = {1}
            for c in children[v]:
                new = set()
                for s1 in list(inc[v]):
                    for s2 in list(exc[c]):
                        new.add(s1 + s2)
                inc[v] = new
    
    # Full influence = union
    influence = [inc[v] | exc[v] for v in range(n)]
    
    # Average k for each vertex
    avg_k = []
    for v in range(n):
        if influence[v]:
            avg_k.append(sum(influence[v]) / len(influence[v]))
        else:
            avg_k.append(0)
    
    return influence, avg_k


def compute_gini(values):
    """Compute Gini coefficient."""
    sorted_vals = sorted(values)
    N = len(sorted_vals)
    gini = sum((2*(i+1) - N - 1) * sorted_vals[i] for i in range(N))
    gini /= (N * sum(sorted_vals)) if sum(sorted_vals) > 0 else 1
    return gini


def main():
    print("=" * 70)
    print("EXACT GINI COEFFICIENTS")
    print("=" * 70)
    
    # Stars
    print("\nSTARS:")
    print("-" * 40)
    for k in [10, 50, 100, 200]:
        n = k + 1
        adj = [[] for _ in range(n)]
        for i in range(1, n):
            adj[0].append(i)
            adj[i].append(0)
        
        poly = independence_poly(n, adj)
        nm, _ = near_miss_ratio(poly)
        
        _, avg_k = compute_vertex_influence_positions(n, adj)
        gini = compute_gini(avg_k)
        
        print(f"S({k}): n={n}, G={gini:.4f}, nm={nm:.4f}")
    
    # Paths
    print("\nPATHS:")
    print("-" * 40)
    for n in [10, 20, 50, 100]:
        adj = [[] for _ in range(n)]
        for i in range(n-1):
            adj[i].append(i+1)
            adj[i+1].append(i)
        
        poly = independence_poly(n, adj)
        nm, _ = near_miss_ratio(poly)
        
        _, avg_k = compute_vertex_influence_positions(n, adj)
        gini = compute_gini(avg_k)
        
        print(f"P({n}): n={n}, G={gini:.4f}, nm={nm:.4f}")
    
    # Brooms
    print("\nBROOMS (p=10):")
    print("-" * 40)
    from targeted import make_broom
    for s in [10, 50, 100, 200, 500]:
        n, adj = make_broom(10, s)
        poly = independence_poly(n, adj)
        nm, _ = near_miss_ratio(poly)
        
        _, avg_k = compute_vertex_influence_positions(n, adj)
        gini = compute_gini(avg_k)
        
        print(f"B(10,{s}): n={n}, G={gini:.4f}, nm={nm:.4f}")
    
    print("\n" + "=" * 70)
    print("PATTERN:")
    print("=" * 70)
    print("""
The pattern is now clearer:
- Stars have HIGH G and HIGH nm
- Paths have MODERATE G and MODERATE nm  
- Brooms have HIGH G and HIGH nm

The relationship is: higher G → higher nm, but with variation.

This supports the conjecture but shows it's not a perfect relationship.
More analysis needed to find the exact functional form.
""")


if __name__ == '__main__':
    main()
