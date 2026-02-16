#!/usr/bin/env python3
"""
Search for unimodality counterexamples among Homeomorphically Irreducible (HI) trees.
HI trees (no degree-2 vertices) are the candidates for minimal counterexamples if ECMS holds.
They are rare and structurally dense.
"""

import sys
import random
import networkx as nx
from indpoly import independence_poly, is_unimodal, is_log_concave, near_miss_ratio

def random_hi_tree(n):
    """Generate a random HI tree on n vertices."""
    # Internal vertices count k.
    # Constraints: 
    # 1. Internal vertices induce a tree T_I.
    # 2. Every internal vertex must have degree >= 3 in T.
    # 3. degree_T(v) = degree_T_I(v) + leaves(v).
    # 4. So leaves(v) >= 3 - degree_T_I(v).
    # 5. Sum of leaves = n - k.
    # 6. Sum of required leaves = Sum max(0, 3 - deg_I(v)).
    # 7. We need n - k >= Sum required leaves.
    
    # Max k: min leaves is when T_I is a path (max degree 2, leaves 1).
    # For a path of k nodes:
    # 2 ends need 2 leaves each (total 4).
    # k-2 internal need 1 leaf each (total k-2).
    # Total required = k+2.
    # So n - k >= k + 2 => 2k <= n - 2 => k <= (n-2)/2.
    
    # Min k: 1 (star). degree 0 in T_I. leaves >= 3. n-1 >= 3 => n >= 4.
    
    if n < 4: return None # HI trees start at n=4? No, n=2 (K2) is HI. 
                          # But my algo assumes internal nodes.
                          # K2 has 0 internal nodes? Or 2? 
                          # Vertices of K2 have degree 1. So 0 internal.
                          # HI def: no degree 2.
    
    max_k = (n - 2) // 2
    if max_k < 1: return None
    
    # Pick k
    k = random.randint(1, max_k)
    
    # Generate random tree on k vertices
    # Vertices 0..k-1 are internal.
    T_I = nx.random_labeled_tree(k)
    
    degrees = dict(T_I.degree())
    required_leaves = {}
    total_req = 0
    
    for v in range(k):
        d = degrees[v]
        req = max(0, 3 - d)
        required_leaves[v] = req
        total_req += req
        
    if n - k < total_req:
        return None # Should not happen with k <= (n-2)/2 bound, but just in case
        
    # Distribute leaves
    leaves_count = {v: required_leaves[v] for v in range(k)}
    remaining = (n - k) - total_req
    
    # Distribute remaining leaves randomly among internal vertices
    for _ in range(remaining):
        v = random.randint(0, k-1)
        leaves_count[v] += 1
        
    # Construct adjacency list
    # 0..k-1 are internal.
    # k..n-1 are leaves.
    
    adj = [[] for _ in range(n)]
    
    # Add internal edges
    for u, v in T_I.edges():
        adj[u].append(v)
        adj[v].append(u)
        
    # Add leaf edges
    leaf_idx = k
    for v in range(k):
        count = leaves_count[v]
        for _ in range(count):
            leaf = leaf_idx
            adj[v].append(leaf)
            adj[leaf].append(v)
            leaf_idx += 1
            
    return adj

def main():
    print("Searching HI trees for unimodality counterexamples...")
    print(f"{'n':>4} {'count':>8} {'uni_fail':>8} {'lc_fail':>8} {'max_nm':>8}")
    print("-" * 50)
    
    sizes = [25, 30, 35, 40, 50, 60, 70, 80, 90, 100]
    samples = 1000
    
    for n in sizes:
        uni_fails = 0
        lc_fails = 0
        max_nm = 0.0
        
        for _ in range(samples):
            adj = None
            while adj is None:
                adj = random_hi_tree(n)
            
            poly = independence_poly(n, adj)
            
            if not is_unimodal(poly):
                uni_fails += 1
                # print(f"FOUND COUNTEREXAMPLE at n={n}!")
                # print(adj)
                
            if not is_log_concave(poly):
                lc_fails += 1
                
            nm, _ = near_miss_ratio(poly)
            if nm > max_nm:
                max_nm = nm
                
        print(f"{n:4d} {samples:8d} {uni_fails:8d} {lc_fails:8d} {max_nm:8.4f}")

if __name__ == "__main__":
    main()
