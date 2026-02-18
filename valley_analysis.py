#!/usr/bin/env python3
"""
Why don't boundary violations create valleys?

When ΔA_d + ΔI_d > 0, what's preventing a valley?
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly
import networkx as nx


def analyze_k13_case():
    """Analyze the K_{1,3} case in detail."""
    print("=" * 70)
    print("K_{1,3} CASE STUDY")
    print("=" * 70)
    
    # Original tree
    G = nx.star_graph(3)
    n = 4
    adj = [list(G.neighbors(v)) for v in range(n)]
    poly = independence_poly(n, adj)
    
    print(f"\nOriginal tree: {poly}")
    print(f"Coefficients: i_0={poly[0]}, i_1={poly[1]}, i_2={poly[2]}, i_3={poly[3]}")
    
    # Find d
    d = 1  # First descent
    for i in range(1, len(poly)):
        if poly[i] < poly[i-1]:
            d = i
            break
    
    print(f"\nd(I) = {d}")
    print(f"ΔI_{d} = I_{d} - I_{d-1} = {poly[d]} - {poly[d-1]} = {poly[d] - poly[d-1]}")
    
    # Subdivide edge 0-1
    G2 = G.copy()
    w = 4
    G2.remove_edge(0, 1)
    G2.add_edge(0, w)
    G2.add_edge(w, 1)
    
    adj2 = [list(G2.neighbors(v)) for v in range(5)]
    poly2 = independence_poly(5, adj2)
    
    print(f"\nAfter subdivision: {poly2}")
    
    # Compute A
    max_len = max(len(poly), len(poly2))
    poly = poly + [0] * (max_len - len(poly))
    poly2 = poly2 + [0] * (max_len - len(poly2))
    A = [poly2[i] - poly[i] for i in range(max_len)]
    
    print(f"A(x) = I(T') - I(T): {A}")
    
    # The sums
    S = [poly[i] + (A[i] if i < len(A) else 0) for i in range(max_len)]
    print(f"\nI(T) + A: {S[:5]}")
    
    # Check the critical indices
    print(f"\n--- Analysis at d = {d} ---")
    print(f"I_{d-1} = {poly[d-1]}, A_{d-1} = {A[d-1] if d-1 < len(A) else 0}, sum = {poly[d-1] + (A[d-1] if d-1 < len(A) else 0)}")
    print(f"I_{d} = {poly[d]}, A_{d} = {A[d] if d < len(A) else 0}, sum = {poly[d] + (A[d] if d < len(A) else 0)}")
    
    if d+1 < len(poly):
        print(f"I_{d+1} = {poly[d+1]}, A_{d+1} = {A[d+1] if d+1 < len(A) else 0}, sum = {poly[d+1] + (A[d+1] if d+1 < len(A) else 0)}")
    
    print(f"\n--- Differences ---")
    print(f"ΔI_d = {poly[d] - poly[d-1]}")
    print(f"ΔA_d = {(A[d] if d < len(A) else 0) - (A[d-1] if d-1 < len(A) else 0)}")
    print(f"Δ(I+A)_d = {(poly[d] + (A[d] if d < len(A) else 0)) - (poly[d-1] + (A[d-1] if d-1 < len(A) else 0))}")
    
    print(f"\n--- Is there a valley? ---")
    # Valley would be: S_{d-1} > S_d < S_{d+1}
    S_d_minus_1 = poly[d-1] + (A[d-1] if d-1 < len(A) else 0)
    S_d = poly[d] + (A[d] if d < len(A) else 0)
    S_d_plus_1 = poly[d+1] + (A[d+1] if d+1 < len(A) else 0) if d+1 < len(poly) else 0
    
    print(f"S_{d-1} = {S_d_minus_1}")
    print(f"S_d = {S_d}")
    print(f"S_{d+1} = {S_d_plus_1}")
    
    is_valley = S_d_minus_1 > S_d and S_d < S_d_plus_1
    print(f"Valley at d: {is_valley}")
    
    # Even though ΔA_d + ΔI_d > 0, there's no valley
    # Because the tail is so strongly decreasing


if __name__ == '__main__':
    analyze_k13_case()
