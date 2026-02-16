#!/usr/bin/env python3
"""
Finding the actual condition for unimodality preservation.

The boundary condition ΔA_d <= -ΔI_d is sufficient but NOT NECESSARY.
What IS necessary and sufficient?

For I' = I + A to be unimodal:
- No valley means: S_{k-1} <= S_k >= S_{k+1} for all k

At the boundary d:
- We need S_{d-1} <= S_d  (no valley going up at d)
- Then S_d >= S_{d+1} (tail monotone decreasing)
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly
import networkx as nx


def check_condition(T, u, v):
    """Check if subdivision preserves unimodality and why."""
    n = T.number_of_nodes()
    adj = [list(T.neighbors(i)) for i in range(n)]
    
    poly = independence_poly(n, adj)
    
    # Find d
    d = -1
    for i in range(1, len(poly)):
        if poly[i] < poly[i-1]:
            d = i
            break
    
    if d == -1:
        return True, "no descent"
    
    # Subdivide
    T2 = T.copy()
    w = T2.number_of_nodes()
    T2.remove_edge(u, v)
    T2.add_edge(u, w)
    T2.add_edge(w, v)
    
    n2 = T2.number_of_nodes()
    adj2 = [list(T2.neighbors(i)) for i in range(n2)]
    poly2 = independence_poly(n2, adj2)
    
    # A
    max_len = max(len(poly), len(poly2))
    poly = poly + [0] * (max_len - len(poly))
    poly2 = poly2 + [0] * (max_len - len(poly2))
    A = [poly2[i] - poly[i] for i in range(max_len)]
    
    # S = I + A
    S = [poly[i] + A[i] for i in range(max_len)]
    
    # Check unimodality
    is_uni = True
    for i in range(1, len(S)-1):
        if S[i-1] > S[i] < S[i+1]:
            is_uni = False
            break
    
    # The actual condition at boundary:
    # Need S_{d-1} <= S_d (no valley at boundary)
    # This is: I_{d-1} + A_{d-1} <= I_d + A_d
    # Or: A_d - A_{d-1} <= I_{d-1} - I_d = -ΔI_d
    # Which is exactly: ΔA_d <= -ΔI_d
    
    # So the boundary condition IS necessary!
    # But it can be violated and still have no valley if...?
    
    # Wait - let me check more carefully
    if d > 0 and d < len(A) and d < len(poly):
        delta_A = A[d] - A[d-1]
        delta_I = poly[d] - poly[d-1]
        
        cond_violated = delta_A > -delta_I
        
        # But unimodality holds because S is still decreasing after d
        # For this to fail, we'd need S_d < S_{d+1}
        if d+1 < len(S):
            tail_decreasing = S[d] >= S[d+1]
            
            return is_uni, {
                'd': d,
                'ΔA': delta_A,
                'ΔI': delta_I,
                'boundary_violated': cond_violated,
                'tail_decreasing': tail_decreasing,
                'S': S[:5]
            }
    
    return is_uni, {'d': d, 'S': S[:5]}


def analyze_many_cases():
    """Analyze many edge subdivisions."""
    print("=" * 70)
    print("ANALYZING BOUNDARY CONDITION VS ACTUAL UNIMODALITY")
    print("=" * 70)
    
    import random
    rng = random.Random(42)
    
    boundary_violations = 0
    unimodality_preserved = 0
    
    for _ in range(500):
        n = rng.randint(5, 15)
        G = nx.random_labeled_tree(n, seed=rng.randint(0, 10000))
        
        for u, v in G.edges():
            is_uni, info = check_condition(G, u, v)
            
            if 'boundary_violated' in info:
                boundary_violations += 1
                if is_uni:
                    unimodality_preserved += 1
    
    print(f"\nBoundary violations: {boundary_violations}")
    print(f"With unimodality preserved: {unimodality_preserved}")
    print(f"Success rate despite violation: {100*unimodality_preserved/max(boundary_violations,1):.1f}%")


if __name__ == '__main__':
    analyze_many_cases()
