#!/usr/bin/env python3
"""
Analyze the BOUNDARY condition at k = d(I).

The tail condition is STRONGLY satisfied (by huge margins).
The only potential issue is at k = d(I).
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly, is_unimodal
from targeted import make_broom
import networkx as nx


def analyze_boundary():
    """Analyze the boundary condition."""
    print("=" * 70)
    print("BOUNDARY CONDITION ANALYSIS")
    print("=" * 70)
    
    # The condition: ΔA_d <= -ΔI_d
    # or equivalently: A_d - A_{d-1} <= I_{d-1} - I_d
    
    # For all edges in various trees, check this
    
    test_trees = []
    
    # Brooms
    for p in [5, 10,15]:
        for s in [10, 30, 50, 100]:
            n, adj = make_broom(p, s)
            G = nx.Graph()
            for i in range(n-1):
                G.add_edge(i, i+1)
            for i in range(1, p+1):
                G.add_edge(0, i)
            test_trees.append((f"broom({p},{s})", G))
    
    # Paths
    for n in [10, 20, 30]:
        G = nx.path_graph(n)
        test_trees.append((f"path({n})", G))
    
    print("\nChecking boundary condition at k = d(I):")
    print("-" * 70)
    
    violations = []
    
    for name, G in test_trees:
        n = G.number_of_nodes()
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
        
        # Test each edge
        for u, v in G.edges():
            # Compute A(x) for this edge
            G2 = G.copy()
            w = G2.number_of_nodes()
            G2.remove_edge(u, v)
            G2.add_edge(u, w)
            G2.add_edge(w, v)
            
            n2 = G2.number_of_nodes()
            adj2 = [list(G2.neighbors(v)) for v in range(n2)]
            poly2 = independence_poly(n2, adj2)
            
            # A = poly2 - poly
            max_len = max(len(poly), len(poly2))
            poly = poly + [0] * (max_len - len(poly))
            poly2 = poly2 + [0] * (max_len - len(poly2))
            A = [poly2[i] - poly[i] for i in range(max_len)]
            
            # Check boundary at k=d
            if d > 0 and d < len(A) and d < len(poly):
                delta_A = A[d] - A[d-1]
                delta_I = poly[d] - poly[d-1]
                
                # Condition: delta_A <= -delta_I
                # i.e., delta_A + delta_I <= 0
                if delta_A + delta_I > 0:
                    violations.append((name, u, v, d, delta_A, delta_I, delta_A + delta_I))
    
    if violations:
        print(f"\nFOUND {len(violations)} BOUNDARY VIOLATIONS:")
        for name, u, v, d, dA, dI, sum in violations[:10]:
            print(f"  {name}, edge {u}-{v}, d={d}: ΔA={dA}, ΔI={dI}, sum={sum}")
    else:
        print("\nNo boundary violations found in these tests!")


def find_minimal_boundary_case():
    """Find the minimal case where boundary fails."""
    print("\n" + "=" * 70)
    print("SEARCHING FOR MINIMAL BOUNDARY VIOLATION")
    print("=" * 70)
    
    # We know K_{1,3} (star with 3 leaves) fails
    # Let's verify
    
    G = nx.star_graph(3)  # K_{1,3}
    n = 4
    adj = [list(G.neighbors(v)) for v in range(n)]
    poly = independence_poly(n, adj)
    
    print(f"\nStar K{{1,3}}: {poly}")
    
    # Find d
    d = -1
    for i in range(1, len(poly)):
        if poly[i] < poly[i-1]:
            d = i
            break
    
    print(f"First descent at d = {d}")
    
    # Subdivide edge 0-1
    G2 = G.copy()
    w = 4
    G2.remove_edge(0, 1)
    G2.add_edge(0, w)
    G2.add_edge(w, 1)
    
    adj2 = [list(G2.neighbors(v)) for v in range(5)]
    poly2 = independence_poly(5, adj2)
    
    A = [poly2[i] - (poly[i] if i < len(poly) else 0) for i in range(len(poly2))]
    
    print(f"I(T): {poly}")
    print(f"I(T'): {poly2}")
    print(f"A(x): {A}")
    
    if d > 0 and d < len(A) and d < len(poly):
        delta_A = A[d] - A[d-1]
        delta_I = poly[d] - poly[d-1]
        print(f"\nAt k={d}:")
        print(f"  ΔA = {delta_A}")
        print(f"  ΔI = {delta_I}")
        print(f"  ΔA + ΔI = {delta_A + delta_I}")
        print(f"  Condition violated: {delta_A + delta_I > 0}")


if __name__ == '__main__':
    analyze_boundary()
    find_minimal_boundary_case()
