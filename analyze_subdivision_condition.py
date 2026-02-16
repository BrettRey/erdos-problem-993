#!/usr/bin/env python3
"""
Analyze the algebraic condition for subdivision unimodality.

We need to understand: what would make ΔA_k > -ΔI_k for some k >= d+1?

The key is that for unimodality to fail, we'd need:
i_k + a_k to have a valley after d(I).
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly, is_unimodal
from targeted import make_broom, make_spider
import networkx as nx


def compute_subdivision_terms(T, u, v):
    """Compute the terms in the subdivision identity.
    
    I(T') = I(T) + Q_u Q_v + x P_u P_v = I(T) + A(x)
    """
    n = len(T.nodes())
    adj = [list(T.neighbors(i)) for i in range(n)]
    
    # We need the rooted polynomials at u and v
    # This is expensive to compute directly, so let's use the identity
    # A(x) = Q_u Q_v + x P_u P_v
    
    # For a rough analysis, let's just compute I(T) and I(T')
    poly_T = independence_poly(n, adj)
    
    # Subdivide
    T2 = T.copy()
    w = T2.number_of_nodes()
    T2.remove_edge(u, v)
    T2.add_edge(u, w)
    T2.add_edge(w, v)
    
    n2 = T2.number_of_nodes()
    adj2 = [list(T2.neighbors(i)) for i in range(n2)]
    poly_T2 = independence_poly(n2, adj2)
    
    # A(x) = I(T') - I(T)
    max_len = max(len(poly_T), len(poly_T2))
    poly_T = poly_T + [0] * (max_len - len(poly_T))
    a = [poly_T2[i] - poly_T[i] for i in range(max_len)]
    
    return poly_T, a


def analyze_valley_condition():
    """Analyze when a valley could form."""
    print("=" * 70)
    print("ANALYZING VALLEY FORMATION CONDITION")
    print("=" * 70)
    
    # Let's test on specific trees
    test_trees = []
    
    # Brooms (known to be close to violation)
    for p in [5, 10, 15]:
        for s in [10, 30, 50]:
            n, adj = make_broom(p, s)
            G = nx.Graph()
            G.add_nodes_from(range(n))
            for i in range(n-1):
                G.add_edge(i, i+1)
            for i in range(1, 1+p):
                G.add_edge(0, i)
            test_trees.append((f"broom({p},{s})", G, 0, 1))  # subdivide edge 0-1
    
    # Stars
    for k in [5, 10, 20]:
        G = nx.star_graph(k-1)
        test_trees.append((f"star({k})", G, 0, 1))
    
    print(f"\nTesting {len(test_trees)} tree/edge combinations...")
    
    for name, G, u, v in test_trees:
        poly_T, a = compute_subdivision_terms(G, u, v)
        
        # Find first descent in I(T)
        d = -1
        for i in range(1, len(poly_T)):
            if poly_T[i] < poly_T[i-1]:
                d = i
                break
        
        if d == -1:
            continue
            
        # Check tail differences
        print(f"\n{name}: d(I) = {d}")
        
        # For k >= d+1, check ΔA_k vs -ΔI_k
        violations = []
        for k in range(d+1, min(len(a), len(poly_T))-1):
            if k+1 < len(poly_T) and k < len(a):
                delta_A = a[k+1] - a[k]
                delta_I = poly_T[k+1] - poly_T[k]
                
                if delta_A > -delta_I:  # This would break unimodality!
                    violations.append((k, delta_A, delta_I, delta_A + delta_I))
        
        if violations:
            print(f"  POTENTIAL VIOLATIONS at k={violations}")
        else:
            # Find max violation
            max_viol = 0
            for k in range(d+1, min(len(a), len(poly_T))-1):
                if k+1 < len(poly_T) and k < len(a):
                    delta_A = a[k+1] - a[k]
                    delta_I = poly_T[k+1] - poly_T[k]
                    viol = -delta_I - delta_A
                    if viol > max_viol:
                        max_viol = viol
            print(f"  Max margin: {max_viol:.1f}")


def construct_extreme_case():
    """Try to construct an extreme case."""
    print("\n" + "=" * 70)
    print("CONSTRUCTING EXTREME CASES")
    print("=" * 70)
    
    # What if we start with a tree that's ALREADY non-unimodal?
    # But we need to start with a unimodal tree...
    
    # What tree gives the largest A(x)?
    # That would maximize the perturbation
    
    # Try very unbalanced trees
    for n in [10, 15, 20]:
        # Create very unbalanced tree
        G = nx.path_graph(n)
        
        # Test various edge subdivisions
        for u in range(n-1):
            for v in [u+1]:
                poly_T, a = compute_subdivision_terms(G, u, v)
                
                # Find d
                d = -1
                for i in range(1, len(poly_T)):
                    if poly_T[i] < poly_T[i-1]:
                        d = i
                        break
                
                if d == -1:
                    continue
                
                # Check magnitude of A
                a_sum = sum(a)
                if a_sum > 5:
                    print(f"Path n={n}, edge {u}-{v}: d={d}, |A|={a_sum}")


if __name__ == '__main__':
    analyze_valley_condition()
    construct_extreme_case()
