#!/usr/bin/env python3
"""
Search for counterexample to subdivision lemma.

The lemma: If T is unimodal, subdividing any edge preserves unimodality.

We'll search for trees where this might fail.
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly, is_unimodal, near_miss_ratio
import networkx as nx
import random


def subdivide_edge(G, u, v):
    """Subdivide edge uv with new vertex w."""
    G = G.copy()
    w = G.number_of_nodes()
    G.remove_edge(u, v)
    G.add_edge(u, w)
    G.add_edge(w, v)
    return G


def test_subdivision(G, verbose=False):
    """Test if subdivision preserves unimodality for all edges."""
    n = G.number_of_nodes()
    adj = [list(G.neighbors(v)) for v in range(n)]
    
    # Check original tree
    poly = independence_poly(n, adj)
    orig_uni = is_unimodal(poly)
    
    if not orig_uni:
        return True, "original not unimodal"
    
    # Try subdividing each edge
    for u, v in G.edges():
        G_sub = subdivide_edge(G, u, v)
        n2 = G_sub.number_of_nodes()
        adj2 = [list(G_sub.neighbors(v)) for v in range(n2)]
        poly2 = independence_poly(n2, adj2)
        
        if not is_unimodal(poly2):
            if verbose:
                print(f"FOUND: Edge {u}-{v} subdivision creates non-unimodal!")
                print(f"  Original: {poly}")
                print(f"  Subdivided: {poly2}")
            return False, f"edge {u}-{v}"
    
    return True, "all subdivisions unimodal"


def search_for_counterexample(max_n=20, num_trees=1000, seed=42):
    """Search for counterexample."""
    rng = random.Random(seed)
    
    print(f"Searching for counterexample (n <= {max_n}, {num_trees} trees)...")
    
    found = []
    for i in range(num_trees):
        n = rng.randint(3, max_n)
        G = nx.random_labeled_tree(n, seed=rng.randint(0, 100000))
        
        ok, reason = test_subdivision(G)
        if not ok:
            found.append((n, reason, G))
            print(f"Found potential counterexample: n={n}, {reason}")
    
    if found:
        print(f"\nFound {len(found)} potential counterexamples!")
    else:
        print(f"\nNo counterexamples found in {num_trees} trees")
    
    return found


def focused_search():
    """Focus on structures that might be problematic."""
    
    # What might break it? 
    # Trees with already high near-miss ratio might be close to violation
    
    print("Focused search on high-nm trees...")
    
    # Generate many trees and pick those close to violating
    candidates = []
    rng = random.Random(12345)
    
    for _ in range(5000):
        n = rng.randint(10, 25)
        G = nx.random_labeled_tree(n, seed=rng.randint(0, 100000))
        adj = [list(G.neighbors(v)) for v in range(n)]
        poly = independence_poly(n, adj)
        nm, _ = near_miss_ratio(poly)
        
        if nm > 0.8:  # Close to violation
            candidates.append((nm, G, poly, n))
    
    candidates.sort(reverse=True)
    print(f"Found {len(candidates)} high-nm trees")
    
    # Test subdivisions on top candidates
    for nm, G, poly, n in candidates[:50]:
        ok, reason = test_subdivision(G)
        if not ok:
            print(f"FOUND: n={n}, nm={nm:.4f}, {reason}")
            return
    
    print("No counterexample in top candidates")


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, default=15)
    parser.add_argument('--num', type=int, default=500)
    args = parser.parse_args()
    
    search_for_counterexample(args.n, args.num)
    focused_search()
