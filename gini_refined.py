#!/usr/bin/env python3
"""
Rethinking the Gini-Unimodality connection.

Key insight: For a vertex v, we want to know the RANGE of coefficient positions
that v can affect. A leaf can only contribute to early coefficients (small k).
A vertex deep in a "hub" can contribute to many different k values.

Let's measure: for each vertex, what's the SPREAD of its influence?
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly, near_miss_ratio
from targeted import make_broom
import json


def compute_influence_span(n, adj):
    """
    For each vertex v, compute the range of coefficient positions 
    that v can affect.
    
    More precisely: what's the set of k values where v can appear in 
    an independent set of size k?
    
    We compute this via DP: track the min and max possible IS sizes
    containing each vertex, considering the full tree (not just subtree).
    """
    # This is expensive, so let's approximate with a simpler metric:
    # Distance from each vertex to leaves vs distance to center
    
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for v in range(n):
        for u in adj[v]:
            if u > v:
                G.add_edge(v, u)
    
    # Compute eccentricity: max distance to any other vertex
    # Also compute average distance to all other vertices
    
    # For each vertex: what's the furthest it can "reach" in the tree?
    # This approximates how late in the polynomial it can contribute
    
    from collections import deque
    
    avg_reach = []
    max_reach = []
    for start in range(n):
        dist = [-1] * n
        dist[start] = 0
        q = deque([start])
        while q:
            v = q.popleft()
            for u in adj[v]:
                if dist[u] == -1:
                    dist[u] = dist[v] + 1
                    q.append(u)
        
        valid_dists = [d for d in dist if d >= 0]
        avg_reach.append(sum(valid_dists) / len(valid_dists))
        max_reach.append(max(valid_dists))
    
    return avg_reach, max_reach


def compute_gini_of_reach(n, adj):
    """Compute Gini coefficient of reach (average distance to all vertices)."""
    avg_reach, max_reach = compute_influence_span(n, adj)
    
    # Gini of avg_reach
    sorted_reach = sorted(avg_reach)
    N = len(sorted_reach)
    gini = sum((2*(i+1) - N - 1) * sorted_reach[i] for i in range(N))
    gini /= (N * sum(sorted_reach)) if sum(sorted_reach) > 0 else 1
    
    # Also compute "concentration": max_reach / avg_reach
    concentration = [max_reach[i] / max(avg_reach[i], 1) for i in range(N)]
    
    return {
        'gini_reach': gini,
        'avg_reach': sum(avg_reach) / N,
        'max_reach': max(max_reach),
        'avg_concentration': sum(concentration) / N,
    }


def main():
    print("=" * 70)
    print("REFINED GINI-UNIMODALITY ANALYSIS")
    print("=" * 70)
    
    # Test on many brooms
    print("\nBROOMS: gini of average reach vs near-miss")
    results = []
    for p in [3, 5, 10, 15, 20, 30, 50]:
        for s in [10, 50, 100, 200, 500, 1000]:
            n, adj = make_broom(p, s)
            poly = independence_poly(n, adj)
            nm, _ = near_miss_ratio(poly)
            reach = compute_gini_of_reach(n, adj)
            
            results.append({
                'p': p, 's': s, 'n': n,
                'nm': nm,
                'gini': reach['gini_reach'],
                'avg_reach': reach['avg_reach'],
            })
    
    # Sort by gini
    results.sort(key=lambda x: x['gini'])
    
    print("\nSorted by Gini (reach):")
    for r in results[:15]:
        print(f"   p={r['p']:2d}, s={r['s']:4d}, n={r['n']:4d}, "
              f"gini={r['gini']:.3f}, nm={r['nm']:.4f}")
    
    # Look for pattern: do higher Gini values correlate with higher nm?
    high_gini = [r for r in results if r['gini'] > 0.25]
    low_gini = [r for r in results if r['gini'] <= 0.20]
    
    if high_gini and low_gini:
        avg_nm_high = sum(r['nm'] for r in high_gini) / len(high_gini)
        avg_nm_low = sum(r['nm'] for r in low_gini) / len(low_gini)
        print(f"\nHigh Gini (>0.25): avg nm = {avg_nm_high:.4f}")
        print(f"Low Gini (<=0.20): avg nm = {avg_nm_low:.4f}")
    
    # Now try the reverse: does nm correlate with any reach metric?
    print("\n\nTrying: correlation of nm with various metrics")
    for r in results[:20]:
        print(f"   nm={r['nm']:.4f}, avg_reach={r['avg_reach']:.1f}, gini={r['gini']:.3f}")


if __name__ == '__main__':
    main()
