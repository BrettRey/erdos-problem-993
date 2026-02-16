#!/usr/bin/env python3
"""
New approach: What coefficient positions can EACH VERTEX actually contribute to?

For each vertex v, we compute the SET of k values such that v can be 
in an independent set of size k. This is the "support" of v's influence.

Then we measure: how much do these supports overlap or separate?
This might capture the "timing" better than distance metrics.
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly, near_miss_ratio
from targeted import make_broom
from collections import defaultdict


def compute_vertex_contributions(n, adj):
    """
    For each vertex v, compute the set of k values where v can appear
    in an independent set of size k.
    
    This is expensive - we do it by DP.
    """
    # Root at 0
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
    
    # DP: for each vertex, compute sets of possible sizes
    # When included, when excluded
    inc_sets = [set() for _ in range(n)]  # sizes when v is included
    exc_sets = [set() for _ in range(n)]  # sizes when v is excluded
    
    post_order = queue[::-1]
    
    for v in post_order:
        if not children[v]:
            inc_sets[v] = {1}  # can be in IS of size 1
            exc_sets[v] = {0}   # excluded = empty set contributes 0
        else:
            # exc_sets[v]: combine children (each can be in or out)
            # This is the set of all sums from children
            exc_sets[v] = {0}
            for c in children[v]:
                new_set = set()
                for s1 in exc_sets[v]:
                    for s2 in exc_sets[c]:
                        new_set.add(s1 + s2)
                    for s2 in inc_sets[c]:
                        new_set.add(s1 + s2)
                exc_sets[v] = new_set
            
            # inc_sets[v]: include v, children must be excluded
            inc_sets[v] = {1}
            for c in children[v]:
                new_set = set()
                for s1 in inc_sets[v]:
                    for s2 in exc_sets[c]:
                        new_set.add(s1 + s2)
                inc_sets[v] = new_set
    
    # Full influence set for each vertex = union of when included or excluded
    full_influence = [inc_sets[v] | exc_sets[v] for v in range(n)]
    
    return full_influence


def compute_support_metrics(n, adj):
    """Compute metrics about vertex support overlap/separation."""
    influence = compute_vertex_contributions(n, adj)
    
    # For each k, how many vertices can contribute?
    k_coverage = defaultdict(int)
    for v in range(n):
        for k in influence[v]:
            k_coverage[k] += 1
    
    # Metrics:
    # 1. Support size per vertex
    support_sizes = [len(influence[v]) for v in range(n)]
    avg_support = sum(support_sizes) / n
    
    # 2. Variance in support position
    # For each vertex, what's the average k it can contribute to?
    avg_k_per_vertex = []
    for v in range(n):
        if influence[v]:
            avg_k_per_vertex.append(sum(influence[v]) / len(influence[v]))
        else:
            avg_k_per_vertex.append(0)
    
    # Gini of avg_k
    sorted_avg = sorted(avg_k_per_vertex)
    N = len(sorted_avg)
    gini = sum((2*(i+1) - N - 1) * sorted_avg[i] for i in range(N))
    gini /= (N * sum(sorted_avg)) if sum(sorted_avg) > 0 else 1
    
    # 3. Support span per vertex
    spans = [max(influence[v]) - min(influence[v]) if influence[v] else 0 for v in range(n)]
    avg_span = sum(spans) / n
    
    return {
        'gini_avg_k': gini,
        'avg_support': avg_support,
        'avg_span': avg_span,
        'max_span': max(spans),
    }


def main():
    print("=" * 70)
    print("SUPPORT-BASED GINI ANALYSIS")
    print("=" * 70)
    
    results = []
    for p in [3, 5, 10, 15, 20, 30]:
        for s in [10, 50, 100, 200, 500]:
            n, adj = make_broom(p, s)
            poly = independence_poly(n, adj)
            nm, _ = near_miss_ratio(poly)
            metrics = compute_support_metrics(n, adj)
            
            results.append({
                'p': p, 's': s, 'n': n,
                'nm': nm,
                'gini_avg_k': metrics['gini_avg_k'],
                'avg_span': metrics['avg_span'],
            })
    
    # Sort by gini
    results.sort(key=lambda x: x['gini_avg_k'])
    
    print("\nSorted by Gini of average k:")
    for r in results[:20]:
        print(f"   p={r['p']:2d}, s={r['s']:4d}, n={r['n']:4d}, "
              f"gini={r['gini_avg_k']:.3f}, span={r['avg_span']:.1f}, nm={r['nm']:.4f}")
    
    # Correlation?
    print("\n\nLooking for pattern:")
    print("nm vs gini_avg_k:")
    for r in results:
        print(f"   nm={r['nm']:.4f}, gini={r['gini_avg_k']:.3f}")


if __name__ == '__main__':
    main()
