#!/usr/bin/env python3
"""
Simplified Influence Timing Analysis

For each vertex v, compute:
- min_coefficient[v]: the minimum k such that v can be in an independent set of size k
- max_coefficient[v]: the maximum k such that v can be in an independent set of size k
- influence_span[v]: max - min

Then analyze: do LC failures have different timing patterns than LC passes?
"""

import json
import random
import sys
from collections import defaultdict
from typing import Dict, List, Set

sys.path.insert(0, '.')

from indpoly import independence_poly, is_log_concave, log_concavity_ratio


def compute_min_max_influence(n: int, adj: List[List[int]]) -> Dict[int, Dict]:
    """
    For each vertex v, compute the min and max coefficient positions
    that v can influence.
    
    Approach: Root the tree, then for each vertex compute:
    - min_k: smallest independent set size including v
    - max_k: largest independent set size including v
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
    
    # Post-order processing
    post_order = queue[::-1]
    
    # For each vertex v, we need:
    # - min_k when v is included: 1 + sum of min_k from excluded children
    # - max_k when v is included: 1 + sum of max_k from excluded children
    
    min_k = [0] * n  # min independent set size including v
    max_k = [0] * n  # max independent set size including v
    
    for v in post_order:
        if not children[v]:
            min_k[v] = 1
            max_k[v] = 1
        else:
            # When v is included, children must be excluded
            # min: sum of min_k from each child's subtree (they're excluded)
            # max: sum of max_k from each child's subtree (they're excluded)
            min_sum = 0
            max_sum = 0
            for c in children[v]:
                min_sum += min_k[c]
                max_sum += max_k[c]
            min_k[v] = 1 + min_sum
            max_k[v] = 1 + max_sum
    
    return {
        v: {'min': min_k[v], 'max': max_k[v], 'span': max_k[v] - min_k[v]}
        for v in range(n)
    }


def compute_timing_imbalance(n: int, adj: List[List[int]]) -> Dict:
    """
    Compute various timing imbalance metrics.
    """
    influence = compute_min_max_influence(n, adj)
    
    # Distribution of min_k values
    min_vals = [influence[v]['min'] for v in range(n)]
    max_vals = [influence[v]['max'] for v in range(n)]
    span_vals = [influence[v]['span'] for v in range(n)]
    
    # How spread out are the min values?
    # High spread means some vertices always affect early coefficients,
    # others always affect late coefficients
    min_spread = max(min_vals) - min(min_vals)
    max_spread = max(max_vals) - min(max_vals)
    
    # Variance in spans
    avg_span = sum(span_vals) / len(span_vals)
    span_variance = sum((s - avg_span)**2 for s in span_vals) / len(span_vals)
    
    # Gini coefficient of min values (measure of inequality)
    sorted_mins = sorted(min_vals)
    n_vertices = len(sorted_mins)
    gini = sum((2*(i+1) - n_vertices - 1) * sorted_mins[i] for i in range(n_vertices))
    gini /= (n_vertices * sum(sorted_mins)) if sum(sorted_mins) > 0 else 1
    
    return {
        'min_spread': min_spread,
        'max_spread': max_spread,
        'avg_span': avg_span,
        'span_variance': span_variance,
        'gini': gini,
        'n': n,
    }


def analyze_tree_timing(n: int, adj: List[List[int]], label: str = "") -> Dict:
    """Analyze timing patterns for a tree."""
    poly = independence_poly(n, adj)
    is_lc = is_log_concave(poly)
    lc_viol, lc_viol_pos = log_concavity_ratio(poly)
    
    timing = compute_timing_imbalance(n, adj)
    influence = compute_min_max_influence(n, adj)
    
    # Find vertices at extremes
    early_vertices = [v for v in range(n) if influence[v]['min'] <= 2]
    late_vertices = [v for v in range(n) if influence[v]['max'] >= n // 2]
    
    return {
        'label': label,
        'n': n,
        'log_concave': is_lc,
        'lc_violation': lc_viol,
        'lc_viol_pos': lc_viol_pos,
        **timing,
        'num_early': len(early_vertices),
        'num_late': len(late_vertices),
    }


def run_survey(min_n: int = 10, max_n: int = 30, samples_per_n: int = 50):
    """Survey random trees for timing patterns."""
    import networkx as nx
    
    results = []
    rng = random.Random(12345)
    
    for n in range(min_n, max_n + 1):
        for _ in range(samples_per_n):
            G = nx.random_labeled_tree(n, seed=rng.randint(0, 100000))
            adj = [list(G.neighbors(v)) for v in range(n)]
            
            result = analyze_tree_timing(n, adj, f"random_{n}")
            results.append(result)
    
    # Compare LC failures vs passes
    failures = [r for r in results if not r['log_concave']]
    passes = [r for r in results if r['log_concave']]
    
    print(f"\nSurvey: n={min_n}-{max_n}, {samples_per_n} samples each")
    print(f"Total: {len(results)}, LC failures: {len(failures)}, passes: {len(passes)}")
    
    if failures:
        print("\nLC Failure patterns:")
        for key in ['min_spread', 'max_spread', 'avg_span', 'span_variance', 'gini']:
            fail_vals = [f[key] for f in failures]
            pass_vals = [p[key] for p in passes]
            print(f"  {key}: fail avg={sum(fail_vals)/len(fail_vals):.3f}, pass avg={sum(pass_vals)/len(pass_vals):.3f}")
    
    return results


def analyze_brooms():
    """Analyze timing patterns for brooms (known near-misses)."""
    from targeted import make_broom
    
    print("\nBroom Analysis:")
    results = []
    for p in [5, 10, 15, 20, 30]:
        for s in [10, 50, 100, 200]:
            n, adj = make_broom(p, s)
            result = analyze_tree_timing(n, adj, f"broom({p},{s})")
            results.append(result)
    
    # Show patterns
    for r in results[:10]:
        print(f"  {r['label']}: n={r['n']}, gini={r['gini']:.3f}, min_spread={r['min_spread']}, lc={r['log_concave']}, viol={r['lc_violation']:.4f}")
    
    return results


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Timing imbalance analysis')
    parser.add_argument('--survey', action='store_true', help='Run survey')
    parser.add_argument('--brooms', action='store_true', help='Analyze brooms')
    
    args = parser.parse_args()
    
    if args.survey:
        run_survey()
    elif args.brooms:
        analyze_brooms()
    else:
        # Default: quick test
        import networkx as nx
        G = nx.balanced_tree(2, 3)
        n = len(G.nodes())
        adj = [list(G.neighbors(v)) for v in range(n)]
        result = analyze_tree_timing(n, adj)
        print(json.dumps(result, indent=2))
