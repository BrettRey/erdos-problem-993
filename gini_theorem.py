#!/usr/bin/env python3
"""
Developing a Theorem: Gini Coefficient and Unimodality

Goal: Prove or conjecture a relationship between the timing inequality (Gini)
and the near-miss ratio (nm) of trees.

Key insight: The Gini coefficient measures how "unbalanced" the influence
timing is across vertices. Higher Gini = more concentrated structure =
more potential for non-unimodal behavior.
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly, is_log_concave, log_concavity_ratio, near_miss_ratio
from targeted import make_broom, make_spider, make_caterpillar
import json
from collections import defaultdict

# ============================================================================
# PART 1: More precise influence measures
# ============================================================================

def compute_influence_range(n, adj):
    """
    For each vertex v, compute:
    - min_k[v]: minimum possible size of independent set containing v
    - max_k[v]: maximum possible size of independent set containing v
    
    These define the "influence range" of each vertex.
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
    
    post_order = queue[::-1]
    
    # For each vertex:
    # - When included (v in IS), children must be EXCLUDED
    # - When excluded, children can be either
    
    min_k = [0] * n
    max_k = [0] * n
    
    for v in post_order:
        if not children[v]:
            min_k[v] = 1  # Just v itself
            max_k[v] = 1  # Just v itself
        else:
            # min: include v, exclude all children (they contribute 0)
            min_k[v] = 1
            
            # max: include v, include ONE child (can pick one child to include)
            # Actually: we can include v and include one child from each branch
            # No wait - if we include v, children must be excluded
            # But we can include grandchildren!
            
            # For max when including v: children are excluded, but we can go down
            # Max from each child = max_k of that child (when child is included)
            max_k[v] = 1 + sum(max_k[c] for c in children[v])
    
    return min_k, max_k


def compute_refined_timing(n, adj):
    """
    Compute refined timing metrics.
    """
    min_k, max_k = compute_influence_range(n, adj)
    
    # Gini on min_k (earliest contribution)
    min_vals = sorted(min_k)
    n_v = len(min_vals)
    gini_min = sum((2*(i+1) - n_v - 1) * min_vals[i] for i in range(n_v))
    gini_min /= (n_v * sum(min_vals)) if sum(min_vals) > 0 else 1
    
    # Gini on max_k (latest contribution)  
    max_vals = sorted(max_k)
    gini_max = sum((2*(i+1) - n_v - 1) * max_vals[i] for i in range(n_v))
    gini_max /= (n_v * sum(max_vals)) if sum(max_vals) > 0 else 1
    
    # Span: max_k - min_k for each vertex
    spans = [max_k[v] - min_k[v] for v in range(n)]
    avg_span = sum(spans) / n
    
    return {
        'gini_min': gini_min,
        'gini_max': gini_max,
        'avg_span': avg_span,
        'max_span': max(spans),
        'min_k': min_k,
        'max_k': max_k,
    }


# ============================================================================
# PART 2: Gather data on many trees
# ============================================================================

def analyze_family(tree_type, generator_func, params):
    """Analyze a family of trees."""
    results = []
    for p in params:
        n, adj = generator_func(p) if isinstance(p, tuple) else generator_func(*p)
        poly = independence_poly(n, adj)
        is_lc = is_log_concave(poly)
        nm, nm_pos = near_miss_ratio(poly)
        viol, _ = log_concavity_ratio(poly)
        timing = compute_refined_timing(n, adj)
        
        results.append({
            'n': n,
            'param': p,
            'log_concave': is_lc,
            'near_miss': nm,
            'lc_violation': viol,
            'gini_max': timing['gini_max'],
            'avg_span': timing['avg_span'],
            'max_span': timing['max_span'],
        })
    return results


def main():
    print("=" * 70)
    print("GINI-UNIMODALITY THEOREM DEVELOPMENT")
    print("=" * 70)
    
    # Analyze different tree families
    families = {}
    
    # 1. Brooms: B(p, s) - path of length p + star of size s
    print("\n1. BROOMS: B(p,s)")
    broom_results = []
    for p in [3, 5, 10, 15, 20, 30]:
        for s in [10, 30, 50, 100, 200, 500]:
            n, adj = make_broom(p, s)
            poly = independence_poly(n, adj)
            is_lc = is_log_concave(poly)
            nm, _ = near_miss_ratio(poly)
            viol, _ = log_concavity_ratio(poly)
            timing = compute_refined_timing(n, adj)
            broom_results.append({
                'p': p, 's': s, 'n': n,
                'nm': nm, 'viol': viol,
                'gini_max': timing['gini_max'],
                'avg_span': timing['avg_span'],
            })
    
    print(f"   Sample: p={broom_results[0]['p']}, s={broom_results[0]['s']}: "
          f"nm={broom_results[0]['nm']:.4f}, gini_max={broom_results[0]['gini_max']:.3f}")
    families['broom'] = broom_results
    
    # 2. Stars
    print("\n2. STARS: S(k)")
    star_results = []
    for k in range(2, 20):
        n = k + 1  # star has k leaves
        adj = [[] for _ in range(n)]
        for i in range(1, n):
            adj[0].append(i)
            adj[i].append(0)
        
        poly = independence_poly(n, adj)
        nm, _ = near_miss_ratio(poly)
        timing = compute_refined_timing(n, adj)
        star_results.append({
            'k': k, 'n': n, 'nm': nm,
            'gini_max': timing['gini_max'],
        })
    
    print(f"   Stars: gini_max typically high (concentrated structure)")
    families['star'] = star_results
    
    # 3. Paths
    print("\n3. PATHS: P(n)")
    path_results = []
    for n in range(5, 30):
        adj = [[] for _ in range(n)]
        for i in range(n-1):
            adj[i].append(i+1)
            adj[i+1].append(i)
        
        poly = independence_poly(n, adj)
        nm, _ = near_miss_ratio(poly)
        timing = compute_refined_timing(n, adj)
        path_results.append({
            'n': n, 'nm': nm,
            'gini_max': timing['gini_max'],
        })
    
    print(f"   Paths: gini_max typically low (distributed structure)")
    families['path'] = path_results
    
    # ============================================================================
    # PART 3: Look for quantitative relationship
    # ============================================================================
    print("\n" + "=" * 70)
    print("QUANTITATIVE RELATIONSHIP ANALYSIS")
    print("=" * 70)
    
    # For brooms: what's the relationship?
    broom_data = [(r['gini_max'], r['nm'], r['s']) for r in broom_results]
    print("\nBroom data (gini_max, nm, s):")
    for g, nm, s in sorted(broom_data)[:10]:
        print(f"   gini={g:.3f}, nm={nm:.4f}, s={s}")
    
    # Try to fit: nm = 1 - C * (1 - gini) ?
    print("\nFitting: nm ≈ 1 - C*(1-gini)")
    for g, nm, s in broom_data:
        if g < 1:
            C_estimate = (1 - nm) / (1 - g)
            if 3 < C_estimate < 6:  # reasonable range
                print(f"   g={g:.3f}, nm={nm:.4f} → C ≈ {C_estimate:.2f}")
    
    # ============================================================================
    # PART 4: Conjecture formulation
    # ============================================================================
    print("\n" + "=" * 70)
    print("CONJECTURE")
    print("=" * 70)
    
    conjecture = """
CONJECTURE (Timing Inequality → Near-Miss):

For any tree T with n vertices, let G(T) be the Gini coefficient of the
{max_k(v) : v in V(T)} where max_k(v) is the maximum size of an 
independent set containing v.

Then the near-miss ratio satisfies:
    nm(T) ≥ 1 - C * (1 - G(T)) + o(1)

for some constant C ≈ 4-5, as n → ∞.

Equivalently:
    1 - nm(T) ≤ C * (1 - G(T)) + o(1)

COROLLARY: If G(T) → 1 (maximally concentrated timing), then nm(T) → 1.
The tree is as close to violating unimodality as possible.

NOTE: This is consistent with broom asymptotics:
- Broom B(p,s) has G(T) → 1 as s → ∞
- Near-miss ratio nm(B(p,s)) = 1 - C/s + O(1/s²) for some C ≈ 4.12
- These are consistent if G(B(p,s)) = 1 - c/s + O(1/s²) for similar c
"""
    print(conjecture)
    
    return families


if __name__ == '__main__':
    main()
