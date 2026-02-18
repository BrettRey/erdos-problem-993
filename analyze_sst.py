"""Analyze subdivided stars - the pattern that produces LC failures."""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly, is_log_concave, log_concavity_ratio
from timing_analysis import analyze_tree_timing, compute_timing_imbalance
import json

def make_subdivided_star(arms: int, subdivisions: int = 1):
    """
    Subdivided star: center vertex with 'subdivisions' intermediate vertices on each arm.
    """
    n = 1 + arms * (subdivisions + 1)  # center + arms * (subdivisions + leaf)
    adj = [[] for _ in range(n)]
    
    center = 0
    for i in range(arms):
        prev = center
        for j in range(subdivisions):
            node = 1 + i * (subdivisions + 1) + j
            adj[prev].append(node)
            adj[node].append(prev)
            prev = node
        # Add leaf
        leaf = 1 + i * (subdivisions + 1) + subdivisions
        adj[prev].append(leaf)
        adj[leaf].append(prev)
    
    return n, adj

print("Subdivided Stars Analysis:")
print("=" * 60)

for subdivisions in [1, 2, 3, 4, 5]:
    for arms in [3, 4, 5, 6, 8, 10]:
        n, adj = make_subdivided_star(arms, subdivisions)
        poly = independence_poly(n, adj)
        is_lc = is_log_concave(poly)
        viol, pos = log_concavity_ratio(poly)
        timing = compute_timing_imbalance(n, adj)
        
        status = "LC FAIL" if not is_lc else "OK"
        print(f"SST({arms},{subdivisions}) n={n:3d}: {status:8s} viol={viol:.4f} gini={timing['gini']:.3f}")
