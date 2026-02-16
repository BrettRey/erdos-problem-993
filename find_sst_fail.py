"""Find subdivided stars that actually fail LC."""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly, is_log_concave, log_concavity_ratio
from timing_analysis import compute_timing_imbalance

def make_subdivided_star(arms: int, subdivisions: int = 1):
    n = 1 + arms * (subdivisions + 1)
    adj = [[] for _ in range(n)]
    
    center = 0
    for i in range(arms):
        prev = center
        for j in range(subdivisions):
            node = 1 + i * (subdivisions + 1) + j
            adj[prev].append(node)
            adj[node].append(prev)
            prev = node
        leaf = 1 + i * (subdivisions + 1) + subdivisions
        adj[prev].append(leaf)
        adj[leaf].append(prev)
    
    return n, adj

print("Finding LC failures in subdivided stars:")
print("=" * 60)

found_failure = False
for subdivisions in range(1, 15):
    for arms in range(2, 15):
        n, adj = make_subdivided_star(arms, subdivisions)
        if n > 30:
            continue
        
        poly = independence_poly(n, adj)
        is_lc = is_log_concave(poly)
        
        if not is_lc:
            viol, pos = log_concavity_ratio(poly)
            timing = compute_timing_imbalance(n, adj)
            print(f"SST({arms},{subdivisions}) n={n:3d}: LC FAIL! viol={viol:.4f} gini={timing['gini']:.3f}")
            found_failure = True

if not found_failure:
    print("No LC failures found in subdivided stars up to n=30")
    print("The known failures at n=26 must be different subdivided patterns")
