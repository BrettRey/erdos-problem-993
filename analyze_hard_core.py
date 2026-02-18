#!/usr/bin/env python3
"""Analyze hard-core occupation probabilities P(v) for trees with d_leaf <= 1.

Calculates:
1. Exact P(v) = Number of IS containing v / Total IS.
2. Identifies "High-P" core vertices (P > 1/3).
3. Checks if High-P vertices form a "Cluster" with their neighbors that satisfies Mean(Cluster) <= n_cluster/3.
"""

import subprocess
import sys
from collections import defaultdict

# Use local modules
sys.path.insert(0, ".")
from indpoly import independence_poly
from graph6 import parse_graph6

def get_occupation_probs(n, adj):
    poly_T = independence_poly(n, adj)
    total_IS = sum(poly_T)
    probs = []
    # Optimization: If total_IS is large, use recursive formula?
    # For n <= 20, exact calculation is fast enough.
    
    for v in range(n):
        neighbors = set(adj[v])
        neighbors.add(v) 
        remaining = [u for u in range(n) if u not in neighbors]
        
        if not remaining:
            count_inc = 1 
        else:
            mapping = {old: new for new, old in enumerate(remaining)}
            new_adj = [[] for _ in range(len(remaining))]
            for old_u in remaining:
                for old_w in adj[old_u]:
                    if old_w in mapping: 
                         new_adj[mapping[old_u]].append(mapping[old_w])
            
            poly_sub = independence_poly(len(remaining), new_adj)
            count_inc = sum(poly_sub)
        probs.append(count_inc / total_IS)
    return probs, total_IS

def analyze_cluster(n, adj, probs):
    # Identify Violators
    violators = [v for v in range(n) if probs[v] > 1/3 + 1e-9]
    if not violators:
        return True, 0.0, None

    max_p = max(probs[v] for v in violators)
    
    # Check if every violator has a "safe neighborhood"
    # Strategy: Form a cluster C_v = {v} U {neighbors of v}.
    # Does Sum(P(u) for u in C_v) <= |C_v|/3 ?
    # Note: This is a heuristic. Overlapping clusters make simple summing invalid globally.
    # But if every local cluster is "light", it suggests stability.
    
    worst_cluster_excess = -1.0
    
    for v in violators:
        cluster = [v] + adj[v]
        cluster_load = sum(probs[u] for u in cluster)
        cluster_cap = len(cluster) / 3.0
        
        excess = cluster_load - cluster_cap
        if excess > worst_cluster_excess:
            worst_cluster_excess = excess
            
    if worst_cluster_excess > 1e-9:
        return False, max_p, worst_cluster_excess
        
    return True, max_p, worst_cluster_excess

def main():
    print("ANALYSIS: Cluster Load Verification (d_leaf <= 1)", flush=True)
    print("Checking: For every v with P(v) > 1/3, is Sum(N[v]) <= |N[v]|/3 ?", flush=True)
    
    global_max_p = 0.0
    worst_excess_global = -999.0
    
    for n in range(3, 19): 
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True) 
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        n_high_p = 0
        n_cluster_fail = 0
        
        for line in lines:
            tn, adj = parse_graph6(line)
            
            # Filter d_leaf <= 1
            leaves = [v for v in range(tn) if len(adj[v]) == 1]
            # Must check if any node connects to >= 2 leaves
            is_valid = True
            for v in range(tn):
                if v in leaves: continue
                leaf_neighbors = sum(1 for u in adj[v] if u in leaves)
                if leaf_neighbors > 1:
                    is_valid = False
                    break
            if not is_valid: continue

            probs, total_is = get_occupation_probs(tn, adj)
            ok, max_p, excess = analyze_cluster(tn, adj, probs)
            
            if max_p > global_max_p:
                global_max_p = max_p
                
            if not ok:
                n_cluster_fail += 1
                if excess > worst_excess_global:
                    worst_excess_global = excess
                    print(f"CLUSTER FAIL n={tn}: MaxP={max_p:.3f}, Excess={excess:.3f}", flush=True)
                    print(f"Adj: {adj}", flush=True)
                    print(f"Probs: {['%.3f'%p for p in probs]}", flush=True)
            elif max_p > 1/3:
                n_high_p += 1
                
        if n_cluster_fail == 0 and n_high_p > 0:
            print(f"n={n:2d}: {n_high_p} trees with high-P vertices, all cleared by cluster check.", flush=True)
        elif n_cluster_fail == 0:
            print(f"n={n:2d}: No high-P vertices.", flush=True)

    print("\nFINAL RESULTS", flush=True)
    print(f"Global Max P = {global_max_p:.5f}", flush=True)
    print(f"Worst Cluster Excess = {worst_excess_global:.5f}", flush=True)

if __name__ == "__main__":
    main()
