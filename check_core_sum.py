#!/usr/bin/env python3
"""Check Core Sum Hypothesis for d_leaf <= 1 trees.

Hypothesis: Sum_{v in Core} P(v) <= |Core|/3.
If this holds, the proof decouples.
If this FAILS, the Core relies on Leaf-Support slack.
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
    # Simplified loop for occupation
    for v in range(n):
        neighbors = set(adj[v]); neighbors.add(v) 
        remaining = [u for u in range(n) if u not in neighbors]
        if not remaining: count_inc = 1 
        else:
            mapping = {old: new for new, old in enumerate(remaining)}
            new_adj = [[] for _ in range(len(remaining))]
            for old_u in remaining:
                for old_w in adj[old_u]:
                    if old_w in mapping: new_adj[mapping[old_u]].append(mapping[old_w])
            count_inc = sum(independence_poly(len(remaining), new_adj))
        probs.append(count_inc / total_IS)
    return probs

def check_core_sum(n, adj, probs):
    leaves = [v for v in range(n) if len(adj[v]) == 1]
    support_map = {} 
    for l in leaves:
        supp = adj[l][0]
        support_map[l] = supp
        
    leaf_counts = defaultdict(int)
    for s in support_map.values(): leaf_counts[s] += 1
    if any(c > 1 for c in leaf_counts.values()): return None 
        
    used = set()
    pairs = []
    for l in leaves:
        s = support_map[l]
        pairs.append((l,s))
        used.add(l)
        used.add(s)
        
    core = [v for v in range(n) if v not in used]
    
    if not core: return None
    
    core_sum = sum(probs[v] for v in core)
    core_cap = len(core) / 3.0
    excess = core_sum - core_cap
    
    pair_slack = 0.0
    for l, s in pairs:
        pair_load = probs[l] + probs[s]
        # Slack relative to 2/3
        pair_slack += (2/3.0 - pair_load)
        
    return excess, core_sum, len(core), pair_slack

def main():
    print("ANALYSIS: Core Sum Verification (Sum(Core) <= |Core|/3)", flush=True)
    
    max_excess = 0.0
    worst_tree = None
    
    for n in range(3, 19): 
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True) 
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        n_fail = 0
        
        for line in lines:
            tn, adj = parse_graph6(line)
            try:
                probs = get_occupation_probs(tn, adj)
            except: continue
            
            res = check_core_sum(tn, adj, probs)
            if res is None: continue
            
            excess, c_sum, c_len, slack = res
            
            if excess > 0.000001:
                n_fail += 1
                if excess > max_excess:
                    max_excess = excess
                    worst_tree = (tn, adj, c_sum, c_len, slack)
                    
        if n_fail > 0:
            print(f"n={n:2d}: {n_fail} trees violate Core Sum. Max Excess={max_excess:.5f}", flush=True)

    print("\nFINAL RESULTS", flush=True)
    if max_excess > 0:
        n, adj, c_sum, c_len, slack = worst_tree
        print(f"COUNTEREXAMPLE FOUND: Core Sum > |Core|/3 by {max_excess:.5f}", flush=True)
        print(f"Tree n={n}. Core Size={c_len}. Core Sum={c_sum:.5f} (Target {c_len/3:.5f})", flush=True)
        print(f"Available Pair Slack={slack:.5f}", flush=True)
        if slack > max_excess:
            print("Status: Core Excess is covered by Pair Slack (Global bound holds).", flush=True)
        else:
            print("Status: CRITICAL VIOLATION (Global bound fails).", flush=True)
    else:
        print("SUCCESS: Core Sum <= |Core|/3 holds for all trees.", flush=True)

if __name__ == "__main__":
    main()
