#!/usr/bin/env python3
"""Experiment 3: Does the mean control the mode for tree IS distributions?

Galvin's key question: the paper observes mu < n/3 for d_leaf <= 1 trees,
and wants mode <= floor(n/3) + 1. But the mean only controls the mode if
the distribution has a specific shape (e.g., log-concave distributions
satisfy mode <= ceil(mean)).

Check: for how many trees does mode = ceil(mean)? How far apart can they get?
Is the IS distribution always log-concave (already known: no at n=26)?
For the trees where LC fails, what is mode vs mean?
"""

import subprocess
from indpoly import independence_poly
from graph6 import parse_graph6
import math

def get_mode_and_mean(poly):
    """Return (mode_index, mean_size, is_unimodal)."""
    mode_idx = max(range(len(poly)), key=lambda k: poly[k])
    total = sum(poly)
    mean = sum(k * poly[k] for k in range(len(poly))) / total
    
    # Check unimodality
    is_unimodal = True
    increasing = True
    for k in range(1, len(poly)):
        if poly[k] < poly[k-1]:
            increasing = False
        elif poly[k] > poly[k-1] and not increasing:
            is_unimodal = False
            break
    
    return mode_idx, mean, is_unimodal

def main():
    print("=" * 70)
    print("EXPERIMENT: Mean vs Mode in Tree IS Distributions")
    print("=" * 70)
    print()
    
    print("Question: For trees, is mode ≤ ceil(mean)?")
    print("          If so, then mu < n/3 would imply mode ≤ floor(n/3) + 1.")
    print()
    
    max_gap = 0
    max_gap_tree = None
    mode_vs_mean = {'>': 0, '=': 0, '<': 0}
    mode_exceeds_ceil_mean = 0
    total_trees = 0
    
    for n in range(3, 21):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        n_mode_gt_ceil_mean = 0
        worst_gap_n = 0
        
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            mode_idx, mean, _ = get_mode_and_mean(poly)
            total_trees += 1
            
            ceil_mean = math.ceil(mean)
            gap = mode_idx - ceil_mean
            
            if mode_idx > ceil_mean:
                mode_vs_mean['>'] += 1
                n_mode_gt_ceil_mean += 1
                if gap > worst_gap_n:
                    worst_gap_n = gap
                if gap > max_gap:
                    max_gap = gap
                    max_gap_tree = (n, adj, poly, mode_idx, mean)
            elif mode_idx == ceil_mean:
                mode_vs_mean['='] += 1
            else:
                mode_vs_mean['<'] += 1
        
        if n_mode_gt_ceil_mean > 0:
            print(f"  n={n}: {n_mode_gt_ceil_mean}/{len(lines)} trees have mode > ceil(mean) (worst gap: {worst_gap_n})")
        else:
            print(f"  n={n}: {len(lines)} trees, ALL have mode ≤ ceil(mean) ✓")
    
    print()
    print(f"Summary over {total_trees} trees (n=3..20):")
    print(f"  mode > ceil(mean): {mode_vs_mean['>']}")
    print(f"  mode = ceil(mean): {mode_vs_mean['=']}")
    print(f"  mode < ceil(mean): {mode_vs_mean['<']}")
    print()
    
    if max_gap_tree:
        n, adj, poly, mode_idx, mean = max_gap_tree
        print(f"Worst violation: n={n}, mode={mode_idx}, mean={mean:.4f}, ceil(mean)={math.ceil(mean)}")
        print(f"  Gap = {mode_idx - math.ceil(mean)}")
        print(f"  Poly = {poly[:10]}...")
    else:
        print("NO violations found: mode ≤ ceil(mean) for ALL trees up to n=20.")
        print("⟹ If you can prove mu < n/3, then mode ≤ ceil(n/3) = floor(n/3) + 1.")
        print("   (This would close Galvin's gap!)")
    
    print()
    
    # Also check: is mode <= floor(mean) + 1 (slightly different)?
    print("--- Checking mode ≤ floor(mean) + 1 ---")
    violations = 0
    for n in range(3, 21):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            mode_idx, mean, _ = get_mode_and_mean(poly)
            
            if mode_idx > math.floor(mean) + 1:
                violations += 1
                print(f"  VIOLATION: n={tn}, mode={mode_idx}, mean={mean:.4f}, floor(mean)+1={math.floor(mean)+1}")
    
    if violations == 0:
        print("  No violations: mode ≤ floor(mean) + 1 for all trees up to n=20.")
    
    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
