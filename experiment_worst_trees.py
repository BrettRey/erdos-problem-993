#!/usr/bin/env python3
"""Identify the worst-case tree (max mode-mean) at each n and study its structure.

Key question: are these paths? stars? caterpillars? something else?
"""

import subprocess
import math
from indpoly import independence_poly
from graph6 import parse_graph6

def tree_signature(adj, n):
    """Return structural fingerprint of a tree."""
    degs = sorted([len(adj[v]) for v in range(n)], reverse=True)
    num_leaves = sum(1 for d in degs if d == 1)
    max_deg = degs[0]
    
    # Is it a path?
    if max_deg == 2 and num_leaves == 2:
        return "path"
    
    # Is it a star?
    if max_deg == n - 1:
        return "star"
    
    # Is it a caterpillar? (removing leaves gives a path)
    internal = [v for v in range(n) if len(adj[v]) > 1]
    if len(internal) > 0:
        # Check if internal vertices form a path
        int_degs = [sum(1 for u in adj[v] if len(adj[u]) > 1) for v in internal]
        int_degs_sorted = sorted(int_degs)
        if all(d <= 2 for d in int_degs):
            # Count endpoints (degree 1 in internal subgraph)
            endpoints = sum(1 for d in int_degs if d <= 1)
            if endpoints == 2 or len(internal) == 1:
                return f"caterpillar(spine={len(internal)})"
    
    # Degree sequence
    return f"deg_seq={degs[:5]}..."

def main():
    print("=" * 70)
    print("  WORST-CASE TREE IDENTIFICATION")
    print("=" * 70)
    print()
    
    for n in range(5, 23):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        max_excess = -float('inf')
        best = None
        
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            total = sum(poly)
            mean = sum(k * poly[k] for k in range(len(poly))) / total
            mode = max(range(len(poly)), key=lambda k: poly[k])
            excess = mode - mean
            
            if excess > max_excess:
                max_excess = excess
                degs = sorted([len(adj[v]) for v in range(tn)], reverse=True)
                sig = tree_signature(adj, tn)
                best = (excess, mode, mean, poly, sig, degs)
        
        excess, mode, mean, poly, sig, degs = best
        
        # Look at the polynomial near the mode
        mode_k = mode
        if mode_k > 0 and mode_k < len(poly) - 1:
            left = poly[mode_k - 1]  
            center = poly[mode_k]
            right = poly[mode_k + 1]
            ratio_left = left / center
            ratio_right = right / center
        else:
            ratio_left = ratio_right = 0
        
        print(f"n={n:2d}: excess={excess:.4f}, type={sig}")
        print(f"       degs={degs}")
        print(f"       poly={poly}")
        print(f"       mode={mode}, a[m-1]/a[m]={ratio_left:.4f}, a[m+1]/a[m]={ratio_right:.4f}")
        
        # Near-palindromic check around mode
        alpha = len(poly) - 1
        diffs = []
        for j in range(1, min(mode, alpha - mode) + 1):
            if mode - j >= 0 and mode + j < len(poly):
                diffs.append(poly[mode + j] - poly[mode - j])
        print(f"       palindrome diffs from mode: {diffs[:5]}")
        print()

if __name__ == "__main__":
    main()
