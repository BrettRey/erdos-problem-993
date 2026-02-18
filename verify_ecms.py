#!/usr/bin/env python3
"""
Verify the Edge Contraction Mode Stability (ECMS) conjecture.
And collect statistics on mode shifts.

Conjecture: |mode(I(T)) - mode(I(T/e))| <= 1 for all edges e in T.
"""

import sys
from collections import Counter
from indpoly import independence_poly
from graph6 import parse_graph6

def get_mode(poly):
    m = max(poly)
    indices = [i for i, x in enumerate(poly) if x == m]
    return min(indices), max(indices)

def check_ecms(n, adj, stats):
    poly_T = independence_poly(n, adj)
    min_mode_T, max_mode_T = get_mode(poly_T)
    
    # Iterate edges
    edges = []
    for u in range(n):
        for v in adj[u]:
            if u < v:
                edges.append((u, v))
                
    for u, v in edges:
        # Contract uv
        mapping = {}
        idx = 1
        for i in range(n):
            if i == u or i == v:
                mapping[i] = 0
            else:
                mapping[i] = idx
                idx += 1
                
        new_n = n - 1
        new_adj = [set() for _ in range(new_n)]
        
        for i in range(n):
            src = mapping[i]
            for j in adj[i]:
                dst = mapping[j]
                if src != dst:
                    new_adj[src].add(dst)
                    new_adj[dst].add(src)
        
        new_adj_list = [list(s) for s in new_adj]
        poly_Te = independence_poly(new_n, new_adj_list)
        
        min_mode_Te, max_mode_Te = get_mode(poly_Te)
        
        # Calculate shift
        # We compare the intervals.
        # If intervals overlap, shift is 0.
        # If Te is strictly left, shift is negative.
        # If Te is strictly right, shift is positive.
        
        shift = 0
        if max_mode_Te < min_mode_T:
            shift = max_mode_Te - min_mode_T # Negative
        elif min_mode_Te > max_mode_T:
            shift = min_mode_Te - max_mode_T # Positive
            
        stats[shift] += 1
            
        if abs(shift) > 1:
            print(f"FAIL ECMS: n={n} edge={u}-{v} shift={shift}")
            print(f"  Mode T: [{min_mode_T}, {max_mode_T}]")
            print(f"  Mode T/e: [{min_mode_Te}, {max_mode_Te}]")
            return False
            
    return True

def main():
    count = 0
    stats = Counter()
    for line in sys.stdin.buffer:
        line = line.strip()
        if not line: continue
        n, adj = parse_graph6(line)
        if not check_ecms(n, adj, stats):
            sys.exit(1)
        count += 1
        if count % 1000 == 0:
            print(f"Verified {count} trees...", file=sys.stderr)
            
    print(f"Success! Verified {count} trees.")
    print("Mode shift distribution (T/e - T):")
    total_edges = sum(stats.values())
    for shift in sorted(stats.keys()):
        pct = 100 * stats[shift] / total_edges if total_edges > 0 else 0
        print(f"  Shift {shift:+d}: {stats[shift]:10d} ({pct:6.2f}%)")

if __name__ == "__main__":
    main()
