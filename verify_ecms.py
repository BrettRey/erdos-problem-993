#!/usr/bin/env python3
"""
Verify the Edge Contraction Mode Stability (ECMS) conjecture.

Conjecture: |mode(I(T)) - mode(I(T/e))| <= 1 for all edges e in T.
"""

import sys
from indpoly import independence_poly
from graph6 import parse_graph6

def get_mode(poly):
    m = max(poly)
    # Return the LAST index of the max to be conservative about shifts?
    # Or the first? The conjecture usually applies to the set of modes.
    # Let's check both min_mode and max_mode.
    indices = [i for i, x in enumerate(poly) if x == m]
    return min(indices), max(indices)

def check_ecms(n, adj):
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
        
        # Check shift
        # The mode of T/e should be within [min_mode_T - 1, max_mode_T + 1]
        # or something similar.
        # Strongest version: |mode_Te - mode_T| <= 1 for any choice of mode.
        
        # Let's check the distance between the intervals
        dist = 0
        if max_mode_Te < min_mode_T - 1:
            dist = (min_mode_T - 1) - max_mode_Te
        elif min_mode_Te > max_mode_T + 1:
            dist = min_mode_Te - (max_mode_T + 1)
            
        if dist > 0:
            print(f"FAIL ECMS: n={n} edge={u}-{v}")
            print(f"  Mode T: [{min_mode_T}, {max_mode_T}]")
            print(f"  Mode T/e: [{min_mode_Te}, {max_mode_Te}]")
            return False
            
    return True

def main():
    count = 0
    for line in sys.stdin.buffer:
        line = line.strip()
        if not line: continue
        n, adj = parse_graph6(line)
        if not check_ecms(n, adj):
            sys.exit(1)
        count += 1
        if count % 1000 == 0:
            print(f"Verified {count} trees...", file=sys.stderr)
    print(f"Success! Verified {count} trees.")

if __name__ == "__main__":
    main()
