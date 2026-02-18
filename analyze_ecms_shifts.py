#!/usr/bin/env python3
"""
Analyze ECMS positive shifts.
Find trees T and edges e where mode(I(T/e)) > mode(I(T)).
These are the "dangerous" cases where contraction shifts the mode to the right.
We want to see if shift >= 2 is possible.
"""

import sys
from indpoly import independence_poly
from graph6 import parse_graph6

def get_mode(poly):
    m = max(poly)
    indices = [i for i, x in enumerate(poly) if x == m]
    return min(indices), max(indices)

def analyze_shifts(n, adj):
    poly_T = independence_poly(n, adj)
    min_mode_T, max_mode_T = get_mode(poly_T)
    
    results = []
    
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
        
        # Check for positive shift
        # We define positive shift if the interval moves strictly to the right
        if min_mode_Te > max_mode_T:
            shift = min_mode_Te - max_mode_T
            results.append({
                'edge': (u, v),
                'mode_T': (min_mode_T, max_mode_T),
                'mode_Te': (min_mode_Te, max_mode_Te),
                'shift': shift,
                'poly_T': poly_T,
                'poly_Te': poly_Te
            })
            
    return results

def main():
    count = 0
    positive_shifts = 0
    
    print(f"Scanning for positive mode shifts (T/e > T)...")
    
    for line in sys.stdin.buffer:
        line = line.strip()
        if not line: continue
        n, adj = parse_graph6(line)
        
        shifts = analyze_shifts(n, adj)
        if shifts:
            for res in shifts:
                positive_shifts += 1
                print(f"FOUND +{res['shift']} SHIFT: n={n} edge={res['edge']}")
                print(f"  T  modes: {res['mode_T']}  poly: {res['poly_T']}")
                print(f"  Te modes: {res['mode_Te']}  poly: {res['poly_Te']}")
                print("-" * 40)
                
        count += 1
        
    print(f"Scanned {count} trees. Found {positive_shifts} positive shifts.")

if __name__ == "__main__":
    main()
