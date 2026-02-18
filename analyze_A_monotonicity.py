#!/usr/bin/env python3
"""
Check if A(x) = x * I(T/e) is non-decreasing up to the mode of I(T).
This condition is equivalent to mode(I(T/e)) >= mode(I(T)) - 1,
assuming I(T/e) is unimodal.

We check this for all trees of size n (default 19).
"""

import sys
import subprocess
from indpoly import independence_poly
from graph6 import parse_graph6

def get_mode_indices(poly):
    m = max(poly)
    return [i for i, x in enumerate(poly) if x == m]

def get_mode_range(poly):
    indices = get_mode_indices(poly)
    return min(indices), max(indices)

def check_A_monotonicity(n, adj):
    poly_T = independence_poly(n, adj)
    min_mode_T, max_mode_T = get_mode_range(poly_T)
    
    # We want A(x) to be non-decreasing up to max_mode_T.
    # A(x) = x * I(T/e).
    # Coeffs of A are [0, c0, c1, ...].
    # So we check 0 <= c0 <= c1 <= ... <= c_{max_mode_T - 1}.
    
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
        poly_Te = independence_poly(new_n, new_adj_list) # coefficients c0, c1, ...
        
        # Check monotonicity of poly_Te up to index max_mode_T - 1
        # If max_mode_T - 1 < 0, nothing to check (mode is 0)
        limit = max_mode_T - 1
        if limit < 0:
            continue
            
        # Ensure we don't go out of bounds of poly_Te
        # poly_Te has size new_alpha + 1.
        # If limit >= len(poly_Te), we treat missing coeffs as 0.
        
        prev = 0
        for k in range(limit + 1):
            if k >= len(poly_Te):
                curr = 0
            else:
                curr = poly_Te[k]
                
            if curr < prev:
                print(f"VIOLATION: n={n} edge={u}-{v}")
                print(f"  Mode T: [{min_mode_T}, {max_mode_T}]")
                print(f"  Poly T/e: {poly_Te}")
                print(f"  Decrease at index {k} (A_{k+1}): {curr} < {prev}")
                return False
            prev = curr
            
    return True

def main():
    n = 19
    if len(sys.argv) > 1:
        n = int(sys.argv[1])
        
    print(f"Checking A(x) monotonicity for n={n}...", flush=True)
    cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
    try:
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    except Exception as e:
        print(f"Error running geng: {e}")
        return

    count = 0
    for line in proc.stdout:
        line = line.strip()
        if not line: continue
        
        tn, adj = parse_graph6(line)
        if not check_A_monotonicity(tn, adj):
            print("Found violation! Stopping.")
            return
        
        count += 1
        if count % 10000 == 0:
            print(f"Checked {count} trees...", flush=True)
            
    print(f"Checked {count} trees. All OK.")

if __name__ == "__main__":
    main()
