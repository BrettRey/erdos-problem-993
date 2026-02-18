#!/usr/bin/env python3
"""
Verify the analytic bounds for the Subdivision Lemma proof.

Bound 1: i_k(T/e) <= i_k(T) for all k.
Bound 2: I(T/e) is non-increasing for k >= mode(I(T)) + 1.

If these hold, the Subdivision Lemma follows analytically (assuming T/e is unimodal).
"""

import sys
import argparse
from indpoly import independence_poly
from graph6 import parse_graph6

def get_mode(poly):
    # Find the LAST maximal element to be safe, or check for strict unimodality
    m = max(poly)
    # Return the first index of the max
    return poly.index(m)

def check_tree(n, adj):
    # 1. Compute I(T)
    poly_T = independence_poly(n, adj)
    mode_T = get_mode(poly_T)
    
    # 2. Iterate over all edges
    edges = []
    for u in range(n):
        for v in adj[u]:
            if u < v:
                edges.append((u, v))
                
    for u, v in edges:
        # Contract edge uv -> T/e
        # We can implement contraction or compute poly directly
        # Contraction: merge u,v. 
        # But wait, indpoly.py doesn't have a public contract function?
        # Let's implement a simple one for trees.
        
        # New adjacency for T/e
        # Map u, v to 0. Map others to 1..n-2.
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
                if src != dst: # ignore self loops (the contracted edge)
                    new_adj[src].add(dst)
                    new_adj[dst].add(src)
                    
        # Convert back to list
        new_adj_list = [list(s) for s in new_adj]
        
        poly_Te = independence_poly(new_n, new_adj_list)
        
        # Check Bound 1: i_k(T/e) <= i_k(T)
        # Pad with zeros if lengths differ
        len_T = len(poly_T)
        len_Te = len(poly_Te)
        max_len = max(len_T, len_Te)
        
        pT = poly_T + [0]*(max_len - len_T)
        pTe = poly_Te + [0]*(max_len - len_Te)
        
        for k in range(max_len):
            if pTe[k] > pT[k]:
                print(f"FAIL Bound 1: n={n} edge={u}-{v} k={k} Te={pTe[k]} T={pT[k]}")
                return False
                
        # Check Bound 2: I(T/e) decreasing in tail of I(T)
        # Tail starts at mode_T + 1
        start_k = mode_T + 1
        
        # We need I(T/e) to be non-increasing from start_k onwards
        # i.e., pTe[k] >= pTe[k+1] for all k >= start_k
        
        # Note: if start_k >= len_Te, it's vacuously true (all 0)
        
        current_val = pTe[start_k] if start_k < len_Te else 0
        
        for k in range(start_k, max_len - 1):
            next_val = pTe[k+1]
            if next_val > current_val:
                print(f"FAIL Bound 2: n={n} edge={u}-{v} mode(T)={mode_T} k={k} next={next_val} curr={current_val}")
                print(f"  Poly T: {poly_T}")
                print(f"  Poly Te: {poly_Te}")
                return False
            current_val = next_val

    return True

def main():
    # Read graph6 lines from stdin
    count = 0
    for line in sys.stdin.buffer:
        line = line.strip()
        if not line: continue
        
        n, adj = parse_graph6(line)
        if not check_tree(n, adj):
            print("Verification failed.")
            sys.exit(1)
        count += 1
        if count % 1000 == 0:
            print(f"Verified {count} trees...", file=sys.stderr)
            
    print(f"Success! Verified {count} trees.")

if __name__ == "__main__":
    main()
