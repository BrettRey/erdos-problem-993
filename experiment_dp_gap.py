#!/usr/bin/env python3
"""Experiment: Does the DP decomposition bound the mode-mean gap?

A tree's IS poly decomposes as:
  I(T) = I(T-v) + x·I(T-v-N(v))

when rooting at vertex v. Both factors are IS polys of smaller forests.

If p = a * b (convolution of two polys), then:
  mode(p) ≤ mode(a) + mode(b)
  mean(p) = mean(a) + mean(b)

So: excess(p) = mode(p) - mean(p) ≤ excess(a) + excess(b)

But the tree IS poly is a SUM, not a product:
  I(T) = I(T-v) + x·I(T-v-N(v))

So the convolution bound doesn't directly apply. Instead, the mode
of a sum of two polys can differ from individual modes in complex ways.

Let's check: for the DP decomposition, how does the mode-mean gap
compose through the tree? Track it through the leaf-to-root DP.
"""

import subprocess
import math
from indpoly import independence_poly
from graph6 import parse_graph6

def poly_stats(poly):
    """Return (mode, mean, excess) for a polynomial coefficient sequence."""
    if not poly or sum(poly) == 0:
        return 0, 0, 0
    total = sum(poly)
    mean = sum(k * poly[k] for k in range(len(poly))) / total
    mode = max(range(len(poly)), key=lambda k: poly[k])
    return mode, mean, mode - mean

def dp_decompose(adj, n):
    """Run the DP and track mode-mean gap at each step."""
    # Root at vertex 0, BFS order
    visited = [False] * n
    parent = [-1] * n
    children = [[] for _ in range(n)]
    order = []
    
    queue = [0]
    visited[0] = True
    while queue:
        v = queue.pop(0)
        order.append(v)
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                parent[u] = v
                children[v].append(u)
                queue.append(u)
    
    # Process leaves first
    order.reverse()
    
    # dp0[v] = IS poly when v is NOT in the set
    # dp1[v] = IS poly when v IS in the set (multiply by x)
    dp0 = [None] * n
    dp1 = [None] * n
    
    # Track gap at each step
    gap_info = []
    
    for v in order:
        if not children[v]:
            # Leaf
            dp0[v] = [1]
            dp1[v] = [0, 1]
        else:
            # dp0[v] = product of (dp0[c] + dp1[c]) over children c
            prod = [1]
            for c in children[v]:
                summand = polyadd(dp0[c], dp1[c])
                prod = polymul(prod, summand)
            dp0[v] = prod
            
            # dp1[v] = x * product of dp0[c] over children c
            prod = [1]
            for c in children[v]:
                prod = polymul(prod, dp0[c])
            dp1[v] = [0] + prod  # multiply by x
        
        # Total IS poly if this were the root
        total_poly = polyadd(dp0[v], dp1[v])
        m, mu, ex = poly_stats(total_poly)
        gap_info.append((v, len(children[v]), m, mu, ex))
    
    return gap_info

def polyadd(a, b):
    n = max(len(a), len(b))
    result = [0] * n
    for i in range(len(a)):
        result[i] += a[i]
    for i in range(len(b)):
        result[i] += b[i]
    return result

def polymul(a, b):
    if not a or not b:
        return [0]
    n = len(a) + len(b) - 1
    result = [0] * n
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[i + j] += ai * bj
    return result

def main():
    print("=" * 70)
    print("  DP DECOMPOSITION: Mode-Mean Gap Through Tree DP")
    print("=" * 70)
    print()
    
    # Check the worst-case tree from n=18
    # poly = [1, 18, 136, 569, 1453, 2337, 2338, 1364, 392, 32]
    # That specific tree has excess 0.6188
    
    # For now, check some medium trees
    for n in range(5, 16):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        # Find worst-case tree
        max_excess = -float('inf')
        worst_line = None
        
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            _, _, ex = poly_stats(poly)
            if ex > max_excess:
                max_excess = ex
                worst_line = line
        
        tn, adj = parse_graph6(worst_line)
        poly = independence_poly(tn, adj)
        mode, mean, excess = poly_stats(poly)
        print(f"n={n}: excess={excess:.4f}, poly={poly}")
        
        gap_info = dp_decompose(adj, tn)
        print(f"  DP gaps (vertex, #children, mode, mean, excess):")
        for v, nc, m, mu, ex in gap_info:
            print(f"    v={v}, children={nc}: mode={m}, mean={mu:.4f}, excess={ex:+.4f}")
        
        # What's the max intermediate excess?
        max_intermediate = max(ex for _, _, _, _, ex in gap_info)
        print(f"  Max intermediate excess: {max_intermediate:.4f}")
        print()
    
    print("=" * 70)

if __name__ == "__main__":
    main()
