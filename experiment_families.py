#!/usr/bin/env python3
"""Test mode-mean excess for specific tree families at large n.

The exhaustive search found max excess ≈ 0.62 through n=22.
All worst-case trees are caterpillar-type.

Let's check:
- Paths P_n (exactly computable)
- Stars K_{1,n-1}
- Caterpillars with specific leaf distributions
- The specific worst-case family: "one high-deg vertex on a long spine"
"""

import math
from indpoly import independence_poly

def build_path(n):
    adj = [[] for _ in range(n)]
    for i in range(n-1):
        adj[i].append(i+1)
        adj[i+1].append(i)
    return n, adj

def build_star(n):
    adj = [[] for _ in range(n)]
    for i in range(1, n):
        adj[0].append(i)
        adj[i].append(0)
    return n, adj

def build_caterpillar(spine_len, leaves_per_vertex):
    """Build a caterpillar: path of spine_len, each internal vertex gets leaves."""
    n = spine_len
    adj = [[] for _ in range(1000)]  # will trim
    # Build spine
    for i in range(spine_len - 1):
        adj[i].append(i+1)
        adj[i+1].append(i)
    # Add leaves to internal vertices
    v = spine_len
    for i in range(spine_len):
        num = leaves_per_vertex[i] if i < len(leaves_per_vertex) else 0
        for _ in range(num):
            adj[i].append(v)
            adj[v].append(i)
            n += 1
            v += 1
    adj = adj[:n]
    return n, adj

def build_broom(s, p):
    """Broom: star with s leaves, one arm of length p."""
    n = 1 + s + p
    adj = [[] for _ in range(n)]
    # Hub is vertex 0
    for i in range(1, s + 1):
        adj[0].append(i)
        adj[i].append(0)
    # Path from hub
    prev = 0
    for i in range(s + 1, n):
        adj[prev].append(i)
        adj[i].append(prev)
        prev = i
    return n, adj

def analyze(n, adj, name):
    poly = independence_poly(n, adj)
    total = sum(poly)
    mean = sum(k * poly[k] for k in range(len(poly))) / total
    mode = max(range(len(poly)), key=lambda k: poly[k])
    excess = mode - mean
    alpha = len(poly) - 1
    
    # a[m-1]/a[m]
    if mode > 0:
        near_flat = poly[mode-1] / poly[mode]
    else:
        near_flat = 0
    
    print(f"  {name:40s}  n={n:4d}  excess={excess:+.6f}  mode={mode}  mean={mean:.4f}  "
          f"alpha={alpha}  a[m-1]/a[m]={near_flat:.6f}")
    return excess

def main():
    print("=" * 100)
    print("  FAMILY ANALYSIS: Mode - Mean Excess at Large n")
    print("=" * 100)
    print()
    
    # Paths 
    print("─── PATHS ───")
    for n in [10, 20, 30, 50, 75, 100, 150, 200]:
        _, adj = build_path(n)
        analyze(n, adj, f"P_{n}")
    print()
    
    # Stars
    print("─── STARS ───")
    for n in [10, 20, 50, 100, 200]:
        _, adj = build_star(n)
        analyze(n, adj, f"K_1,{n-1}")
    print()
    
    # Brooms — found in worst-case analysis
    print("─── BROOMS ───")
    for n in [20, 30, 50, 100, 200]:
        # Try different hub sizes
        best_excess = -float('inf')
        best_name = ""
        for s in range(2, n - 2):
            p = n - 1 - s
            if p < 1:
                continue
            _, adj = build_broom(s, p)
            poly = independence_poly(n, adj)
            total = sum(poly)
            mean = sum(k * pk for k, pk in enumerate(poly)) / total
            mode = max(range(len(poly)), key=lambda k: poly[k])
            excess = mode - mean
            if excess > best_excess:
                best_excess = excess
                best_s = s
                best_p = p
                best_poly = poly
                best_mean = mean
                best_mode = mode
        
        analyze_result = best_excess
        name = f"Broom({best_s},{best_p}) n={n}"
        alpha = len(best_poly) - 1
        near_flat = best_poly[best_mode-1] / best_poly[best_mode] if best_mode > 0 else 0
        print(f"  {'BEST '+name:40s}  n={n:4d}  excess={best_excess:+.6f}  mode={best_mode}  mean={best_mean:.4f}  "
              f"alpha={alpha}  a[m-1]/a[m]={near_flat:.6f}")
    print()
    
    # Caterpillars with one high-degree vertex (like the worst-case trees found)
    print("─── CATERPILLARS WITH ONE HIGH-DEG VERTEX ───")
    for n in [20, 30, 50, 75, 100]:
        best_excess = -float('inf')
        # Try: spine of length L, one vertex gets k extra leaves
        for L in range(3, min(n, 20)):
            max_leaves = n - L
            for hub_pos in range(L):
                leaves = [0] * L
                leaves[hub_pos] = max_leaves
                _, adj = build_caterpillar(L, leaves)
                nn = len(adj)
                if nn != n:
                    continue
                poly = independence_poly(nn, adj)
                total = sum(poly)
                mean = sum(k * pk for k, pk in enumerate(poly)) / total
                mode = max(range(len(poly)), key=lambda k: poly[k])
                excess = mode - mean
                if excess > best_excess:
                    best_excess = excess
                    best_config = (L, hub_pos, max_leaves)
        
        L, hp, nl = best_config
        print(f"  Best caterpillar n={n}: spine={L}, hub@{hp} with {nl} leaves, excess={best_excess:.6f}")
    print()
    
    print("=" * 100)

if __name__ == "__main__":
    main()
