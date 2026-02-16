#!/usr/bin/env python3
"""Attempt to prove: excess(total) < 1 at every vertex in the tree DP.

THEOREM (mode-mean bound via DP):
  For every tree T, the IS poly I_T(x) satisfies mode(I_T) < mean(I_T) + 1.

PROOF ATTEMPT:
  Root T at any vertex. Process vertices bottom-up.
  
  At each vertex v with children c_1, ..., c_k:
    dp0[v] = prod_i total[c_i]  (v excluded)
    dp1[v] = x * prod_i dp0[c_i]  (v included)
    total[v] = dp0[v] + dp1[v]
  
  The excess at v is:
    excess(total[v]) = mode(total[v]) - mean(total[v])
  
  CASE 1: mode(total[v]) = mode(dp0[v])
    Then excess = mode(dp0[v]) - mean(total[v])
               < mode(dp0[v]) - mean(dp0[v])  [since mean(total) > mean(dp0)]
               = excess(dp0[v])
    
    Since dp0[v] is a product of subtree totals, and the mode(product) ≤ sum(modes),
    while the mean(product) = sum(means), we need excess(dp0[v]) < 1.
    
    KEY FACT: excess of a product of polys with nonneg coefficients ≤ sum of excesses.
    If each subtree total has excess < 1, then by induction, dp0 has excess < k 
    (where k = number of children). But that's too weak!
    
    Actually: since dp0 is a product of independent distributions, 
    the CLT makes dp0 approach a normal distribution (mode ≈ mean) as k grows.
    For k=1, dp0 = total[c_1], which has excess < 1 by induction.
    For k≥2, the convolution concentrates the distribution.
  
  CASE 2: mode(total[v]) = mode(dp0[v]) + 1 = mode(dp1[v])
    This is the "critical merge" case. Here:
    excess = 1 + excess(dp0) - w1*(mean(dp1) - mean(dp0))
    
    Since excess(dp0) < 0 at the critical merge (empirically always),
    and w1*(mean(dp1)-mean(dp0)) ≥ 0, we get excess < 1.
    
    BUT: can we PROVE excess(dp0) < 0?
    
    dp0[v] is a product of subtree totals. Each subtree total has some excess.
    Product of polys: the excess of a convolution concentrates toward 0.
    
    Hmm, this needs more work. Let me check: is excess(dp0) < 0 ALWAYS,
    not just at critical merges?

Let me verify: is excess(dp0[v]) ≤ 0 at every vertex, for every tree?
"""

import subprocess
from indpoly import independence_poly
from graph6 import parse_graph6

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

def poly_stats(poly):
    if not poly or sum(poly) == 0:
        return 0, 0, 0, 0
    total = sum(poly)
    mean = sum(k * poly[k] for k in range(len(poly))) / total
    mode = max(range(len(poly)), key=lambda k: poly[k])
    return mode, mean, mode - mean, total

def check_dp0_excess(adj, n, root=0):
    """Check if excess(dp0[v]) ≤ 0 at every internal vertex."""
    visited = [False] * n
    children = [[] for _ in range(n)]
    order = []
    
    queue = [root]
    visited[root] = True
    while queue:
        v = queue.pop(0)
        order.append(v)
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                children[v].append(u)
                queue.append(u)
    
    order.reverse()
    dp0 = [None] * n
    dp1 = [None] * n
    
    max_dp0_excess = -100
    violation = False
    
    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
        else:
            prod = [1]
            for c in children[v]:
                prod = polymul(prod, polyadd(dp0[c], dp1[c]))
            dp0[v] = prod
            
            prod_out = [1]
            for c in children[v]:
                prod_out = polymul(prod_out, dp0[c])
            dp1[v] = [0] + prod_out
        
        if children[v]:
            m, mu, ex, _ = poly_stats(dp0[v])
            if ex > max_dp0_excess:
                max_dp0_excess = ex
            if ex > 0:
                violation = True
    
    return max_dp0_excess, violation

def main():
    print("=" * 70)
    print("  IS excess(dp0[v]) ≤ 0 AT EVERY INTERNAL VERTEX?")
    print("=" * 70)
    print()
    
    total_violations = 0
    total_trees = 0
    global_max_dp0_excess = -100
    
    for n in range(3, 18):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        n_violations = 0
        n_max_excess = -100
        
        for line in lines:
            tn, adj = parse_graph6(line)
            max_ex, viol = check_dp0_excess(adj, tn)
            total_trees += 1
            if viol:
                n_violations += 1
                total_violations += 1
            if max_ex > n_max_excess:
                n_max_excess = max_ex
            if max_ex > global_max_dp0_excess:
                global_max_dp0_excess = max_ex
        
        print(f"  n={n:2d}: trees={len(lines):6d}, violations={n_violations}, "
              f"max excess(dp0)={n_max_excess:+.6f}")
    
    print()
    if total_violations == 0:
        print("  ═══════════════════════════════════════════════════════════════")
        print(f"  RESULT: excess(dp0) ≤ 0 at EVERY internal vertex")
        print(f"  in EVERY tree through n=17 ({total_trees} trees)")
        print(f"  Global max excess(dp0) = {global_max_dp0_excess:+.6f}")
        print("  ═══════════════════════════════════════════════════════════════")
    else:
        print(f"  VIOLATIONS FOUND: {total_violations} trees have excess(dp0) > 0")
        print(f"  Global max excess(dp0) = {global_max_dp0_excess:+.6f}")
    
    print()
    
    # Now check: is excess(dp0) ≤ 0 because dp0 is a product?
    # Products of unimodal distributions: the mode of a product ≤ sum of modes
    # The mean of a product = sum of means
    # So excess(product) ≤ sum(excess_i)
    # For this to give excess ≤ 0, we need sum(excess_i) ≤ 0
    # But each subtree total has excess POSITIVE (for worst-case trees)!
    
    # So the story must be subtler. Let me check:
    # For 2+ children, does the product of unimodal distributions 
    # always have mode ≤ mean?
    
    print("  Check: mode ≤ mean for products of IS polys?")
    print("  (Products of subtree total polynomials)")
    print()
    
    # Check specifically: if we take 2 or more IS polys and multiply them,
    # does the product always have mode ≤ mean?
    
    # Test: multiply random pairs of IS polys and check
    from itertools import combinations
    
    # Collect all IS polys for small trees
    all_polys = {}
    for n in range(3, 9):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        polys = []
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            polys.append(poly)
        all_polys[n] = polys
    
    product_violations = 0
    product_tests = 0
    max_product_excess = -100
    
    for n1 in range(3, 8):
        for n2 in range(n1, 8):
            for p1 in all_polys[n1]:
                for p2 in all_polys[n2]:
                    prod = polymul(p1, p2)
                    m, mu, ex, _ = poly_stats(prod)
                    product_tests += 1
                    if ex > max_product_excess:
                        max_product_excess = ex
                    if ex > 0:
                        product_violations += 1
    
    print(f"  Products tested: {product_tests}")
    print(f"  Products with mode > mean: {product_violations}")
    print(f"  Max excess of product: {max_product_excess:+.6f}")
    
    if product_violations > 0:
        print("  → Products of IS polys CAN have mode > mean!")
        print("  → The dp0 excess ≤ 0 isn't just from being a product")
        print("  → Something specific about the TREE structure enforces it")
    else:
        print("  → All products have mode ≤ mean! The product structure suffices.")
    
    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
