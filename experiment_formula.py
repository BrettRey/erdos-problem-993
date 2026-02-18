#!/usr/bin/env python3
"""Understand the mechanism: why does mode(dp1) ≈ mode(dp0) + 1?

Theory:
  dp0[v] = prod_i total[c_i]    -- v excluded, children unconstrained
  dp1[v] = x * prod_i dp0[c_i]  -- v included, children forced out

  mode(dp0[v]) ≈ sum_i mode(total[c_i])
  mode(dp1[v]) ≈ 1 + sum_i mode(dp0[c_i])

  gap = mode(dp1) - mode(dp0) ≈ 1 + sum_i (mode(dp0[c_i]) - mode(total[c_i]))

  Since mode(total[c_i]) = mode(dp0[c_i] + dp1[c_i]) ≥ mode(dp0[c_i]) typically
  (because dp1 adds mass above dp0's mode), we get:

  gap ≈ 1 - sum_i (mode(total[c_i]) - mode(dp0[c_i]))

  If mode(total) = mode(dp0) for every child (common case), then gap = 1.
  If mode(total) > mode(dp0) for some child (i.e. dp1's mass shifted the mode up),
  then gap < 1, so gap = 0.

This explains the observation that gap is almost always 0 or 1!

Now for the excess bound: if total = dp0 + dp1 and:
  - dp0 has weight w0 ≈ 0.75, mode m
  - dp1 has weight w1 ≈ 0.25, mode m+1 (or m)
  - mean(total) = w0*mean(dp0) + w1*mean(dp1)
  - mode(total) = m or m+1

Then excess = mode(total) - mean(total).

If mode(total) = m+1 (same as dp1's mode):
  excess = (m+1) - (w0*μ0 + w1*μ1)
         = (m+1) - μ0 - w1*(μ1 - μ0)
         = (m+1 - μ0) - w1*(μ1 - μ0)

Since m = mode(dp0) and excess(dp0) = m - μ0:
  excess = 1 + excess(dp0) - w1*(μ1 - μ0)

Since dp1 peaks higher, μ1 > μ0, so w1*(μ1-μ0) > 0.
The excess is controlled by: 1 + excess(dp0) - (positive correction).

The question: can 1 + excess(dp0) ever be large enough to overcome the correction?

Let's check this empirically.
"""

import subprocess
import math

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

from indpoly import independence_poly
from graph6 import parse_graph6

def dp_excess_decomposition(adj, n, root=0):
    """Decompose the excess at the root into the formula components."""
    visited = [False] * n
    parent = [-1] * n
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
                parent[u] = v
                children[v].append(u)
                queue.append(u)
    
    order.reverse()
    dp0 = [None] * n
    dp1 = [None] * n
    
    results = []
    
    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
        else:
            prod_total = [1]
            for c in children[v]:
                prod_total = polymul(prod_total, polyadd(dp0[c], dp1[c]))
            dp0[v] = prod_total
            
            prod_out = [1]
            for c in children[v]:
                prod_out = polymul(prod_out, dp0[c])
            dp1[v] = [0] + prod_out
        
        total = polyadd(dp0[v], dp1[v])
        m0, mu0, ex0, w0 = poly_stats(dp0[v])
        m1, mu1, ex1, w1 = poly_stats(dp1[v])
        mt, mut, ext, wt = poly_stats(total)
        
        if w0 + w1 > 0:
            weight0 = w0 / (w0 + w1)
            weight1 = w1 / (w0 + w1)
        else:
            weight0 = weight1 = 0
        
        # The formula: excess(total) = 1 + excess(dp0) - weight1*(mean1 - mean0)
        # (when mode(total) = mode(dp1) = mode(dp0)+1)
        if m1 == m0 + 1 and mt == m1:
            formula_val = 1 + ex0 - weight1 * (mu1 - mu0)
            results.append({
                'v': v, 'children': len(children[v]),
                'excess': ext, 'formula': formula_val,
                'ex_dp0': ex0, 'w1': weight1,
                'mu_diff': mu1 - mu0,
                'correction': weight1 * (mu1 - mu0),
            })
    
    return results

def main():
    print("=" * 80)
    print("  EXCESS FORMULA VERIFICATION")
    print("  excess(total) = 1 + excess(dp0) - w1*(μ1 - μ0)")
    print("=" * 80)
    print()
    
    match_count = 0
    total_count = 0
    max_correction = 0
    min_correction = float('inf')
    
    for n in range(5, 18):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        # Find worst-case tree
        max_excess = -float('inf')
        worst_line = None
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            m, mu, ex, _ = poly_stats(poly)
            if ex > max_excess:
                max_excess = ex
                worst_line = line
        
        tn, adj = parse_graph6(worst_line)
        results = dp_excess_decomposition(adj, tn)
        
        # Find the critical merge (max excess)
        if results:
            critical = max(results, key=lambda r: r['excess'])
            err = abs(critical['excess'] - critical['formula'])
            print(f"n={n:2d}: excess={critical['excess']:.6f}, formula={critical['formula']:.6f}, "
                  f"err={err:.1e}, 1+ex_dp0={1+critical['ex_dp0']:.4f}, "
                  f"correction={critical['correction']:.4f}, w1={critical['w1']:.4f}")
            
            if err < 1e-6:
                match_count += 1
            total_count += 1
            max_correction = max(max_correction, critical['correction'])
            min_correction = min(min_correction, critical['correction'])
    
    print()
    print(f"Formula matches: {match_count}/{total_count}")
    print(f"Correction range: [{min_correction:.4f}, {max_correction:.4f}]")
    print()
    
    # So excess < 1 iff correction > excess(dp0)
    # i.e., w1*(μ1 - μ0) > excess(dp0)
    # i.e., the weight-adjusted mean shift exceeds dp0's own excess
    
    print("KEY INSIGHT:")
    print("  excess(total) = 1 + excess(dp0) - w1·(μ1 - μ0)")
    print("  excess < 1  ⟺  w1·(μ1 - μ0) > excess(dp0)")
    print()
    print("  Since dp0 is a PRODUCT of subtree totals, and products of")
    print("  bell-shaped distributions are nearly symmetric, excess(dp0)")
    print("  is typically < 0 (right-skewed by the product structure).")
    print("  This makes the bound easier to satisfy.")
    print()
    
    # Verify: is excess(dp0) < 0 at the critical merge?
    print("  Check: is excess(dp0) < 0 at critical merges?")
    for n in range(5, 18):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        max_excess = -float('inf')
        worst_line = None
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            m, mu, ex, _ = poly_stats(poly)
            if ex > max_excess:
                max_excess = ex
                worst_line = line
        
        tn, adj = parse_graph6(worst_line)
        results = dp_excess_decomposition(adj, tn)
        if results:
            critical = max(results, key=lambda r: r['excess'])
            print(f"    n={n}: excess(dp0) = {critical['ex_dp0']:+.4f}")
    
    print()
    print("=" * 80)

if __name__ == "__main__":
    main()
