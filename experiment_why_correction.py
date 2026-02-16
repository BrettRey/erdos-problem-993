#!/usr/bin/env python3
"""Investigate: WHY does w1*(μ1-μ0) > excess(dp0) always hold at gap=+1 merges?

At a vertex v with children c_1,...,c_k:
  dp0[v] = prod_i total[c_i]
  dp1[v] = x * prod_i dp0[c_i]

  w0 = |dp0| / |total|, w1 = |dp1| / |total|  (unnormalized totals)
  
  μ0 = mean of dp0 (= sum of means of subtree totals)
  μ1 = mean of dp1 (= 1 + sum of means of subtree dp0's)

So:
  μ1 - μ0 = 1 + sum_i mean(dp0[c_i]) - sum_i mean(total[c_i])
           = 1 - sum_i [mean(total[c_i]) - mean(dp0[c_i])]
           = 1 - sum_i [w1_i * (mean(dp1[c_i]) - mean(dp0[c_i])) / 1]
           ... this is getting circular. Let me just measure.

Key insight: μ1 - μ0 = 1 - Δ where Δ = sum_i [mean(total_i) - mean(dp0_i)]
           = 1 - sum_i [contribution of dp1_i to the mean of total_i]

So the condition becomes:
  w1*(1-Δ) > excess(dp0)

Since dp0 is a product, its excess is ≤ sum of individual subtree excesses.
But the correction w1*(1-Δ) depends on the weight balance.

Let me track these quantities precisely.
"""

import subprocess
import math
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

def main():
    print("=" * 90)
    print("  WHY THE CORRECTION ALWAYS WINS")
    print("=" * 90)
    print()
    
    # For each gap=+1 merge, decompose:
    # excess(total) = 1 + excess(dp0) - w1*(μ1-μ0)
    # w1 = |dp1|/|total|
    # μ1 - μ0 = ?
    
    print(f"{'n':>3} {'v':>3} {'nc':>2} {'ex_dp0':>8} {'w1':>6} {'μ1-μ0':>8} "
          f"{'correction':>10} {'slack':>7} {'excess':>7}")
    print("-" * 70)
    
    # Collect all gap=+1 merges from worst-case trees
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
            _, _, ex, _ = poly_stats(poly)
            if ex > max_excess:
                max_excess = ex
                worst_line = line
        
        tn, adj = parse_graph6(worst_line)
        
        # Run DP, find the critical merge
        visited = [False] * tn
        children = [[] for _ in range(tn)]
        order = []
        queue = [0]
        visited[0] = True
        while queue:
            v = queue.pop(0)
            order.append(v)
            for u in adj[v]:
                if not visited[u]:
                    visited[u] = True
                    children[v].append(u)
                    queue.append(u)
        order.reverse()
        
        dp0 = [None] * tn
        dp1 = [None] * tn
        
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
            
            total = polyadd(dp0[v], dp1[v])
            mt, mut, ext, _ = poly_stats(total)
            m0, mu0, ex0, w0_raw = poly_stats(dp0[v])
            m1, mu1, ex1, w1_raw = poly_stats(dp1[v])
            
            gap = mt - m0
            
            if gap == 1 and (w0_raw + w1_raw) > 0:
                w1 = w1_raw / (w0_raw + w1_raw)
                mu_diff = mu1 - mu0
                correction = w1 * mu_diff
                slack = correction - ex0
                
                if v == 0 or ext == max_excess:
                    print(f"{n:3d} {v:3d} {len(children[v]):2d} {ex0:+8.4f} {w1:6.4f} "
                          f"{mu_diff:8.4f} {correction:10.4f} {slack:7.4f} {ext:7.4f}")
    
    print()
    
    # KEY OBSERVATION: what is w1 and μ1-μ0 for the critical merge?
    # w1 ≈ 0.2-0.3 (dp1 is the minority)
    # μ1-μ0 ≈ 1 - small_number ≈ 0.5-1.0
    # So correction ≈ 0.2 * 0.8 ≈ 0.16
    # And excess(dp0) ≈ -0.3 to -0.4
    # So slack ≈ 0.16 + 0.3 = 0.46 -- HUGE
    
    # Actually, let me look at the tightest cases across ALL trees
    print("TIGHTEST CASES (smallest slack across all gap=+1 merges):")
    print()
    
    tight_cases = []
    
    for n in range(3, 17):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        for line in lines:
            tn, adj = parse_graph6(line)
            
            visited = [False] * tn
            children = [[] for _ in range(tn)]
            order = []
            queue = [0]
            visited[0] = True
            while queue:
                v = queue.pop(0)
                order.append(v)
                for u in adj[v]:
                    if not visited[u]:
                        visited[u] = True
                        children[v].append(u)
                        queue.append(u)
            order.reverse()
            
            dp0 = [None] * tn
            dp1 = [None] * tn
            
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
                
                total = polyadd(dp0[v], dp1[v])
                mt, mut, ext, _ = poly_stats(total)
                m0, mu0, ex0, w0_raw = poly_stats(dp0[v])
                m1, mu1, ex1, w1_raw = poly_stats(dp1[v])
                
                gap = mt - m0
                
                if gap == 1 and (w0_raw + w1_raw) > 0:
                    w1 = w1_raw / (w0_raw + w1_raw)
                    mu_diff = mu1 - mu0
                    correction = w1 * mu_diff
                    slack = correction - ex0
                    
                    tight_cases.append((slack, n, v, len(children[v]), ex0, w1, mu_diff, correction, ext))
    
    tight_cases.sort()
    print(f"{'slack':>7} {'n':>3} {'v':>3} {'nc':>2} {'ex_dp0':>8} {'w1':>6} {'μ1-μ0':>8} "
          f"{'correct':>8} {'excess':>7}")
    for slack, n, v, nc, ex0, w1, mu_diff, correction, ext in tight_cases[:15]:
        print(f"{slack:7.4f} {n:3d} {v:3d} {nc:2d} {ex0:+8.4f} {w1:6.4f} "
              f"{mu_diff:8.4f} {correction:8.4f} {ext:7.4f}")
    
    print()
    
    # Pattern in the tightest cases
    if tight_cases:
        min_slack_case = tight_cases[0]
        slack, n, v, nc, ex0, w1, mu_diff, correction, ext = min_slack_case
        print(f"TIGHTEST: n={n}, v={v}, children={nc}")
        print(f"  excess(dp0)={ex0:+.6f}, w1={w1:.6f}, μ1-μ0={mu_diff:.6f}")
        print(f"  correction={correction:.6f}, slack={slack:.6f}")
        print(f"  excess(total)={ext:.6f}")
        print()
        print(f"  For the bound to fail, we'd need slack < 0")
        print(f"  i.e., w1*(μ1-μ0) < excess(dp0)")
        print(f"  Current tightest: {correction:.4f} vs {ex0:.4f}, ratio = {correction/ex0:.4f}" if ex0 > 0 else "")
    
    print()
    print("=" * 90)

if __name__ == "__main__":
    main()
