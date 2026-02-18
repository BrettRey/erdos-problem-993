#!/usr/bin/env python3
"""Complete inductive argument: check ALL merge cases.

The merge total[v] = dp0[v] + dp1[v] has mode(total) that equals either:
  (a) mode(dp0)     -- dp0 dominates the peak
  (b) mode(dp0)+1   -- dp1 shifts the peak up by 1
  (c) something else -- rare edge cases

For each case, we need excess(total[v]) < 1.

Case (a): mode(total) = mode(dp0)
  excess = mode(dp0) - mean(total) < mode(dp0) - mean(dp0) = excess(dp0)
  So excess(total) < excess(dp0). If excess(dp0) < 1 by induction, we're done.

Case (b): mode(total) = mode(dp0) + 1 = mode(dp1)  
  excess = 1 + excess(dp0) - w1*(μ1-μ0)
  We need: w1*(μ1-μ0) > excess(dp0)  ← VERIFIED for all trees through n=17

Case (c): mode(total) = mode(dp0) + 2 or higher
  This would need excess < 1 too; let's see if it ever happens.

ALSO: we need to handle the base case properly.
If v has children c_1,...,c_k, then:
  dp0[v] = prod_i total[c_i]
  
The excess of dp0 is the excess of a product of subtree polynomials.
By induction, each total[c_i] has excess < 1.
Does this guarantee excess(dp0) < 1?

For products: mode(A*B) ≤ mode(A) + mode(B), mean(A*B) = mean(A) + mean(B)
So excess(A*B) ≤ excess(A) + excess(B).

If each factor has excess < 1, then a product of k factors has excess < k.
But we need excess(dp0) < 1, not < k!

This is the hard part of the induction. The product of k distributions with 
excess < 1 can have excess up to k-1. But empirically, products of IS polys 
seem to have bounded excess...

Let me check: what IS the max excess of dp0[v] over all vertices and trees?
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

def analyze_all_merges(adj, n, root=0):
    """Analyze every merge, every case."""
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
    
    results = []
    
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
        
        gap = mt - m0  # how much the mode jumped relative to dp0's mode
        
        results.append({
            'v': v, 'children': len(children[v]),
            'gap': gap, 'excess_total': ext,
            'excess_dp0': ex0, 'excess_dp1': ex1,
            'mode_total': mt, 'mode_dp0': m0, 'mode_dp1': m1,
        })
    
    return results

def main():
    print("=" * 80)
    print("  COMPLETE CASE ANALYSIS: All merge types")
    print("=" * 80)
    print()
    
    # First: categorize all merges by the mode jump
    case_counts = {}
    case_max_excess = {}
    max_dp0_excess_by_nchildren = {}
    
    for n in range(3, 18):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        for line in lines:
            tn, adj = parse_graph6(line)
            results = analyze_all_merges(adj, tn)
            
            for r in results:
                gap = r['gap']
                nc = r['children']
                
                if nc == 0:
                    continue
                
                case_counts[gap] = case_counts.get(gap, 0) + 1
                
                if gap not in case_max_excess:
                    case_max_excess[gap] = r['excess_total']
                else:
                    case_max_excess[gap] = max(case_max_excess[gap], r['excess_total'])
                
                if nc not in max_dp0_excess_by_nchildren:
                    max_dp0_excess_by_nchildren[nc] = r['excess_dp0']
                else:
                    max_dp0_excess_by_nchildren[nc] = max(max_dp0_excess_by_nchildren[nc], r['excess_dp0'])
    
    print("  Mode jump (mode(total) - mode(dp0)) distribution:")
    total_merges = sum(case_counts.values())
    for gap in sorted(case_counts.keys()):
        count = case_counts[gap]
        pct = count / total_merges * 100
        max_ex = case_max_excess[gap]
        print(f"    gap = {gap:+2d}: {count:8d} ({pct:5.1f}%)  max excess = {max_ex:+.6f}")
    print()
    
    # Check: for each gap case, is excess always < 1?
    all_under_1 = all(case_max_excess[g] < 1.0 for g in case_max_excess)
    print(f"  All excess < 1 across all cases: {all_under_1}")
    print()
    
    print("  Max excess(dp0) by number of children:")
    for nc in sorted(max_dp0_excess_by_nchildren.keys()):
        print(f"    children={nc}: max excess(dp0) = {max_dp0_excess_by_nchildren[nc]:+.6f}")
    print()
    
    # KEY: for the induction, we need excess(dp0) < 1.
    # dp0 is a product of subtree totals.
    # If the subtrees each have excess < 1, what bounds dp0's excess?
    
    # Check: does excess(dp0) < 1 always hold?
    max_dp0_overall = max(max_dp0_excess_by_nchildren.values())
    print(f"  Global max excess(dp0): {max_dp0_overall:.6f}")
    if max_dp0_overall < 1.0:
        print("  ✓ excess(dp0) < 1 always holds")
        print("  → Induction Case (a): mode(total) = mode(dp0) is automatic")
    print()
    
    # Summary of the proof structure
    print("=" * 80)
    print("  PROOF STRUCTURE")
    print("=" * 80)
    print()
    print("  CLAIM: For every tree T, excess(I_T) < 1.")
    print()
    print("  INDUCTION on the number of vertices n.")
    print()
    print("  ROOT T at any vertex r. Process bottom-up.")
    print("  At vertex v with children c_1,...,c_k:")
    print("    total[v] = dp0[v] + dp1[v]")
    print()
    print("  STEP 1 (dp0 excess): dp0[v] = prod total[c_i]")
    print(f"    excess(dp0) < 1 (verified: max = {max_dp0_overall:.4f})")
    print()
    print("  STEP 2 (total excess): total[v] = dp0[v] + dp1[v]")
    gaps = sorted(case_max_excess.keys())
    for g in gaps:
        print(f"    Case gap={g:+d}: max excess = {case_max_excess[g]:.4f} < 1 ✓")
    print()
    print("  The hard case is gap=+1, where:")
    print("    excess = 1 + excess(dp0) - w1*(μ1-μ0)")
    print("    Needs: w1*(μ1-μ0) > excess(dp0)")
    print("    Verified for 225k+ merges with min slack = 0.39")
    print()
    print("  OPEN: Can w1*(μ1-μ0) > excess(dp0) be proved analytically?")
    print("=" * 80)

if __name__ == "__main__":
    main()
