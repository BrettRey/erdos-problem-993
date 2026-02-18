#!/usr/bin/env python3
"""Deep analysis of the critical merge step.

Observation: the excess at the root equals the max intermediate excess.
The jump happens at a specific vertex. Let's understand what happens there.

At vertex v with children c_1, ..., c_k:
  dp0[v] = prod_i (dp0[c_i] + dp1[c_i])     -- "v out": children free
  dp1[v] = x * prod_i dp0[c_i]               -- "v in":  children forced out
  total[v] = dp0[v] + dp1[v]

The critical observation is that dp0[v] is a product of subtree totals,
and dp1[v] is x times a product of subtree "out" polynomials.

For the mode-mean excess of dp0[v] + dp1[v]:
  - dp0[v] has mode m0, mean μ0
  - dp1[v] has mode m1, mean μ1 (shifted by 1 due to x factor)
  - The combined mode depends on where the sum peaks

Key question: can we bound mode(A+B) - mean(A+B) in terms of 
properties of A and B?
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
        return {'mode': 0, 'mean': 0, 'excess': 0, 'total': 0}
    total = sum(poly)
    mean = sum(k * poly[k] for k in range(len(poly))) / total
    mode = max(range(len(poly)), key=lambda k: poly[k])
    var = sum(k**2 * poly[k] for k in range(len(poly))) / total - mean**2
    return {
        'mode': mode, 'mean': mean, 'excess': mode - mean,
        'total': total, 'var': var, 'sd': math.sqrt(max(0, var)),
        'alpha': len(poly) - 1,
    }

def dp_with_merge_analysis(adj, n, root=0):
    """Run DP and analyze the critical merge at each vertex."""
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
    
    merges = []
    
    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
        else:
            prod_total = [1]
            for c in children[v]:
                subtotal = polyadd(dp0[c], dp1[c])
                prod_total = polymul(prod_total, subtotal)
            dp0[v] = prod_total
            
            prod_out = [1]
            for c in children[v]:
                prod_out = polymul(prod_out, dp0[c])
            dp1[v] = [0] + prod_out
        
        total = polyadd(dp0[v], dp1[v])
        s_total = poly_stats(total)
        s_dp0 = poly_stats(dp0[v])
        s_dp1 = poly_stats(dp1[v])
        
        # The merge: total = dp0 + dp1
        # dp0 comes from "v excluded" (product of subtree totals)
        # dp1 comes from "v included" (x * product of subtree outs)
        
        # Weight of dp0 vs dp1 in the sum
        w0 = s_dp0['total'] / (s_dp0['total'] + s_dp1['total']) if s_dp0['total'] + s_dp1['total'] > 0 else 0
        w1 = 1 - w0
        
        merges.append({
            'v': v,
            'children': len(children[v]),
            'total_excess': s_total['excess'],
            'dp0_excess': s_dp0['excess'],
            'dp1_excess': s_dp1['excess'],
            'dp0_mode': s_dp0['mode'],
            'dp1_mode': s_dp1['mode'],
            'dp0_mean': s_dp0['mean'],
            'dp1_mean': s_dp1['mean'],
            'w0': w0,
            'w1': w1,
            'combined_mean': w0 * s_dp0['mean'] + w1 * s_dp1['mean'],
            'mode_of_sum': s_total['mode'],
            'total_poly': total,
        })
    
    return merges

def main():
    print("=" * 80)
    print("  CRITICAL MERGE ANALYSIS")
    print("=" * 80)
    print()
    
    # Analyze worst-case trees at n=14 and n=18
    for n in [8, 14, 18, 21]:
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        max_excess = -float('inf')
        worst_line = None
        
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            s = poly_stats(poly)
            if s['excess'] > max_excess:
                max_excess = s['excess']
                worst_line = line
        
        tn, adj = parse_graph6(worst_line)
        poly = independence_poly(tn, adj)
        
        print(f"─── Worst-case tree at n={n}, excess={max_excess:.4f} ───")
        print(f"    poly = {poly}")
        print()
        
        # Try multiple roots to find the one that makes the merge structure clearest
        best_root = 0
        best_critical_excess = float('inf')
        
        for root in range(n):
            merges = dp_with_merge_analysis(adj, tn, root)
            critical = max(merges, key=lambda m: m['total_excess'])
            if critical['total_excess'] < best_critical_excess:
                best_critical_excess = critical['total_excess']
                best_root = root
        
        # But for now, just analyze with root 0
        merges = dp_with_merge_analysis(adj, tn, root=0)
        critical = max(merges, key=lambda m: m['total_excess'])
        
        print(f"    Critical merge at v={critical['v']} (children={critical['children']})")
        print(f"    dp0: mode={critical['dp0_mode']}, mean={critical['dp0_mean']:.4f}, "
              f"excess={critical['dp0_excess']:.4f}, weight={critical['w0']:.4f}")
        print(f"    dp1: mode={critical['dp1_mode']}, mean={critical['dp1_mean']:.4f}, "
              f"excess={critical['dp1_excess']:.4f}, weight={critical['w1']:.4f}")
        print(f"    SUM: mode={critical['mode_of_sum']}, mean={critical['combined_mean']:.4f}, "
              f"excess={critical['total_excess']:.4f}")
        print()
        
        # Key analysis: WHY does mode(dp0+dp1) - mean(dp0+dp1) stay < 1?
        # Is it because dp0 and dp1 have overlapping supports near the mode?
        dp0 = merges[-1]  # root stats
        
        print("    All merge steps (sorted by excess):")
        sorted_merges = sorted(merges, key=lambda m: m['total_excess'], reverse=True)
        for m in sorted_merges[:5]:
            print(f"      v={m['v']:2d} children={m['children']}: "
                  f"excess={m['total_excess']:+.4f}  "
                  f"[dp0: mode={m['dp0_mode']}, dp1: mode={m['dp1_mode']}, "
                  f"gap={m['dp1_mode']-m['dp0_mode']:+d}]")
        print()
        
        # CRUCIAL: What is mode(dp1) - mode(dp0) at the critical merge?
        # If dp1 = x * prod(dp0[c_i]), its mode is roughly 1 + sum(mode(dp0[c_i]))
        # If dp0 = prod(total[c_i]), its mode is roughly sum(mode(total[c_i]))
        # The difference = 1 + sum(mode(dp0[c_i])) - sum(mode(total[c_i]))
        # Since mode(total[c_i]) >= mode(dp0[c_i]) typically, this difference is ≤ 1
        
        print(f"    *** mode(dp1) - mode(dp0) = {critical['dp1_mode'] - critical['dp0_mode']:+d} "
              f"at critical merge ***")
        print()
    
    # Summary: check mode(dp1) - mode(dp0) across ALL merges in ALL trees
    print("=" * 80)
    print("  GLOBAL CHECK: mode(dp1) - mode(dp0) at every DP merge")
    print("=" * 80)
    print()
    
    max_mode_gap = -100
    total_merges = 0
    gap_counts = {}
    
    for n in range(3, 17):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        for line in lines:
            tn, adj = parse_graph6(line)
            merges = dp_with_merge_analysis(adj, tn, root=0)
            for m in merges:
                if m['children'] > 0:
                    gap = m['dp1_mode'] - m['dp0_mode']
                    total_merges += 1
                    gap_counts[gap] = gap_counts.get(gap, 0) + 1
                    if gap > max_mode_gap:
                        max_mode_gap = gap
    
    print(f"  Total non-leaf merges analyzed: {total_merges}")
    print(f"  Max mode(dp1) - mode(dp0) = {max_mode_gap}")
    print()
    print("  Distribution of mode(dp1) - mode(dp0):")
    for gap in sorted(gap_counts.keys()):
        pct = gap_counts[gap] / total_merges * 100
        print(f"    gap = {gap:+2d}: {gap_counts[gap]:8d} ({pct:5.1f}%)")
    
    print()
    
    if max_mode_gap <= 1:
        print("  ═══════════════════════════════════════════════════════════════")
        print("  KEY FINDING: mode(dp1) - mode(dp0) ≤ 1 at EVERY merge step.")
        print("  This means dp1 (v-in branch) peaks at most 1 position above")
        print("  dp0 (v-out branch). When combined, the sum's mode is pulled")
        print("  toward the heavier component, limiting the excess.")
        print("  ═══════════════════════════════════════════════════════════════")
    
    print()
    print("=" * 80)

if __name__ == "__main__":
    main()
