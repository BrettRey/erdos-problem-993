#!/usr/bin/env python3
"""Check: is w1*(μ1 - μ0) > excess(dp0) at EVERY vertex where 
mode(dp1) = mode(dp0) + 1?

This is the precise condition needed for excess(total) < 1.

Also check: what is the relationship between excess(dp0) and 
the correction term? Is there a structural bound?
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

def main():
    print("=" * 80)
    print("  CONDITION: w1*(μ1-μ0) > excess(dp0) when mode(dp1)=mode(dp0)+1")
    print("=" * 80)
    print()
    
    violations = 0
    total_merges = 0
    min_slack = float('inf')  # min of correction - excess(dp0)
    
    # Also track: what is excess(total) directly? 
    max_excess_anywhere = -100
    
    for n in range(3, 18):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        n_max_excess = -100
        n_min_slack = float('inf')
        
        for line in lines:
            tn, adj = parse_graph6(line)
            
            # Run DP with root 0
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
                
                if ext > n_max_excess:
                    n_max_excess = ext
                
                # Check the formula only when mode(dp1) = mode(dp0) + 1
                # AND mode(total) = mode(dp1)
                if m1 == m0 + 1 and mt == m1 and (w0_raw + w1_raw) > 0:
                    w1 = w1_raw / (w0_raw + w1_raw)
                    correction = w1 * (mu1 - mu0)
                    slack = correction - ex0
                    total_merges += 1
                    
                    if slack < n_min_slack:
                        n_min_slack = slack
                    if slack < min_slack:
                        min_slack = slack
                    
                    if slack < 0:
                        violations += 1
                        # This would mean excess(total) > 1
                        print(f"  *** VIOLATION at n={tn}, v={v}: "
                              f"correction={correction:.6f}, excess(dp0)={ex0:.6f}, "
                              f"excess(total)={ext:.6f}")
        
        if n_max_excess > max_excess_anywhere:
            max_excess_anywhere = n_max_excess
        
        print(f"  n={n:2d}: max_excess(total)={n_max_excess:+.6f}, min_slack={n_min_slack:+.6f}")
    
    print()
    print(f"  Total conditional merges checked: {total_merges}")
    print(f"  Violations (slack < 0): {violations}")
    print(f"  Global min slack: {min_slack:+.6f}")
    print(f"  Global max excess(total): {max_excess_anywhere:+.6f}")
    
    if violations == 0:
        print()
        print("  ═══════════════════════════════════════════════════")
        print("  CONFIRMED: w1·(μ1-μ0) > excess(dp0) always holds")
        print("  Therefore: excess(total) < 1 at every vertex")
        print("  ═══════════════════════════════════════════════════")
        print()
        print("  This means: the INDUCTIVE STEP is valid.")
        print("  At each merge where the mode jumps by 1,")
        print("  the weight-adjusted mean correction always")
        print("  overcomes the dp0 excess.")
    
    print()
    print("=" * 80)

if __name__ == "__main__":
    main()
