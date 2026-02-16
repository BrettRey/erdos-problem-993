#!/usr/bin/env python3
"""
Single-child gap=+1 analysis using the existing indpoly DP.
Key question: what is excess(dp₀) at gap=+1 merges?
"""

import subprocess
from graph6 import parse_graph6

def polymul(a, b):
    n = len(a) + len(b) - 1
    r = [0]*n
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            r[i+j] += ai*bj
    return r

def polyadd(a, b):
    n = max(len(a), len(b))
    r = [0]*n
    for i in range(len(a)): r[i] += a[i]
    for i in range(len(b)): r[i] += b[i]
    return r

def poly_mode(p):
    return max(range(len(p)), key=lambda k: p[k])

def poly_mean(p):
    s = sum(p)
    if s == 0: return 0
    return sum(k*p[k] for k in range(len(p))) / s

def poly_excess(p):
    return poly_mode(p) - poly_mean(p)

def dp_tree(adj_list, root):
    n = len(adj_list)
    dp0 = [None]*n
    dp1 = [None]*n
    total = [None]*n
    
    parent = [-1]*n
    order = []
    visited = [False]*n
    stack = [root]
    visited[root] = True
    while stack:
        v = stack.pop()
        order.append(v)
        for u in adj_list[v]:
            if not visited[u]:
                visited[u] = True
                parent[u] = v
                stack.append(u)
    
    for v in reversed(order):
        children = [u for u in adj_list[v] if parent[u] == v]
        if not children:
            dp0[v] = [1]
            dp1[v] = [0, 1]
            total[v] = [1, 1]
        else:
            d0 = [1]
            d1 = [1]
            for c in children:
                d0 = polymul(d0, total[c])
                d1 = polymul(d1, dp0[c])
            dp0[v] = d0
            dp1[v] = [0] + d1
            total[v] = polyadd(d0, [0] + d1)
    
    return dp0, dp1, total, parent

def main():
    print("=" * 70)
    print("  GAP=+1 MERGE ANALYSIS: Excess of dp₀ by #children")
    print("=" * 70)
    print()
    
    all_single = []
    all_multi = {}
    
    for n in range(3, 19):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        for line in lines:
            tn, adj = parse_graph6(line)
            
            for root in range(tn):
                dp0_v, dp1_v, total_v, par = dp_tree(adj, root)
                
                children = [u for u in adj[root] if par[u] == root]
                m_dp0 = poly_mode(dp0_v[root])
                m_total = poly_mode(total_v[root])
                gap = m_total - m_dp0
                
                if gap == 1:
                    ex_dp0 = poly_excess(dp0_v[root])
                    nc = len(children)
                    
                    # Compute correction and slack
                    S_A = sum(dp0_v[root])
                    S_B = sum(dp1_v[root])
                    w = S_B / (S_A + S_B)
                    mu_A = poly_mean(dp0_v[root])
                    mu_B = poly_mean(dp1_v[root])
                    correction = w * (mu_B - mu_A)
                    slack = correction - ex_dp0
                    
                    if nc == 1:
                        all_single.append((n, ex_dp0, correction, slack))
                    else:
                        all_multi.setdefault(nc, []).append((n, ex_dp0, correction, slack))
        
    # Report
    print(f"  SINGLE-CHILD gap=+1 merges: {len(all_single)}")
    if all_single:
        exs = [x[1] for x in all_single]
        slacks = [x[3] for x in all_single]
        print(f"    excess(dp₀): range [{min(exs):+.4f}, {max(exs):+.4f}]")
        print(f"    correction slack: range [{min(slacks):+.4f}, {max(slacks):+.4f}]")
        neg = sum(1 for e in exs if e <= 0)
        print(f"    excess(dp₀) ≤ 0: {neg}/{len(exs)} ({100*neg/len(exs):.1f}%)")
    print()
    
    for nc in sorted(all_multi):
        data = all_multi[nc]
        exs = [x[1] for x in data]
        slacks = [x[3] for x in data]
        print(f"  {nc}-CHILDREN gap=+1: {len(data)} merges")
        print(f"    excess(dp₀): range [{min(exs):+.4f}, {max(exs):+.4f}]")
        print(f"    slack: range [{min(slacks):+.4f}, {max(slacks):+.4f}]")
        neg = sum(1 for e in exs if e <= 0)
        print(f"    excess(dp₀) ≤ 0: {neg}/{len(exs)} ({100*neg/len(exs):.1f}%)")
    print()
    
    # KEY OBSERVATION: For single-child, excess(dp₀) = excess(total[c]).
    # The gap=+1 condition at the parent SELECTS for specific excess values
    # at the child. What IS this selection?
    
    # At gap=+1 with single child c:
    # total[v] = total[c] + x·dp₀[c]
    # mode(total[v]) = mode(total[c]) + 1
    # This requires dp₀[c] to be increasing at position mode(total[c]).
    # i.e., dp₀[c][mode(total[c])] > dp₀[c][mode(total[c])-1]}
    # i.e., mode(dp₀[c]) > mode(total[c])
    
    print("  STRUCTURAL CONSTRAINT at single-child gap=+1:")
    print("  mode(dp₀[c]) vs mode(total[c]):")
    
    mode_diffs = []
    for n_val in range(3, 16):
        cmd = f"/opt/homebrew/bin/geng {n_val} {n_val-1}:{n_val-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        for line in lines:
            tn, adj_l = parse_graph6(line)
            for root in range(tn):
                dp0_v, dp1_v, total_v, par = dp_tree(adj_l, root)
                children = [u for u in adj_l[root] if par[u] == root]
                
                if len(children) == 1:
                    m_dp0 = poly_mode(dp0_v[root])
                    m_total = poly_mode(total_v[root])
                    if m_total - m_dp0 == 1:
                        c = children[0]
                        m_dp0c = poly_mode(dp0_v[c])
                        m_totalc = poly_mode(total_v[c])
                        mode_diffs.append(m_dp0c - m_totalc)
    
    from collections import Counter
    ct = Counter(mode_diffs)
    print(f"    mode(dp₀[c]) - mode(total[c]) distribution:")
    for d in sorted(ct):
        print(f"      diff = {d:+d}: {ct[d]} ({100*ct[d]/len(mode_diffs):.1f}%)")
    
    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
