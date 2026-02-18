#!/usr/bin/env python3
"""Gap analysis for mode(A) = d+2 cases with CORRECT formula.

A = P_u*P_v + x*R_u*R_v (CORRECT)

For the 0.09% of cases where first_descent(A) = d+2:
- Sign pattern of Δ(I+A) at d and d+1?
- How tight is the combined tail at d+1?
- Structural patterns?
"""
import subprocess
import sys
import time
from collections import Counter, deque

from indpoly import _polyadd, _polymul, independence_poly

GENG = "/opt/homebrew/bin/geng"


def parse_graph6(s):
    s = s.strip()
    data = [c - 63 for c in s.encode("ascii")]
    n = data[0]
    idx = 1
    adj = [[] for _ in range(n)]
    bit = 5
    word = data[idx] if idx < len(data) else 0
    for j in range(n):
        for i in range(j):
            if word & (1 << bit):
                adj[i].append(j)
                adj[j].append(i)
            if bit == 0:
                bit = 5
                idx += 1
                word = data[idx] if idx < len(data) else 0
            else:
                bit -= 1
    return n, adj


def first_descent(seq):
    for k in range(len(seq) - 1):
        if seq[k] > seq[k + 1]:
            return k
    return len(seq) - 1


def split_at_edge(n, adj, u, v):
    A = set()
    queue = deque([u])
    A.add(u)
    while queue:
        x = queue.popleft()
        for y in adj[x]:
            if y not in A and y != v:
                A.add(y)
                queue.append(y)
    return A, set(range(n)) - A


def rooted_is_poly(adj, vertices, root):
    vset = set(vertices)
    parent = {root: -1}
    order = [root]
    queue = deque([root])
    while queue:
        x = queue.popleft()
        for y in adj[x]:
            if y in vset and y not in parent:
                parent[y] = x
                queue.append(y)
                order.append(y)
    dp_in = {}
    dp_out = {}
    for v in reversed(order):
        children = [y for y in adj[v] if y in vset and parent.get(y) == v]
        if not children:
            dp_in[v] = [0, 1]
            dp_out[v] = [1]
        else:
            prod_out = [1]
            for c in children:
                prod_out = _polymul(prod_out, dp_out[c])
            dp_in[v] = [0] + prod_out
            prod_both = [1]
            for c in children:
                prod_both = _polymul(prod_both, _polyadd(dp_in[c], dp_out[c]))
            dp_out[v] = prod_both
    return dp_in[root], dp_out[root]


def main():
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 19
    print(f"CORRECTED gap analysis, n up to {max_n}", flush=True)
    print("=" * 78, flush=True)

    t0 = time.time()
    gap_cases = []
    total_edges = 0
    mode_dist = Counter()

    for n in range(3, max_n + 1):
        tn = time.time()
        n_edges = 0
        n_gap = 0

        cmd = [GENG, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        for line in proc.stdout:
            g6 = line.decode().strip()
            nn, adj = parse_graph6(g6)
            I_T = independence_poly(nn, adj)
            d = first_descent(I_T)

            edges = [(u, v) for u in range(nn) for v in adj[u] if u < v]
            for u, v in edges:
                n_edges += 1
                side_u, side_v = split_at_edge(nn, adj, u, v)
                P_u, R_u = rooted_is_poly(adj, side_u, u)
                P_v, R_v = rooted_is_poly(adj, side_v, v)

                # CORRECT formula
                A = _polyadd(_polymul(P_u, P_v), [0] + _polymul(R_u, R_v))
                d_A = first_descent(A)
                diff = d_A - d
                mode_dist[diff] += 1

                if diff >= 2:
                    n_gap += 1
                    I_Tp = _polyadd(I_T, A)

                    # Signs at d and d+1
                    delta_sum_d = I_Tp[d + 1] - I_Tp[d] if d + 1 < len(I_Tp) else -I_Tp[d]
                    delta_sum_d1 = (I_Tp[d + 2] - I_Tp[d + 1]) if d + 2 < len(I_Tp) else -I_Tp[d + 1]

                    # Margins: positive means I(T') is nonincreasing
                    margin_d1 = I_Tp[d + 1] - I_Tp[d + 2] if d + 2 < len(I_Tp) else I_Tp[d + 1]

                    # Drops and rises
                    I_drop_d1 = I_T[d + 1] - I_T[d + 2] if d + 2 < len(I_T) else I_T[d + 1]
                    A_rise_d1 = (A[d + 2] if d + 2 < len(A) else 0) - (A[d + 1] if d + 1 < len(A) else 0)

                    deg_u = len(adj[u])
                    deg_v = len(adj[v])

                    gap_cases.append({
                        'g6': g6, 'n': nn, 'u': u, 'v': v,
                        'd': d, 'd_A': d_A, 'diff': diff,
                        'deg_u': deg_u, 'deg_v': deg_v,
                        'su': len(side_u), 'sv': len(side_v),
                        'delta_sum_d': delta_sum_d,
                        'delta_sum_d1': delta_sum_d1,
                        'margin_d1': margin_d1,
                        'I_drop_d1': I_drop_d1,
                        'A_rise_d1': A_rise_d1,
                        'I_d': I_T[d], 'I_d1': I_T[d + 1] if d + 1 < len(I_T) else 0,
                        'I_d2': I_T[d + 2] if d + 2 < len(I_T) else 0,
                        'A_d': A[d] if d < len(A) else 0,
                        'A_d1': A[d + 1] if d + 1 < len(A) else 0,
                        'A_d2': A[d + 2] if d + 2 < len(A) else 0,
                    })

        proc.wait()
        elapsed = time.time() - tn
        total_edges += n_edges
        print(f"n={n:2d}: edges={n_edges:>10,}  gap_cases(d_A>=d+2)={n_gap:4d}  ({elapsed:.1f}s)",
              flush=True)

    total_time = time.time() - t0
    print("=" * 78, flush=True)
    print(f"Total: {total_edges:,} edges, {len(gap_cases)} gap cases, {total_time:.1f}s", flush=True)
    print()

    if not gap_cases:
        print("No gap cases!", flush=True)
        return

    # Sign patterns
    sign_patterns = Counter()
    for c in gap_cases:
        s1 = '+' if c['delta_sum_d'] > 0 else ('-' if c['delta_sum_d'] < 0 else '0')
        s2 = '+' if c['delta_sum_d1'] > 0 else ('-' if c['delta_sum_d1'] < 0 else '0')
        sign_patterns[(s1, s2)] += 1

    print(f"Sign patterns (?₁, ?₂) at (d, d+1) [{len(gap_cases)} gap cases]:", flush=True)
    for (s1, s2), count in sign_patterns.most_common():
        pct = 100.0 * count / len(gap_cases)
        print(f"  ({s1}, {s2}): {count:>6} ({pct:.1f}%)", flush=True)
    print()

    # Valley risk
    valley_risk = sum(1 for c in gap_cases
                      if c['delta_sum_d'] < 0 and c['delta_sum_d1'] > 0)
    print(f"VALLEY RISK (?₁<0 and ?₂>0): {valley_risk}/{len(gap_cases)}", flush=True)
    print()

    # Margin analysis
    margins = [c['margin_d1'] for c in gap_cases]
    margins.sort()
    print(f"Combined tail margin at d+1 (positive = safe):", flush=True)
    print(f"  min: {margins[0]}", flush=True)
    print(f"  p01: {margins[int(0.01 * len(margins))]}", flush=True)
    print(f"  p10: {margins[int(0.10 * len(margins))]}", flush=True)
    print(f"  median: {margins[len(margins) // 2]}", flush=True)
    print(f"  max: {margins[-1]}", flush=True)
    print()

    # Ratio of A_rise to I_drop at d+1
    ratios = []
    for c in gap_cases:
        if c['I_drop_d1'] > 0 and c['A_rise_d1'] > 0:
            ratios.append(c['A_rise_d1'] / c['I_drop_d1'])
    if ratios:
        ratios.sort()
        print(f"A_rise/I_drop ratio at d+1 (must be < 1 for safety):", flush=True)
        print(f"  max: {ratios[-1]:.6f}", flush=True)
        print(f"  p99: {ratios[int(0.99 * len(ratios))]:.6f}", flush=True)
        print(f"  p90: {ratios[int(0.90 * len(ratios))]:.6f}", flush=True)
        print(f"  median: {ratios[len(ratios) // 2]:.6f}", flush=True)
        print()

    # Tightest cases
    tightest = sorted(gap_cases, key=lambda c: c['margin_d1'])[:10]
    print("10 tightest cases:", flush=True)
    for c in tightest:
        ratio = c['A_rise_d1'] / c['I_drop_d1'] if c['I_drop_d1'] > 0 else float('inf')
        print(f"  n={c['n']} e=({c['u']},{c['v']}) d={c['d']} d_A={c['d_A']} "
              f"deg=({c['deg_u']},{c['deg_v']}) sides=({c['su']},{c['sv']}) "
              f"margin={c['margin_d1']} ratio={ratio:.4f}", flush=True)
        print(f"    I[d..d+2]={c['I_d']},{c['I_d1']},{c['I_d2']}  "
              f"A[d..d+2]={c['A_d']},{c['A_d1']},{c['A_d2']}  "
              f"I_drop={c['I_drop_d1']} A_rise={c['A_rise_d1']}", flush=True)
    print()

    # Structural patterns: degree pairs
    print("Edge endpoint degrees in gap cases:", flush=True)
    deg_pairs = Counter()
    for c in gap_cases:
        du, dv = sorted([c['deg_u'], c['deg_v']])
        deg_pairs[(du, dv)] += 1
    for (du, dv), count in deg_pairs.most_common(10):
        print(f"  ({du},{dv}): {count}", flush=True)
    print()

    # Per-diff breakdown
    for d_val in sorted(set(c['diff'] for c in gap_cases)):
        subset = [c for c in gap_cases if c['diff'] == d_val]
        margins_sub = [c['margin_d1'] for c in subset]
        print(f"diff={d_val}: {len(subset)} cases, "
              f"min_margin={min(margins_sub)}, max_margin={max(margins_sub)}", flush=True)


if __name__ == "__main__":
    main()
