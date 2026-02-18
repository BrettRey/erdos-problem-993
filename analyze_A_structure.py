#!/usr/bin/env python3
"""Analyze the structure of A = P_u*P_v + x*R_u*R_v.

Key questions:
1. What are the modes of P_u*P_v and x*R_u*R_v individually?
2. Is the sum always LC because the modes are "aligned"?
3. Does some interlacing or dominance relationship hold between the two terms?
4. Can we express A as the IS polynomial of some auxiliary graph?

A = (IS with u,v both in S) + x*(IS with w in S, u,v out)
The first term counts IS of T-uv that include both u,v.
The second counts IS of T-{u,v} scaled by x (for w).

Note: P_u*P_v has minimum degree 2 (since P_u, P_v start at x).
And x*R_u*R_v has minimum degree 1 (since R_u, R_v start at constant 1).
So the low-degree coefficients come from x*R_u*R_v and the high-degree from P_u*P_v.
"""
import subprocess
import sys
from collections import Counter, deque

from indpoly import _polyadd, _polymul, independence_poly, is_log_concave, is_unimodal

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
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    print(f"Analyzing A = PP + xRR structure, n up to {max_n}", flush=True)
    print("=" * 78, flush=True)

    total = 0

    # Track: mode(PP) vs mode(xRR) vs mode(A) vs d
    mode_pp_vs_xrr = Counter()  # (mode(PP) - d, mode(xRR) - d)
    mode_pp_vs_a = Counter()  # mode(PP) - mode(A)
    # Track: is PP always LC? is xRR always LC?
    pp_lc_fail = 0
    xrr_lc_fail = 0
    # Track: does PP dominate xRR for large k?
    pp_dom_xrr_fail = 0
    # Track: ratio PP/xRR at mode(A)
    ratios_at_mode = []
    # Key: does PP[k] / xRR[k] increase or decrease?
    ratio_mono_fail = 0

    # Track: A as IS polynomial of a graph?
    # A counts IS of T-uv with {u,v} in S, plus IS of T-{u,v} with new vertex
    # The first type is IS of the forest T-edge(uv) restricted to include both u,v
    # Can we combine these into a single graph's IS polynomial?

    for n in range(3, max_n + 1):
        cmd = [GENG, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        n_edges = 0

        for line in proc.stdout:
            g6 = line.decode().strip()
            nn, adj = parse_graph6(g6)
            I_T = independence_poly(nn, adj)
            d = first_descent(I_T)

            edges = [(u, v) for u in range(nn) for v in adj[u] if u < v]
            for u, v in edges:
                n_edges += 1
                total += 1
                side_u, side_v = split_at_edge(nn, adj, u, v)
                P_u, R_u = rooted_is_poly(adj, side_u, u)
                P_v, R_v = rooted_is_poly(adj, side_v, v)

                PP = _polymul(P_u, P_v)
                RR = _polymul(R_u, R_v)
                xRR = [0] + RR
                A = _polyadd(PP, xRR)

                d_pp = first_descent(PP)
                d_xrr = first_descent(xRR)
                d_A = first_descent(A)

                mode_pp_vs_xrr[(d_pp - d, d_xrr - d)] += 1
                mode_pp_vs_a[d_pp - d_A] += 1

                if not is_log_concave(PP):
                    pp_lc_fail += 1
                if not is_log_concave(xRR):
                    xrr_lc_fail += 1

                # Check if PP dominates xRR in tail (k >= mode(A))
                for k in range(d_A, max(len(PP), len(xRR))):
                    ppk = PP[k] if k < len(PP) else 0
                    xrrk = xRR[k] if k < len(xRR) else 0
                    if xrrk > ppk:
                        pp_dom_xrr_fail += 1
                        break

                # Ratio at mode
                ppk = PP[d_A] if d_A < len(PP) else 0
                xrrk = xRR[d_A] if d_A < len(xRR) else 0
                if xrrk > 0:
                    ratios_at_mode.append(ppk / xrrk)

        proc.wait()
        print(f"n={n:2d}: edges={n_edges:>8,}", flush=True)

    print()
    print(f"Total: {total:,} edges", flush=True)
    print()
    print(f"PP (P_u*P_v) always LC: {'YES' if pp_lc_fail == 0 else f'NO ({pp_lc_fail} fails)'}",
          flush=True)
    print(f"xRR (x*R_u*R_v) always LC: {'YES' if xrr_lc_fail == 0 else f'NO ({xrr_lc_fail} fails)'}",
          flush=True)
    print()

    # Mode alignment
    print("mode(PP) - d and mode(xRR) - d distributions:", flush=True)
    pp_diffs = Counter()
    xrr_diffs = Counter()
    for (pp_d, xrr_d), count in mode_pp_vs_xrr.items():
        pp_diffs[pp_d] += count
        xrr_diffs[xrr_d] += count

    print("  mode(PP) - d:", flush=True)
    for diff in sorted(pp_diffs.keys()):
        print(f"    {diff:+d}: {pp_diffs[diff]:>10,} ({100*pp_diffs[diff]/total:.1f}%)", flush=True)
    print("  mode(xRR) - d:", flush=True)
    for diff in sorted(xrr_diffs.keys()):
        print(f"    {diff:+d}: {xrr_diffs[diff]:>10,} ({100*xrr_diffs[diff]/total:.1f}%)", flush=True)
    print()

    print("mode(PP) - mode(A):", flush=True)
    for diff in sorted(mode_pp_vs_a.keys()):
        count = mode_pp_vs_a[diff]
        print(f"    {diff:+d}: {count:>10,} ({100*count/total:.1f}%)", flush=True)
    print()

    print(f"PP dominates xRR in tail (k >= mode(A)): "
          f"{'ALWAYS' if pp_dom_xrr_fail == 0 else f'FAILS ({pp_dom_xrr_fail})'}",
          flush=True)
    print()

    if ratios_at_mode:
        ratios_at_mode.sort()
        print(f"PP/xRR ratio at mode(A):", flush=True)
        print(f"  min: {ratios_at_mode[0]:.4f}", flush=True)
        print(f"  median: {ratios_at_mode[len(ratios_at_mode)//2]:.4f}", flush=True)
        print(f"  max: {ratios_at_mode[-1]:.4f}", flush=True)


if __name__ == "__main__":
    main()
