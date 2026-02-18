#!/usr/bin/env python3
"""Test whether mode(I_u·I_v) and mode(x·R_u·R_v) differ by at most 1.

If they do, then I(T') = I_u·I_v + x·R_u·R_v is automatically unimodal
(sum of two unimodal sequences whose modes differ by ≤ 1).

This avoids the C1+C2 conditional framework entirely.

Also tracks mode(R_u·R_v) vs d(I) to understand the relationship.
"""
import subprocess
import sys
import time
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
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 19
    print(f"Testing mode gap between I_u·I_v and x·R_u·R_v, n up to {max_n}", flush=True)
    print("If |mode(I_u·I_v) - mode(x·RR)| <= 1 always, unimodality is automatic.", flush=True)
    print("=" * 90, flush=True)

    t0 = time.time()
    total = 0
    gap_dist = Counter()  # mode(I_u·I_v) - mode(xRR)
    gap_too_big = 0  # |gap| >= 2
    max_gap = 0

    # Also track: does I_u·I_v ascend at k < d?
    ii_not_ascending = 0
    # And: does xRR ascend at k < d?
    xrr_not_ascending = 0

    # Track mode relationships
    d_vs_modeII = Counter()  # d(I) - mode(I_u·I_v)
    d_vs_modeRR = Counter()  # d(I) - mode(RR)

    for n in range(3, max_n + 1):
        tn = time.time()
        n_edges = 0
        n_gap_big = 0
        n_ii_notasc = 0
        n_xrr_notasc = 0

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
                total += 1
                side_u, side_v = split_at_edge(nn, adj, u, v)
                P_u, R_u = rooted_is_poly(adj, side_u, u)
                P_v, R_v = rooted_is_poly(adj, side_v, v)

                I_u = _polyadd(P_u, R_u)
                I_v = _polyadd(P_v, R_v)

                II = _polymul(I_u, I_v)
                RR = _polymul(R_u, R_v)
                xRR = [0] + RR

                m_II = first_descent(II)
                m_xRR = first_descent(xRR)
                m_RR = first_descent(RR)

                gap = m_II - m_xRR
                gap_dist[gap] += 1
                if abs(gap) > max_gap:
                    max_gap = abs(gap)
                if abs(gap) >= 2:
                    n_gap_big += 1
                    gap_too_big += 1

                d_vs_modeII[d - m_II] += 1
                d_vs_modeRR[d - m_RR] += 1

                # Check ascending at k < d
                for k in range(d):
                    if k + 1 < len(II) and II[k + 1] < II[k]:
                        n_ii_notasc += 1
                        ii_not_ascending += 1
                        break

                for k in range(d):
                    if k + 1 < len(xRR) and xRR[k + 1] < xRR[k]:
                        n_xrr_notasc += 1
                        xrr_not_ascending += 1
                        break

        proc.wait()
        elapsed = time.time() - tn
        print(f"n={n:2d}: edges={n_edges:>10,}  |gap|>=2={n_gap_big}  "
              f"II_desc<d={n_ii_notasc}  xRR_desc<d={n_xrr_notasc}  ({elapsed:.1f}s)",
              flush=True)

    total_time = time.time() - t0
    print("=" * 90, flush=True)
    print(f"Total: {total:,} edges in {total_time:.1f}s", flush=True)
    print(flush=True)

    print(f"|mode(I_u·I_v) - mode(x·RR)| >= 2:  "
          f"{'NONE' if gap_too_big == 0 else f'{gap_too_big} ({100*gap_too_big/total:.3f}%)'}",
          flush=True)
    print(f"Max |gap|: {max_gap}", flush=True)
    print(flush=True)

    print("mode(I_u·I_v) - mode(x·RR) distribution:", flush=True)
    for g in sorted(gap_dist.keys()):
        count = gap_dist[g]
        pct = 100 * count / total
        print(f"  {g:+d}: {count:>10,} ({pct:.2f}%)", flush=True)

    print(flush=True)
    print(f"I_u·I_v descends before d(I):  "
          f"{'NEVER' if ii_not_ascending == 0 else f'{ii_not_ascending} ({100*ii_not_ascending/total:.3f}%)'}",
          flush=True)
    print(f"x·RR descends before d(I):     "
          f"{'NEVER' if xrr_not_ascending == 0 else f'{xrr_not_ascending} ({100*xrr_not_ascending/total:.3f}%)'}",
          flush=True)

    print(flush=True)
    print("d(I) - mode(I_u·I_v) distribution:", flush=True)
    for g in sorted(d_vs_modeII.keys()):
        count = d_vs_modeII[g]
        pct = 100 * count / total
        print(f"  {g:+d}: {count:>10,} ({pct:.2f}%)", flush=True)

    print(flush=True)
    print("d(I) - mode(RR) distribution:", flush=True)
    for g in sorted(d_vs_modeRR.keys()):
        count = d_vs_modeRR[g]
        pct = 100 * count / total
        print(f"  {g:+d}: {count:>10,} ({pct:.2f}%)", flush=True)


if __name__ == "__main__":
    main()
