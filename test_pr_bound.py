#!/usr/bin/env python3
"""Test the P_u[i] <= R_u[i-1] bound and PP[k] <= RR[k-2] bound.

Also test whether a stronger bound holds:
  PP[k] - PP[k+1] <= RR[k] - RR[k-1]  (for k >= mode(PP))
  (this would prove C2 if true)

And test the condition:
  RR[k] >= RR[k-1] + PP[k]  (sufficient for A ascending at k)
"""
import subprocess
import sys
import time
from collections import deque

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
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    print(f"Testing P<=R bound and derived conditions, n up to {max_n}", flush=True)
    print("=" * 80, flush=True)

    t0 = time.time()
    total = 0
    pr_fail = 0  # P_u[i] > R_u[i-1]
    pprr_fail = 0  # PP[k] > RR[k-2]
    c2_sufficient_fail = 0  # RR[k] < RR[k-1] + PP[k] for k < d
    delta_bound_fail = 0  # PP[k]-PP[k+1] > RR[k]-RR[k-1] for k >= mode(PP) and k < d

    # Track: in how many cases is PP descending at some k < d?
    pp_desc_before_d = 0
    # Track: max PP/RR ratio at k-2
    max_pp_rr_ratio = 0

    for n in range(3, max_n + 1):
        tn = time.time()
        n_edges = 0
        n_c2_fail = 0
        n_delta_fail = 0

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

                # Test P_u[i] <= R_u[i-1]
                for i in range(1, len(P_u)):
                    pu_i = P_u[i]
                    ru_im1 = R_u[i - 1] if i - 1 < len(R_u) else 0
                    if pu_i > ru_im1:
                        pr_fail += 1
                        break
                for j in range(1, len(P_v)):
                    pv_j = P_v[j]
                    rv_jm1 = R_v[j - 1] if j - 1 < len(R_v) else 0
                    if pv_j > rv_jm1:
                        pr_fail += 1
                        break

                PP = _polymul(P_u, P_v)
                RR = _polymul(R_u, R_v)

                # Test PP[k] <= RR[k-2]
                for k in range(2, len(PP)):
                    ppk = PP[k]
                    rrk2 = RR[k - 2] if k - 2 < len(RR) else 0
                    if ppk > rrk2:
                        pprr_fail += 1
                        break
                    if rrk2 > 0:
                        ratio = ppk / rrk2
                        if ratio > max_pp_rr_ratio:
                            max_pp_rr_ratio = ratio

                # Check if PP descends before d
                d_PP = first_descent(PP)
                if d_PP < d:
                    pp_desc_before_d += 1

                    # For positions k where PP is descending and k < d:
                    for k in range(d_PP, d):
                        ppk = PP[k] if k < len(PP) else 0
                        ppk1 = PP[k + 1] if k + 1 < len(PP) else 0
                        pp_drop = ppk - ppk1  # > 0 since PP descending

                        rrk = RR[k] if k < len(RR) else 0
                        rrkm1 = RR[k - 1] if k - 1 >= 0 and k - 1 < len(RR) else 0
                        rr_rise = rrk - rrkm1  # xRR rise at this position

                        if pp_drop > rr_rise:
                            n_delta_fail += 1
                            delta_bound_fail += 1

                # Test sufficient condition: RR[k] >= RR[k-1] + PP[k] for k < d
                for k in range(1, d):
                    rrk = RR[k] if k < len(RR) else 0
                    rrkm1 = RR[k - 1] if k - 1 >= 0 and k - 1 < len(RR) else 0
                    ppk = PP[k] if k < len(PP) else 0
                    if rrk < rrkm1 + ppk:
                        n_c2_fail += 1
                        c2_sufficient_fail += 1
                        break

        proc.wait()
        elapsed = time.time() - tn
        print(f"n={n:2d}: edges={n_edges:>10,}  "
              f"suff_fail={n_c2_fail}  delta_fail={n_delta_fail}  ({elapsed:.1f}s)",
              flush=True)

    total_time = time.time() - t0
    print("=" * 80, flush=True)
    print(f"Total: {total:,} edges in {total_time:.1f}s", flush=True)
    print(flush=True)
    print(f"P_u[i] <= R_u[i-1]:              {'HOLDS' if pr_fail == 0 else f'FAILS ({pr_fail})'}",
          flush=True)
    print(f"PP[k] <= RR[k-2]:                {'HOLDS' if pprr_fail == 0 else f'FAILS ({pprr_fail})'}",
          flush=True)
    print(f"Max PP[k]/RR[k-2] ratio:         {max_pp_rr_ratio:.6f}", flush=True)
    print(f"PP descends before d(I):         {pp_desc_before_d}/{total} "
          f"({100*pp_desc_before_d/max(total,1):.1f}%)", flush=True)
    print(f"Delta bound (PP drop <= RR rise): {'HOLDS' if delta_bound_fail == 0 else f'FAILS ({delta_bound_fail})'}",
          flush=True)
    print(f"Sufficient (RR[k]>=RR[k-1]+PP):  {'HOLDS' if c2_sufficient_fail == 0 else f'FAILS ({c2_sufficient_fail})'}",
          flush=True)


if __name__ == "__main__":
    main()
