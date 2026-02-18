#!/usr/bin/env python3
"""Test the mixed LC condition for I(T') = a + b where a = I_u*I_v, b = x*R_u*R_v.

The sum of two LC sequences a, b is LC iff for all k:
  2*a_k*b_k >= a_{k-1}*b_{k+1} + b_{k-1}*a_{k+1}

This is equivalent to: the "cross-Turán" inequality holds for the pair (a, b).

If this holds, then I(T') being LC follows from I_u, I_v, R_u, R_v being LC
(which they are as IS polynomials of trees/forests).

Also test the stronger condition where we decompose differently:
I(T') = I(T) + A where A = PP + xRR.
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
    print(f"Testing mixed LC condition for I(T') = IuIv + xRuRv, n up to {max_n}",
          flush=True)
    print("=" * 78, flush=True)

    t0 = time.time()
    total = 0
    mixed_lc_fail = 0
    mixed_lc_fail_details = []
    # Also track the margin: how much slack does the inequality have?
    min_margin_ratio = float('inf')

    for n in range(3, max_n + 1):
        tn = time.time()
        n_edges = 0
        n_fail = 0

        cmd = [GENG, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        for line in proc.stdout:
            g6 = line.decode().strip()
            nn, adj = parse_graph6(g6)

            edges = [(u, v) for u in range(nn) for v in adj[u] if u < v]
            for u, v in edges:
                n_edges += 1
                total += 1
                side_u, side_v = split_at_edge(nn, adj, u, v)
                P_u, R_u = rooted_is_poly(adj, side_u, u)
                P_v, R_v = rooted_is_poly(adj, side_v, v)

                # a = I_u * I_v (forest IS poly)
                I_u = _polyadd(P_u, R_u)
                I_v = _polyadd(P_v, R_v)
                a = _polymul(I_u, I_v)

                # b = x * R_u * R_v
                RR = _polymul(R_u, R_v)
                b = [0] + RR

                # Check mixed LC: 2*a[k]*b[k] >= a[k-1]*b[k+1] + b[k-1]*a[k+1]
                max_len = max(len(a), len(b))
                failed = False
                for k in range(1, max_len - 1):
                    ak = a[k] if k < len(a) else 0
                    bk = b[k] if k < len(b) else 0
                    ak1 = a[k-1] if k-1 < len(a) else 0
                    bk1 = b[k+1] if k+1 < len(b) else 0
                    bkm1 = b[k-1] if k-1 < len(b) else 0
                    akp1 = a[k+1] if k+1 < len(a) else 0

                    lhs = 2 * ak * bk
                    rhs = ak1 * bk1 + bkm1 * akp1

                    if lhs < rhs:
                        if not failed:
                            failed = True
                            n_fail += 1
                            mixed_lc_fail += 1
                            if mixed_lc_fail <= 5:
                                mixed_lc_fail_details.append(
                                    (nn, u, v, k, lhs, rhs, g6))
                    elif rhs > 0:
                        ratio = lhs / rhs
                        if ratio < min_margin_ratio:
                            min_margin_ratio = ratio

        proc.wait()
        elapsed = time.time() - tn
        print(f"n={n:2d}: edges={n_edges:>10,}  mixed_LC_fail={n_fail}  ({elapsed:.1f}s)",
              flush=True)

    total_time = time.time() - t0
    print("=" * 78, flush=True)
    print(f"Total: {total:,} edges, {total_time:.1f}s", flush=True)
    print(flush=True)

    if mixed_lc_fail == 0:
        print("MIXED LC CONDITION HOLDS for all edges!", flush=True)
        print(f"Min margin ratio (LHS/RHS): {min_margin_ratio:.6f}", flush=True)
        print(flush=True)
        print("This PROVES: I(T') is LC whenever I_u, I_v are LC", flush=True)
        print("(since products of LC polys are LC, and the mixed condition", flush=True)
        print("guarantees the sum preserves LC).", flush=True)
    else:
        print(f"MIXED LC FAILS in {mixed_lc_fail} cases!", flush=True)
        for nn, u, v, k, lhs, rhs, g6 in mixed_lc_fail_details:
            print(f"  n={nn} edge=({u},{v}) k={k}: 2ab={lhs} < cross={rhs}  tree={g6}",
                  flush=True)


if __name__ == "__main__":
    main()
