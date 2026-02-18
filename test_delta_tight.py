#!/usr/bin/env python3
"""Analyze tight cases of the delta bound: edges where A[d] = A[d-1] (margin 0).

For each edge uv, the delta bound says A[k+1] >= A[k] for k < d.
The tightest case is at k = d-1: A[d] >= A[d-1].

We characterize edges where A[d] - A[d-1] is smallest relative to A[d-1],
to understand what makes the bound tight and find structural patterns.

Also compute: for each edge, the minimum ratio A[k+1]/A[k] over k < d,
and track which position k achieves the minimum.
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
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    print(f"Analyzing tight delta bound cases, n up to {max_n}", flush=True)
    print("=" * 90, flush=True)

    t0 = time.time()
    total = 0
    tight_cases = []  # (min_ratio, min_k, n, g6, u, v, deg_u, deg_v)
    min_ratio_dist = Counter()  # floor(min_ratio * 100) / 100

    # Track: at which position k is the bound tightest?
    tightest_k_from_d = Counter()  # k relative to d: tightest at d-1, d-2, etc.

    for n in range(3, max_n + 1):
        tn = time.time()
        n_edges = 0
        n_min_ratio = float('inf')

        cmd = [GENG, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        for line in proc.stdout:
            g6 = line.decode().strip()
            nn, adj = parse_graph6(g6)
            I_T = independence_poly(nn, adj)
            d = first_descent(I_T)

            if d <= 0:
                continue

            edges = [(u, v) for u in range(nn) for v in adj[u] if u < v]
            for u, v in edges:
                n_edges += 1
                total += 1
                side_u, side_v = split_at_edge(nn, adj, u, v)
                P_u, R_u = rooted_is_poly(adj, side_u, u)
                P_v, R_v = rooted_is_poly(adj, side_v, v)

                PP = _polymul(P_u, P_v)
                RR = _polymul(R_u, R_v)

                # A[k] = PP[k] + RR[k-1]
                deg_A = max(len(PP), len(RR) + 1) - 1
                A = [0] * (deg_A + 1)
                for k in range(deg_A + 1):
                    ppk = PP[k] if k < len(PP) else 0
                    rrk1 = RR[k - 1] if 0 < k <= len(RR) else 0
                    A[k] = ppk + rrk1

                # Find minimum ratio A[k+1]/A[k] for k < d
                min_ratio = float('inf')
                min_k = -1
                for k in range(min(d, len(A) - 1)):
                    if A[k] > 0:
                        ratio = A[k + 1] / A[k]
                        if ratio < min_ratio:
                            min_ratio = ratio
                            min_k = k

                if min_ratio < n_min_ratio:
                    n_min_ratio = min_ratio

                if min_k >= 0:
                    tightest_k_from_d[d - 1 - min_k] += 1

                # Bucket the ratio
                if min_ratio < float('inf'):
                    bucket = int(min_ratio * 100) // 10  # 10% buckets
                    min_ratio_dist[bucket] += 1

                # Keep track of very tight cases
                if min_ratio < 1.05 and len(tight_cases) < 100:
                    tight_cases.append((
                        min_ratio, min_k, d, nn, g6, u, v,
                        len(adj[u]), len(adj[v])
                    ))

        proc.wait()
        elapsed = time.time() - tn
        print(f"n={n:2d}: edges={n_edges:>10,}  min_ratio={n_min_ratio:.6f}  ({elapsed:.1f}s)",
              flush=True)

    total_time = time.time() - t0
    print("=" * 90, flush=True)
    print(f"Total: {total:,} edges in {total_time:.1f}s", flush=True)
    print(flush=True)

    # Sort tight cases by ratio
    tight_cases.sort()
    print(f"Tightest {min(20, len(tight_cases))} cases:", flush=True)
    for ratio, k, d, nn, g6, u, v, du, dv in tight_cases[:20]:
        print(f"  ratio={ratio:.6f} at k={k}, d={d}, n={nn}, "
              f"edge=({u},{v}), degs=({du},{dv}), g6={g6}", flush=True)

    print(flush=True)
    print("Position of tightest ratio (distance from d-1):", flush=True)
    for dist in sorted(tightest_k_from_d.keys()):
        count = tightest_k_from_d[dist]
        pct = 100 * count / total
        print(f"  d-1-{dist} (k={{}}):" f" {count:>10,} ({pct:.2f}%)", flush=True)

    print(flush=True)
    print("Min ratio distribution (10% buckets):", flush=True)
    for bucket in sorted(min_ratio_dist.keys()):
        count = min_ratio_dist[bucket]
        pct = 100 * count / total
        lo = bucket * 10
        hi = lo + 10
        print(f"  [{lo:3d}%,{hi:3d}%): {count:>10,} ({pct:.2f}%)", flush=True)


if __name__ == "__main__":
    main()
