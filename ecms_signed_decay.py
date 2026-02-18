#!/usr/bin/env python3
"""Test whether signed per-distance sums decay geometrically.

The unsigned |ΔP| sums at each distance decay at ~0.38 per step, but
the total unsigned sum is ~0.6 > 1/2. So proving |remote| < 1/2
requires using sign structure.

Key question: for each edge, does the signed sum at distance d
(i.e., Σ_{w at dist d} ΔP(w)) decay faster than the unsigned sum?
If so, the signed total might be provably < 1/2.

Also check: do signed sums at consecutive distances alternate in sign?
If so, the partial sums oscillate and stay bounded.
"""

import json
import os
import time
from collections import defaultdict

from trees import trees

MAX_N = 18


def cavity_messages(n, adj):
    if n <= 1:
        return {}
    parent = [-1] * n
    children = [[] for _ in range(n)]
    visited = [False] * n
    order = []
    visited[0] = True
    queue = [0]
    head = 0
    while head < len(queue):
        v = queue[head]; head += 1
        order.append(v)
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                parent[u] = v
                children[v].append(u)
                queue.append(u)
    msgs = {}
    for v in reversed(order):
        p = parent[v]
        if p == -1:
            continue
        prod = 1.0
        for c in children[v]:
            prod *= 1.0 / (1.0 + msgs[(c, v)])
        msgs[(v, p)] = prod
    for v in order:
        for c in children[v]:
            prod = 1.0
            if parent[v] != -1:
                prod *= 1.0 / (1.0 + msgs[(parent[v], v)])
            for c2 in children[v]:
                if c2 != c:
                    prod *= 1.0 / (1.0 + msgs[(c2, v)])
            msgs[(v, c)] = prod
    return msgs


def occupation_probs(n, adj, msgs=None):
    if msgs is None:
        msgs = cavity_messages(n, adj)
    P = [0.0] * n
    for v in range(n):
        Rv = 1.0
        for u in adj[v]:
            Rv *= 1.0 / (1.0 + msgs.get((u, v), 0.0))
        P[v] = Rv / (1.0 + Rv)
    return P


def bfs_distances(n, adj, sources):
    dist = [-1] * n
    queue = []
    for s in sources:
        dist[s] = 0
        queue.append(s)
    head = 0
    while head < len(queue):
        v = queue[head]; head += 1
        for u in adj[v]:
            if dist[u] == -1:
                dist[u] = dist[v] + 1
                queue.append(u)
    return dist


def contract_edge(n, adj, u, v):
    merged_neighbors = set()
    for w in adj[u]:
        if w != v:
            merged_neighbors.add(w)
    for w in adj[v]:
        if w != u:
            merged_neighbors.add(w)
    old_to_new = {}
    for i in range(n):
        if i == v:
            continue
        old_to_new[i] = i if i < v else i - 1
    n_new = n - 1
    adj_new = [[] for _ in range(n_new)]
    u_new = old_to_new[u]
    for w in sorted(merged_neighbors):
        w_new = old_to_new[w]
        adj_new[u_new].append(w_new)
        adj_new[w_new].append(u_new)
    for i in range(n):
        if i == u or i == v:
            continue
        i_new = old_to_new[i]
        for j in adj[i]:
            if j == u or j == v:
                continue
            j_new = old_to_new[j]
            if j_new not in adj_new[i_new]:
                adj_new[i_new].append(j_new)
    for i in range(n_new):
        adj_new[i].sort()
    return n_new, adj_new, old_to_new


def main():
    t0 = time.time()

    # For each distance d, track:
    # - max signed sum (positive)
    # - min signed sum (negative)
    # - distribution of signed sum magnitudes
    max_signed_by_dist = defaultdict(float)
    min_signed_by_dist = defaultdict(lambda: float('inf'))

    # Track partial sums: cumulative signed sum through distance d
    max_partial_by_dist = defaultdict(float)
    min_partial_by_dist = defaultdict(lambda: float('inf'))

    # Sign alternation: how often does sgn(S_d) ≠ sgn(S_{d+1})?
    alternations = defaultdict(int)
    same_sign_consecutive = defaultdict(int)
    total_consecutive = defaultdict(int)

    # Track worst cases for signed sums
    worst_partial_sums = []

    total_edges = 0

    for n in range(3, MAX_N + 1):
        tn = time.time()
        n_trees = 0
        n_edges = 0

        for _, adj in trees(n):
            n_trees += 1
            msgs = cavity_messages(n, adj)
            P_T = occupation_probs(n, adj, msgs)

            seen = set()
            for u in range(n):
                for v in adj[u]:
                    e_key = (min(u, v), max(u, v))
                    if e_key in seen:
                        continue
                    seen.add(e_key)
                    total_edges += 1
                    n_edges += 1

                    dist = bfs_distances(n, adj, [u, v])

                    nc, adjc, o2n = contract_edge(n, adj, u, v)
                    msgs_c = cavity_messages(nc, adjc)
                    P_Te = occupation_probs(nc, adjc, msgs_c)

                    # Compute signed sum at each distance
                    signed_by_d = defaultdict(float)
                    max_dist = 0
                    for w in range(n):
                        if w == u or w == v:
                            continue
                        d = dist[w]
                        w_new = o2n[w]
                        dP = P_T[w] - P_Te[w_new]
                        signed_by_d[d] += dP
                        if d > max_dist:
                            max_dist = d

                    # Track per-distance extremes
                    for d in range(1, max_dist + 1):
                        s = signed_by_d.get(d, 0)
                        if s > max_signed_by_dist[d]:
                            max_signed_by_dist[d] = s
                        if s < min_signed_by_dist[d]:
                            min_signed_by_dist[d] = s

                    # Partial sums and alternation
                    partial = 0
                    prev_sign = 0
                    for d in range(1, max_dist + 1):
                        s = signed_by_d.get(d, 0)
                        partial += s

                        if partial > max_partial_by_dist[d]:
                            max_partial_by_dist[d] = partial
                        if partial < min_partial_by_dist[d]:
                            min_partial_by_dist[d] = partial

                        curr_sign = 1 if s > 1e-15 else (-1 if s < -1e-15 else 0)
                        if prev_sign != 0 and curr_sign != 0:
                            total_consecutive[d] += 1
                            if curr_sign != prev_sign:
                                alternations[d] += 1
                            else:
                                same_sign_consecutive[d] += 1
                        prev_sign = curr_sign

        elapsed = time.time() - tn
        print(f"n={n}: {n_trees} trees, {n_edges} edges, {elapsed:.1f}s",
              flush=True)

    total_time = time.time() - t0

    print(f"\n=== SIGNED SUM BY DISTANCE ===")
    print(f"{'dist':>4} {'max':>12} {'min':>12} {'alternation%':>14}")
    for d in range(1, 17):
        mx = max_signed_by_dist.get(d, 0)
        mn = min_signed_by_dist.get(d, float('inf'))
        if mn == float('inf'):
            break
        tc = total_consecutive.get(d, 0)
        alt = alternations.get(d, 0)
        alt_pct = 100 * alt / tc if tc > 0 else 0
        print(f"{d:4d} {mx:12.8f} {mn:12.8f} {alt_pct:13.1f}%")

    print(f"\n=== PARTIAL SUM (cumulative through distance d) ===")
    print(f"{'dist':>4} {'max partial':>14} {'min partial':>14}")
    for d in range(1, 17):
        mx = max_partial_by_dist.get(d, 0)
        mn = min_partial_by_dist.get(d, float('inf'))
        if mn == float('inf'):
            break
        print(f"{d:4d} {mx:14.8f} {mn:14.8f}")

    print(f"\n=== KEY BOUND ===")
    # The partial sum at d → ∞ is the remote.
    # max partial at any d is an upper bound on the remote.
    max_partial_ever = max(max_partial_by_dist.values())
    min_partial_ever = min(v for v in min_partial_by_dist.values()
                          if v != float('inf'))
    print(f"Max partial sum over all d: {max_partial_ever:.8f}")
    print(f"Min partial sum over all d: {min_partial_ever:.8f}")
    print(f"=> |remote| ≤ max(|max|, |min|) = "
          f"{max(abs(max_partial_ever), abs(min_partial_ever)):.8f}")

    results = {
        "max_n": MAX_N,
        "total_edges": total_edges,
        "signed_by_dist": {
            str(d): {"max": round(max_signed_by_dist[d], 10),
                     "min": round(min_signed_by_dist[d], 10)}
            for d in sorted(max_signed_by_dist.keys())
        },
        "partial_by_dist": {
            str(d): {"max": round(max_partial_by_dist[d], 10),
                     "min": round(min_partial_by_dist[d], 10)}
            for d in sorted(max_partial_by_dist.keys())
            if min_partial_by_dist[d] != float('inf')
        },
        "alternation_by_dist": {
            str(d): round(100 * alternations[d] / total_consecutive[d], 2)
            if total_consecutive[d] > 0 else 0
            for d in sorted(alternations.keys())
        },
        "total_time_s": round(total_time, 2),
    }

    out_path = "results/ecms_signed_decay.json"
    os.makedirs("results", exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
