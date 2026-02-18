#!/usr/bin/env python3
"""Decompose the remote term of δμ by distance from the contracted edge.

For each tree T and edge e = uv, after contracting e to get T/e:
  remote = Σ_{w ≠ u,v} [P_T(w) - P_{T/e}(w)]

We decompose this by dist(w, {u,v}) in T, where dist is the minimum
distance to either endpoint. Track:
  1. Max |ΔP(w)| at each distance
  2. Sum of |ΔP(w)| at each distance
  3. Whether |ΔP(w)| decays geometrically with distance
  4. The effective decay rate

Also measure the "initial perturbation" at the merged vertex and
how cavity messages change step by step.
"""

import json
import os
import time
from collections import defaultdict

from indpoly import independence_poly
from trees import trees

MAX_N = 18
PROGRESS_INTERVAL = 5000


def cavity_messages(n, adj):
    """Compute cavity messages R_{u→v} for a tree at λ=1."""
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
    """P(v) = R_v/(1+R_v)."""
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
    """BFS distance from the source set."""
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
    """Contract edge (u, v): merge v into u, return (n', adj', old_to_new)."""
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

    # Per-distance statistics
    max_dist_seen = 0
    # dist → max |ΔP(w)| across all edges
    max_abs_dP_by_dist = defaultdict(float)
    # dist → sum of |ΔP(w)| across all edges
    sum_abs_dP_by_dist = defaultdict(float)
    # dist → count of vertices at that distance
    count_by_dist = defaultdict(int)
    # dist → max sum|ΔP| at that distance for a single edge
    max_sum_dP_single_edge_by_dist = defaultdict(float)

    # Track worst-case remote decomposition
    worst_remote_pos = {"remote": 0.0}
    worst_remote_neg = {"remote": 0.0}
    worst_abs_dmu = {"abs_dmu": 0.0}

    # Decay rate: for each edge, fit |ΔP| vs distance
    decay_rates = []  # (rate, n, edge)

    total_edges = 0

    for n in range(3, MAX_N + 1):
        tn = time.time()
        n_trees = 0
        n_edges = 0
        n_max_abs_dmu = 0.0
        n_max_remote = -1e18
        n_min_remote = 1e18

        for _, adj in trees(n):
            n_trees += 1
            msgs_T = cavity_messages(n, adj)
            P_T = occupation_probs(n, adj, msgs_T)

            seen = set()
            for u in range(n):
                for v in adj[u]:
                    e = (min(u, v), max(u, v))
                    if e in seen:
                        continue
                    seen.add(e)
                    n_edges += 1
                    total_edges += 1

                    # Distances from {u, v}
                    dist = bfs_distances(n, adj, [u, v])

                    # Contract edge
                    nc, adjc, o2n = contract_edge(n, adj, u, v)
                    msgs_Te = cavity_messages(nc, adjc)
                    P_Te = occupation_probs(nc, adjc, msgs_Te)

                    # Local term
                    A = msgs_T.get((u, v), 0.0)
                    B = msgs_T.get((v, u), 0.0)
                    num = A + B - A * B
                    den = (1 + A + B) * (1 + A * B)
                    local = num / den if den > 0 else 0.0

                    # Per-distance remote decomposition
                    per_dist_sum = defaultdict(float)
                    per_dist_abs_sum = defaultdict(float)
                    per_dist_max_abs = defaultdict(float)

                    remote = 0.0
                    for w in range(n):
                        if w == u or w == v:
                            continue
                        w_new = o2n[w]
                        dP = P_T[w] - P_Te[w_new]
                        remote += dP
                        d = dist[w]
                        if d > max_dist_seen:
                            max_dist_seen = d
                        per_dist_sum[d] += dP
                        per_dist_abs_sum[d] += abs(dP)
                        per_dist_max_abs[d] = max(per_dist_max_abs[d], abs(dP))

                        # Global stats
                        max_abs_dP_by_dist[d] = max(max_abs_dP_by_dist[d], abs(dP))
                        sum_abs_dP_by_dist[d] += abs(dP)
                        count_by_dist[d] += 1

                    for d in per_dist_abs_sum:
                        max_sum_dP_single_edge_by_dist[d] = max(
                            max_sum_dP_single_edge_by_dist[d],
                            per_dist_abs_sum[d])

                    dmu = local + remote
                    abs_dmu = abs(dmu)

                    if abs_dmu > n_max_abs_dmu:
                        n_max_abs_dmu = abs_dmu
                    if remote > n_max_remote:
                        n_max_remote = remote
                    if remote < n_min_remote:
                        n_min_remote = remote

                    # Track worst cases with full decomposition
                    if remote > worst_remote_pos.get("remote", -1e18):
                        worst_remote_pos = {
                            "remote": remote, "local": local, "dmu": dmu,
                            "n": n, "edge": [u, v],
                            "A": round(A, 6), "B": round(B, 6),
                            "per_dist": {str(d): round(v, 8)
                                         for d, v in sorted(per_dist_sum.items())},
                        }
                    if remote < worst_remote_neg.get("remote", 1e18):
                        worst_remote_neg = {
                            "remote": remote, "local": local, "dmu": dmu,
                            "n": n, "edge": [u, v],
                            "A": round(A, 6), "B": round(B, 6),
                            "per_dist": {str(d): round(v, 8)
                                         for d, v in sorted(per_dist_sum.items())},
                        }
                    if abs_dmu > worst_abs_dmu.get("abs_dmu", 0):
                        worst_abs_dmu = {
                            "abs_dmu": abs_dmu, "local": local, "remote": remote,
                            "dmu": dmu, "n": n, "edge": [u, v],
                            "A": round(A, 6), "B": round(B, 6),
                            "per_dist": {str(d): round(v, 8)
                                         for d, v in sorted(per_dist_sum.items())},
                        }

                    # Estimate decay rate for this edge
                    # Fit: max|ΔP(w)| at dist d ≈ C * rate^d
                    # Use dist 1 and max available dist
                    if per_dist_max_abs.get(1, 0) > 1e-12:
                        max_d = max(d for d in per_dist_max_abs if d >= 1)
                        if max_d >= 2 and per_dist_max_abs.get(max_d, 0) > 1e-15:
                            import math
                            rate = (per_dist_max_abs[max_d] /
                                    per_dist_max_abs[1]) ** (1.0 / (max_d - 1))
                            if rate < 2.0:  # sanity
                                decay_rates.append((rate, n, e))

        elapsed = time.time() - tn
        print(f"n={n}: {n_trees} trees, {n_edges} edges, "
              f"remote=[{n_min_remote:.4f}, {n_max_remote:.4f}], "
              f"|dmu|<={n_max_abs_dmu:.4f}, {elapsed:.1f}s", flush=True)

    total_time = time.time() - t0

    # Compute average decay rate
    if decay_rates:
        rates = [r for r, _, _ in decay_rates]
        avg_rate = sum(rates) / len(rates)
        max_rate = max(rates)
        p95_rate = sorted(rates)[int(0.95 * len(rates))]
        p99_rate = sorted(rates)[int(0.99 * len(rates))]
    else:
        avg_rate = max_rate = p95_rate = p99_rate = 0.0

    # Print distance decomposition
    print(f"\n=== REMOTE CHANGE BY DISTANCE FROM CONTRACTED EDGE ===")
    print(f"{'dist':>4} {'max|ΔP|':>12} {'avg|ΔP|':>12} {'max Σ|ΔP| (1 edge)':>20} {'count':>10}")
    for d in range(1, max_dist_seen + 1):
        if count_by_dist[d] == 0:
            continue
        avg = sum_abs_dP_by_dist[d] / count_by_dist[d]
        print(f"{d:4d} {max_abs_dP_by_dist[d]:12.8f} {avg:12.8f} "
              f"{max_sum_dP_single_edge_by_dist[d]:20.8f} {count_by_dist[d]:10d}")

    print(f"\n=== DECAY RATES (max|ΔP| ~ C·r^d) ===")
    print(f"Average decay rate: {avg_rate:.4f}")
    print(f"Max decay rate: {max_rate:.4f}")
    print(f"95th percentile: {p95_rate:.4f}")
    print(f"99th percentile: {p99_rate:.4f}")

    print(f"\n=== WORST REMOTE (positive) ===")
    for k, v in sorted(worst_remote_pos.items()):
        print(f"  {k}: {v}")

    print(f"\n=== WORST REMOTE (negative) ===")
    for k, v in sorted(worst_remote_neg.items()):
        print(f"  {k}: {v}")

    print(f"\n=== WORST |δμ| ===")
    for k, v in sorted(worst_abs_dmu.items()):
        print(f"  {k}: {v}")

    # Save results
    results = {
        "max_n": MAX_N,
        "total_edges": total_edges,
        "total_time_s": round(total_time, 2),
        "distance_decomposition": {
            str(d): {
                "max_abs_dP": round(max_abs_dP_by_dist[d], 10),
                "avg_abs_dP": round(sum_abs_dP_by_dist[d] / count_by_dist[d], 10)
                    if count_by_dist[d] > 0 else 0,
                "max_sum_abs_dP_single_edge": round(
                    max_sum_dP_single_edge_by_dist[d], 10),
                "count": count_by_dist[d],
            }
            for d in range(1, max_dist_seen + 1)
            if count_by_dist[d] > 0
        },
        "decay_rates": {
            "average": round(avg_rate, 6),
            "max": round(max_rate, 6),
            "p95": round(p95_rate, 6),
            "p99": round(p99_rate, 6),
            "n_samples": len(decay_rates),
        },
        "worst_remote_positive": worst_remote_pos,
        "worst_remote_negative": worst_remote_neg,
        "worst_abs_dmu": worst_abs_dmu,
    }

    out_path = "results/ecms_remote_decay.json"
    os.makedirs("results", exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
