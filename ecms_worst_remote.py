#!/usr/bin/env python3
"""Identify the trees that maximize |remote| at each n.

For the ECMS proof: we need |remote| < 1/2. The exhaustive search through
n=18 shows max|remote| = 0.264, slowly increasing. Is it converging?
What trees achieve the maximum?

For each n, find the edge with max |remote| and print the tree structure.
"""

import json
import os
import time
from collections import defaultdict

from indpoly import independence_poly
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


def tree_description(n, adj):
    """Describe tree structure: degree sequence and type."""
    degs = sorted([len(adj[v]) for v in range(n)], reverse=True)
    leaves = sum(1 for d in degs if d == 1)
    max_deg = degs[0] if degs else 0

    # Classify
    if max_deg == n - 1:
        return f"Star K_{{1,{n-1}}}", degs
    if max_deg == 2:
        return f"Path P_{n}", degs
    if degs.count(1) == n - 2 and n >= 4:
        # Two non-leaf vertices
        hubs = [v for v in range(n) if len(adj[v]) > 1]
        if len(hubs) == 2 and hubs[0] in adj[hubs[1]]:
            a = len(adj[hubs[0]]) - 1
            b = len(adj[hubs[1]]) - 1
            return f"Double star S({max(a,b)},{min(a,b)})", degs

    # Check caterpillar: all non-leaf vertices form a path
    internal = [v for v in range(n) if len(adj[v]) > 1]
    if len(internal) >= 2:
        # Check if internal vertices form a path
        is_path = all(
            sum(1 for u in adj[v] if len(adj[u]) > 1) <= 2
            for v in internal
        )
        if is_path:
            spine_degs = sorted([len(adj[v]) for v in internal])
            return f"Caterpillar spine={len(internal)} degs={spine_degs}", degs

    # Spider: one vertex of high degree, rest form paths
    if sum(1 for d in degs if d >= 3) == 1:
        hub = [v for v in range(n) if len(adj[v]) >= 3][0]
        arms = []
        for c in adj[hub]:
            # Walk to leaf
            arm_len = 1
            prev, curr = hub, c
            while len(adj[curr]) == 2:
                arm_len += 1
                nxt = [x for x in adj[curr] if x != prev][0]
                prev, curr = curr, nxt
            arms.append(arm_len)
        arms.sort(reverse=True)
        return f"Spider({','.join(map(str,arms))})", degs

    return f"deg_seq={degs[:8]}", degs


def main():
    t0 = time.time()
    results_per_n = []

    for n in range(3, MAX_N + 1):
        tn = time.time()
        n_trees = 0

        worst_neg = {"remote": 0, "n": n}
        worst_pos = {"remote": 0, "n": n}
        worst_dmu = {"abs_dmu": 0, "n": n}

        for _, adj in trees(n):
            n_trees += 1
            msgs = cavity_messages(n, adj)
            P_T = occupation_probs(n, adj, msgs)

            seen = set()
            for u in range(n):
                for v in adj[u]:
                    e = (min(u, v), max(u, v))
                    if e in seen:
                        continue
                    seen.add(e)

                    A = msgs.get((u, v), 0.0)
                    B = msgs.get((v, u), 0.0)
                    num = A + B - A * B
                    den = (1 + A + B) * (1 + A * B)
                    loc = num / den if den > 0 else 0.0

                    nc, adjc, o2n = contract_edge(n, adj, u, v)
                    msgs_c = cavity_messages(nc, adjc)
                    P_Te = occupation_probs(nc, adjc, msgs_c)
                    mu_Te = sum(P_Te)
                    mu_T = sum(P_T)

                    dmu = mu_T - mu_Te
                    remote = dmu - loc

                    if remote < worst_neg["remote"]:
                        desc, degs = tree_description(n, adj)
                        worst_neg = {
                            "remote": remote, "local": loc, "dmu": dmu,
                            "n": n, "edge": [u, v],
                            "A": round(A, 6), "B": round(B, 6),
                            "tree": desc,
                            "deg_u": len(adj[u]), "deg_v": len(adj[v]),
                            "adj": [list(adj[i]) for i in range(n)],
                        }
                    if remote > worst_pos["remote"]:
                        desc, _ = tree_description(n, adj)
                        worst_pos = {
                            "remote": remote, "local": loc, "dmu": dmu,
                            "n": n, "edge": [u, v],
                            "A": round(A, 6), "B": round(B, 6),
                            "tree": desc,
                            "deg_u": len(adj[u]), "deg_v": len(adj[v]),
                        }
                    if abs(dmu) > worst_dmu["abs_dmu"]:
                        desc, _ = tree_description(n, adj)
                        worst_dmu = {
                            "abs_dmu": abs(dmu), "dmu": dmu,
                            "local": loc, "remote": remote,
                            "n": n, "edge": [u, v],
                            "A": round(A, 6), "B": round(B, 6),
                            "tree": desc,
                        }

        elapsed = time.time() - tn
        results_per_n.append({
            "n": n, "trees": n_trees,
            "worst_neg_remote": worst_neg,
            "worst_pos_remote": worst_pos,
            "worst_dmu": worst_dmu,
        })

        print(f"\nn={n} ({n_trees} trees, {elapsed:.1f}s):", flush=True)
        print(f"  Worst neg remote: {worst_neg['remote']:.6f} "
              f"(tree: {worst_neg.get('tree','?')}, edge={worst_neg.get('edge','?')}, "
              f"deg_u={worst_neg.get('deg_u','?')}, deg_v={worst_neg.get('deg_v','?')})")
        print(f"  Worst pos remote: {worst_pos['remote']:.6f} "
              f"(tree: {worst_pos.get('tree','?')})")
        print(f"  Worst |δμ|: {worst_dmu['abs_dmu']:.6f} "
              f"(tree: {worst_dmu.get('tree','?')})")

    total_time = time.time() - t0
    print(f"\n=== TREND: worst negative remote by n ===")
    for r in results_per_n:
        print(f"  n={r['n']:2d}: {r['worst_neg_remote']['remote']:.6f}  "
              f"tree={r['worst_neg_remote'].get('tree','?')}")

    print(f"\n=== TREND: worst |δμ| by n ===")
    for r in results_per_n:
        print(f"  n={r['n']:2d}: {r['worst_dmu']['abs_dmu']:.6f}  "
              f"tree={r['worst_dmu'].get('tree','?')}")

    out_path = "results/ecms_worst_remote.json"
    os.makedirs("results", exist_ok=True)
    # Remove adj from output to keep file small
    for r in results_per_n:
        for key in ["worst_neg_remote", "worst_pos_remote", "worst_dmu"]:
            r[key].pop("adj", None)
    with open(out_path, "w") as f:
        json.dump({"per_n": results_per_n, "total_time_s": round(total_time, 2)},
                  f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
