#!/usr/bin/env python3
"""Test mixed caterpillars: part with k leaves per spine, part bare path.

The exhaustive search shows the worst remote trees have degree sequences
like [4,4,4,3,2,2,...] -- caterpillars where some spine vertices have
3 leaves and some have 0. The transition edge creates the worst remote.

Build: caterpillar with m spine vertices having k leaves, followed by
p bare path vertices. Contract the edge at the transition.
"""

import json
import os
import time
from collections import defaultdict


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


def make_mixed_caterpillar(m, k, p):
    """m spine vertices with k leaves each, then p bare path vertices.

    Spine: 0,1,...,m-1 (each with k leaves)
    Path: m, m+1, ..., m+p-1 (no leaves)
    Transition edge: (m-1, m)
    """
    n_leaves = m * k
    n = m + p + n_leaves
    adj = [[] for _ in range(n)]

    # Spine + path: 0-1-2-..-(m+p-1)
    for i in range(m + p - 1):
        adj[i].append(i + 1)
        adj[i + 1].append(i)

    # Leaves on spine vertices
    leaf_idx = m + p
    for i in range(m):
        for _ in range(k):
            adj[i].append(leaf_idx)
            adj[leaf_idx].append(i)
            leaf_idx += 1

    return n, adj


def analyze_edge(n, adj, u, v):
    """Compute local, remote, dmu for contracting edge (u,v)."""
    msgs = cavity_messages(n, adj)
    P_T = occupation_probs(n, adj, msgs)
    mu_T = sum(P_T)

    A = msgs.get((u, v), 0.0)
    B = msgs.get((v, u), 0.0)
    num = A + B - A * B
    den = (1 + A + B) * (1 + A * B)
    loc = num / den if den > 0 else 0.0

    nc, adjc, o2n = contract_edge(n, adj, u, v)
    msgs_c = cavity_messages(nc, adjc)
    P_Te = occupation_probs(nc, adjc, msgs_c)
    mu_Te = sum(P_Te)

    dmu = mu_T - mu_Te
    remote = dmu - loc
    return {"A": A, "B": B, "local": loc, "remote": remote, "dmu": dmu}


def analyze_all_edges(n, adj):
    """Find worst edge for remote."""
    msgs = cavity_messages(n, adj)
    P_T = occupation_probs(n, adj, msgs)
    mu_T = sum(P_T)

    worst = {"remote": 0}
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

            dmu = mu_T - mu_Te
            remote = dmu - loc

            if abs(remote) > abs(worst["remote"]):
                worst = {"edge": (u, v), "A": A, "B": B,
                         "local": loc, "remote": remote, "dmu": dmu}
    return worst


def main():
    print("=== MIXED CATERPILLAR: m vertices with k=3 leaves + p bare path ===")
    print("Testing transition edge (m-1, m) and all edges\n")

    # Vary m and p to find the worst case
    results = []

    for k in [2, 3, 4]:
        print(f"\n--- k={k} leaves per caterpillar vertex ---")
        for p in [0, 1, 2, 3, 5, 10, 20, 50, 100]:
            print(f"\n  p={p} (bare path length):")
            for m in [1, 2, 3, 5, 10, 20, 50, 100]:
                n, adj = make_mixed_caterpillar(m, k, p)
                if n < 3:
                    continue

                # Transition edge
                if m > 0 and p > 0:
                    u, v = m - 1, m
                    res_trans = analyze_edge(n, adj, u, v)
                    # Also check all edges for the worst
                    res_worst = analyze_all_edges(n, adj)
                    results.append({
                        "k": k, "m": m, "p": p, "n": n,
                        "trans_remote": res_trans["remote"],
                        "trans_dmu": res_trans["dmu"],
                        "worst_remote": res_worst["remote"],
                        "worst_dmu": res_worst["dmu"],
                        "worst_edge": res_worst.get("edge"),
                    })
                    print(f"    m={m:3d}, n={n:4d}: "
                          f"trans_remote={res_trans['remote']:.6f}, "
                          f"worst_remote={res_worst['remote']:.6f}, "
                          f"worst_|δμ|={abs(res_worst['dmu']):.6f}")
                else:
                    res_worst = analyze_all_edges(n, adj)
                    results.append({
                        "k": k, "m": m, "p": p, "n": n,
                        "trans_remote": None,
                        "worst_remote": res_worst["remote"],
                        "worst_dmu": res_worst["dmu"],
                    })
                    print(f"    m={m:3d}, n={n:4d}: "
                          f"worst_remote={res_worst['remote']:.6f}, "
                          f"worst_|δμ|={abs(res_worst['dmu']):.6f}")

    # Find the global worst
    print("\n=== GLOBAL WORST (by |remote|) ===")
    results.sort(key=lambda r: abs(r["worst_remote"]), reverse=True)
    for r in results[:10]:
        print(f"  k={r['k']}, m={r['m']}, p={r['p']}, n={r['n']}: "
              f"remote={r['worst_remote']:.8f}, dmu={r['worst_dmu']:.8f}")

    # Now focus on k=3 and find optimal m, p
    print("\n=== OPTIMAL m,p for k=3 (worst |remote|) ===")
    k = 3
    best = {"remote": 0}
    for m in range(1, 201):
        for p in [0, 1, 2, 3, 5, 10, 20, 50, 100, 200]:
            n, adj = make_mixed_caterpillar(m, k, p)
            if n < 3:
                continue
            res = analyze_all_edges(n, adj)
            if abs(res["remote"]) > abs(best.get("remote", 0)):
                best = {"m": m, "p": p, "n": n,
                        "remote": res["remote"], "dmu": res["dmu"]}

    print(f"Best: m={best['m']}, p={best['p']}, n={best['n']}: "
          f"remote={best['remote']:.8f}")

    # Also try "double caterpillar": m1 with k1 leaves, m2 with k2 leaves
    print("\n=== DOUBLE CATERPILLAR: k1 leaves section + k2 leaves section ===")
    double_results = []
    for k1 in [2, 3, 4, 5]:
        for k2 in [0, 1, 2]:
            if k1 == k2:
                continue
            for m1 in [3, 10, 50]:
                for m2 in [3, 10, 50]:
                    # Build: m1 spine with k1 leaves + m2 spine with k2 leaves
                    total_spine = m1 + m2
                    n_leaves = m1 * k1 + m2 * k2
                    n = total_spine + n_leaves
                    adj = [[] for _ in range(n)]

                    # Spine path
                    for i in range(total_spine - 1):
                        adj[i].append(i + 1)
                        adj[i + 1].append(i)

                    # Leaves
                    leaf_idx = total_spine
                    for i in range(m1):
                        for _ in range(k1):
                            adj[i].append(leaf_idx)
                            adj[leaf_idx].append(i)
                            leaf_idx += 1
                    for i in range(m1, total_spine):
                        for _ in range(k2):
                            adj[i].append(leaf_idx)
                            adj[leaf_idx].append(i)
                            leaf_idx += 1

                    res = analyze_all_edges(n, adj)
                    double_results.append({
                        "k1": k1, "k2": k2, "m1": m1, "m2": m2,
                        "n": n, "remote": res["remote"], "dmu": res["dmu"],
                    })

    double_results.sort(key=lambda r: abs(r["remote"]), reverse=True)
    print("Top 10 by |remote|:")
    for r in double_results[:10]:
        print(f"  k1={r['k1']}, k2={r['k2']}, m1={r['m1']}, m2={r['m2']}, "
              f"n={r['n']}: remote={r['remote']:.8f}")

    out_path = "results/ecms_mixed_caterpillar.json"
    os.makedirs("results", exist_ok=True)
    with open(out_path, "w") as f:
        json.dump({"mixed": results[:20], "double": double_results[:20],
                   "best_k3": best}, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
