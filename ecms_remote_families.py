#!/usr/bin/env python3
"""Study remote bound for specific tree families at large n.

Key question: does max|remote| stay bounded as n → ∞?
If it converges to some C < 1/2, we get |δμ| < 1 (since local < 1/2).

Families:
  - Paths P_n (long correlation length, γ → 1)
  - Stars K_{1,s} (maximum branching)
  - Double stars S(a,b) (two hubs)
  - Caterpillars (path of hubs)
  - Spiders S(l1,...,lk) (central branching)
  - Balanced tripods T(a,a,a)
  - Extended stars S(m,1^k) (one long arm + short arms)

For each family, compute local and remote for every edge, track max|remote|
and max|δμ| as n grows. Look for convergence.
"""

import json
import os
import time
from collections import defaultdict

from indpoly import independence_poly


def make_path(n):
    """Path on n vertices: 0-1-2-..-(n-1)."""
    adj = [[] for _ in range(n)]
    for i in range(n - 1):
        adj[i].append(i + 1)
        adj[i + 1].append(i)
    return n, adj


def make_star(s):
    """Star K_{1,s}: center 0, leaves 1..s."""
    n = s + 1
    adj = [[] for _ in range(n)]
    for i in range(1, n):
        adj[0].append(i)
        adj[i].append(0)
    return n, adj


def make_double_star(a, b):
    """Double star S(a,b): centers 0,1 connected, a leaves on 0, b leaves on 1."""
    n = a + b + 2
    adj = [[] for _ in range(n)]
    adj[0].append(1)
    adj[1].append(0)
    for i in range(2, 2 + a):
        adj[0].append(i)
        adj[i].append(0)
    for i in range(2 + a, n):
        adj[1].append(i)
        adj[i].append(1)
    return n, adj


def make_caterpillar(spine_len, leaves_per_vertex):
    """Caterpillar: path of spine_len vertices, each with leaves_per_vertex leaves."""
    n = spine_len + spine_len * leaves_per_vertex
    adj = [[] for _ in range(n)]
    # Spine: 0, 1, ..., spine_len-1
    for i in range(spine_len - 1):
        adj[i].append(i + 1)
        adj[i + 1].append(i)
    # Leaves
    leaf_idx = spine_len
    for i in range(spine_len):
        for _ in range(leaves_per_vertex):
            adj[i].append(leaf_idx)
            adj[leaf_idx].append(i)
            leaf_idx += 1
    return n, adj


def make_spider(arms):
    """Spider: center 0, arms of given lengths."""
    n = 1 + sum(arms)
    adj = [[] for _ in range(n)]
    idx = 1
    for arm_len in arms:
        prev = 0
        for _ in range(arm_len):
            adj[prev].append(idx)
            adj[idx].append(prev)
            prev = idx
            idx += 1
    return n, adj


def make_tripod(a):
    """Balanced tripod T(a,a,a): spider with 3 arms of length a+1 (a leaves on each)."""
    # Actually: center connected to 3 bridges, each bridge has a leaves
    # = spider with arms [a+1, a+1, a+1]? No...
    # T(a,a,a) from the paper: center vertex, 3 paths of length a+1 each
    # Let me use spider with 3 arms of length a
    return make_spider([a, a, a])


def make_broom(handle_len, bristles):
    """Broom: path of handle_len, then star with bristles at the end."""
    n = handle_len + bristles
    adj = [[] for _ in range(n)]
    # Handle: 0-1-..-(handle_len-1)
    for i in range(handle_len - 1):
        adj[i].append(i + 1)
        adj[i + 1].append(i)
    # Bristles on vertex handle_len-1
    hub = handle_len - 1
    for i in range(handle_len, n):
        adj[hub].append(i)
        adj[i].append(hub)
    return n, adj


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


def analyze_edges(n, adj, max_edges=None):
    """For each edge, compute local, remote, δμ. Return summary."""
    msgs_T = cavity_messages(n, adj)
    P_T = occupation_probs(n, adj, msgs_T)
    mu_T = sum(P_T)

    results = []
    seen = set()
    count = 0
    for u in range(n):
        for v in adj[u]:
            e = (min(u, v), max(u, v))
            if e in seen:
                continue
            seen.add(e)
            count += 1
            if max_edges and count > max_edges:
                break

            A = msgs_T.get((u, v), 0.0)
            B = msgs_T.get((v, u), 0.0)
            num = A + B - A * B
            den = (1 + A + B) * (1 + A * B)
            local = num / den if den > 0 else 0.0

            nc, adjc, o2n = contract_edge(n, adj, u, v)
            msgs_Te = cavity_messages(nc, adjc)
            P_Te = occupation_probs(nc, adjc, msgs_Te)
            mu_Te = sum(P_Te)

            dmu = mu_T - mu_Te
            remote = dmu - local

            results.append({
                "edge": (u, v), "A": A, "B": B,
                "local": local, "remote": remote, "dmu": dmu,
            })
        if max_edges and count > max_edges:
            break

    return results


def main():
    t0 = time.time()
    all_results = {}

    # 1. Paths
    print("=== PATHS ===")
    path_data = []
    for n_val in range(4, 201):
        n, adj = make_path(n_val)
        res = analyze_edges(n, adj)
        max_remote = max(abs(r["remote"]) for r in res)
        max_dmu = max(abs(r["dmu"]) for r in res)
        max_local = max(r["local"] for r in res)
        # Find the extremal edge
        ext = max(res, key=lambda r: abs(r["remote"]))
        path_data.append({
            "n": n_val, "max_remote": max_remote,
            "max_dmu": max_dmu, "max_local": max_local,
            "ext_A": ext["A"], "ext_B": ext["B"],
        })
        if n_val <= 20 or n_val % 20 == 0:
            print(f"  P_{n_val}: max|remote|={max_remote:.6f}, "
                  f"|δμ|={max_dmu:.6f}, local={max_local:.6f}", flush=True)
    all_results["paths"] = path_data

    # 2. Stars
    print("\n=== STARS ===")
    star_data = []
    for s in range(2, 201):
        n, adj = make_star(s)
        res = analyze_edges(n, adj)
        max_remote = max(abs(r["remote"]) for r in res)
        max_dmu = max(abs(r["dmu"]) for r in res)
        max_local = max(r["local"] for r in res)
        star_data.append({
            "n": s + 1, "s": s, "max_remote": max_remote,
            "max_dmu": max_dmu, "max_local": max_local,
        })
        if s <= 20 or s % 20 == 0:
            print(f"  K_{{1,{s}}}: max|remote|={max_remote:.6f}, "
                  f"|δμ|={max_dmu:.6f}", flush=True)
    all_results["stars"] = star_data

    # 3. Balanced tripods
    print("\n=== BALANCED TRIPODS T(a,a,a) ===")
    tripod_data = []
    for a in range(2, 80):
        n, adj = make_tripod(a)
        res = analyze_edges(n, adj)
        max_remote = max(abs(r["remote"]) for r in res)
        max_dmu = max(abs(r["dmu"]) for r in res)
        max_local = max(r["local"] for r in res)
        ext = max(res, key=lambda r: abs(r["dmu"]))
        tripod_data.append({
            "n": 3 * a + 1, "a": a, "max_remote": max_remote,
            "max_dmu": max_dmu, "max_local": max_local,
            "ext_edge": list(ext["edge"]),
        })
        if a <= 10 or a % 10 == 0:
            print(f"  T({a},{a},{a}) n={3*a+1}: max|remote|={max_remote:.6f}, "
                  f"|δμ|={max_dmu:.6f}", flush=True)
    all_results["tripods"] = tripod_data

    # 4. Brooms (long handle + star)
    print("\n=== BROOMS (handle + bristles) ===")
    broom_data = []
    for handle in [5, 10, 20, 50]:
        for bristles in [3, 5, 10, 20, 50]:
            n, adj = make_broom(handle, bristles)
            res = analyze_edges(n, adj)
            max_remote = max(abs(r["remote"]) for r in res)
            max_dmu = max(abs(r["dmu"]) for r in res)
            broom_data.append({
                "handle": handle, "bristles": bristles, "n": n,
                "max_remote": max_remote, "max_dmu": max_dmu,
            })
            print(f"  Broom({handle},{bristles}) n={n}: "
                  f"max|remote|={max_remote:.6f}, |δμ|={max_dmu:.6f}", flush=True)
    all_results["brooms"] = broom_data

    # 5. Extended stars S(m, 1^k) -- one long arm + k leaves
    print("\n=== EXTENDED STARS S(m, 1^k) ===")
    ext_star_data = []
    for m in [5, 10, 20, 50]:
        for k in [2, 5, 10, 20]:
            arms = [m] + [1] * k
            n, adj = make_spider(arms)
            res = analyze_edges(n, adj)
            max_remote = max(abs(r["remote"]) for r in res)
            max_dmu = max(abs(r["dmu"]) for r in res)
            ext_star_data.append({
                "m": m, "k": k, "n": n,
                "max_remote": max_remote, "max_dmu": max_dmu,
            })
            print(f"  S({m},1^{k}) n={m+k+1}: "
                  f"max|remote|={max_remote:.6f}, |δμ|={max_dmu:.6f}", flush=True)
    all_results["extended_stars"] = ext_star_data

    # 6. Caterpillars
    print("\n=== CATERPILLARS ===")
    cat_data = []
    for spine in [5, 10, 20, 50]:
        for lpv in [1, 2, 3, 5]:
            n, adj = make_caterpillar(spine, lpv)
            res = analyze_edges(n, adj)
            max_remote = max(abs(r["remote"]) for r in res)
            max_dmu = max(abs(r["dmu"]) for r in res)
            cat_data.append({
                "spine": spine, "leaves_per_vertex": lpv, "n": n,
                "max_remote": max_remote, "max_dmu": max_dmu,
            })
            print(f"  Cat({spine},{lpv}) n={n}: "
                  f"max|remote|={max_remote:.6f}, |δμ|={max_dmu:.6f}", flush=True)
    all_results["caterpillars"] = cat_data

    # Summary
    print(f"\n=== CONVERGENCE SUMMARY ===")

    # Path convergence
    if path_data:
        print(f"Paths: max|remote| converges to {path_data[-1]['max_remote']:.6f}, "
              f"|δμ| to {path_data[-1]['max_dmu']:.6f}")

    # Star convergence
    if star_data:
        print(f"Stars: max|remote| converges to {star_data[-1]['max_remote']:.6f}, "
              f"|δμ| to {star_data[-1]['max_dmu']:.6f}")

    # Tripod convergence
    if tripod_data:
        print(f"Tripods: max|remote| converges to {tripod_data[-1]['max_remote']:.6f}, "
              f"|δμ| to {tripod_data[-1]['max_dmu']:.6f}")

    # Find global max across all families
    all_remotes = []
    all_dmus = []
    for family, data in all_results.items():
        for entry in data:
            all_remotes.append((entry["max_remote"], family, entry.get("n", 0)))
            all_dmus.append((entry["max_dmu"], family, entry.get("n", 0)))

    all_remotes.sort(reverse=True)
    all_dmus.sort(reverse=True)

    print(f"\nTop 5 max|remote|:")
    for val, fam, n_val in all_remotes[:5]:
        print(f"  {val:.6f} ({fam}, n={n_val})")

    print(f"\nTop 5 max|δμ|:")
    for val, fam, n_val in all_dmus[:5]:
        print(f"  {val:.6f} ({fam}, n={n_val})")

    total_time = time.time() - t0
    all_results["total_time_s"] = round(total_time, 2)

    out_path = "results/ecms_remote_families.json"
    os.makedirs("results", exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
