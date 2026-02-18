#!/usr/bin/env python3
"""Attempt to bound |δμ| = |μ(T) - μ(T/e)| using the hard-core cavity method.

Key decomposition (proved in hard-core model at λ=1):
  μ(T) - μ(T/e) = local_change + remote_change

where:
  local_change = P_T(u) + P_T(v) - P_{T/e}(u*)
               = (A+B-AB) / [(1+A+B)(1+AB)]
  where A = R_{u→v}, B = R_{v→u} are cavity messages.

  remote_change = Σ_{w ≠ u,v} [P_T(w) - P_{T/e}(w)]

We compute both terms for all trees through n=18 and check:
  1. Is local_change always in (0, 1/2)?
  2. How large can |remote_change| get?
  3. What fraction of δμ comes from local vs remote?
  4. Is there a universal bound δμ < C for some C < 1?
"""

import json
import math
import os
import time
from fractions import Fraction

from indpoly import independence_poly, is_log_concave
from trees import trees

MAX_N = 16  # Smaller n since we compute cavity messages per vertex
PROGRESS_INTERVAL = 2000


def cavity_messages(n, adj):
    """Compute all cavity messages R_{u→v} for a tree.

    In the hard-core model at λ=1 on a tree:
      R_{v→u} = Π_{c ∈ N(v)\{u}} 1/(1+R_{c→v})

    Returns dict: (v, u) → R_{v→u}
    """
    if n <= 1:
        return {}

    # Root at 0, compute messages bottom-up then top-down
    parent = [-1] * n
    children = [[] for _ in range(n)]
    visited = [False] * n
    order = []  # BFS order

    visited[0] = True
    queue = [0]
    head = 0
    while head < len(queue):
        v = queue[head]
        head += 1
        order.append(v)
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                parent[u] = v
                children[v].append(u)
                queue.append(u)

    msgs = {}  # (from, to) → R

    # Bottom-up: compute R_{v→parent(v)} for all v
    for v in reversed(order):
        p = parent[v]
        if p == -1:
            continue
        # R_{v→p} = Π_{c ∈ children(v)} 1/(1+R_{c→v})
        prod = 1.0
        for c in children[v]:
            prod *= 1.0 / (1.0 + msgs[(c, v)])
        msgs[(v, p)] = prod

    # Top-down: compute R_{parent→v} for all v
    for v in order:
        for c in children[v]:
            # R_{v→c} = Π_{j ∈ N(v)\{c}} 1/(1+R_{j→v})
            # N(v)\{c} = parent(v) (if exists) + other children
            prod = 1.0
            if parent[v] != -1:
                prod *= 1.0 / (1.0 + msgs[(parent[v], v)])
            for c2 in children[v]:
                if c2 != c:
                    prod *= 1.0 / (1.0 + msgs[(c2, v)])
            msgs[(v, c)] = prod

    return msgs


def occupation_probs(n, adj, msgs):
    """Compute P(v) = R_v/(1+R_v) for each vertex.

    R_v = Π_{c ∈ N(v)} 1/(1+R_{c→v})
    """
    P = [0.0] * n
    for v in range(n):
        Rv = 1.0
        for u in adj[v]:
            Rv *= 1.0 / (1.0 + msgs.get((u, v), 0.0))
        P[v] = Rv / (1.0 + Rv)
    return P


def contract_edge(n, adj, u, v):
    """Contract edge (u, v): merge u and v into u, remove v."""
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


def poly_mean(poly):
    s = sum(poly)
    if s == 0:
        return 0.0
    return sum(k * poly[k] for k in range(len(poly))) / s


def main():
    results = {
        "max_n": MAX_N,
        "total_edges": 0,
        "local_min": float("inf"),
        "local_max": float("-inf"),
        "remote_min": float("inf"),
        "remote_max": float("-inf"),
        "local_always_positive": True,
        "max_abs_dmu": 0.0,
        "max_local_fraction": 0.0,  # max |local/dmu|
        "min_local_fraction": float("inf"),
        "per_n": {},
        # Track: local = f(A, B) where A = R_{u→v}, B = R_{v→u}
        "AB_extremal": [],
    }

    t0 = time.time()

    for n in range(3, MAX_N + 1):
        tn = time.time()
        n_trees = 0
        n_edges = 0
        n_local_max = float("-inf")
        n_local_min = float("inf")
        n_remote_max = float("-inf")
        n_remote_min = float("inf")
        n_max_abs_dmu = 0.0

        for _, adj in trees(n):
            n_trees += 1
            msgs_T = cavity_messages(n, adj)
            P_T = occupation_probs(n, adj, msgs_T)
            mu_T = sum(P_T)

            # Cross-check with polynomial mean
            poly_T = independence_poly(n, adj)
            mu_T_poly = poly_mean(poly_T)
            assert abs(mu_T - mu_T_poly) < 1e-10, \
                f"Cavity mean {mu_T} != poly mean {mu_T_poly}"

            seen = set()
            for u in range(n):
                for v in adj[u]:
                    e = (min(u, v), max(u, v))
                    if e in seen:
                        continue
                    seen.add(e)
                    n_edges += 1
                    results["total_edges"] += 1

                    # Get cavity messages A = R_{u→v}, B = R_{v→u}
                    A = msgs_T.get((u, v), 0.0)
                    B = msgs_T.get((v, u), 0.0)

                    # Local change formula
                    num = A + B - A * B
                    den = (1 + A + B) * (1 + A * B)
                    local_change = num / den if den > 0 else 0.0

                    # Contract edge and compute P_{T/e}
                    nc, adjc, o2n = contract_edge(n, adj, u, v)
                    msgs_Te = cavity_messages(nc, adjc)
                    P_Te = occupation_probs(nc, adjc, msgs_Te)
                    mu_Te = sum(P_Te)

                    dmu = mu_T - mu_Te
                    remote_change = dmu - local_change

                    # Track extremes
                    if local_change > n_local_max:
                        n_local_max = local_change
                    if local_change < n_local_min:
                        n_local_min = local_change
                    if remote_change > n_remote_max:
                        n_remote_max = remote_change
                    if remote_change < n_remote_min:
                        n_remote_min = remote_change

                    abs_dmu = abs(dmu)
                    if abs_dmu > n_max_abs_dmu:
                        n_max_abs_dmu = abs_dmu

                    # Update global stats
                    if local_change > results["local_max"]:
                        results["local_max"] = local_change
                    if local_change < results["local_min"]:
                        results["local_min"] = local_change
                    if remote_change > results["remote_max"]:
                        results["remote_max"] = remote_change
                    if remote_change < results["remote_min"]:
                        results["remote_min"] = remote_change
                    if local_change <= 0:
                        results["local_always_positive"] = False
                    if abs_dmu > results["max_abs_dmu"]:
                        results["max_abs_dmu"] = abs_dmu

                    if abs_dmu > 1e-10:
                        frac = abs(local_change / dmu)
                        if frac > results["max_local_fraction"]:
                            results["max_local_fraction"] = frac
                        if frac < results["min_local_fraction"]:
                            results["min_local_fraction"] = frac

                    # Save extremal A,B pairs
                    if abs_dmu > 0.5:
                        results["AB_extremal"].append({
                            "n": n,
                            "edge": [u, v],
                            "A": round(A, 6),
                            "B": round(B, 6),
                            "local": round(local_change, 6),
                            "remote": round(remote_change, 6),
                            "dmu": round(dmu, 6),
                            "P_u": round(P_T[u], 6),
                            "P_v": round(P_T[v], 6),
                        })

            if n_trees % PROGRESS_INTERVAL == 0 and n_trees > 0:
                print(f"  n={n}: {n_trees} trees...", flush=True)

        elapsed = time.time() - tn
        results["per_n"][str(n)] = {
            "trees": n_trees,
            "edges": n_edges,
            "local_max": round(n_local_max, 8),
            "local_min": round(n_local_min, 8),
            "remote_max": round(n_remote_max, 8),
            "remote_min": round(n_remote_min, 8),
            "max_abs_dmu": round(n_max_abs_dmu, 8),
            "time_s": round(elapsed, 2),
        }
        print(f"n={n}: {n_trees} trees, {n_edges} edges, "
              f"local=[{n_local_min:.4f}, {n_local_max:.4f}], "
              f"remote=[{n_remote_min:.4f}, {n_remote_max:.4f}], "
              f"|dmu|<={n_max_abs_dmu:.4f}, {elapsed:.1f}s", flush=True)

    results["total_time_s"] = round(time.time() - t0, 2)

    # Clean up for JSON
    results["local_max"] = round(results["local_max"], 10)
    results["local_min"] = round(results["local_min"], 10)
    results["remote_max"] = round(results["remote_max"], 10)
    results["remote_min"] = round(results["remote_min"], 10)
    results["max_abs_dmu"] = round(results["max_abs_dmu"], 10)
    results["max_local_fraction"] = round(results["max_local_fraction"], 6)
    results["min_local_fraction"] = round(results["min_local_fraction"], 6)

    print(f"\n=== SUMMARY ===")
    print(f"Total edges: {results['total_edges']}")
    print(f"Local change range: [{results['local_min']}, {results['local_max']}]")
    print(f"Remote change range: [{results['remote_min']}, {results['remote_max']}]")
    print(f"Local always positive: {results['local_always_positive']}")
    print(f"max |δμ|: {results['max_abs_dmu']}")
    print(f"Local fraction of δμ: [{results['min_local_fraction']}, "
          f"{results['max_local_fraction']}]")
    print(f"Extremal (A,B) pairs: {len(results['AB_extremal'])}")
    for ex in results["AB_extremal"][:10]:
        print(f"  n={ex['n']}, A={ex['A']}, B={ex['B']}, "
              f"local={ex['local']}, remote={ex['remote']}, dmu={ex['dmu']}")

    out_path = "results/ecms_cavity_bound.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
