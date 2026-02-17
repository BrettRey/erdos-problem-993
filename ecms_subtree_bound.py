#!/usr/bin/env python3
"""Test whether subtree contributions to the remote term are bounded.

Key idea for proving |remote| < 1/2:
  When we contract edge e = uv, the perturbation propagates outward.
  For each neighbor w of {u,v} (at distance 1), define:
    subtree(w) = all vertices reachable from w without going through {u,v}
    subtree_sum(w) = Σ_{z ∈ subtree(w)} [P_T(z) - P_{T/e}(z')]

  If |subtree_sum(w)| ≤ C · |ΔP(w)| for some universal constant C,
  then the total remote is bounded.

  Better: track the ratio subtree_sum / initial_perturbation at w,
  where initial_perturbation = change in cavity message entering w.

  Even better: check if the subtree sums have OPPOSITE signs from the
  initial perturbation, causing cancellation. Then the remote is
  SMALLER than the sum of distance-1 contributions.

  Also test: the "amplification ratio" -- how much does the total
  perturbation in a subtree exceed the perturbation at its root?
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


def get_subtree(n, adj, root, exclude):
    """Get all vertices reachable from root without going through exclude set."""
    visited = set(exclude)
    visited.add(root)
    subtree = [root]
    queue = [root]
    head = 0
    while head < len(queue):
        v = queue[head]; head += 1
        for w in adj[v]:
            if w not in visited:
                visited.add(w)
                subtree.append(w)
                queue.append(w)
    return subtree


def main():
    t0 = time.time()

    # Statistics
    max_amplification = 0.0  # max |subtree_sum| / |ΔP(root)|
    max_amplification_info = {}
    amplification_histogram = defaultdict(int)

    # Track signed ratio: subtree_sum / ΔP(root)
    # If always negative (subtree opposes root), there's cancellation
    sign_same = 0
    sign_opposite = 0
    sign_zero = 0

    # Track total |subtree_sum| across all neighbors at dist 1
    max_total_subtree = 0.0
    max_remote_via_subtree = 0.0

    # Track the "direct bound": Σ |ΔP(w)| at dist 1 vs |remote|
    max_remote_over_dist1_sum = 0.0

    total_edges = 0

    for n in range(3, MAX_N + 1):
        tn = time.time()
        n_trees = 0
        n_max_amp = 0.0
        n_max_remote = 0.0

        for _, adj in trees(n):
            n_trees += 1
            msgs_T = cavity_messages(n, adj)
            P_T = occupation_probs(n, adj, msgs_T)

            seen = set()
            for u in range(n):
                for v in adj[u]:
                    e_key = (min(u, v), max(u, v))
                    if e_key in seen:
                        continue
                    seen.add(e_key)
                    total_edges += 1

                    # Contract
                    nc, adjc, o2n = contract_edge(n, adj, u, v)
                    msgs_Te = cavity_messages(nc, adjc)
                    P_Te = occupation_probs(nc, adjc, msgs_Te)

                    # Compute ΔP for each vertex
                    dP = {}
                    for w in range(n):
                        if w == u or w == v:
                            continue
                        w_new = o2n[w]
                        dP[w] = P_T[w] - P_Te[w_new]

                    # For each neighbor of {u,v} (at distance 1):
                    neighbors_of_edge = set()
                    for w in adj[u]:
                        if w != v:
                            neighbors_of_edge.add(w)
                    for w in adj[v]:
                        if w != u:
                            neighbors_of_edge.add(w)

                    dist1_sum = 0.0
                    for w in neighbors_of_edge:
                        dist1_sum += abs(dP.get(w, 0))

                    total_remote = sum(dP.values())
                    abs_remote = abs(total_remote)

                    if dist1_sum > 1e-15 and abs_remote > 1e-15:
                        ratio = abs_remote / dist1_sum
                        if ratio > max_remote_over_dist1_sum:
                            max_remote_over_dist1_sum = ratio

                    for w in neighbors_of_edge:
                        subtree = get_subtree(n, adj, w, {u, v})
                        subtree_sum = sum(dP.get(z, 0) for z in subtree)
                        root_dP = dP.get(w, 0)

                        if abs(root_dP) > 1e-15:
                            amp = abs(subtree_sum) / abs(root_dP)
                            bin_idx = min(int(amp * 10), 100)
                            amplification_histogram[bin_idx] += 1

                            if amp > max_amplification:
                                max_amplification = amp
                                max_amplification_info = {
                                    "n": n, "edge": [u, v], "w": w,
                                    "subtree_sum": subtree_sum,
                                    "root_dP": root_dP,
                                    "amp": amp,
                                    "subtree_size": len(subtree),
                                }

                            if amp > n_max_amp:
                                n_max_amp = amp

                            # Sign analysis
                            if subtree_sum * root_dP > 0:
                                sign_same += 1
                            elif subtree_sum * root_dP < 0:
                                sign_opposite += 1
                            else:
                                sign_zero += 1

                    if abs_remote > n_max_remote:
                        n_max_remote = abs_remote

        elapsed = time.time() - tn
        print(f"n={n}: {n_trees} trees, max_amp={n_max_amp:.4f}, "
              f"max|remote|={n_max_remote:.4f}, {elapsed:.1f}s", flush=True)

    total_time = time.time() - t0
    total_signs = sign_same + sign_opposite + sign_zero

    print(f"\n=== SUBTREE AMPLIFICATION ===")
    print(f"Max |subtree_sum| / |ΔP(root)|: {max_amplification:.6f}")
    print(f"Info: {max_amplification_info}")

    print(f"\n=== SIGN ANALYSIS ===")
    print(f"Same sign (subtree amplifies root): "
          f"{sign_same} ({100*sign_same/total_signs:.1f}%)")
    print(f"Opposite sign (subtree opposes root): "
          f"{sign_opposite} ({100*sign_opposite/total_signs:.1f}%)")
    print(f"Zero: {sign_zero}")

    print(f"\n=== REMOTE vs DIST-1 SUM ===")
    print(f"Max |remote| / Σ|ΔP(dist=1)|: {max_remote_over_dist1_sum:.6f}")
    print(f"(If < 1, the remote NEVER exceeds the dist-1 perturbation sum)")

    print(f"\n=== AMPLIFICATION HISTOGRAM ===")
    for bin_idx in range(min(20, max(amplification_histogram.keys()) + 1)):
        lo = bin_idx / 10.0
        hi = (bin_idx + 1) / 10.0
        cnt = amplification_histogram.get(bin_idx, 0)
        if cnt > 0:
            print(f"  [{lo:.1f}, {hi:.1f}): {cnt}")

    results = {
        "max_n": MAX_N,
        "total_edges": total_edges,
        "max_amplification": round(max_amplification, 8),
        "max_amplification_info": max_amplification_info,
        "sign_same": sign_same,
        "sign_opposite": sign_opposite,
        "sign_zero": sign_zero,
        "max_remote_over_dist1": round(max_remote_over_dist1_sum, 8),
        "total_time_s": round(total_time, 2),
    }

    out_path = "results/ecms_subtree_bound.json"
    os.makedirs("results", exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
