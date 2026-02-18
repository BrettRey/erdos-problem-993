#!/usr/bin/env python3
"""Recursive bound on |remote| using the S_1 formula at each level.

Key idea: The S_1 formula gives ΔP(w) = F_edge × P_w^{cavity} for dist-1 vertices.
The perturbation entering each subtree is a change in cavity message.
At the next level, the SAME formula structure applies with modified parameters.

Define the "mean response function" M(subtree, δR) = Σ_{z in subtree} ΔP(z)
for a change δR in the cavity message entering the subtree root.

We want: |remote| = |Σ M(subtree_w, δR_w)| < 1/2.

Test whether M(subtree, δR) ≈ c × δR for some bounded constant c,
and whether |remote| ≤ 1/2 follows from the recursive structure.

Also test: define h(T,e) = |remote| / [|local| × correction(A,B)].
Find the tightest correction that makes h < 1 universally.
"""

import json
import math
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


def get_subtree_rooted(n, adj, root, parent_vertex):
    """Get subtree rooted at root, going away from parent_vertex."""
    visited = {parent_vertex, root}
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

    # Track mean response M_w = subtree_sum(w) per unit cavity perturbation
    # δR_w = change in cavity message entering w from the edge side
    max_M = 0.0
    max_M_info = {}

    # Track |remote| / (some bound)
    # Candidate bound: local × K for some constant K
    max_remote_over_local = 0.0

    # Better: |remote| ≤ Σ_w |ΔP(w)| × |M_w / ΔP_w|
    # where M_w / ΔP_w is the "amplification" from root to subtree

    # Check: |remote| ≤ g(A,B) × max_M_ratio
    # where max_M_ratio = max subtree_sum / root_dP
    max_remote_from_gM = 0.0

    # Direct: is |remote| < local always? (would be sufficient since local < 1/2)
    remote_exceeds_local = 0

    # Direct: is |remote| < f(A,B) for f = local × φ(A,B)?
    # where φ captures the subtree amplification
    # Test: |remote| < local + some_bound_on_S1_bound
    max_remote_plus_check = 0.0

    # The KEY test: compute per-subtree mean response δR → M
    # M(w) = subtree_sum(w) / δR(w)
    # where δR(w) = change in cavity message from edge-side to w
    M_values = []

    total_edges = 0

    for n in range(3, MAX_N + 1):
        tn = time.time()
        n_edges = 0

        for _, adj in trees(n):
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

                    A = msgs.get((u, v), 0.0)
                    B = msgs.get((v, u), 0.0)
                    denom = (1 + A + B) * (1 + A * B)
                    local = (A + B - A * B) / denom if denom > 0 else 0

                    nc, adjc, o2n = contract_edge(n, adj, u, v)
                    msgs_c = cavity_messages(nc, adjc)
                    P_Te = occupation_probs(nc, adjc, msgs_c)

                    dmu = sum(P_T) - sum(P_Te)
                    remote = dmu - local
                    abs_remote = abs(remote)

                    if abs(local) > 1e-15:
                        r = abs_remote / abs(local)
                        if r > max_remote_over_local:
                            max_remote_over_local = r
                    if abs_remote > abs(local):
                        remote_exceeds_local += 1

                    # Per-subtree mean response
                    for side, endpoint, other in [(u, v, B), (v, u, A)]:
                        for w in adj[side]:
                            if w == endpoint:
                                continue
                            # Compute δR entering w
                            q_w = msgs.get((w, side), 0.0)
                            # α_w = A_side × (1+q_w) / (1+B_side)
                            # β_w = A_side × (1+q_w) × B_side
                            # δR = β_w - α_w = A_side(1+q_w)(B²+B-1)/(1+B)
                            A_side = msgs.get((side, endpoint), 0.0)
                            B_side = msgs.get((endpoint, side), 0.0)
                            delta_R = (A_side * (1 + q_w) *
                                       (B_side**2 + B_side - 1) /
                                       (1 + B_side)) if (1 + B_side) > 0 else 0

                            if abs(delta_R) < 1e-15:
                                continue

                            sub = get_subtree_rooted(n, adj, w, side)
                            sub_sum = sum(P_T[z] - P_Te[o2n[z]]
                                          for z in sub if z in o2n)

                            M = sub_sum / delta_R
                            M_values.append(M)

                            if abs(M) > max_M:
                                max_M = abs(M)
                                max_M_info = {
                                    "n": n, "edge": (u, v), "w": w,
                                    "M": M, "delta_R": delta_R,
                                    "sub_sum": sub_sum,
                                    "subtree_size": len(sub),
                                }

        elapsed = time.time() - tn
        print(f"n={n}: {n_edges} edges, max|M|={max_M:.4f}, "
              f"max|rem/loc|={max_remote_over_local:.4f}, {elapsed:.1f}s",
              flush=True)

    total_time = time.time() - t0

    print(f"\n=== MEAN RESPONSE ===")
    print(f"Max |M| (subtree_sum / δR): {max_M:.8f}")
    print(f"Info: {max_M_info}")

    # Histogram of M values
    print(f"\n=== M DISTRIBUTION ===")
    M_abs = [abs(m) for m in M_values]
    for threshold in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]:
        frac = sum(1 for m in M_abs if m > threshold) / len(M_abs)
        print(f"  |M| > {threshold:.1f}: {100*frac:.2f}%")
    print(f"  Mean |M|: {sum(M_abs)/len(M_abs):.6f}")
    print(f"  Median |M|: {sorted(M_abs)[len(M_abs)//2]:.6f}")

    print(f"\n=== |remote| vs |local| ===")
    print(f"|remote| > |local| in {remote_exceeds_local}/{total_edges} "
          f"({100*remote_exceeds_local/total_edges:.2f}%)")
    print(f"Max |remote|/|local|: {max_remote_over_local:.6f}")

    # The bound: |remote| = |Σ M_w × δR_w|
    # ≤ Σ |M_w| × |δR_w|
    # ≤ max|M| × Σ |δR_w|
    # But Σ |δR_w| depends on edge structure...
    # Better: from the formula, δR_w ∝ F × (1+q_w) for each w.
    # So |remote| ≤ max|M| × |F| × Σ(1+q_w)
    # And Σ(1+q_w) = deg_side - 1 + Σ q_w
    print(f"\n=== BOUND ANALYSIS ===")
    print(f"If |M| ≤ C universally, then:")
    print(f"  |remote| ≤ C × Σ|δR_w|")
    print(f"  ≤ C × g(A,B) (since Σ|δR_w| relates to |S_1| via P_w^cav)")
    print(f"  This gives |remote| ≤ {max_M:.4f} × 0.355 = "
          f"{max_M * 0.355:.4f}")
    print(f"  {'< 1/2' if max_M * 0.355 < 0.5 else '>= 1/2'}")

    results = {
        "max_n": MAX_N,
        "total_edges": total_edges,
        "max_M": round(max_M, 8),
        "max_M_info": {k: (round(v, 8) if isinstance(v, float) else
                            list(v) if isinstance(v, tuple) else v)
                       for k, v in max_M_info.items()},
        "max_remote_over_local": round(max_remote_over_local, 6),
        "remote_exceeds_local": remote_exceeds_local,
        "M_stats": {
            "mean": round(sum(M_abs)/len(M_abs), 6),
            "median": round(sorted(M_abs)[len(M_abs)//2], 6),
            "max": round(max(M_abs), 6),
            "pct_gt_1": round(100*sum(1 for m in M_abs if m > 1)/len(M_abs), 2),
        },
        "total_time_s": round(total_time, 2),
    }

    out_path = "results/ecms_recursive_bound.json"
    os.makedirs("results", exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
