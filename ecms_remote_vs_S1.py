#!/usr/bin/env python3
"""Test whether |remote| ≤ |S_1| always (the tail after dist 1 always cancels).

If true, combined with the algebraic bound |S_1| ≤ g(A,B) < 0.355 < 1/2,
this would prove |remote| < 1/2 and hence ECMS.

Also test: what is the tightest constant c such that |remote| ≤ c × |S_1|?
And: what is the max ratio of the total subtree sum to the root perturbation?
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


def get_subtree(n, adj, root, exclude):
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

    # Test: |remote| ≤ |S_1|?
    remote_exceeds_S1 = 0
    max_remote_over_S1 = 0.0
    max_remote_over_S1_info = {}

    # Track max |remote| and max |S_1|
    max_abs_remote = 0.0
    max_abs_S1 = 0.0

    # For the full bound: is |remote| < g(A,B)?
    remote_exceeds_g = 0
    max_remote_over_g = 0.0

    # Subtree telescoping: does subtree_sum(w) / ΔP(w) always have |ratio| < C?
    max_subtree_ratio = 0.0

    # Does the tail have opposite sign to S_1?
    tail_opp_sign = 0
    tail_same_sign = 0
    tail_zero = 0

    # Tightest multiplicative bound
    max_abs_remote_value = 0.0
    max_abs_remote_info = {}

    total_edges = 0

    for n in range(3, MAX_N + 1):
        tn = time.time()
        n_edges = 0
        n_exceed = 0

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
                    F_u = A * (B * B + B - 1) / denom if denom > 0 else 0
                    F_v = B * (A * A + A - 1) / denom if denom > 0 else 0

                    nc, adjc, o2n = contract_edge(n, adj, u, v)
                    msgs_c = cavity_messages(nc, adjc)
                    P_Te = occupation_probs(nc, adjc, msgs_c)

                    dmu = sum(P_T) - sum(P_Te)
                    remote = dmu - local

                    # Compute S_1
                    S1 = 0.0
                    dist = bfs_distances(n, adj, [u, v])
                    for w in range(n):
                        if w == u or w == v:
                            continue
                        if dist[w] == 1:
                            w_new = o2n[w]
                            S1 += P_T[w] - P_Te[w_new]

                    abs_remote = abs(remote)
                    abs_S1 = abs(S1)

                    if abs_remote > max_abs_remote_value:
                        max_abs_remote_value = abs_remote
                        max_abs_remote_info = {
                            "n": n, "edge": (u, v),
                            "remote": remote, "S1": S1,
                            "A": A, "B": B, "local": local,
                        }

                    if abs_remote > max_abs_remote:
                        max_abs_remote = abs_remote
                    if abs_S1 > max_abs_S1:
                        max_abs_S1 = abs_S1

                    # Test: |remote| ≤ |S_1|?
                    if abs_S1 > 1e-15:
                        ratio = abs_remote / abs_S1
                        if ratio > max_remote_over_S1:
                            max_remote_over_S1 = ratio
                            max_remote_over_S1_info = {
                                "n": n, "edge": (u, v),
                                "remote": remote, "S1": S1, "ratio": ratio,
                                "A": A, "B": B,
                            }
                        if abs_remote > abs_S1 + 1e-12:
                            remote_exceeds_S1 += 1
                            n_exceed += 1

                    # Tail sign analysis
                    tail = remote - S1
                    if abs(tail) > 1e-15 and abs(S1) > 1e-15:
                        if tail * S1 < 0:
                            tail_opp_sign += 1
                        else:
                            tail_same_sign += 1
                    else:
                        tail_zero += 1

                    # Test: |remote| < g(A,B)?
                    if A > 0 and B > 0:
                        g = (abs(F_u) * (-math.log(A)) +
                             abs(F_v) * (-math.log(B)))
                        if abs_remote > g + 1e-10:
                            remote_exceeds_g += 1
                        if g > 1e-15:
                            r = abs_remote / g
                            if r > max_remote_over_g:
                                max_remote_over_g = r

                    # Subtree ratio analysis
                    for w in adj[u]:
                        if w == v:
                            continue
                        sub = get_subtree(n, adj, w, {u, v})
                        sub_sum = sum(P_T[z] - P_Te[o2n[z]]
                                      for z in sub if z in o2n)
                        dP_w = P_T[w] - P_Te[o2n[w]] if w in o2n else 0
                        if abs(dP_w) > 1e-15:
                            sr = abs(sub_sum) / abs(dP_w)
                            if sr > max_subtree_ratio:
                                max_subtree_ratio = sr
                    for w in adj[v]:
                        if w == u:
                            continue
                        sub = get_subtree(n, adj, w, {u, v})
                        sub_sum = sum(P_T[z] - P_Te[o2n[z]]
                                      for z in sub if z in o2n)
                        dP_w = P_T[w] - P_Te[o2n[w]] if w in o2n else 0
                        if abs(dP_w) > 1e-15:
                            sr = abs(sub_sum) / abs(dP_w)
                            if sr > max_subtree_ratio:
                                max_subtree_ratio = sr

        elapsed = time.time() - tn
        print(f"n={n}: {n_edges} edges, exceed={n_exceed}, "
              f"max_ratio={max_remote_over_S1:.6f}, {elapsed:.1f}s", flush=True)

    total_time = time.time() - t0
    total_sign = tail_opp_sign + tail_same_sign + tail_zero

    print(f"\n=== |remote| vs |S_1| ===")
    print(f"|remote| > |S_1| in {remote_exceeds_S1} / {total_edges} edges "
          f"({100*remote_exceeds_S1/total_edges:.2f}%)")
    print(f"Max |remote|/|S_1|: {max_remote_over_S1:.8f}")
    print(f"Info: {max_remote_over_S1_info}")
    print(f"Max |remote|: {max_abs_remote:.8f}")
    print(f"Max |S_1|: {max_abs_S1:.8f}")

    print(f"\n=== TAIL SIGN ===")
    print(f"Opposite to S_1 (cancels): {tail_opp_sign} ({100*tail_opp_sign/total_sign:.1f}%)")
    print(f"Same sign (adds): {tail_same_sign} ({100*tail_same_sign/total_sign:.1f}%)")

    print(f"\n=== |remote| vs g(A,B) ===")
    print(f"|remote| > g(A,B) in {remote_exceeds_g} edges")
    print(f"Max |remote|/g: {max_remote_over_g:.6f}")
    print(f"{'PASSES' if remote_exceeds_g == 0 else 'FAILS'}: "
          f"|remote| ≤ g(A,B) always")

    print(f"\n=== SUBTREE RATIO ===")
    print(f"Max |subtree_sum|/|ΔP(root)|: {max_subtree_ratio:.4f}")

    print(f"\n=== WORST REMOTE ===")
    print(f"Info: {max_abs_remote_info}")

    results = {
        "max_n": MAX_N,
        "total_edges": total_edges,
        "remote_exceeds_S1": remote_exceeds_S1,
        "max_remote_over_S1": round(max_remote_over_S1, 8),
        "max_abs_remote": round(max_abs_remote, 8),
        "max_abs_S1": round(max_abs_S1, 8),
        "tail_opp_sign": tail_opp_sign,
        "tail_same_sign": tail_same_sign,
        "remote_exceeds_g": remote_exceeds_g,
        "max_remote_over_g": round(max_remote_over_g, 6),
        "max_subtree_ratio": round(max_subtree_ratio, 4),
        "total_time_s": round(total_time, 2),
    }

    out_path = "results/ecms_remote_vs_S1.json"
    os.makedirs("results", exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
