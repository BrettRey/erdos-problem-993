#!/usr/bin/env python3
"""Algebraic analysis of the caterpillar fixed point for remote bound.

On an infinite caterpillar with k leaves per spine vertex:
- Each spine vertex s has 2 spine neighbors + k leaf neighbors
- Cavity fixed point: R_{spine→spine} satisfies
    R = (1/2)^k / (1 + R)
  i.e., R(1+R) = 2^{-k}, giving R = (-1 + sqrt(1 + 2^{2-k})) / 2

- R_{spine→leaf} = (1/2)^{k-1} · (1/(1+R))^2
  (exclude one leaf and both spine neighbors' incoming messages)

Wait, let me be more careful. Spine vertex s has neighbors:
  - 2 spine vertices s_L, s_R
  - k leaf vertices l_1, ..., l_k

R_{s→s_R} = product over N(s)\{s_R}:
  = (1/(1+R_{s_L→s})) · product_{i=1}^k (1/(1+R_{l_i→s}))
  = (1/(1+R)) · (1/(1+1))^k
  = (1/(1+R)) · (1/2)^k

So R = (1/2)^k / (1+R), giving R(1+R) = (1/2)^k.

For contracting a spine-spine edge e = (s_L, s_R):
Local term: A = R_{s_L→s_R} = R (by symmetry), B = R_{s_R→s_L} = R
  local = (2R - R^2) / ((1+2R)(1+R^2))

The remote term involves all changes in P(v) for vertices not on the edge.
We compute this by:
1. Setting up the cavity messages on the infinite caterpillar
2. Contracting one spine edge and computing new messages on a finite
   approximation
3. Tracking the convergence as the caterpillar length grows

This gives exact algebraic expressions for the limiting remote.
"""

import json
import os
from fractions import Fraction
from math import sqrt


def caterpillar_fixed_point(k):
    """Compute the cavity fixed point for infinite caterpillar with k leaves/spine.

    Returns R = R_{spine→spine}, P_spine, P_leaf.
    """
    # R(1+R) = (1/2)^k
    # R^2 + R - 2^{-k} = 0
    # R = (-1 + sqrt(1 + 4·2^{-k})) / 2 = (-1 + sqrt(1 + 2^{2-k})) / 2
    val = 1.0 + 2.0 ** (2 - k)
    R = (-1.0 + sqrt(val)) / 2.0

    # P_spine = R_s / (1 + R_s) where R_s = full product
    # R_s = (1/(1+R))^2 · (1/2)^k  (both spine neighbors + k leaves)
    R_s = (1.0 / (1.0 + R)) ** 2 * (0.5) ** k
    P_spine = R_s / (1.0 + R_s)

    # P_leaf:
    # R_leaf = 1/(1+R_{s→leaf})
    # R_{s→leaf} = (1/(1+R))^2 · (1/2)^{k-1}  (2 spine nbrs + k-1 other leaves)
    R_s_to_leaf = (1.0 / (1.0 + R)) ** 2 * (0.5) ** (k - 1)
    R_leaf = 1.0 / (1.0 + R_s_to_leaf)
    P_leaf = R_leaf / (1.0 + R_leaf)

    return R, P_spine, P_leaf, R_s_to_leaf


def local_change(R):
    """Local change for contracting spine-spine edge with A=B=R."""
    A, B = R, R
    num = A + B - A * B
    den = (1 + A + B) * (1 + A * B)
    return num / den


def make_caterpillar(spine_len, k):
    """Build a caterpillar with spine_len spine vertices, k leaves each."""
    n = spine_len * (1 + k)
    adj = [[] for _ in range(n)]
    # Spine: 0, 1, ..., spine_len-1
    for i in range(spine_len - 1):
        adj[i].append(i + 1)
        adj[i + 1].append(i)
    # Leaves
    leaf_idx = spine_len
    for i in range(spine_len):
        for _ in range(k):
            adj[i].append(leaf_idx)
            adj[leaf_idx].append(i)
            leaf_idx += 1
    return n, adj


def cavity_messages(n, adj):
    """Compute cavity messages."""
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


def analyze_middle_edge(spine_len, k):
    """Contract the middle spine-spine edge and compute remote."""
    n, adj = make_caterpillar(spine_len, k)
    msgs = cavity_messages(n, adj)
    P_T = occupation_probs(n, adj, msgs)
    mu_T = sum(P_T)

    # Contract middle edge
    mid = spine_len // 2
    u, v = mid, mid + 1

    A = msgs.get((u, v), 0.0)
    B = msgs.get((v, u), 0.0)
    loc = local_change(A)  # Using actual A, B (not necessarily equal)

    nc, adjc, o2n = contract_edge(n, adj, u, v)
    msgs_c = cavity_messages(nc, adjc)
    P_Te = occupation_probs(nc, adjc, msgs_c)
    mu_Te = sum(P_Te)

    dmu = mu_T - mu_Te
    remote = dmu - loc

    return {
        "spine": spine_len, "k": k, "n": n,
        "A": A, "B": B, "local": loc, "remote": remote, "dmu": dmu,
        "mu_T": mu_T, "mu_Te": mu_Te,
    }


def main():
    print("=== FIXED POINT VALUES ===")
    for k in range(1, 11):
        R, P_s, P_l, R_sl = caterpillar_fixed_point(k)
        loc = local_change(R)
        mu_per_unit = P_s + k * P_l  # mean per spine vertex + its leaves
        print(f"k={k}: R={R:.6f}, P_spine={P_s:.6f}, P_leaf={P_l:.6f}, "
              f"local={loc:.6f}, μ/unit={mu_per_unit:.6f}")

    print("\n=== CONVERGENCE OF REMOTE FOR MIDDLE EDGE ===")
    for k in [1, 2, 3, 4, 5]:
        print(f"\n--- k={k} leaves per spine ---")
        for spine in [10, 20, 40, 80, 160, 320]:
            res = analyze_middle_edge(spine, k)
            print(f"  spine={spine:3d}, n={res['n']:4d}: "
                  f"A={res['A']:.6f}, B={res['B']:.6f}, "
                  f"local={res['local']:.6f}, remote={res['remote']:.6f}, "
                  f"δμ={res['dmu']:.6f}")

    # Now find the worst k
    print("\n=== WORST k (by remote) FOR SPINE=320 ===")
    results = []
    for k in range(1, 21):
        res = analyze_middle_edge(320, k)
        results.append(res)
        print(f"k={k:2d}: remote={res['remote']:.6f}, |δμ|={abs(res['dmu']):.6f}")

    # Find extremal edge across all spine positions for Cat(320, 3)
    print("\n=== ALL EDGES FOR Cat(320, 3) ===")
    spine = 320
    k = 3
    n, adj = make_caterpillar(spine, k)
    msgs = cavity_messages(n, adj)
    P_T = occupation_probs(n, adj, msgs)
    mu_T = sum(P_T)

    worst_remote_neg = 0
    worst_remote_pos = 0
    worst_dmu = 0
    worst_info = {}

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

            if remote < worst_remote_neg:
                worst_remote_neg = remote
                worst_info = {"edge": e, "A": A, "B": B, "local": loc,
                              "remote": remote, "dmu": dmu,
                              "u_is_spine": u < spine, "v_is_spine": v < spine}
            if remote > worst_remote_pos:
                worst_remote_pos = remote
            if abs(dmu) > abs(worst_dmu):
                worst_dmu = dmu

    print(f"Worst negative remote: {worst_remote_neg:.8f}")
    print(f"Worst positive remote: {worst_remote_pos:.8f}")
    print(f"Worst |δμ|: {abs(worst_dmu):.8f}")
    print(f"Worst edge info: {worst_info}")

    # Algebraic analysis of the limiting remote
    print("\n=== ALGEBRAIC LIMITING VALUES ===")
    for k in range(1, 11):
        R, P_s, P_l, R_sl = caterpillar_fixed_point(k)
        loc = local_change(R)
        # The limiting remote for the middle edge in an infinite caterpillar
        # should be very close to the spine=320 value
        res = analyze_middle_edge(320, k)
        print(f"k={k}: R_spine={R:.8f}, local={loc:.8f}, "
              f"remote_320={res['remote']:.8f}, δμ_320={res['dmu']:.8f}")


if __name__ == "__main__":
    main()
