#!/usr/bin/env python3
"""Test proof-ready inequalities for the alternating remote bound.

Key hypothesis: the signed sums S_d at each distance d from the contracted
edge form an alternating series (S_d and S_{d+1} have opposite sign), with
magnitudes decaying geometrically. If this holds strictly, then:

  |remote| = |Σ S_d| ≤ |S_1|

So we need: (1) alternating dominance (|partial sum| ≤ |S_1| for all d),
and (2) |S_1| < 1/2.

This script checks both conditions exhaustively and decomposes S_1 into
u-side and v-side contributions to find what bounds it.
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

    # Test 1: Alternating dominance -- |partial(d)| ≤ |S_1| for all d?
    alt_dom_failures = 0
    alt_dom_total = 0

    # Test 2: |S_1| < 1/2 always?
    max_abs_S1 = 0.0
    max_S1_info = {}

    # Test 3: Decompose S_1 into contributions from N(u)\v and N(v)\u
    max_S1_uside = 0.0
    max_S1_vside = 0.0

    # Test 4: Per-vertex sign alternation -- does sgn(ΔP(w)) alternate with dist?
    sign_alt_vertex_total = 0
    sign_alt_vertex_match = 0

    # Test 5: |S_1| as function of cavity messages A, B
    # S_1 vs local decomposition
    max_S1_over_local = 0.0

    # Test 6: Tighter bound -- |S_1| ≤ |local|? Or |S_1| ≤ some f(A,B)?
    S1_exceeds_local = 0
    S1_vs_local_max_ratio = 0.0

    # Track S_1 distributions
    S1_histogram = defaultdict(int)

    total_edges = 0

    for n in range(3, MAX_N + 1):
        tn = time.time()
        n_edges = 0
        n_fail = 0

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

                    # Cavity messages for this edge
                    A = msgs.get((u, v), 0.0)
                    B = msgs.get((v, u), 0.0)

                    # Local term
                    num = A + B - A * B
                    den = (1 + A + B) * (1 + A * B)
                    local = num / den if den > 0 else 0.0

                    # Contract and compute ΔP
                    nc, adjc, o2n = contract_edge(n, adj, u, v)
                    msgs_c = cavity_messages(nc, adjc)
                    P_Te = occupation_probs(nc, adjc, msgs_c)

                    dist = bfs_distances(n, adj, [u, v])

                    # Compute per-distance signed sums
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

                    if max_dist == 0:
                        continue

                    S1 = signed_by_d.get(1, 0.0)
                    abs_S1 = abs(S1)

                    # Test 2: |S_1| < 1/2
                    if abs_S1 > max_abs_S1:
                        max_abs_S1 = abs_S1
                        max_S1_info = {"n": n, "edge": (u, v),
                                       "S1": S1, "A": A, "B": B,
                                       "local": local}

                    # S1 histogram
                    bin_idx = min(int(abs_S1 * 100), 100)
                    S1_histogram[bin_idx] += 1

                    # Test 1: Alternating dominance
                    partial = 0.0
                    for d in range(1, max_dist + 1):
                        s = signed_by_d.get(d, 0.0)
                        partial += s
                        alt_dom_total += 1
                        if abs(partial) > abs_S1 + 1e-12:
                            alt_dom_failures += 1
                            n_fail += 1

                    # Test 3: Decompose S_1
                    S1_u = 0.0  # from neighbors of u (excluding v)
                    S1_v = 0.0  # from neighbors of v (excluding u)
                    for w in adj[u]:
                        if w == v:
                            continue
                        if w in o2n:
                            S1_u += P_T[w] - P_Te[o2n[w]]
                    for w in adj[v]:
                        if w == u:
                            continue
                        if w in o2n:
                            S1_v += P_T[w] - P_Te[o2n[w]]
                    if abs(S1_u) > max_S1_uside:
                        max_S1_uside = abs(S1_u)
                    if abs(S1_v) > max_S1_vside:
                        max_S1_vside = abs(S1_v)

                    # Test 5 & 6: S_1 vs local
                    if abs(local) > 1e-15:
                        ratio = abs_S1 / abs(local)
                        if ratio > S1_vs_local_max_ratio:
                            S1_vs_local_max_ratio = ratio
                        if abs_S1 > abs(local):
                            S1_exceeds_local += 1
                        if ratio > max_S1_over_local:
                            max_S1_over_local = ratio

                    # Test 4: Per-vertex sign alternation
                    for w in range(n):
                        if w == u or w == v:
                            continue
                        d = dist[w]
                        if d >= 2:
                            w_new = o2n[w]
                            dP = P_T[w] - P_Te[w_new]
                            # Expected sign: opposite of parity of d
                            # At d=1, ΔP is typically negative (contraction
                            # reduces occupancy of neighbors)
                            sign_alt_vertex_total += 1
                            # Check if sign matches parity
                            expected_neg = (d % 2 == 1)  # odd dist → same as d=1
                            actual_neg = (dP < -1e-15)
                            actual_pos = (dP > 1e-15)
                            if (expected_neg and actual_neg) or \
                               (not expected_neg and actual_pos):
                                sign_alt_vertex_match += 1

        elapsed = time.time() - tn
        print(f"n={n}: {n_edges} edges, alt_dom_fail={n_fail}, "
              f"max|S1|={max_abs_S1:.6f}, {elapsed:.1f}s", flush=True)

    total_time = time.time() - t0

    print(f"\n=== TEST 1: ALTERNATING DOMINANCE ===")
    print(f"|partial(d)| ≤ |S_1| for all d?")
    print(f"Total checks: {alt_dom_total}")
    print(f"Failures: {alt_dom_failures} "
          f"({100*alt_dom_failures/alt_dom_total:.4f}%)" if alt_dom_total > 0 else "")

    print(f"\n=== TEST 2: |S_1| < 1/2 ===")
    print(f"Max |S_1|: {max_abs_S1:.10f}")
    print(f"Info: {max_S1_info}")
    print(f"{'PASSES' if max_abs_S1 < 0.5 else 'FAILS'}: |S_1| < 1/2")

    print(f"\n=== TEST 3: S_1 DECOMPOSITION ===")
    print(f"Max |S_1 from u-side|: {max_S1_uside:.8f}")
    print(f"Max |S_1 from v-side|: {max_S1_vside:.8f}")

    print(f"\n=== TEST 4: PER-VERTEX SIGN ALTERNATION ===")
    if sign_alt_vertex_total > 0:
        print(f"Vertices at dist ≥ 2 matching expected sign: "
              f"{sign_alt_vertex_match}/{sign_alt_vertex_total} "
              f"({100*sign_alt_vertex_match/sign_alt_vertex_total:.1f}%)")

    print(f"\n=== TEST 5-6: S_1 vs LOCAL ===")
    print(f"Max |S_1| / |local|: {S1_vs_local_max_ratio:.6f}")
    print(f"|S_1| > |local| in {S1_exceeds_local} cases")

    print(f"\n=== |S_1| HISTOGRAM ===")
    for b in range(min(50, max(S1_histogram.keys()) + 1)):
        lo = b / 100.0
        hi = (b + 1) / 100.0
        cnt = S1_histogram.get(b, 0)
        if cnt > 0:
            print(f"  [{lo:.2f}, {hi:.2f}): {cnt}")

    results = {
        "max_n": MAX_N,
        "total_edges": total_edges,
        "alt_dom_total": alt_dom_total,
        "alt_dom_failures": alt_dom_failures,
        "max_abs_S1": round(max_abs_S1, 10),
        "max_S1_info": {k: (round(v, 8) if isinstance(v, float) else
                            list(v) if isinstance(v, tuple) else v)
                        for k, v in max_S1_info.items()},
        "max_S1_uside": round(max_S1_uside, 8),
        "max_S1_vside": round(max_S1_vside, 8),
        "S1_vs_local_max_ratio": round(S1_vs_local_max_ratio, 6),
        "S1_exceeds_local": S1_exceeds_local,
        "vertex_sign_alt_match_pct": round(
            100 * sign_alt_vertex_match / sign_alt_vertex_total, 2)
        if sign_alt_vertex_total > 0 else 0,
        "total_time_s": round(total_time, 2),
    }

    out_path = "results/ecms_alternating_proof.json"
    os.makedirs("results", exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
