#!/usr/bin/env python3
"""Test the mode sandwich lemma and its consequences for ECMS.

From I(T_e) = I(T) + x·I(T/e) and unimodality of both summands:
  mode(T_e) ∈ [min(m, m'+1), max(m, m'+1)]
where m = mode(I(T)), m' = mode(I(T/e)).

Also test: for each tree T, subdivide every edge e.
  (a) mode(T_e) ∈ [min(m, m'+1), max(m, m'+1)]?  (must hold)
  (b) |mode(T_e) - mode(T)| ≤ 1?  (subdivision mode stability)
  (c) |mode(T_e) - mode(T/e)| ≤ 1?  (another ECMS instance)

If both (b) and (c) hold, then from (b): m_e ∈ {m-1, m, m+1}.
Combined with (a), this constrains m' (and hence ECMS).
"""

import json
import os
import time

from indpoly import independence_poly
from trees import trees

MAX_N = 18


def poly_mode(poly):
    return max(range(len(poly)), key=lambda k: poly[k])


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
    return n_new, adj_new


def subdivide_edge(n, adj, u, v):
    n_new = n + 1
    w = n
    adj_new = [list(nbrs) for nbrs in adj]
    adj_new.append([])
    adj_new[u] = [x if x != v else w for x in adj_new[u]]
    adj_new[v] = [x if x != u else w for x in adj_new[v]]
    adj_new[w] = [u, v]
    return n_new, adj_new


def main():
    t0 = time.time()

    stats = {
        "sandwich_holds": 0,
        "sandwich_fails": 0,
        "subdiv_mode_shift": {},  # mode(T_e) - mode(T) distribution
        "ecms_Te_Te_div_e": {},   # mode(T_e) - mode(T/e) distribution (= ECMS for T_e)
        "ecms_T_Te": {},          # mode(T) - mode(T/e) distribution (= ECMS)
        # Key test: if |mode(T_e) - mode(T)| ≤ 1, does ECMS follow?
        "subdiv_ecms_implies_ecms": True,
        "total_edges": 0,
    }

    for n in range(3, MAX_N + 1):
        tn = time.time()
        n_trees = 0
        n_edges = 0

        for _, adj in trees(n):
            n_trees += 1
            poly_T = independence_poly(n, adj)
            m = poly_mode(poly_T)

            seen = set()
            for u in range(n):
                for v in adj[u]:
                    e = (min(u, v), max(u, v))
                    if e in seen:
                        continue
                    seen.add(e)
                    n_edges += 1
                    stats["total_edges"] += 1

                    # Contraction
                    nc, adjc = contract_edge(n, adj, u, v)
                    poly_Te = independence_poly(nc, adjc)
                    m_prime = poly_mode(poly_Te)

                    # Subdivision
                    ns, adjs = subdivide_edge(n, adj, u, v)
                    poly_sub = independence_poly(ns, adjs)
                    m_e = poly_mode(poly_sub)

                    # Mode sandwich: mode(T_e) ∈ [min(m, m'+1), max(m, m'+1)]
                    lo = min(m, m_prime + 1)
                    hi = max(m, m_prime + 1)
                    if lo <= m_e <= hi:
                        stats["sandwich_holds"] += 1
                    else:
                        stats["sandwich_fails"] += 1

                    # Subdivision mode shift
                    d1 = m_e - m
                    stats["subdiv_mode_shift"][d1] = \
                        stats["subdiv_mode_shift"].get(d1, 0) + 1

                    # ECMS for T_e: |mode(T_e) - mode(T/e)|
                    d2 = m_e - m_prime
                    stats["ecms_Te_Te_div_e"][d2] = \
                        stats["ecms_Te_Te_div_e"].get(d2, 0) + 1

                    # Standard ECMS
                    d3 = m - m_prime
                    stats["ecms_T_Te"][d3] = \
                        stats["ecms_T_Te"].get(d3, 0) + 1

                    # Does subdivision mode stability + sandwich → ECMS?
                    if abs(d1) <= 1:  # |mode(T_e) - mode(T)| ≤ 1
                        # From sandwich: m_e ∈ [min(m,m'+1), max(m,m'+1)]
                        # From subdiv stability: m_e ∈ {m-1, m, m+1}
                        # Check if this constrains |m - m'| ≤ 1
                        if abs(m - m_prime) > 1:
                            stats["subdiv_ecms_implies_ecms"] = False

        elapsed = time.time() - tn
        print(f"n={n}: {n_trees} trees, {n_edges} edges, {elapsed:.1f}s",
              flush=True)

    stats["total_time_s"] = round(time.time() - t0, 2)

    print(f"\n=== MODE SANDWICH ===")
    print(f"Holds: {stats['sandwich_holds']}")
    print(f"Fails: {stats['sandwich_fails']}")

    print(f"\n=== SUBDIVISION MODE SHIFT: mode(T_e) - mode(T) ===")
    for k in sorted(stats["subdiv_mode_shift"].keys()):
        v = stats["subdiv_mode_shift"][k]
        print(f"  {k:+d}: {v} ({100*v/stats['total_edges']:.1f}%)")

    print(f"\n=== ECMS: mode(T) - mode(T/e) ===")
    for k in sorted(stats["ecms_T_Te"].keys()):
        v = stats["ecms_T_Te"][k]
        print(f"  {k:+d}: {v} ({100*v/stats['total_edges']:.1f}%)")

    print(f"\n=== ECMS for T_e: mode(T_e) - mode(T/e) ===")
    for k in sorted(stats["ecms_Te_Te_div_e"].keys()):
        v = stats["ecms_Te_Te_div_e"][k]
        print(f"  {k:+d}: {v} ({100*v/stats['total_edges']:.1f}%)")

    print(f"\nSubdiv stability + sandwich → ECMS: "
          f"{stats['subdiv_ecms_implies_ecms']}")
    print(f"Total: {stats['total_edges']} edges, {stats['total_time_s']}s")

    out_path = "results/ecms_mode_sandwich.json"
    os.makedirs("results", exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"Saved to {out_path}")


if __name__ == "__main__":
    main()
