#!/usr/bin/env python3
"""Analyze ECMS straddling cases: when μ(T) and μ(T/e) straddle an integer.

Key algebraic identity (from I(T_e) = I(T) + x·I(T/e)):
    α - β = 1 + (μ_e - 1 - β) / p
where α=μ(T), β=μ(T/e), μ_e=μ(T_e), p=I(T;1)/I(T_e;1).

This means α - β < 1  ⟺  α < 1 + β  ⟺  μ_e < 1 + β.

For LC polynomials, mode ∈ {⌊μ⌋, ⌈μ⌉}. If |δμ| < 1 but α,β straddle
integer k, then mode(T) ∈ {k, k+1} and mode(T/e) ∈ {k-1, k}. The gap
case |Δmode|=2 requires mode(T)=k+1 AND mode(T/e)=k-1 simultaneously.

This script checks: does this gap case ever occur?
"""

import json
import math
import sys
import time
from fractions import Fraction

from indpoly import independence_poly, is_log_concave
from trees import trees

MAX_N = 18
PROGRESS_INTERVAL = 5000


def poly_mean(poly):
    """Compute μ = I'(1)/I(1)."""
    s = sum(poly)
    if s == 0:
        return 0.0
    return sum(k * poly[k] for k in range(len(poly))) / s


def poly_mean_exact(poly):
    """Compute μ exactly as a Fraction."""
    s = sum(poly)
    if s == 0:
        return Fraction(0)
    return Fraction(sum(k * poly[k] for k in range(len(poly))), s)


def poly_mode(poly):
    """Return the mode (index of max coefficient)."""
    return max(range(len(poly)), key=lambda k: poly[k])


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

    return n_new, adj_new


def subdivide_edge(n, adj, u, v):
    """Subdivide edge (u,v): add new vertex w between u and v."""
    n_new = n + 1
    w = n  # new vertex index
    adj_new = [list(nbrs) for nbrs in adj]
    adj_new.append([])  # vertex w

    # Remove u-v edge, add u-w and w-v
    adj_new[u] = [x if x != v else w for x in adj_new[u]]
    adj_new[v] = [x if x != u else w for x in adj_new[v]]
    adj_new[w] = [u, v]

    return n_new, adj_new


def main():
    results = {
        "max_n": MAX_N,
        "identity_verified": True,
        "max_identity_err": 0.0,
        "straddling_cases": 0,
        "straddling_mode_gap_2": 0,  # KEY: does |Δmode|=2 ever occur in straddling?
        "straddling_examples": [],
        "alpha_lt_1_plus_beta": True,
        "beta_lt_1_plus_alpha": True,
        "max_alpha_minus_beta": 0.0,
        "max_beta_minus_alpha": 0.0,
        "total_edges": 0,
        "per_n": {},
    }

    # Detailed straddling analysis
    straddling_mode_stats = {
        "both_at_k": 0,      # mode(T)=k, mode(T/e)=k
        "T_up": 0,           # mode(T)=k+1, mode(T/e)=k  (or k, k-1)
        "Te_up": 0,          # mode(T)=k, mode(T/e)=k-1  (reversed)
        "gap_2": 0,          # mode(T)=k+1, mode(T/e)=k-1 -- THE GAP CASE
    }

    t0 = time.time()
    total_edges = 0

    for n in range(2, MAX_N + 1):
        tn = time.time()
        n_trees = 0
        n_edges = 0
        n_straddle = 0
        n_gap2 = 0
        max_dmu_pos = 0.0
        max_dmu_neg = 0.0

        for _, adj in trees(n):
            n_trees += 1
            poly_T = independence_poly(n, adj)
            mu_T = poly_mean(poly_T)
            mode_T = poly_mode(poly_T)
            lc_T = is_log_concave(poly_T)
            I_T = sum(poly_T)

            # Process each edge (avoid duplicates)
            seen_edges = set()
            for u in range(n):
                for v in adj[u]:
                    if (min(u, v), max(u, v)) in seen_edges:
                        continue
                    seen_edges.add((min(u, v), max(u, v)))
                    n_edges += 1
                    total_edges += 1

                    # Contract edge
                    nc, adjc = contract_edge(n, adj, u, v)
                    poly_Te = independence_poly(nc, adjc)
                    mu_Te = poly_mean(poly_Te)
                    mode_Te = poly_mode(poly_Te)
                    lc_Te = is_log_concave(poly_Te)
                    I_Te = sum(poly_Te)

                    # Subdivide edge (to verify identity)
                    ns, adjs = subdivide_edge(n, adj, u, v)
                    poly_sub = independence_poly(ns, adjs)
                    mu_sub = poly_mean(poly_sub)
                    I_sub = sum(poly_sub)

                    alpha = mu_T
                    beta = mu_Te

                    # Verify algebraic identity: α - β = 1 + (μ_e - 1 - β)/p
                    # where p = I(T;1)/I(T_e;1)
                    p = I_T / I_sub
                    identity_lhs = alpha - beta
                    identity_rhs = 1 + (mu_sub - 1 - beta) / p
                    identity_err = abs(identity_lhs - identity_rhs)
                    if identity_err > results["max_identity_err"]:
                        results["max_identity_err"] = identity_err
                    if identity_err > 1e-10:
                        results["identity_verified"] = False

                    # Check |δμ| < 1 (both directions)
                    dmu = alpha - beta
                    if dmu > max_dmu_pos:
                        max_dmu_pos = dmu
                    if -dmu > max_dmu_neg:
                        max_dmu_neg = -dmu

                    if alpha >= 1 + beta:
                        results["alpha_lt_1_plus_beta"] = False
                    if beta >= 1 + alpha:
                        results["beta_lt_1_plus_alpha"] = False

                    if dmu > results["max_alpha_minus_beta"]:
                        results["max_alpha_minus_beta"] = dmu
                    if -dmu > results["max_beta_minus_alpha"]:
                        results["max_beta_minus_alpha"] = -dmu

                    # Check for straddling: ⌊α⌋ ≠ ⌊β⌋
                    floor_alpha = math.floor(alpha)
                    floor_beta = math.floor(beta)

                    if floor_alpha != floor_beta and lc_T and lc_Te:
                        n_straddle += 1
                        results["straddling_cases"] += 1

                        # The straddled integer
                        k = max(floor_alpha, floor_beta)

                        # Classify mode positions relative to straddled integer
                        # If α > β (α straddles above k, β below):
                        #   mode(T) ∈ {k, k+1}, mode(T/e) ∈ {k-1, k}
                        # Gap case: mode(T)=k+1, mode(T/e)=k-1
                        delta_mode = mode_T - mode_Te
                        abs_delta_mode = abs(delta_mode)

                        if abs_delta_mode >= 2:
                            n_gap2 += 1
                            results["straddling_mode_gap_2"] += 1

                        # Classify straddling case
                        if alpha > beta:
                            # α above k, β below k
                            if mode_T == k and mode_Te == k:
                                straddling_mode_stats["both_at_k"] += 1
                            elif mode_T >= k and mode_Te >= k:
                                straddling_mode_stats["both_at_k"] += 1
                            elif mode_T > mode_Te:
                                straddling_mode_stats["T_up"] += 1
                            elif mode_Te > mode_T:
                                straddling_mode_stats["Te_up"] += 1

                            if mode_T == k + 1 and mode_Te == k - 1:
                                straddling_mode_stats["gap_2"] += 1
                        else:
                            # β above k, α below k
                            if mode_T == k and mode_Te == k:
                                straddling_mode_stats["both_at_k"] += 1
                            elif mode_Te > mode_T:
                                straddling_mode_stats["Te_up"] += 1
                            elif mode_T > mode_Te:
                                straddling_mode_stats["T_up"] += 1

                            if mode_Te == k + 1 and mode_T == k - 1:
                                straddling_mode_stats["gap_2"] += 1

                        # Save examples (first 50 straddling cases total)
                        if len(results["straddling_examples"]) < 50:
                            results["straddling_examples"].append({
                                "n": n,
                                "edge": [u, v],
                                "alpha": round(alpha, 6),
                                "beta": round(beta, 6),
                                "floor_alpha": floor_alpha,
                                "floor_beta": floor_beta,
                                "mode_T": mode_T,
                                "mode_Te": mode_Te,
                                "delta_mode": delta_mode,
                                "lc_T": lc_T,
                                "lc_Te": lc_Te,
                            })

            if n_trees % PROGRESS_INTERVAL == 0 and n_trees > 0:
                print(f"  n={n}: {n_trees} trees, {n_edges} edges...",
                      flush=True)

        elapsed = time.time() - tn
        results["per_n"][str(n)] = {
            "trees": n_trees,
            "edges": n_edges,
            "straddling": n_straddle,
            "gap_2": n_gap2,
            "max_dmu_pos": round(max_dmu_pos, 10),
            "max_dmu_neg": round(max_dmu_neg, 10),
            "time_s": round(elapsed, 2),
        }
        pct_straddle = 100 * n_straddle / max(n_edges, 1)
        print(f"n={n}: {n_trees} trees, {n_edges} edges, "
              f"{n_straddle} straddling ({pct_straddle:.1f}%), "
              f"{n_gap2} gap-2, {elapsed:.1f}s", flush=True)

    results["total_edges"] = total_edges
    results["total_time_s"] = round(time.time() - t0, 2)
    results["straddling_mode_stats"] = straddling_mode_stats

    # Summary
    print(f"\n=== SUMMARY ===")
    print(f"Total edges: {total_edges}")
    print(f"Identity verified: {results['identity_verified']} "
          f"(max err: {results['max_identity_err']:.2e})")
    print(f"α < 1+β always: {results['alpha_lt_1_plus_beta']}")
    print(f"β < 1+α always: {results['beta_lt_1_plus_alpha']}")
    print(f"max(α-β) = {results['max_alpha_minus_beta']:.6f}")
    print(f"max(β-α) = {results['max_beta_minus_alpha']:.6f}")
    print(f"Straddling cases: {results['straddling_cases']}")
    print(f"Straddling with |Δmode|=2: {results['straddling_mode_gap_2']}")
    print(f"Straddling mode stats: {straddling_mode_stats}")

    out_path = "results/ecms_straddling.json"
    os.makedirs("results", exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


import os
if __name__ == "__main__":
    main()
