#!/usr/bin/env python3
"""Extend ECMS mean shift test to n=17 and n=18.

Previous result (ecms_mean_shift.py, n<=16): max |delta_mu| = 0.5263.
This script tests n=17 (48,629 trees) and n=18 (123,867 trees).

For each tree T and each edge e = {u,v}:
  1. Compute I(T) and I(T/e) using independence_poly
  2. Compute mu(T) = sum(k*i_k) / sum(i_k) using exact Fraction arithmetic
  3. Track max |delta_mu| and ECMS violations (|Delta_mode| > 1)
"""

import json
import os
import sys
import time
from collections import defaultdict
from fractions import Fraction

from indpoly import independence_poly, is_log_concave
from trees import trees


PROGRESS_INTERVAL = 5000


def poly_mean_exact(poly):
    """Compute mu = I'(1)/I(1) using exact Fraction arithmetic."""
    i_sum = sum(poly)
    if i_sum == 0:
        return Fraction(0)
    i_prime = sum(k * poly[k] for k in range(len(poly)))
    return Fraction(i_prime, i_sum)


def poly_mode(poly):
    """Return the index of the maximum coefficient (mode)."""
    return max(range(len(poly)), key=lambda k: poly[k])


def contract_edge(n, adj, u, v):
    """Contract edge (u, v): merge v into u, redirect v's neighbors, remove v.

    Returns (n_new, adj_new) with vertices renumbered 0..n-2.
    """
    # Collect neighbors of merged vertex
    merged_neighbors = set()
    for w in adj[u]:
        if w != v:
            merged_neighbors.add(w)
    for w in adj[v]:
        if w != u:
            merged_neighbors.add(w)

    # Mapping: old -> new (v removed, vertices > v shift down)
    old_to_new = {}
    for i in range(n):
        if i == v:
            continue
        old_to_new[i] = i if i < v else i - 1

    n_new = n - 1
    adj_new = [[] for _ in range(n_new)]

    # Merged vertex edges
    u_new = old_to_new[u]
    for w in sorted(merged_neighbors):
        w_new = old_to_new[w]
        adj_new[u_new].append(w_new)
        adj_new[w_new].append(u_new)

    # All other edges (not involving u or v)
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


def get_edges(n, adj):
    """Extract edges as list of (u, v) with u < v."""
    edges = []
    for u in range(n):
        for v in adj[u]:
            if u < v:
                edges.append((u, v))
    return edges


def main():
    print("ECMS Mean Shift Test: extending to n=17 and n=18")
    print(f"Using exact Fraction arithmetic for mean computation")
    print(f"Progress every {PROGRESS_INTERVAL} trees")
    print()

    total_trees = 0
    total_edges = 0

    max_abs_dmu = Fraction(0)
    max_abs_dmu_info = None

    # |delta_mu| >= 1 tracker
    dmu_ge_1_count = 0
    dmu_ge_1_examples = []

    # ECMS: |Delta_mode| > 1
    ecms_violations = 0
    ecms_violation_examples = []
    max_mode_shift = 0
    max_mode_shift_info = None

    # Histogram of delta_mu in 0.05 buckets
    hist = defaultdict(int)

    per_n_stats = {}

    t_start = time.time()

    for n in [17, 18]:
        t_n = time.time()
        n_trees = 0
        n_edges = 0
        n_max_abs_dmu = Fraction(0)
        n_max_mode_shift = 0
        n_dmu_ge_1 = 0
        n_ecms_viol = 0

        for tree_n, adj in trees(n, backend="geng"):
            n_trees += 1
            total_trees += 1

            edges = get_edges(tree_n, adj)
            poly_T = independence_poly(tree_n, adj)
            mu_T = poly_mean_exact(poly_T)
            mode_T = poly_mode(poly_T)

            for u, v in edges:
                n_edges += 1
                total_edges += 1

                # Contract edge
                n_c, adj_c = contract_edge(tree_n, adj, u, v)
                poly_Te = independence_poly(n_c, adj_c)
                mu_Te = poly_mean_exact(poly_Te)
                mode_Te = poly_mode(poly_Te)

                dmu = mu_T - mu_Te
                abs_dmu = abs(dmu)

                # Histogram bucket (0.05 resolution)
                dmu_float = float(dmu)
                bucket = round(dmu_float * 20) / 20  # round to nearest 0.05
                hist[bucket] += 1

                # Track max |delta_mu|
                if abs_dmu > n_max_abs_dmu:
                    n_max_abs_dmu = abs_dmu

                if abs_dmu > max_abs_dmu:
                    max_abs_dmu = abs_dmu
                    max_abs_dmu_info = {
                        "n": tree_n,
                        "edge": [u, v],
                        "dmu": str(dmu),
                        "dmu_float": float(dmu),
                        "mu_T": str(mu_T),
                        "mu_T_float": float(mu_T),
                        "mu_Te": str(mu_Te),
                        "mu_Te_float": float(mu_Te),
                        "mode_T": mode_T,
                        "mode_Te": mode_Te,
                    }

                # Check |delta_mu| >= 1
                if abs_dmu >= 1:
                    dmu_ge_1_count += 1
                    n_dmu_ge_1 += 1
                    if len(dmu_ge_1_examples) < 10:
                        dmu_ge_1_examples.append({
                            "n": tree_n,
                            "edge": [u, v],
                            "dmu": float(dmu),
                            "mu_T": float(mu_T),
                            "mu_Te": float(mu_Te),
                        })

                # ECMS cross-check
                mode_shift = abs(mode_T - mode_Te)
                if mode_shift > n_max_mode_shift:
                    n_max_mode_shift = mode_shift

                if mode_shift > max_mode_shift:
                    max_mode_shift = mode_shift
                    max_mode_shift_info = {
                        "n": tree_n,
                        "edge": [u, v],
                        "mode_T": mode_T,
                        "mode_Te": mode_Te,
                        "shift": mode_shift,
                    }

                if mode_shift > 1:
                    ecms_violations += 1
                    n_ecms_viol += 1
                    if len(ecms_violation_examples) < 10:
                        ecms_violation_examples.append({
                            "n": tree_n,
                            "edge": [u, v],
                            "mode_T": mode_T,
                            "mode_Te": mode_Te,
                        })

            # Progress
            if n_trees % PROGRESS_INTERVAL == 0:
                elapsed = time.time() - t_start
                print(f"  [{elapsed:7.1f}s] n={n}: {n_trees:>10,} trees, "
                      f"{n_edges:>12,} edges, max |dmu| = {float(n_max_abs_dmu):.6f}",
                      flush=True)

        elapsed_n = time.time() - t_n
        per_n_stats[n] = {
            "trees": n_trees,
            "edges": n_edges,
            "max_abs_dmu": float(n_max_abs_dmu),
            "max_abs_dmu_exact": str(n_max_abs_dmu),
            "max_mode_shift": n_max_mode_shift,
            "dmu_ge_1": n_dmu_ge_1,
            "ecms_violations": n_ecms_viol,
            "time_s": round(elapsed_n, 2),
        }

        print(f"n={n:>3}: {n_trees:>8,} trees, {n_edges:>10,} edges, "
              f"max |dmu| = {float(n_max_abs_dmu):.6f}, max |Dmode| = {n_max_mode_shift}, "
              f"|dmu|>=1: {n_dmu_ge_1}, ECMS viol: {n_ecms_viol} "
              f"[{elapsed_n:.1f}s]", flush=True)

    total_time = time.time() - t_start

    # Summary
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Trees examined:    {total_trees:,}")
    print(f"Edges examined:    {total_edges:,}")
    print(f"Total time:        {total_time:.1f}s")
    print()
    print(f"Max |delta_mu|:    {float(max_abs_dmu):.8f}")
    print(f"  Exact:           {max_abs_dmu}")
    if max_abs_dmu_info:
        info = max_abs_dmu_info
        print(f"  Achieved at:     n={info['n']}, edge={info['edge']}")
        print(f"  mu(T) = {info['mu_T_float']:.8f}, mu(T/e) = {info['mu_Te_float']:.8f}")
        print(f"  delta_mu = {info['dmu_float']:.8f}")
        print(f"  mode(T) = {info['mode_T']}, mode(T/e) = {info['mode_Te']}")
    print()
    print(f"|delta_mu| < 1:    {'HOLDS' if dmu_ge_1_count == 0 else 'FAILS'}")
    print(f"  Violations:      {dmu_ge_1_count:,}")
    if dmu_ge_1_examples:
        print(f"  First examples:")
        for ex in dmu_ge_1_examples[:5]:
            print(f"    n={ex['n']}, edge={ex['edge']}, dmu={ex['dmu']:.6f}")
    print()
    print(f"Max |Delta_mode|:  {max_mode_shift}")
    print(f"ECMS (|Dmode|<=1): {'HOLDS' if ecms_violations == 0 else 'FAILS'}")
    print(f"  Violations:      {ecms_violations:,}")
    if ecms_violation_examples:
        print(f"  First examples:")
        for ex in ecms_violation_examples[:5]:
            print(f"    n={ex['n']}, edge={ex['edge']}, "
                  f"mode(T)={ex['mode_T']}, mode(T/e)={ex['mode_Te']}")

    # Histogram
    print()
    print("delta_mu histogram (bucket: count):")
    for bucket in sorted(hist.keys()):
        count = hist[bucket]
        bar_len = min(count // max(1, total_edges // 500), 60)
        bar = "#" * bar_len
        print(f"  {bucket:+6.2f}: {count:>10,}  {bar}")

    # Save results
    hist_sorted = sorted(hist.items())
    results = {
        "description": "ECMS mean shift test for n=17 and n=18",
        "prior_max_abs_dmu": 0.5263,
        "prior_max_n": 16,
        "tested_n": [17, 18],
        "total_trees": total_trees,
        "total_edges": total_edges,
        "total_time_s": round(total_time, 2),
        "max_abs_dmu": float(max_abs_dmu),
        "max_abs_dmu_exact": str(max_abs_dmu),
        "max_abs_dmu_info": max_abs_dmu_info,
        "dmu_lt_1_holds": dmu_ge_1_count == 0,
        "dmu_ge_1_count": dmu_ge_1_count,
        "dmu_ge_1_examples": dmu_ge_1_examples,
        "max_mode_shift": max_mode_shift,
        "ecms_holds": ecms_violations == 0,
        "ecms_violations": ecms_violations,
        "ecms_violation_examples": ecms_violation_examples,
        "max_mode_shift_info": max_mode_shift_info,
        "histogram": {f"{k:+.2f}": v for k, v in hist_sorted},
        "per_n": {str(k): v for k, v in per_n_stats.items()},
    }

    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "results",
        "ecms_mean_n18.json",
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
