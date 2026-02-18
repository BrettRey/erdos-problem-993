#!/usr/bin/env python3
"""Test whether |μ(T) - μ(T/e)| < 1 for all tree edges.

For all trees through n=16 and all edges e in each tree:
  1. Compute μ(T) = I'(T;1)/I(T;1)
  2. Construct T/e (edge contraction)
  3. Compute μ(T/e)
  4. Record δμ = μ(T) - μ(T/e) and check |δμ| < 1

Also cross-checks against integer ECMS: |mode(I(T)) - mode(I(T/e))| <= 1.
"""

import json
import os
import sys
import time
from collections import defaultdict

from indpoly import independence_poly, is_log_concave
from trees import trees

MAX_N = 16
PROGRESS_INTERVAL = 5000


def poly_mean(poly):
    """Compute μ = I'(1)/I(1) for independence polynomial."""
    # I(x) = sum_k i_k x^k
    # I'(x) = sum_k k * i_k x^{k-1}
    # I(1) = sum_k i_k
    # I'(1) = sum_k k * i_k
    i_sum = sum(poly)
    if i_sum == 0:
        return 0.0
    i_prime = sum(k * poly[k] for k in range(len(poly)))
    return i_prime / i_sum


def poly_mode(poly):
    """Return the index of the maximum coefficient (mode)."""
    return max(range(len(poly)), key=lambda k: poly[k])


def contract_edge(n, adj, u, v):
    """Contract edge (u, v) in graph with n vertices and adjacency list adj.

    Merges u and v into a single vertex (we keep u, remove v).
    The merged vertex w=u has neighbors = (N(u) | N(v)) - {u, v}.
    Returns (n_new, adj_new).
    """
    # Collect neighbors of merged vertex (everything adjacent to u or v,
    # except u and v themselves)
    merged_neighbors = set()
    for w in adj[u]:
        if w != v:
            merged_neighbors.add(w)
    for w in adj[v]:
        if w != u:
            merged_neighbors.add(w)

    # Build new vertex set: remove v, relabel
    # New vertices: 0..n-2
    # Mapping: old -> new
    # v is removed; vertices > v shift down by 1
    old_to_new = {}
    for i in range(n):
        if i == v:
            continue
        if i < v:
            old_to_new[i] = i
        else:
            old_to_new[i] = i - 1

    n_new = n - 1
    adj_new = [[] for _ in range(n_new)]

    # For the merged vertex (u in old labeling, old_to_new[u] in new):
    u_new = old_to_new[u]
    for w in sorted(merged_neighbors):
        w_new = old_to_new[w]
        adj_new[u_new].append(w_new)
        adj_new[w_new].append(u_new)

    # For all other edges not involving u or v:
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

    # Sort adjacency lists
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
    print(f"ECMS Mean Shift Test: |μ(T) - μ(T/e)| for all trees through n={MAX_N}")
    print(f"Progress every {PROGRESS_INTERVAL} trees")
    print()

    # Accumulators
    total_trees = 0
    total_edges = 0

    max_abs_dmu = 0.0
    max_abs_dmu_info = None

    # Track whether |δμ| < 1 always holds
    dmu_ge_1_count = 0
    dmu_ge_1_examples = []

    # Track whether |δμ| < 1 conditional on both LC
    dmu_ge_1_both_lc_count = 0
    dmu_ge_1_both_lc_examples = []

    # Histogram of δμ in buckets of 0.1 (from -2.0 to +2.0)
    hist = defaultdict(int)

    # ECMS cross-check: |mode(T) - mode(T/e)| <= 1
    ecms_violations = 0
    ecms_violation_examples = []
    max_mode_shift = 0
    max_mode_shift_info = None

    # Per-n statistics
    per_n_stats = {}

    t_start = time.time()

    for n in range(1, MAX_N + 1):
        t_n = time.time()
        n_trees = 0
        n_edges = 0
        n_max_abs_dmu = 0.0
        n_max_mode_shift = 0
        n_dmu_ge_1 = 0
        n_ecms_viol = 0

        for tree_n, adj in trees(n, backend="auto"):
            n_trees += 1
            total_trees += 1

            if n <= 2:
                # Trees with 0 or 1 edges -- trivial
                edges = get_edges(tree_n, adj)
                if not edges:
                    continue

            edges = get_edges(tree_n, adj)
            poly_T = independence_poly(tree_n, adj)
            mu_T = poly_mean(poly_T)
            mode_T = poly_mode(poly_T)
            lc_T = is_log_concave(poly_T)

            for u, v in edges:
                n_edges += 1
                total_edges += 1

                # Contract edge (u, v)
                n_c, adj_c = contract_edge(tree_n, adj, u, v)
                poly_Te = independence_poly(n_c, adj_c)
                mu_Te = poly_mean(poly_Te)
                mode_Te = poly_mode(poly_Te)
                lc_Te = is_log_concave(poly_Te)

                dmu = mu_T - mu_Te
                abs_dmu = abs(dmu)

                # Histogram bucket
                bucket = round(dmu, 1)
                # Clamp to range for histogram
                if bucket < -2.0:
                    bucket = -2.0
                elif bucket > 2.0:
                    bucket = 2.0
                hist[bucket] += 1

                # Track max |δμ|
                if abs_dmu > n_max_abs_dmu:
                    n_max_abs_dmu = abs_dmu

                if abs_dmu > max_abs_dmu:
                    max_abs_dmu = abs_dmu
                    max_abs_dmu_info = {
                        "n": tree_n,
                        "edge": [u, v],
                        "adj": adj,
                        "dmu": dmu,
                        "mu_T": mu_T,
                        "mu_Te": mu_Te,
                        "mode_T": mode_T,
                        "mode_Te": mode_Te,
                        "lc_T": lc_T,
                        "lc_Te": lc_Te,
                    }

                # Check |δμ| >= 1
                if abs_dmu >= 1.0:
                    dmu_ge_1_count += 1
                    n_dmu_ge_1 += 1
                    if len(dmu_ge_1_examples) < 10:
                        dmu_ge_1_examples.append({
                            "n": tree_n,
                            "edge": [u, v],
                            "dmu": dmu,
                            "mu_T": mu_T,
                            "mu_Te": mu_Te,
                        })
                    # Conditional on both LC
                    if lc_T and lc_Te:
                        dmu_ge_1_both_lc_count += 1
                        if len(dmu_ge_1_both_lc_examples) < 10:
                            dmu_ge_1_both_lc_examples.append({
                                "n": tree_n,
                                "edge": [u, v],
                                "dmu": dmu,
                                "mu_T": mu_T,
                                "mu_Te": mu_Te,
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
                            "poly_T": poly_T,
                            "poly_Te": poly_Te,
                        })

            # Progress
            if total_trees % PROGRESS_INTERVAL == 0:
                elapsed = time.time() - t_start
                print(f"  [{elapsed:7.1f}s] {total_trees:>10,} trees, "
                      f"{total_edges:>12,} edges, max |δμ| = {max_abs_dmu:.6f}",
                      flush=True)

        elapsed_n = time.time() - t_n
        per_n_stats[n] = {
            "trees": n_trees,
            "edges": n_edges,
            "max_abs_dmu": n_max_abs_dmu,
            "max_mode_shift": n_max_mode_shift,
            "dmu_ge_1": n_dmu_ge_1,
            "ecms_violations": n_ecms_viol,
            "time_s": round(elapsed_n, 2),
        }

        print(f"n={n:>3}: {n_trees:>8,} trees, {n_edges:>10,} edges, "
              f"max |δμ| = {n_max_abs_dmu:.6f}, max |Δmode| = {n_max_mode_shift}, "
              f"|δμ|≥1: {n_dmu_ge_1}, ECMS viol: {n_ecms_viol} "
              f"[{elapsed_n:.2f}s]", flush=True)

    total_time = time.time() - t_start

    # Build histogram as sorted list
    hist_sorted = sorted(hist.items())

    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Trees examined:    {total_trees:,}")
    print(f"Edges examined:    {total_edges:,}")
    print(f"Total time:        {total_time:.2f}s")
    print()
    print(f"Max |δμ|:          {max_abs_dmu:.8f}")
    if max_abs_dmu_info:
        info = max_abs_dmu_info
        print(f"  Achieved at:     n={info['n']}, edge={info['edge']}")
        print(f"  μ(T) = {info['mu_T']:.8f}, μ(T/e) = {info['mu_Te']:.8f}")
        print(f"  δμ = {info['dmu']:.8f}")
        print(f"  mode(T) = {info['mode_T']}, mode(T/e) = {info['mode_Te']}")
        print(f"  LC(T) = {info['lc_T']}, LC(T/e) = {info['lc_Te']}")
    print()
    print(f"|δμ| < 1 holds:    {'YES' if dmu_ge_1_count == 0 else 'NO'}")
    print(f"  Violations:      {dmu_ge_1_count:,}")
    if dmu_ge_1_examples:
        print(f"  First examples:")
        for ex in dmu_ge_1_examples[:5]:
            print(f"    n={ex['n']}, edge={ex['edge']}, "
                  f"δμ={ex['dmu']:.6f}")
    print()
    print(f"|δμ| < 1 (both LC): {'YES' if dmu_ge_1_both_lc_count == 0 else 'NO'}")
    print(f"  Violations:      {dmu_ge_1_both_lc_count:,}")
    if dmu_ge_1_both_lc_examples:
        print(f"  First examples:")
        for ex in dmu_ge_1_both_lc_examples[:5]:
            print(f"    n={ex['n']}, edge={ex['edge']}, "
                  f"δμ={ex['dmu']:.6f}")
    print()
    print(f"Max |Δmode|:       {max_mode_shift}")
    print(f"ECMS (|Δmode|≤1):  {'HOLDS' if ecms_violations == 0 else 'FAILS'}")
    print(f"  Violations:      {ecms_violations:,}")
    if ecms_violation_examples:
        print(f"  First examples:")
        for ex in ecms_violation_examples[:5]:
            print(f"    n={ex['n']}, edge={ex['edge']}, "
                  f"mode(T)={ex['mode_T']}, mode(T/e)={ex['mode_Te']}")
    print()
    print("δμ histogram (bucket: count):")
    for bucket, count in hist_sorted:
        bar = "#" * min(count // max(1, total_edges // 1000), 60)
        print(f"  {bucket:+5.1f}: {count:>10,}  {bar}")

    # Save results
    results = {
        "max_n": MAX_N,
        "total_trees": total_trees,
        "total_edges": total_edges,
        "total_time_s": round(total_time, 2),
        "max_abs_dmu": max_abs_dmu,
        "max_abs_dmu_info": max_abs_dmu_info,
        "dmu_lt_1_holds": dmu_ge_1_count == 0,
        "dmu_ge_1_count": dmu_ge_1_count,
        "dmu_ge_1_examples": dmu_ge_1_examples,
        "dmu_lt_1_both_lc_holds": dmu_ge_1_both_lc_count == 0,
        "dmu_ge_1_both_lc_count": dmu_ge_1_both_lc_count,
        "dmu_ge_1_both_lc_examples": dmu_ge_1_both_lc_examples,
        "max_mode_shift": max_mode_shift,
        "ecms_holds": ecms_violations == 0,
        "ecms_violations": ecms_violations,
        "ecms_violation_examples": [
            {k: v for k, v in ex.items() if k != "poly_T" and k != "poly_Te"}
            for ex in ecms_violation_examples
        ],
        "max_mode_shift_info": max_mode_shift_info,
        "histogram": {f"{k:+.1f}": v for k, v in hist_sorted},
        "per_n": per_n_stats,
    }

    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "results",
        "ecms_mean_shift.json",
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
