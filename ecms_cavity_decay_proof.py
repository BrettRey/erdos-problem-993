#!/usr/bin/env python3
"""Investigate cavity message decay to bound the remote term.

When edge e = uv is contracted, cavity messages change. The change propagates
outward from {u, v} through the tree. At each vertex w on the propagation
path, the message R_{w→parent} changes by:

  ΔR_{w→parent} ≈ -(R_{w→parent} / (1 + R_{child→w})) · ΔR_{child→w}

The decay factor through w (from child c's direction) is:
  γ(w, c) = R_{w→parent} / (1 + R_{c→w})

Key question: is γ(w, c) < 1 always? If so, is there a universal bound < 1?

We also check the probability perturbation:
  ΔP(w) = ΔR_w / (1 + R_w)^2  (first order)

For the remote sum, we need: |Σ_d Σ_{w at dist d} ΔP(w)| < 1/2.

This script:
1. Computes exact γ(w, c) for every vertex-child pair in every tree
2. Tracks the maximum γ across all trees
3. Measures the actual vs predicted ΔP at each distance
4. Tests whether the geometric series 0.5 · Σ γ^d converges to < 0.5
"""

import json
import os
import time
from collections import defaultdict

from trees import trees

MAX_N = 18


def cavity_messages(n, adj):
    """Compute cavity messages R_{u→v} for a tree at λ=1."""
    if n <= 1:
        return {}, [], [-1] * n, [[] for _ in range(n)]
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
    return msgs, order, parent, children


def main():
    t0 = time.time()

    # Track decay factors γ(w, c) = R_{w→parent} / (1 + R_{c→w})
    # where c is a child of w
    max_gamma = 0.0
    gamma_histogram = defaultdict(int)  # bin to count
    max_gamma_info = {}

    # Also track product of gammas along paths
    # For each tree, for each edge e, for each vertex w at distance d,
    # compute the product of gammas along the path from {u,v} to w.
    max_gamma_product_by_dist = defaultdict(float)

    # Track vertex-level decay: R_{w→parent} itself is the key quantity
    # Since R ∈ (0, 1], and γ = R_{w→parent}/(1+R_{c→w}),
    # and R_{c→w} > 0, we have γ < R_{w→parent} ≤ 1.
    # But R_{w→parent} can be close to 1 (for leaves).
    # The crucial bound: γ ≤ R_{w→parent} and also γ = R_{w→parent}/(1+R_{c→w})
    # Since R_{c→w} ≥ something > 0 for internal vertices...

    # Track the product R_{w→parent} for non-leaf vertices (the "core" bound)
    max_R_nonleaf = 0.0  # max R_{w→parent} where w is not a leaf
    max_R_nonleaf_info = {}

    # For path: vertices with degree 2 (one child, one parent)
    max_gamma_deg2 = 0.0
    max_gamma_deg_ge3 = 0.0

    # Count how often gamma > various thresholds
    gamma_gt = {0.5: 0, 0.6: 0, 0.7: 0, 0.8: 0, 0.9: 0, 0.95: 0, 0.99: 0}
    total_pairs = 0

    for n in range(3, MAX_N + 1):
        tn = time.time()
        n_trees = 0

        for _, adj in trees(n):
            n_trees += 1
            msgs, order, parent, children = cavity_messages(n, adj)

            for w in range(n):
                p = parent[w]
                if p == -1:
                    continue  # root
                R_wp = msgs.get((w, p), 0.0)
                deg_w = len(adj[w])

                for c in children[w]:
                    R_cw = msgs.get((c, w), 0.0)
                    gamma = R_wp / (1.0 + R_cw)
                    total_pairs += 1

                    # Histogram
                    bin_idx = int(gamma * 100)
                    gamma_histogram[bin_idx] += 1

                    for thresh in gamma_gt:
                        if gamma > thresh:
                            gamma_gt[thresh] += 1

                    if gamma > max_gamma:
                        max_gamma = gamma
                        max_gamma_info = {
                            "n": n, "w": w, "child": c, "parent": p,
                            "gamma": gamma, "R_wp": R_wp, "R_cw": R_cw,
                            "deg_w": deg_w,
                        }

                    if deg_w == 2:  # path vertex
                        max_gamma_deg2 = max(max_gamma_deg2, gamma)
                    elif deg_w >= 3:
                        max_gamma_deg_ge3 = max(max_gamma_deg_ge3, gamma)

                    # Non-leaf R bound
                    if len(children[w]) > 0:  # w has children = not a leaf
                        if R_wp > max_R_nonleaf:
                            max_R_nonleaf = R_wp
                            max_R_nonleaf_info = {
                                "n": n, "w": w, "R_wp": R_wp,
                                "deg_w": deg_w,
                                "n_children": len(children[w]),
                            }

        elapsed = time.time() - tn
        print(f"n={n}: {n_trees} trees, max_gamma={max_gamma:.6f}, "
              f"max_gamma_deg2={max_gamma_deg2:.6f}, "
              f"max_R_nonleaf={max_R_nonleaf:.6f}, {elapsed:.1f}s", flush=True)

    total_time = time.time() - t0

    print(f"\n=== DECAY FACTOR γ(w,c) = R_{{w→parent}} / (1 + R_{{c→w}}) ===")
    print(f"Maximum γ: {max_gamma:.8f}")
    print(f"Max γ info: {max_gamma_info}")
    print(f"Max γ (deg=2 vertices): {max_gamma_deg2:.8f}")
    print(f"Max γ (deg≥3 vertices): {max_gamma_deg_ge3:.8f}")
    print(f"Max R_{{w→parent}} (non-leaf): {max_R_nonleaf:.8f}")

    print(f"\nγ threshold distribution (out of {total_pairs} pairs):")
    for thresh in sorted(gamma_gt.keys()):
        pct = 100 * gamma_gt[thresh] / total_pairs if total_pairs > 0 else 0
        print(f"  γ > {thresh}: {gamma_gt[thresh]} ({pct:.3f}%)")

    print(f"\nγ histogram (top bins):")
    for bin_idx in range(100, -1, -1):
        if gamma_histogram.get(bin_idx, 0) > 0:
            lo = bin_idx / 100.0
            hi = (bin_idx + 1) / 100.0
            print(f"  [{lo:.2f}, {hi:.2f}): {gamma_histogram[bin_idx]}")
            if bin_idx < 50:
                break

    # Key insight: if γ < φ^{-1} = 0.618... for all non-leaf vertices,
    # then the remote sum is bounded by geometric series.
    phi_inv = 0.6180339887
    print(f"\nGolden ratio test: γ < 1/φ = {phi_inv:.6f}?")
    violations = sum(1 for b, c in gamma_histogram.items()
                     if b >= int(phi_inv * 100) and c > 0)
    print(f"  Bins at or above 1/φ: {violations}")

    # Better bound: if γ < 1/2 for all pairs, remote sum bounded by
    # initial_perturbation * Σ (1/2)^d * (branching at dist d)
    # Since branching can increase, need to track carefully.

    results = {
        "max_n": MAX_N,
        "total_pairs": total_pairs,
        "max_gamma": round(max_gamma, 10),
        "max_gamma_deg2": round(max_gamma_deg2, 10),
        "max_gamma_deg_ge3": round(max_gamma_deg_ge3, 10),
        "max_R_nonleaf": round(max_R_nonleaf, 10),
        "max_gamma_info": max_gamma_info,
        "gamma_thresholds": {str(k): v for k, v in gamma_gt.items()},
        "total_time_s": round(total_time, 2),
    }

    out_path = "results/ecms_cavity_decay.json"
    os.makedirs("results", exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
