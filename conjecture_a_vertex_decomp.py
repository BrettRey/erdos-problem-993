#!/usr/bin/env python3
"""Analyze Conjecture A via hard-core occupation probabilities.

For all d_leaf <= 1 trees through n=18, compute per-vertex occupation
probabilities P(v) at fugacity lambda=1 using the cavity method on trees.

Classifies vertices as leaf / support / interior and decomposes mu = sum P(v)
into contributions from leaf-support pairs and interior vertices.

Cross-checks mu against I'(1)/I(1) from the independence polynomial.
"""

import json
import os
import subprocess
import sys
import time

from graph6 import parse_graph6
from indpoly import independence_poly


# ---------------------------------------------------------------------------
# Cavity method for hard-core model on trees at lambda = 1
# ---------------------------------------------------------------------------

def compute_occupation_probs(n, adj):
    """Compute P(v) for each vertex of a tree at fugacity lambda=1.

    Uses the cavity method:
      - Root the tree at vertex 0.
      - Compute cavity ratios R_{u->parent(u)} bottom-up.
        R_{leaf} = 1 (lambda = 1).
        R_{u->v} = prod_{c in children(u)} 1/(1 + R_{c->u})
                 = lambda * prod ... but lambda=1 so just the product.
      - Then P(v) = lambda / (lambda + prod_{neighbors u} (1 + R_{u->v}))
                  = 1 / (1 + prod_{neighbors u} (1 + R_{u->v}))
        where R_{u->v} is the cavity ratio from u toward v (in the subtree
        on u's side when v is removed).

    For the root, all neighbors are children, so:
      P(root) = 1 / (1 + prod_{c} (1 + R_{c->root}))

    For non-root vertex v with parent p:
      We need R_{p->v} (cavity ratio from parent toward v).
      R_{p->v} = prod_{c in neighbors(p) \\ {v}} 1/(1 + R_{c->p})
      Then P(v) = 1 / (1 + (1 + R_{p->v}) * prod_{c in children(v)} (1 + R_{c->v}))

    Returns list of floats P[v] for v = 0..n-1.
    """
    if n == 1:
        return [1.0]  # single vertex: always in the IS
    if n == 2:
        return [0.5, 0.5]

    # Root at vertex 0
    parent = [-1] * n
    children = [[] for _ in range(n)]
    visited = [False] * n
    order = []  # BFS order

    visited[0] = True
    queue = [0]
    head = 0
    while head < len(queue):
        v = queue[head]
        head += 1
        order.append(v)
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                parent[u] = v
                children[v].append(u)
                queue.append(u)

    # Bottom-up: compute R[v] = cavity ratio R_{v -> parent(v)}
    # For a leaf: R = lambda = 1
    # For internal: R[v] = prod_{c in children(v)} 1/(1 + R[c])
    # (This is R_{v->parent(v)} = lambda * prod_{c} 1/(1+R_{c->v}))
    # At lambda=1: R[v] = prod_{c in children(v)} 1/(1 + R[c])

    R_up = [0.0] * n  # R_up[v] = R_{v -> parent(v)}
    for v in reversed(order):
        if not children[v]:
            R_up[v] = 1.0  # leaf in rooted tree: lambda = 1
        else:
            prod = 1.0
            for c in children[v]:
                prod *= 1.0 / (1.0 + R_up[c])
            R_up[v] = prod  # lambda * prod at lambda=1

    # Now compute P(v) for each vertex.
    # We need the cavity ratios from ALL neighbors, not just children.
    # For the root, all neighbors are children.
    # For non-root v, we also need R_{parent(v) -> v}.
    #
    # R_{parent -> v} = lambda * prod_{u in neighbors(parent) \ {v}} 1/(1 + R_{u->parent})
    #
    # For the root's children: R_{root -> c} = lambda * prod_{c' != c} 1/(1 + R_{c'->root})
    # = (R_up[root] / (1/(1+R_up[c])))  ... but R_up[root] uses all children
    # = R_up[root] * (1 + R_up[c])
    #
    # For general v with parent p:
    # R_{p -> v} is computed from R_down[p] (cavity from p's parent toward p)
    # and R_up[c'] for all siblings c' of v.
    # R_{p -> v} = lambda * (1/(1 + R_down[p])) * prod_{c' sibling of v} 1/(1 + R_up[c'])
    # = (1/(1 + R_down[p])) * prod_{c' sibling} 1/(1 + R_up[c'])
    # = (1/(1 + R_down[p])) * (R_up[p] * (1 + R_up[v]))
    #   because R_up[p] = prod_{all children c of p} 1/(1+R_up[c])
    #   so prod_{c' != v} 1/(1+R_up[c']) = R_up[p] / (1/(1+R_up[v])) = R_up[p] * (1 + R_up[v])

    # R_down[v] = cavity ratio from parent(v) toward v
    R_down = [0.0] * n  # not defined for root

    # Compute R_down top-down
    for v in order:
        for c in children[v]:
            if v == 0:
                # Root: R_{root -> c} = prod_{c' != c in children(root)} 1/(1+R_up[c'])
                # = R_up[root] * (1 + R_up[c])
                R_down[c] = R_up[v] * (1.0 + R_up[c])
            else:
                # v is non-root with parent. R_down[v] is known.
                # R_{v -> c} = (1/(1 + R_down[v])) * prod_{c' != c in children(v)} 1/(1+R_up[c'])
                # = (1/(1+R_down[v])) * R_up[v] * (1 + R_up[c])
                R_down[c] = (1.0 / (1.0 + R_down[v])) * R_up[v] * (1.0 + R_up[c])

    # Now P(v) = 1 / (1 + prod_{neighbors u of v} (1 + R_{u->v}))
    # For root: all neighbors are children, R_{c->root} = R_up[c]
    #   P(root) = 1 / (1 + prod_c (1 + R_up[c]))
    # For non-root v: neighbors = children + parent
    #   prod = (1 + R_down[v]) * prod_{c in children(v)} (1 + R_up[c])
    #   But prod_{c} (1 + R_up[c]) = 1 / R_up[v]  (since R_up[v] = prod_c 1/(1+R_up[c]))
    #   ... wait, only if v has children. If v is a leaf in the rooted tree, R_up[v] = 1.
    #   Actually for leaves, children(v) is empty, so prod is 1.
    #   For non-leaves: prod_{c in children(v)} (1+R_up[c]) = 1/R_up[v]
    #   P(v) = 1 / (1 + (1+R_down[v]) / R_up[v])  [non-root, non-leaf-in-rooted-tree]
    #   P(v) = 1 / (1 + (1+R_down[v]))             [non-root leaf-in-rooted-tree]
    #   P(root) = 1 / (1 + 1/R_up[root])           [root with children]

    P = [0.0] * n

    for v in range(n):
        if v == 0:
            # Root
            if not children[v]:
                P[v] = 1.0  # isolated (shouldn't happen for n>=2 tree)
            else:
                # prod_{c} (1 + R_up[c]) = 1/R_up[0]
                P[v] = 1.0 / (1.0 + 1.0 / R_up[v])
        else:
            if not children[v]:
                # Leaf in rooted tree
                P[v] = 1.0 / (1.0 + (1.0 + R_down[v]))
            else:
                # Non-root, has children
                # prod over all neighbors = (1 + R_down[v]) * prod_c (1 + R_up[c])
                # = (1 + R_down[v]) / R_up[v]
                P[v] = 1.0 / (1.0 + (1.0 + R_down[v]) / R_up[v])

    return P


# ---------------------------------------------------------------------------
# d_leaf classification
# ---------------------------------------------------------------------------

def classify_vertices(n, adj):
    """Classify vertices as 'leaf', 'support', or 'interior'.

    - leaf: degree 1
    - support: not a leaf, has exactly 1 leaf-neighbor
    - interior: not a leaf, has 0 leaf-neighbors

    Also returns d_leaf[v] = number of leaf-neighbors for each v.
    """
    degree = [len(adj[v]) for v in range(n)]
    is_leaf = [degree[v] == 1 for v in range(n)]

    d_leaf = [0] * n
    for v in range(n):
        if is_leaf[v]:
            continue
        for u in adj[v]:
            if is_leaf[u]:
                d_leaf[v] += 1

    label = [''] * n
    for v in range(n):
        if is_leaf[v]:
            label[v] = 'leaf'
        elif d_leaf[v] == 1:
            label[v] = 'support'
        elif d_leaf[v] == 0:
            label[v] = 'interior'
        else:
            label[v] = 'high_dleaf'  # d_leaf >= 2 (should not appear in d_leaf<=1 trees)

    return label, d_leaf


def is_dleaf_le1(n, adj):
    """Check whether every vertex has at most 1 leaf-child (d_leaf <= 1)."""
    degree = [len(adj[v]) for v in range(n)]
    is_leaf = [degree[v] == 1 for v in range(n)]

    for v in range(n):
        if is_leaf[v]:
            continue
        leaf_count = 0
        for u in adj[v]:
            if is_leaf[u]:
                leaf_count += 1
                if leaf_count >= 2:
                    return False
    return True


# ---------------------------------------------------------------------------
# Cross-check mu via independence polynomial
# ---------------------------------------------------------------------------

def mu_from_poly(n, adj):
    """Compute mu = I'(1)/I(1) from the independence polynomial."""
    poly = independence_poly(n, adj)
    # I(1) = sum of coefficients
    I_1 = sum(poly)
    # I'(1) = sum of k * poly[k]
    Iprime_1 = sum(k * poly[k] for k in range(len(poly)))
    return Iprime_1 / I_1


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def trees_geng(n, res=None, mod=None):
    """Enumerate trees on n vertices using geng."""
    edges = n - 1
    cmd = ["/opt/homebrew/bin/geng", "-cq", str(n), f"{edges}:{edges}"]
    if res is not None and mod is not None:
        cmd.append(f"{res}/{mod}")

    with subprocess.Popen(cmd, stdout=subprocess.PIPE) as proc:
        for line in proc.stdout:
            yield parse_graph6(line)


def main():
    MAX_N = 18
    MIN_N = 3  # d_leaf<=1 trees don't exist below n=3 in interesting ways

    print("Conjecture A: Hard-core occupation probability decomposition")
    print(f"Analyzing d_leaf <= 1 trees for n = {MIN_N}..{MAX_N}")
    print(f"Fugacity lambda = 1")
    print()

    # Global tracking
    global_max_mu_ratio = 0.0
    global_max_mu_ratio_tree = None
    global_max_mu_ratio_n = 0

    global_max_P_interior = 0.0
    global_max_P_interior_tree = None
    global_max_P_interior_n = 0

    global_max_P_support = 0.0
    global_max_P_support_tree = None
    global_max_P_support_n = 0

    global_max_avg_P_interior = 0.0
    global_max_avg_P_interior_tree = None
    global_max_avg_P_interior_n = 0

    global_max_crosscheck_err = 0.0

    decomp_always_holds = True  # sum_pairs + sum_interior < n/3
    decomp_failure_count = 0

    total_dleaf1_trees = 0
    per_n_results = {}

    t_start = time.time()

    for n in range(MIN_N, MAX_N + 1):
        t0 = time.time()
        count_all = 0
        count_dleaf1 = 0

        n_max_mu_ratio = 0.0
        n_max_mu_ratio_adj = None
        n_max_P_interior = 0.0
        n_max_P_interior_adj = None
        n_max_P_support = 0.0
        n_max_avg_P_interior = 0.0
        n_max_crosscheck_err = 0.0

        n_decomp_failures = 0

        # Accumulators for per-n stats
        n_sum_mu_ratio = 0.0
        n_count_interior_gt_third = 0  # trees with max P(interior) > 1/3
        n_count_support_gt_third = 0

        for tree_n, adj in trees_geng(n):
            count_all += 1

            if not is_dleaf_le1(tree_n, adj):
                continue

            count_dleaf1 += 1
            total_dleaf1_trees += 1

            # Compute occupation probabilities
            P = compute_occupation_probs(tree_n, adj)

            # Cross-check: mu from P vs mu from polynomial
            mu_P = sum(P)
            mu_poly = mu_from_poly(tree_n, adj)
            err = abs(mu_P - mu_poly)
            if err > n_max_crosscheck_err:
                n_max_crosscheck_err = err
            if err > global_max_crosscheck_err:
                global_max_crosscheck_err = err
            if err > 1e-8:
                print(f"  WARNING: crosscheck error {err:.2e} at n={n}, tree #{count_dleaf1}",
                      flush=True)

            # Classify vertices
            label, d_leaf = classify_vertices(tree_n, adj)

            # Compute decomposition
            mu = mu_P
            n_third = tree_n / 3.0
            mu_ratio = mu / n_third

            n_sum_mu_ratio += mu_ratio

            # Sum over leaf-support pairs
            # Each support vertex has exactly 1 leaf-child. Pair them.
            sum_pairs = 0.0
            sum_interior = 0.0
            max_P_int = 0.0
            max_P_sup = 0.0
            count_interior = 0

            for v in range(tree_n):
                if label[v] == 'support':
                    # Find the leaf-child
                    degree_v = len(adj[v])
                    is_leaf_v = [len(adj[u]) == 1 for u in adj[v]]
                    for idx_u, u in enumerate(adj[v]):
                        if is_leaf_v[idx_u]:
                            sum_pairs += P[v] + P[u]
                            break
                    if P[v] > max_P_sup:
                        max_P_sup = P[v]
                elif label[v] == 'interior':
                    sum_interior += P[v]
                    count_interior += 1
                    if P[v] > max_P_int:
                        max_P_int = P[v]

            # Check decomposition: sum_pairs + sum_interior < n/3
            # Note: mu = sum_pairs + sum_interior (leaves covered by pairs,
            # plus standalone interior vertices). Actually, mu = sum over ALL v of P(v).
            # We need to verify: sum_pairs + sum_interior covers all of mu?
            # Leaves not paired with support are standalone. But in d_leaf<=1 trees,
            # every leaf has a unique support vertex... wait, not necessarily.
            # A leaf's parent might have d_leaf >= 2 -- but we filtered for d_leaf <= 1.
            # So every leaf's parent has d_leaf = 1 (it's a support vertex), OR
            # the leaf's parent has d_leaf = 0 -- impossible since the leaf IS a leaf-child.
            # So: in d_leaf <= 1 trees, every leaf has exactly one parent, and that parent
            # has d_leaf = 1 (it's a support vertex). Thus every leaf is paired.
            # mu = sum_pairs + sum_interior.

            decomp_total = sum_pairs + sum_interior
            decomp_check = decomp_total  # should equal mu

            if decomp_total >= n_third:
                decomp_always_holds = False
                decomp_failure_count += 1
                n_decomp_failures += 1

            avg_P_int = sum_interior / count_interior if count_interior > 0 else 0.0

            # Update per-n maxima
            if mu_ratio > n_max_mu_ratio:
                n_max_mu_ratio = mu_ratio
                n_max_mu_ratio_adj = [list(nbrs) for nbrs in adj]
            if max_P_int > n_max_P_interior:
                n_max_P_interior = max_P_int
                n_max_P_interior_adj = [list(nbrs) for nbrs in adj]
            if max_P_sup > n_max_P_support:
                n_max_P_support = max_P_sup
            if avg_P_int > n_max_avg_P_interior:
                n_max_avg_P_interior = avg_P_int

            if max_P_int > 1.0 / 3.0:
                n_count_interior_gt_third += 1
            if max_P_sup > 1.0 / 3.0:
                n_count_support_gt_third += 1

            # Update global maxima
            if mu_ratio > global_max_mu_ratio:
                global_max_mu_ratio = mu_ratio
                global_max_mu_ratio_tree = [list(nbrs) for nbrs in adj]
                global_max_mu_ratio_n = tree_n
            if max_P_int > global_max_P_interior:
                global_max_P_interior = max_P_int
                global_max_P_interior_tree = [list(nbrs) for nbrs in adj]
                global_max_P_interior_n = tree_n
            if max_P_sup > global_max_P_support:
                global_max_P_support = max_P_sup
                global_max_P_support_tree = [list(nbrs) for nbrs in adj]
                global_max_P_support_n = tree_n
            if avg_P_int > global_max_avg_P_interior:
                global_max_avg_P_interior = avg_P_int
                global_max_avg_P_interior_tree = [list(nbrs) for nbrs in adj]
                global_max_avg_P_interior_n = tree_n

            if total_dleaf1_trees % 10000 == 0:
                print(f"  [{total_dleaf1_trees:>8,} d_leaf<=1 trees processed, "
                      f"n={n}, max mu/(n/3)={global_max_mu_ratio:.6f}]", flush=True)

        elapsed = time.time() - t0
        avg_mu_ratio = n_sum_mu_ratio / count_dleaf1 if count_dleaf1 > 0 else 0.0

        per_n_results[n] = {
            "total_trees": count_all,
            "dleaf1_trees": count_dleaf1,
            "max_mu_ratio": n_max_mu_ratio,
            "avg_mu_ratio": avg_mu_ratio,
            "max_P_interior": n_max_P_interior,
            "max_P_support": n_max_P_support,
            "max_avg_P_interior": n_max_avg_P_interior,
            "max_crosscheck_err": n_max_crosscheck_err,
            "decomp_failures": n_decomp_failures,
            "trees_with_interior_gt_third": n_count_interior_gt_third,
            "trees_with_support_gt_third": n_count_support_gt_third,
            "time_s": round(elapsed, 2),
        }

        print(f"n={n:>3}: {count_all:>9,} total, {count_dleaf1:>8,} d_leaf<=1 | "
              f"max mu/(n/3)={n_max_mu_ratio:.6f}, "
              f"max P(int)={n_max_P_interior:.6f}, "
              f"max P(sup)={n_max_P_support:.6f}, "
              f"int>1/3: {n_count_interior_gt_third}, "
              f"decomp_fail: {n_decomp_failures} | "
              f"{elapsed:.1f}s", flush=True)

    total_elapsed = time.time() - t_start

    # Print detailed summary
    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total d_leaf<=1 trees analyzed: {total_dleaf1_trees:,}")
    print(f"Total time: {total_elapsed:.1f}s")
    print()
    print(f"Max mu/(n/3): {global_max_mu_ratio:.10f} (n={global_max_mu_ratio_n})")
    print(f"  Tree adj: {global_max_mu_ratio_tree}")
    print()
    print(f"Max P(interior vertex): {global_max_P_interior:.10f} (n={global_max_P_interior_n})")
    print(f"  Tree adj: {global_max_P_interior_tree}")
    print(f"  Is > 1/3? {'YES' if global_max_P_interior > 1.0/3.0 else 'NO'}")
    print()
    print(f"Max P(support vertex): {global_max_P_support:.10f} (n={global_max_P_support_n})")
    print(f"  Tree adj: {global_max_P_support_tree}")
    print(f"  Is > 1/3? {'YES' if global_max_P_support > 1.0/3.0 else 'NO'}")
    print()
    print(f"Max avg P(interior): {global_max_avg_P_interior:.10f} (n={global_max_avg_P_interior_n})")
    print(f"  Tree adj: {global_max_avg_P_interior_tree}")
    print()
    print(f"Decomposition sum_pairs + sum_interior < n/3:")
    print(f"  Always holds? {'YES' if decomp_always_holds else 'NO'}")
    print(f"  Failures: {decomp_failure_count}")
    print()
    print(f"Max cross-check error (mu_P vs mu_poly): {global_max_crosscheck_err:.2e}")
    print()

    # Print the extremal tree details
    if global_max_mu_ratio_tree is not None:
        n_ext = global_max_mu_ratio_n
        adj_ext = global_max_mu_ratio_tree
        P_ext = compute_occupation_probs(n_ext, adj_ext)
        label_ext, dleaf_ext = classify_vertices(n_ext, adj_ext)
        print("Extremal tree (highest mu/(n/3)):")
        print(f"  n = {n_ext}, mu = {sum(P_ext):.10f}, n/3 = {n_ext/3:.10f}")
        print(f"  Vertex details:")
        for v in range(n_ext):
            deg = len(adj_ext[v])
            print(f"    v={v}: P={P_ext[v]:.8f}, deg={deg}, "
                  f"label={label_ext[v]}, d_leaf={dleaf_ext[v]}")

    # Save results
    results = {
        "max_n": MAX_N,
        "total_dleaf1_trees": total_dleaf1_trees,
        "total_time_s": round(total_elapsed, 2),
        "global": {
            "max_mu_ratio": global_max_mu_ratio,
            "max_mu_ratio_n": global_max_mu_ratio_n,
            "max_mu_ratio_tree": global_max_mu_ratio_tree,
            "max_P_interior": global_max_P_interior,
            "max_P_interior_n": global_max_P_interior_n,
            "max_P_interior_tree": global_max_P_interior_tree,
            "P_interior_exceeds_third": global_max_P_interior > 1.0 / 3.0,
            "max_P_support": global_max_P_support,
            "max_P_support_n": global_max_P_support_n,
            "max_P_support_tree": global_max_P_support_tree,
            "P_support_exceeds_third": global_max_P_support > 1.0 / 3.0,
            "max_avg_P_interior": global_max_avg_P_interior,
            "max_avg_P_interior_n": global_max_avg_P_interior_n,
            "max_avg_P_interior_tree": global_max_avg_P_interior_tree,
            "decomp_always_holds": decomp_always_holds,
            "decomp_failure_count": decomp_failure_count,
            "max_crosscheck_err": global_max_crosscheck_err,
        },
        "per_n": per_n_results,
    }

    os.makedirs("results", exist_ok=True)
    outpath = "results/conjecture_a_vertex_decomp.json"
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {outpath}")


if __name__ == "__main__":
    main()
