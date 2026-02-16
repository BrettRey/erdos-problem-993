#!/usr/bin/env python3
"""Verify that core vertices have average P(v) < 1/3 in d_leaf <= 1 trees.

For trees with d_leaf(v) <= 1 (every vertex has at most one leaf-child):
- L = leaves (degree 1)
- S = support vertices (exactly one leaf-child each, so |S| = |L|)
- C = core vertices (no leaf children, not leaves)

We know from edge bounds: Σ_{L∪S} P(v) <= |L∪S|/3.
So μ < n/3 reduces to: Σ_{core} P(v) < |C|/3 (core average below 1/3).

This script verifies:
1. Is core_sum < |C|/3 for all d_leaf <= 1 trees?
2. For each heavy core vertex h (P(h) > 1/3), does h have a light neighbor
   in the core? (I.e., can we match each heavy core vertex to a distinct
   light core neighbor via edges?)
"""

import subprocess
import sys
import time

sys.path.insert(0, ".")
from indpoly import independence_poly


def parse_graph6(s):
    """Parse graph6 format (string version)."""
    s = s.strip()
    if not s:
        return 0, []
    idx = 0
    if ord(s[0]) - 63 < 63:
        n = ord(s[0]) - 63
        idx = 1
    else:
        idx = 1
        n = 0
        for _ in range(3):
            n = n * 64 + (ord(s[idx]) - 63)
            idx += 1
    adj = [[] for _ in range(n)]
    bits = []
    for c in s[idx:]:
        val = ord(c) - 63
        for b in range(5, -1, -1):
            bits.append((val >> b) & 1)
    bit_idx = 0
    for j in range(n):
        for i in range(j):
            if bit_idx < len(bits) and bits[bit_idx]:
                adj[i].append(j)
                adj[j].append(i)
            bit_idx += 1
    return n, adj


def vertex_occupation_probs(n, adj):
    """Compute P(v in S) for each vertex v.

    P(v) = (# IS containing v) / (# total IS)
    """
    # Total IS count
    total_poly = independence_poly(n, adj)
    total = sum(total_poly)
    if total == 0:
        return [0.0] * n

    probs = []
    nbr = [set(adj[v]) for v in range(n)]

    for v in range(n):
        # IS containing v: must exclude all neighbors of v
        remaining = [u for u in range(n) if u != v and u not in nbr[v]]
        if not remaining:
            # Only v can be in the IS
            probs.append(1.0 / total)
            continue

        # Build subgraph on remaining vertices
        old_to_new = {u: i for i, u in enumerate(remaining)}
        m = len(remaining)
        sub_adj = [[] for _ in range(m)]
        for u in remaining:
            for w in adj[u]:
                if w in old_to_new:
                    sub_adj[old_to_new[u]].append(old_to_new[w])

        sub_poly = independence_poly(m, sub_adj)
        count_v = sum(sub_poly)
        probs.append(count_v / total)

    return probs


def classify_vertices(n, adj):
    """Classify vertices into L (leaves), S (support), C (core).

    Returns (L, S, C) as sets of vertex indices.
    """
    leaves = set()
    for v in range(n):
        if len(adj[v]) == 1:
            leaves.add(v)

    support = set()
    core = set()
    for v in range(n):
        if v in leaves:
            continue
        leaf_children = sum(1 for u in adj[v] if u in leaves)
        if leaf_children == 1:
            support.add(v)
        elif leaf_children == 0:
            core.add(v)
        # If leaf_children >= 2, this tree doesn't satisfy d_leaf <= 1

    return leaves, support, core


def check_heavy_light_matching(core_vertices, probs, adj):
    """Check if each heavy core vertex has a light core neighbor.

    Returns (ok, problem_cases) where:
    - ok: True if every heavy vertex has a light core neighbor
    - problem_cases: list of (heavy_v, P(heavy_v)) for violations
    """
    heavy = [v for v in core_vertices if probs[v] > 1/3 + 1e-12]
    if not heavy:
        return True, []

    core_set = set(core_vertices)
    problems = []

    for h in heavy:
        # Check neighbors of h in the core
        core_neighbors = [u for u in adj[h] if u in core_set]
        # Is any neighbor light (P <= 1/3)?
        has_light_neighbor = any(probs[u] <= 1/3 + 1e-12 for u in core_neighbors)
        if not has_light_neighbor:
            problems.append((h, probs[h]))

    return len(problems) == 0, problems


def main():
    print("CORE VERTEX AVERAGE VERIFICATION FOR d_leaf <= 1 TREES", flush=True)
    print("=" * 70, flush=True)
    print(flush=True)
    print("Checking:", flush=True)
    print("  1. Is Σ_{v∈C} P(v) < |C|/3 for all d_leaf <= 1 trees?", flush=True)
    print("  2. Does each heavy core vertex (P > 1/3) have a light core", flush=True)
    print("     neighbor (P <= 1/3)?", flush=True)
    print(flush=True)

    for n in range(3, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_total = 0
        n_with_core = 0
        n_core_violations = 0
        n_matching_violations = 0
        max_core_avg = 0.0
        max_core_info = None

        for line in lines:
            tn, adj_data = parse_graph6(line)

            # Check d_leaf <= 1 condition
            leaves, support, core = classify_vertices(tn, adj_data)

            # Skip if any vertex has d_leaf >= 2
            leaves_set = leaves
            ok = True
            for v in range(tn):
                if v in leaves_set:
                    continue
                d_leaf_v = sum(1 for u in adj_data[v] if u in leaves_set)
                if d_leaf_v >= 2:
                    ok = False
                    break
            if not ok:
                continue

            n_total += 1

            if not core:
                # No core vertices (all leaves or support)
                continue

            n_with_core += 1

            # Compute occupation probabilities
            probs = vertex_occupation_probs(tn, adj_data)

            # Check core average
            core_sum = sum(probs[v] for v in core)
            core_size = len(core)
            core_avg = core_sum / core_size if core_size > 0 else 0.0

            if core_avg > max_core_avg:
                max_core_avg = core_avg
                # Store info about this tree
                deg_seq = sorted([len(adj_data[v]) for v in range(tn)], reverse=True)
                max_core_info = (tn, core_size, deg_seq[:8])

            # Check if core_avg < 1/3
            if core_avg >= 1/3 - 1e-12:
                n_core_violations += 1
                if n_core_violations <= 3:
                    print(f"  CORE AVG VIOLATION n={tn}: core_avg={core_avg:.6f}, "
                          f"|C|={core_size}, core_sum={core_sum:.6f}", flush=True)

            # Check heavy-light matching
            matching_ok, problems = check_heavy_light_matching(core, probs, adj_data)
            if not matching_ok:
                n_matching_violations += 1
                if n_matching_violations <= 3:
                    print(f"  MATCHING VIOLATION n={tn}: {len(problems)} heavy "
                          f"vertices without light core neighbors", flush=True)
                    for h, p_h in problems[:2]:
                        print(f"    v={h}, P(v)={p_h:.6f}, "
                              f"deg={len(adj_data[h])}", flush=True)

        elapsed = time.time() - t0
        info_str = ""
        if max_core_info:
            tn, cs, ds = max_core_info
            info_str = f" max_tree: n={tn}, |C|={cs}, deg={ds}"

        print(f"n={n:2d}: {n_total:6d} trees ({n_with_core:6d} with core), "
              f"max_core_avg={max_core_avg:.6f}, "
              f"core_viol={n_core_violations}, match_viol={n_matching_violations} "
              f"({elapsed:.1f}s){info_str}", flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
