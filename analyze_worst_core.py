#!/usr/bin/env python3
"""Analyze the structure of trees with worst-case core averages.

For each n, find the tree with maximum core average and report:
- Detailed vertex classification (L, S, C)
- Per-vertex occupation probabilities
- Core subgraph structure (is it a path? star? other?)
- Degree sequence
"""

import subprocess
import sys

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
    """Compute P(v in S) for each vertex v."""
    total_poly = independence_poly(n, adj)
    total = sum(total_poly)
    if total == 0:
        return [0.0] * n

    probs = []
    nbr = [set(adj[v]) for v in range(n)]

    for v in range(n):
        remaining = [u for u in range(n) if u != v and u not in nbr[v]]
        if not remaining:
            probs.append(1.0 / total)
            continue

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
    """Classify vertices into L (leaves), S (support), C (core)."""
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

    return leaves, support, core


def analyze_core_structure(core, adj):
    """Analyze the structure of the core subgraph.

    Returns a description string.
    """
    core_list = sorted(core)
    if len(core_list) == 0:
        return "empty"
    if len(core_list) == 1:
        return "singleton"

    # Build core subgraph
    core_set = set(core)
    core_adj = {v: [u for u in adj[v] if u in core_set] for v in core}
    core_degs = [len(core_adj[v]) for v in core_list]

    # Check if it's a path
    if all(d <= 2 for d in core_degs):
        if core_degs.count(1) == 2 and all(d in [1, 2] for d in core_degs):
            return f"path (length {len(core_list)})"

    # Check if it's a star
    if max(core_degs) == len(core_list) - 1 and core_degs.count(1) == len(core_list) - 1:
        return f"star (K_{1,{len(core_list)-1}})"

    # Count leaves and internal nodes in core subgraph
    deg_dist = {}
    for d in core_degs:
        deg_dist[d] = deg_dist.get(d, 0) + 1

    return f"general (deg dist: {deg_dist})"


def main():
    print("WORST-CASE CORE STRUCTURE ANALYSIS", flush=True)
    print("=" * 70, flush=True)
    print(flush=True)

    for n in range(7, 19):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        best_core_avg = 0.0
        best_tree = None
        best_probs = None
        best_cls = None

        for line in lines:
            tn, adj_data = parse_graph6(line)

            leaves, support, core = classify_vertices(tn, adj_data)

            # Check d_leaf <= 1
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

            if not core:
                continue

            probs = vertex_occupation_probs(tn, adj_data)
            core_sum = sum(probs[v] for v in core)
            core_avg = core_sum / len(core)

            if core_avg > best_core_avg:
                best_core_avg = core_avg
                best_tree = (tn, adj_data)
                best_probs = probs
                best_cls = (leaves, support, core)

        if best_tree is None:
            continue

        tn, adj = best_tree
        leaves, support, core = best_cls

        print(f"n={n}: core_avg = {best_core_avg:.6f}", flush=True)
        print(f"  |L|={len(leaves)}, |S|={len(support)}, |C|={len(core)}", flush=True)

        deg_seq = sorted([len(adj[v]) for v in range(tn)], reverse=True)
        print(f"  deg: {deg_seq}", flush=True)

        core_structure = analyze_core_structure(core, adj)
        print(f"  core structure: {core_structure}", flush=True)

        # Print per-vertex details
        core_list = sorted(core)
        print(f"  core vertices and P(v):", flush=True)
        for v in core_list:
            p = best_probs[v]
            d = len(adj[v])
            d_core = sum(1 for u in adj[v] if u in core)
            d_supp = sum(1 for u in adj[v] if u in support)
            print(f"    v={v}: P={p:.6f}, deg={d} ({d_core} core, {d_supp} support)",
                  flush=True)

        print(flush=True)

    print("DONE", flush=True)


if __name__ == "__main__":
    main()
