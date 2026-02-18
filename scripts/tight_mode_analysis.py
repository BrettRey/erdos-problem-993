#!/usr/bin/env python3
"""Analyze tight cases where mode(I(T-w)) = d(I(T)).

For each tree T and vertex w, the conjecture (verified for n <= 18) is:
    mode(I(T-w)) <= d(I(T))
where d(I(T)) is the first descent index of I(T).

This script finds all (T, w) pairs where equality holds, and classifies
the structural patterns.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from typing import Any

sys.path.insert(0, "/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993")
from indpoly import independence_poly, _polymul, _polyadd
from trees import trees


def mode_index(seq: list[int]) -> int:
    """Return the index of the last maximum in seq."""
    if not seq:
        return -1
    maxv = max(seq)
    idx = -1
    for i, v in enumerate(seq):
        if v == maxv:
            idx = i
    return idx


def first_mode_index(seq: list[int]) -> int:
    """Return the index of the first maximum in seq."""
    if not seq:
        return -1
    maxv = max(seq)
    for i, v in enumerate(seq):
        if v == maxv:
            return i
    return -1


def first_descent(seq: list[int]) -> int:
    """Return the first index i where seq[i] < seq[i-1], or -1 if non-decreasing."""
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            return i
    return -1


def encode_graph6_small(adj: list[list[int]]) -> str:
    """Encode a small graph as graph6."""
    n = len(adj)
    if n >= 63:
        raise ValueError("n too large for small graph6 encoder")
    aset = [set(nei) for nei in adj]
    bits: list[int] = []
    for j in range(1, n):
        sj = aset[j]
        for i in range(j):
            bits.append(1 if i in sj else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [chr(n + 63)]
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k : k + 6]:
            v = (v << 1) | b
        out.append(chr(v + 63))
    return "".join(out)


def degree_sequence(adj: list[list[int]]) -> list[int]:
    """Return sorted degree sequence (descending)."""
    return sorted([len(adj[v]) for v in range(len(adj))], reverse=True)


def classify_tree(adj: list[list[int]]) -> str:
    """Classify a tree into structural types."""
    n = len(adj)
    degs = [len(adj[v]) for v in range(n)]
    max_deg = max(degs) if degs else 0
    leaves = sum(1 for d in degs if d == 1)
    internal = sum(1 for d in degs if d >= 2)
    high_deg = sum(1 for d in degs if d >= 3)

    if n <= 2:
        return "path"

    # Path: all internal vertices have degree 2
    if max_deg <= 2:
        return "path"

    # Star: one vertex has degree n-1
    if max_deg == n - 1:
        return "star"

    # Caterpillar: all vertices are within distance 1 of a central path
    # (i.e., removing all leaves gives a path)
    spine = [v for v in range(n) if degs[v] >= 2]
    if spine:
        spine_degs = [sum(1 for u in adj[v] if degs[u] >= 2) for v in spine]
        is_caterpillar = all(d <= 2 for d in spine_degs)
    else:
        is_caterpillar = True

    # Spider: exactly one vertex with degree >= 3
    if high_deg == 1:
        return "spider"

    # Broom: caterpillar with one high-degree endpoint
    # (path with leaves attached to one end)
    if is_caterpillar and high_deg == 1:
        return "caterpillar-broom"

    if is_caterpillar:
        return "caterpillar"

    # Double star: two adjacent high-degree vertices, all others are leaves
    if high_deg == 2:
        hv = [v for v in range(n) if degs[v] >= 3]
        if len(hv) == 2 and hv[1] in adj[hv[0]]:
            other_degs = [degs[v] for v in range(n) if v not in hv]
            if all(d == 1 for d in other_degs):
                return "double-star"

    return "other"


def component_sizes(adj: list[list[int]], removed: set[int]) -> list[int]:
    """Return sorted component sizes of the forest after removing vertices."""
    n = len(adj)
    remaining = [v for v in range(n) if v not in removed]
    if not remaining:
        return []
    rem_set = set(remaining)
    seen = set()
    sizes = []
    for start in remaining:
        if start in seen:
            continue
        comp_size = 0
        stack = [start]
        seen.add(start)
        while stack:
            u = stack.pop()
            comp_size += 1
            for v in adj[u]:
                if v in rem_set and v not in seen:
                    seen.add(v)
                    stack.append(v)
        sizes.append(comp_size)
    return sorted(sizes, reverse=True)


def forest_poly_from_adj(adj: list[list[int]], removed: set[int]) -> list[int]:
    """Compute independence polynomial of forest T - removed."""
    n = len(adj)
    remaining = [v for v in range(n) if v not in removed]
    if not remaining:
        return [1]
    rem_set = set(remaining)
    seen = set()
    out = [1]
    for start in remaining:
        if start in seen:
            continue
        comp = []
        stack = [start]
        seen.add(start)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v in rem_set and v not in seen:
                    seen.add(v)
                    stack.append(v)
        # Build adjacency for component
        mapping = {old: i for i, old in enumerate(comp)}
        cadj: list[list[int]] = [[] for _ in range(len(comp))]
        for old in comp:
            ni = mapping[old]
            for v in adj[old]:
                j = mapping.get(v)
                if j is not None:
                    cadj[ni].append(j)
        poly_c = independence_poly(len(comp), cadj)
        out = _polymul(out, poly_c)
    return out


def eccentricity(adj: list[list[int]], v: int) -> int:
    """BFS eccentricity of vertex v."""
    n = len(adj)
    dist = [-1] * n
    dist[v] = 0
    queue = [v]
    head = 0
    while head < len(queue):
        u = queue[head]
        head += 1
        for w in adj[u]:
            if dist[w] == -1:
                dist[w] = dist[u] + 1
                queue.append(w)
    return max(dist)


def center_vertices(adj: list[list[int]]) -> set[int]:
    """Return the center vertices (minimum eccentricity)."""
    n = len(adj)
    eccs = [eccentricity(adj, v) for v in range(n)]
    min_ecc = min(eccs)
    return {v for v in range(n) if eccs[v] == min_ecc}


def is_center(adj: list[list[int]], v: int) -> bool:
    """Check if v is a center vertex."""
    return v in center_vertices(adj)


def main() -> None:
    ap = argparse.ArgumentParser(description="Analyze tight mode(I(T-w)) = d(I(T)) cases")
    ap.add_argument("--min-n", type=int, default=5)
    ap.add_argument("--max-n", type=int, default=16)
    ap.add_argument("--backend", default="auto", choices=["auto", "geng", "networkx"])
    ap.add_argument("--out", default="/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/tight_mode_cases.json")
    args = ap.parse_args()

    print(f"Analyzing tight mode cases for n = {args.min_n} to {args.max_n}")
    print("=" * 60)

    all_tight_cases: list[dict[str, Any]] = []
    summary_by_n: list[dict[str, Any]] = []

    # Global counters
    total_trees = 0
    total_vertex_cases = 0
    total_tight_cases = 0

    # Vertex-degree distribution in tight cases
    tight_w_degree_counts: Counter = Counter()
    # Tree type distribution
    tight_tree_type_counts: Counter = Counter()
    # Whether w is a leaf
    tight_w_is_leaf_count = 0
    tight_w_is_center_count = 0
    # mode(I(T-N[w])) stats
    mode_h_values: list[int] = []
    mode_g_values: list[int] = []
    d_f_values: list[int] = []

    for n in range(args.min_n, args.max_n + 1):
        tree_count = 0
        vertex_cases_n = 0
        tight_cases_n = 0
        non_tight_cases_n = 0
        tight_cases_list: list[dict[str, Any]] = []

        # Per-n counters
        w_deg_counts_n: Counter = Counter()
        tree_type_counts_n: Counter = Counter()
        w_is_leaf_n = 0
        w_is_center_n = 0

        for _, adj in trees(n, backend=args.backend):
            tree_count += 1
            f = independence_poly(n, adj)
            d_f = first_descent(f)
            mode_f = mode_index(f)
            first_mode_f = first_mode_index(f)
            deg_seq = degree_sequence(adj)
            tree_type = classify_tree(adj)
            g6 = encode_graph6_small(adj)

            # If the sequence is non-decreasing, d_f = -1, skip
            if d_f == -1:
                vertex_cases_n += n
                continue

            centers = center_vertices(adj)
            degs = [len(adj[v]) for v in range(n)]

            for w in range(n):
                vertex_cases_n += 1

                # Compute I(T-w)
                g = forest_poly_from_adj(adj, {w})
                mode_g = mode_index(g)

                if mode_g == d_f:
                    # TIGHT CASE
                    tight_cases_n += 1

                    # Compute I(T - N[w]) (closed neighborhood)
                    closed_neighborhood = {w} | set(adj[w])
                    h = forest_poly_from_adj(adj, closed_neighborhood)
                    mode_h = mode_index(h)

                    # Component sizes of T-w
                    comp_sizes_tw = component_sizes(adj, {w})
                    # Component sizes of T-N[w]
                    comp_sizes_tnw = component_sizes(adj, closed_neighborhood)

                    deg_w = degs[w]
                    w_is_leaf = (deg_w == 1)
                    w_in_center = (w in centers)

                    case_info = {
                        "n": n,
                        "graph6": g6,
                        "deg_seq": deg_seq,
                        "tree_type": tree_type,
                        "w": w,
                        "deg_w": deg_w,
                        "w_is_leaf": w_is_leaf,
                        "w_in_center": w_in_center,
                        "d_f": d_f,
                        "mode_f": mode_f,
                        "first_mode_f": first_mode_f,
                        "mode_g": mode_g,
                        "mode_h": mode_h,
                        "mode_g_minus_mode_h": mode_g - mode_h,
                        "comp_sizes_tw": comp_sizes_tw,
                        "comp_sizes_tnw": comp_sizes_tnw,
                        "coeffs_f": f,
                        "coeffs_g": g,
                        "coeffs_h": h,
                    }
                    tight_cases_list.append(case_info)
                    all_tight_cases.append(case_info)

                    # Update global counters
                    tight_w_degree_counts[deg_w] += 1
                    tight_tree_type_counts[tree_type] += 1
                    if w_is_leaf:
                        tight_w_is_leaf_count += 1
                    if w_in_center:
                        tight_w_is_center_count += 1
                    mode_h_values.append(mode_h)
                    mode_g_values.append(mode_g)
                    d_f_values.append(d_f)
                else:
                    non_tight_cases_n += 1

        total_trees += tree_count
        total_vertex_cases += vertex_cases_n
        total_tight_cases += tight_cases_n

        frac = tight_cases_n / vertex_cases_n if vertex_cases_n > 0 else 0

        n_summary = {
            "n": n,
            "trees": tree_count,
            "vertex_cases": vertex_cases_n,
            "tight_cases": tight_cases_n,
            "tight_fraction": round(frac, 6),
            "tight_trees": len(set(c["graph6"] for c in tight_cases_list)),
        }

        # Count how many tight cases per vertex degree for this n
        deg_dist_n: Counter = Counter()
        type_dist_n: Counter = Counter()
        leaf_count_n = 0
        center_count_n = 0
        for c in tight_cases_list:
            deg_dist_n[c["deg_w"]] += 1
            type_dist_n[c["tree_type"]] += 1
            if c["w_is_leaf"]:
                leaf_count_n += 1
            if c["w_in_center"]:
                center_count_n += 1

        n_summary["tight_w_deg_dist"] = dict(sorted(deg_dist_n.items()))
        n_summary["tight_tree_type_dist"] = dict(sorted(type_dist_n.items()))
        n_summary["tight_w_is_leaf"] = leaf_count_n
        n_summary["tight_w_in_center"] = center_count_n

        summary_by_n.append(n_summary)
        print(f"n={n:2d}: {tree_count:6d} trees, {vertex_cases_n:8d} vertex-cases, "
              f"{tight_cases_n:6d} tight ({frac:.4%}), "
              f"tight-trees={n_summary['tight_trees']}")

    print("=" * 60)
    tight_frac_total = total_tight_cases / total_vertex_cases if total_vertex_cases > 0 else 0
    print(f"Total: {total_trees} trees, {total_vertex_cases} vertex-cases, "
          f"{total_tight_cases} tight ({tight_frac_total:.4%})")

    # Analyze patterns
    print("\n--- Vertex Degree Distribution in Tight Cases ---")
    for deg, count in sorted(tight_w_degree_counts.items()):
        pct = count / total_tight_cases * 100 if total_tight_cases > 0 else 0
        print(f"  deg(w) = {deg}: {count} ({pct:.1f}%)")

    print(f"\n  w is leaf: {tight_w_is_leaf_count} ({tight_w_is_leaf_count / total_tight_cases * 100:.1f}%)" if total_tight_cases > 0 else "")
    print(f"  w is center: {tight_w_is_center_count} ({tight_w_is_center_count / total_tight_cases * 100:.1f}%)" if total_tight_cases > 0 else "")

    print("\n--- Tree Type Distribution in Tight Cases ---")
    for tt, count in sorted(tight_tree_type_counts.items(), key=lambda x: -x[1]):
        pct = count / total_tight_cases * 100 if total_tight_cases > 0 else 0
        print(f"  {tt}: {count} ({pct:.1f}%)")

    # mode(I(T-N[w])) analysis
    print("\n--- mode(I(T-N[w])) vs mode(I(T-w)) in Tight Cases ---")
    diff_counter: Counter = Counter()
    for mg, mh in zip(mode_g_values, mode_h_values):
        diff_counter[mg - mh] += 1
    for diff, count in sorted(diff_counter.items()):
        pct = count / total_tight_cases * 100 if total_tight_cases > 0 else 0
        print(f"  mode(g) - mode(h) = {diff}: {count} ({pct:.1f}%)")

    # d_f vs mode_f relationship in tight cases
    print("\n--- d(I(T)) vs mode(I(T)) in Tight Cases ---")
    df_mf_counter: Counter = Counter()
    for df, mf in zip(d_f_values, [c["mode_f"] for c in all_tight_cases]):
        df_mf_counter[df - mf] += 1
    for diff, count in sorted(df_mf_counter.items()):
        pct = count / total_tight_cases * 100 if total_tight_cases > 0 else 0
        print(f"  d(f) - mode(f) = {diff}: {count} ({pct:.1f}%)")

    # Analyze: is there a "gap" pattern?  Do tight cases cluster
    # around particular mode values relative to n?
    print("\n--- d(I(T)) / n ratio in Tight Cases ---")
    if all_tight_cases:
        ratios = [c["d_f"] / c["n"] for c in all_tight_cases]
        avg_ratio = sum(ratios) / len(ratios)
        min_ratio = min(ratios)
        max_ratio = max(ratios)
        print(f"  avg d(f)/n = {avg_ratio:.4f}, range [{min_ratio:.4f}, {max_ratio:.4f}]")

    # Which specific trees produce the most tight vertices?
    print("\n--- Trees with Most Tight Vertices (top 20) ---")
    tree_tight_count: Counter = Counter()
    tree_info: dict[str, dict] = {}
    for c in all_tight_cases:
        key = (c["n"], c["graph6"])
        tree_tight_count[key] += 1
        if key not in tree_info:
            tree_info[key] = {"deg_seq": c["deg_seq"], "tree_type": c["tree_type"], "n": c["n"]}
    for (nn, g6), count in tree_tight_count.most_common(20):
        info = tree_info[(nn, g6)]
        print(f"  n={nn} {g6} type={info['tree_type']} deg_seq={info['deg_seq']}: {count} tight vertices")

    # Save results
    results = {
        "min_n": args.min_n,
        "max_n": args.max_n,
        "total_trees": total_trees,
        "total_vertex_cases": total_vertex_cases,
        "total_tight_cases": total_tight_cases,
        "tight_fraction": round(tight_frac_total, 8),
        "tight_w_degree_dist": {str(k): v for k, v in sorted(tight_w_degree_counts.items())},
        "tight_tree_type_dist": dict(sorted(tight_tree_type_counts.items(), key=lambda x: -x[1])),
        "tight_w_is_leaf_count": tight_w_is_leaf_count,
        "tight_w_is_leaf_fraction": round(tight_w_is_leaf_count / total_tight_cases, 6) if total_tight_cases > 0 else 0,
        "tight_w_in_center_count": tight_w_is_center_count,
        "tight_w_in_center_fraction": round(tight_w_is_center_count / total_tight_cases, 6) if total_tight_cases > 0 else 0,
        "mode_g_minus_mode_h_dist": {str(k): v for k, v in sorted(diff_counter.items())},
        "by_n": summary_by_n,
        # Store individual tight cases (only for n <= 12 to keep file size reasonable)
        "tight_cases_detail": [c for c in all_tight_cases if c["n"] <= 12],
    }

    with open(args.out, "w", encoding="utf-8") as fout:
        json.dump(results, fout, indent=2)
    print(f"\nResults saved to {args.out}")


if __name__ == "__main__":
    main()
