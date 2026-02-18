#!/usr/bin/env python3
"""Check mode alignment: mode(I(T-w)) vs mode(I(T)) for all vertices w."""

from __future__ import annotations

import argparse
import json
from typing import Dict, List

from indpoly import independence_poly, _polymul
from trees import trees


def poly_add(a: List[int], b: List[int]) -> List[int]:
    out = [0] * max(len(a), len(b))
    for i, v in enumerate(a):
        out[i] += v
    for i, v in enumerate(b):
        out[i] += v
    return out


def mode_index(seq: List[int]) -> int:
    maxv = max(seq) if seq else 0
    idx = -1
    for i, v in enumerate(seq):
        if v == maxv:
            idx = i
    return idx


def first_descent(seq: List[int]) -> int:
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            return i
    return -1


def encode_graph6_small(adj: List[List[int]]) -> str:
    n = len(adj)
    if n >= 63:
        raise ValueError("n too large for small graph6 encoder")
    aset = [set(nei) for nei in adj]
    bits: List[int] = []
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


def compute_messages(adj: List[List[int]]) -> tuple[List[Dict[int, List[int]]], List[Dict[int, List[int]]]]:
    """Compute dp0, dp1 messages for all directed edges."""
    n = len(adj)
    parent = [-1] * n
    order: List[int] = []
    parent[0] = 0
    stack = [0]
    while stack:
        u = stack.pop()
        order.append(u)
        for v in adj[u]:
            if parent[v] == -1:
                parent[v] = u
                stack.append(v)

    children = [[] for _ in range(n)]
    for v in range(1, n):
        children[parent[v]].append(v)

    msg0: List[Dict[int, List[int]]] = [dict() for _ in range(n)]
    msg1: List[Dict[int, List[int]]] = [dict() for _ in range(n)]

    # Downward messages: child -> parent
    for u in reversed(order):
        prodA = [1]
        prodB = [1]
        for c in children[u]:
            a = poly_add(msg0[c][u], msg1[c][u])
            prodA = _polymul(prodA, a)
            prodB = _polymul(prodB, msg0[c][u])
        if u != 0:
            p = parent[u]
            msg0[u][p] = prodA
            msg1[u][p] = [0] + prodB
        else:
            # root has no parent message
            pass

    # Upward messages via prefix/suffix products
    for u in order:
        neighs = adj[u]
        m = len(neighs)
        if m == 0:
            continue
        A = []
        B = []
        for v in neighs:
            # message from v to u must already exist (child->u or parent->u)
            a = poly_add(msg0[v][u], msg1[v][u])
            A.append(a)
            B.append(msg0[v][u])

        prefixA = [[1]]
        prefixB = [[1]]
        for i in range(m):
            prefixA.append(_polymul(prefixA[-1], A[i]))
            prefixB.append(_polymul(prefixB[-1], B[i]))
        suffixA = [[1]] * (m + 1)
        suffixB = [[1]] * (m + 1)
        suffixA[m] = [1]
        suffixB[m] = [1]
        for i in range(m - 1, -1, -1):
            suffixA[i] = _polymul(A[i], suffixA[i + 1])
            suffixB[i] = _polymul(B[i], suffixB[i + 1])

        for i, v in enumerate(neighs):
            prodA_excl = _polymul(prefixA[i], suffixA[i + 1])
            prodB_excl = _polymul(prefixB[i], suffixB[i + 1])
            msg0[u][v] = prodA_excl
            msg1[u][v] = [0] + prodB_excl

    return msg0, msg1


def forest_poly_at_vertex(adj: List[List[int]], w: int, msg0, msg1) -> List[int]:
    """Compute I(T-w) using messages for each neighbor component."""
    polys = []
    for u in adj[w]:
        polys.append(poly_add(msg0[u][w], msg1[u][w]))
    out = [1]
    for p in polys:
        out = _polymul(out, p)
    return out


def forest_poly_removed(adj: List[List[int]], removed: List[bool]) -> List[int]:
    n = len(adj)
    remaining = [i for i in range(n) if not removed[i]]
    if not remaining:
        return [1]
    rem_set = set(remaining)
    seen = set()
    out = [1]
    for start in remaining:
        if start in seen:
            continue
        # BFS component
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
        # build adjacency for component
        mapping = {old: i for i, old in enumerate(comp)}
        cadj: List[List[int]] = [[] for _ in range(len(comp))]
        for old in comp:
            i = mapping[old]
            for v in adj[old]:
                j = mapping.get(v)
                if j is not None:
                    cadj[i].append(j)
        poly_c = independence_poly(len(comp), cadj)
        out = _polymul(out, poly_c)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=2)
    ap.add_argument("--max-n", type=int, default=18)
    ap.add_argument("--backend", default="auto", choices=["auto", "geng", "networkx"])
    ap.add_argument("--out", default="")
    ap.add_argument("--max-witnesses", type=int, default=5)
    args = ap.parse_args()

    stats = {
        "min_n": args.min_n,
        "max_n": args.max_n,
        "backend": args.backend,
        "trees": 0,
        "vertex_cases": 0,
        "mode_gt_failures": 0,
        "mode_abs_gt1_failures": 0,
        "mode_drop_failures": 0,
        "mode_gt_d_f_failures": 0,
        "mode_g_h_failures": 0,
        "max_mode_g_minus_h": 0,
        "mode_gt_d_f_failures": 0,
        "max_mode_diff": 0,
        "witnesses": [],
        "by_n": [],
    }

    for n in range(args.min_n, args.max_n + 1):
        tree_count = 0
        vertex_cases = 0
        mode_gt_failures = 0
        mode_abs_gt1_failures = 0
        mode_drop_failures = 0
        mode_gt_d_f_failures = 0
        mode_g_h_failures = 0
        max_mode_g_minus_h = 0
        max_mode_diff = 0

        for _, adj in trees(n, backend=args.backend):
            tree_count += 1
            f = independence_poly(n, adj)
            mode_f = mode_index(f)
            d_f = first_descent(f)
            msg0, msg1 = compute_messages(adj)

            for w in range(n):
                vertex_cases += 1
                g = forest_poly_at_vertex(adj, w, msg0, msg1)
                mode_g = mode_index(g)
                removed = [False] * n
                removed[w] = True
                for nb in adj[w]:
                    removed[nb] = True
                h0 = forest_poly_removed(adj, removed)
                mode_h0 = mode_index(h0)
                diff = mode_g - mode_f
                if diff > 0:
                    mode_gt_failures += 1
                if abs(diff) > 1:
                    mode_abs_gt1_failures += 1
                if abs(diff) > max_mode_diff:
                    max_mode_diff = abs(diff)
                if mode_g - mode_h0 > 1:
                    mode_g_h_failures += 1
                    if len(stats["witnesses"]) < args.max_witnesses:
                        stats["witnesses"].append(
                            {
                                "n": n,
                                "graph6": encode_graph6_small(adj),
                                "w": w,
                                "mode_f": mode_f,
                                "mode_g": mode_g,
                                "mode_h0": mode_h0,
                                "diff_g_h0": mode_g - mode_h0,
                                "f": f,
                                "g": g,
                                "h0": h0,
                            }
                        )
                if mode_g - mode_h0 > max_mode_g_minus_h:
                    max_mode_g_minus_h = mode_g - mode_h0
                # boundary drop check at p = mode_f:
                # f_p - f_{p+1} >= g_p - g_{p-1}
                p = mode_f
                f_p = f[p]
                f_p1 = f[p + 1] if p + 1 < len(f) else 0
                g_p = g[p] if p < len(g) else 0
                g_pm1 = g[p - 1] if p - 1 >= 0 and p - 1 < len(g) else 0
                if (f_p - f_p1) < (g_p - g_pm1):
                    mode_drop_failures += 1
                    if len(stats["witnesses"]) < args.max_witnesses:
                        stats["witnesses"].append(
                            {
                                "n": n,
                                "graph6": encode_graph6_small(adj),
                                "w": w,
                                "mode_f": mode_f,
                                "mode_g": mode_g,
                                "diff": diff,
                                "f": f,
                                "g": g,
                                "drop_lhs": f_p - f_p1,
                                "drop_rhs": g_p - g_pm1,
                            }
                        )
                if d_f != -1 and mode_g > d_f:
                    mode_gt_d_f_failures += 1
                if (abs(diff) > 1 or diff > 0) and len(stats["witnesses"]) < args.max_witnesses:
                    stats["witnesses"].append(
                        {
                            "n": n,
                            "graph6": encode_graph6_small(adj),
                            "w": w,
                            "mode_f": mode_f,
                            "mode_g": mode_g,
                            "diff": diff,
                            "f": f,
                            "g": g,
                        }
                    )

        stats["trees"] += tree_count
        stats["vertex_cases"] += vertex_cases
        stats["mode_gt_failures"] += mode_gt_failures
        stats["mode_abs_gt1_failures"] += mode_abs_gt1_failures
        stats["mode_drop_failures"] += mode_drop_failures
        stats["mode_gt_d_f_failures"] += mode_gt_d_f_failures
        stats["mode_g_h_failures"] += mode_g_h_failures
        stats["max_mode_g_minus_h"] = max(stats["max_mode_g_minus_h"], max_mode_g_minus_h)
        stats["max_mode_diff"] = max(stats["max_mode_diff"], max_mode_diff)
        stats["by_n"].append(
            {
                "n": n,
                "trees": tree_count,
                "vertex_cases": vertex_cases,
                "mode_gt_failures": mode_gt_failures,
                "mode_abs_gt1_failures": mode_abs_gt1_failures,
                "mode_drop_failures": mode_drop_failures,
                "mode_gt_d_f_failures": mode_gt_d_f_failures,
                "mode_g_h_failures": mode_g_h_failures,
                "max_mode_g_minus_h": max_mode_g_minus_h,
                "max_mode_diff": max_mode_diff,
            }
        )

    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(stats, f, indent=2)
    else:
        print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
