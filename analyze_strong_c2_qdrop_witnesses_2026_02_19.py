#!/usr/bin/env python3
"""Deep analysis of Q-drop witnesses for STRONG C2.

Given witness graph6 strings, this script:
1) Reconstructs canonical (leaf, support, hub u), B, P, Q.
2) Roots B at u and builds child factors F_i=dp0[child], G_i=dp1[child].
3) Forms R=prod(F_i), P=prod(F_i+G_i), Q=xR.
4) Decomposes P by number of G factors:
     P = sum_{t=0}^d H_t,  H_t = sum_{|S|=t} prod_{i in S} G_i prod_{i notin S} F_i.
5) Reports slope contributions around mode-related windows.
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Any

from diagnose_bridge_decomposition import compute_hub_polys, mode_index_leftmost, remove_vertices
from graph6 import parse_graph6
from indpoly import _polyadd, _polymul, independence_poly


def coeff(poly: list[int], k: int) -> int:
    return poly[k] if 0 <= k < len(poly) else 0


def rooted_children(adj: list[list[int]], root: int) -> tuple[list[list[int]], list[int], list[int]]:
    """Return (children, parent, bfs_order) for tree rooted at root."""
    n = len(adj)
    parent = [-1] * n
    children = [[] for _ in range(n)]
    seen = [False] * n
    seen[root] = True
    q = [root]
    order = []
    head = 0
    while head < len(q):
        v = q[head]
        head += 1
        order.append(v)
        for w in adj[v]:
            if not seen[w]:
                seen[w] = True
                parent[w] = v
                children[v].append(w)
                q.append(w)
    return children, parent, order


def rooted_dp(adj: list[list[int]], root: int) -> tuple[list[list[int]], list[list[int]], list[list[int]]]:
    children, _, _ = rooted_children(adj, root)
    n = len(adj)
    order = []
    stack = [(root, False)]
    while stack:
        v, done = stack.pop()
        if done:
            order.append(v)
            continue
        stack.append((v, True))
        for c in children[v]:
            stack.append((c, False))

    dp0 = [[] for _ in range(n)]
    dp1 = [[] for _ in range(n)]
    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
            continue

        p0 = [1]
        for c in children[v]:
            p0 = _polymul(p0, _polyadd(dp0[c], dp1[c]))
        dp0[v] = p0

        p1 = [1]
        for c in children[v]:
            p1 = _polymul(p1, dp0[c])
        dp1[v] = [0] + p1

    return children, dp0, dp1


def subtree_sizes(children: list[list[int]], root: int) -> list[int]:
    n = len(children)
    size = [1] * n
    order = []
    stack = [(root, False)]
    while stack:
        v, done = stack.pop()
        if done:
            order.append(v)
            continue
        stack.append((v, True))
        for c in children[v]:
            stack.append((c, False))
    for v in order:
        for c in children[v]:
            size[v] += size[c]
    return size


def poly_sum(polys: list[list[int]]) -> list[int]:
    out: list[int] = []
    for p in polys:
        out = _polyadd(out, p)
    return out


def poly_sub(a: list[int], b: list[int]) -> list[int]:
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = coeff(a, i) - coeff(b, i)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def analyze_one_g6(g6: str) -> dict[str, Any]:
    nn, adj = parse_graph6((g6.strip() + "\n").encode("ascii"))

    poly_t = independence_poly(nn, adj)
    m = mode_index_leftmost(poly_t)

    deg = [len(nb) for nb in adj]
    leaves = [v for v in range(nn) if deg[v] == 1]
    if not leaves:
        raise ValueError("Input graph has no leaves; expected a tree")
    min_parent_deg = min(deg[adj[l][0]] for l in leaves)
    leaf = min(l for l in leaves if deg[adj[l][0]] == min_parent_deg)
    support = adj[leaf][0]
    if deg[support] != 2:
        raise ValueError("Canonical support is not degree-2; this is outside bridge regime")
    u = [x for x in adj[support] if x != leaf][0]

    b_adj = remove_vertices(adj, {leaf, support})
    keep = [v for v in range(nn) if v not in {leaf, support}]
    idx_map = {v: i for i, v in enumerate(keep)}
    u_in_b = idx_map[u]

    b_poly = independence_poly(len(b_adj), b_adj)
    p_poly, q_poly = compute_hub_polys(b_adj, u_in_b)

    children, dp0, dp1 = rooted_dp(b_adj, u_in_b)
    size = subtree_sizes(children, u_in_b)
    child_nodes = children[u_in_b]
    d = len(child_nodes)

    F = [dp0[c] for c in child_nodes]
    G = [dp1[c] for c in child_nodes]
    I = [_polyadd(dp0[c], dp1[c]) for c in child_nodes]

    # Build R and rebuild P from factors.
    r_poly = [1]
    for f in F:
        r_poly = _polymul(r_poly, f)
    p_rebuild = [1]
    for i_poly in I:
        p_rebuild = _polymul(p_rebuild, i_poly)
    if p_rebuild != p_poly:
        raise RuntimeError("P reconstruction mismatch")

    # Decompose P by number of G-factors: H_t.
    H = [[1]] + [[] for _ in range(d)]
    for i in range(d):
        newH = [[] for _ in range(d + 1)]
        for t in range(i + 1):
            if H[t]:
                term_f = _polymul(H[t], F[i])
                newH[t] = _polyadd(newH[t], term_f)
                term_g = _polymul(H[t], G[i])
                newH[t + 1] = _polyadd(newH[t + 1], term_g)
        H = newH

    if H[0] != r_poly:
        raise RuntimeError("H_0 mismatch with R")
    if poly_sum(H) != p_poly:
        raise RuntimeError("Sum_t H_t mismatch with P")

    s_poly = poly_sub(p_poly, r_poly)
    s_from_H = poly_sum(H[1:])
    if s_poly != s_from_H:
        raise RuntimeError("S mismatch: P-R != sum_{t>=1} H_t")

    k_list = [m - 3, m - 2, m - 1, m]
    k_windows = [
        {"name": "q_window", "k0": m - 3, "k1": m - 2},
        {"name": "same_window", "k0": m - 2, "k1": m - 1},
    ]

    h_levels = []
    for t in range(d + 1):
        ht = H[t]
        info = {"t": t, "coeffs": {}, "deltas": {}}
        for k in k_list:
            info["coeffs"][f"k{k}"] = coeff(ht, k)
        for w in k_windows:
            info["deltas"][w["name"]] = coeff(ht, w["k1"]) - coeff(ht, w["k0"])
        h_levels.append(info)

    child_info = []
    for c in child_nodes:
        fc = dp0[c]
        gc = dp1[c]
        ic = _polyadd(fc, gc)
        child_info.append(
            {
                "child_index_in_B": c,
                "subtree_size": size[c],
                "deg_in_B": len(b_adj[c]),
                "mode_F": mode_index_leftmost(fc),
                "mode_G": mode_index_leftmost(gc) if gc else 0,
                "mode_I": mode_index_leftmost(ic),
                "F": fc,
                "G": gc,
                "I": ic,
                "delta_F_same_window": coeff(fc, m - 1) - coeff(fc, m - 2),
                "delta_I_same_window": coeff(ic, m - 1) - coeff(ic, m - 2),
            }
        )

    p0 = coeff(p_poly, m - 2)
    p1 = coeff(p_poly, m - 1)
    q0 = coeff(q_poly, m - 2)
    q1 = coeff(q_poly, m - 1)
    b0 = coeff(b_poly, m - 2)
    b1 = coeff(b_poly, m - 1)

    return {
        "n": nn,
        "g6": g6,
        "mode_T": m,
        "mode_B": mode_index_leftmost(b_poly),
        "shift": mode_index_leftmost(b_poly) - (m - 1),
        "leaf": leaf,
        "support": support,
        "hub_u_in_T": u,
        "hub_u_in_B": u_in_b,
        "deg_u_in_B": len(b_adj[u_in_b]),
        "I_T": poly_t,
        "I_B": b_poly,
        "P": p_poly,
        "Q": q_poly,
        "R": r_poly,
        "S": s_poly,
        "q_drop": q1 < q0,
        "delta_q": q1 - q0,
        "delta_b": b1 - b0,
        "delta_p": p1 - p0,
        "p1_over_b1": (p1 / b1) if b1 > 0 else None,
        "q_index_check": {
            "q_m2": q0,
            "R_m3": coeff(r_poly, m - 3),
            "q_m1": q1,
            "R_m2": coeff(r_poly, m - 2),
            "ok": (q0 == coeff(r_poly, m - 3) and q1 == coeff(r_poly, m - 2)),
        },
        "windows": {
            "q_window": {
                "k0": m - 3,
                "k1": m - 2,
                "delta_R": coeff(r_poly, m - 2) - coeff(r_poly, m - 3),
                "delta_Q": q1 - q0,
                "delta_P": coeff(p_poly, m - 2) - coeff(p_poly, m - 3),
                "delta_S": coeff(s_poly, m - 2) - coeff(s_poly, m - 3),
            },
            "same_window": {
                "k0": m - 2,
                "k1": m - 1,
                "delta_R": coeff(r_poly, m - 1) - coeff(r_poly, m - 2),
                "delta_P": coeff(p_poly, m - 1) - coeff(p_poly, m - 2),
                "delta_S": coeff(s_poly, m - 1) - coeff(s_poly, m - 2),
            },
        },
        "H_levels": h_levels,
        "children": child_info,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Analyze STRONG C2 Q-drop witnesses in depth.")
    ap.add_argument(
        "--g6",
        action="append",
        default=[
            "S???????C?G?G?C?@??G??_?@??@?F~_?",
            "V???????????_?O?C??_?A??C??C??A???_?{E??^g??",
        ],
        help="Witness graph6 string (may be passed multiple times).",
    )
    ap.add_argument("--out", default="results/analyze_strong_c2_qdrop_witnesses_2026_02_19.json")
    args = ap.parse_args()

    # argparse append + list default means defaults stay when user also passes --g6.
    # Deduplicate while preserving order.
    seen = set()
    g6_list = []
    for s in args.g6:
        if s not in seen:
            seen.add(s)
            g6_list.append(s)

    out = {
        "params": {"witness_count": len(g6_list)},
        "witnesses": [],
    }
    for g6 in g6_list:
        out["witnesses"].append(analyze_one_g6(g6))

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)

    print(f"Wrote {args.out}")
    for w in out["witnesses"]:
        print(
            f"n={w['n']}: q_drop={w['q_drop']} shift={w['shift']} "
            f"deg_u={w['deg_u_in_B']} delta_q={w['delta_q']} "
            f"delta_P_same={w['windows']['same_window']['delta_P']} "
            f"delta_R_same={w['windows']['same_window']['delta_R']} "
            f"delta_S_same={w['windows']['same_window']['delta_S']}",
            flush=True,
        )


if __name__ == "__main__":
    main()
