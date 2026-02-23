#!/usr/bin/env python3
"""Subset-contribution diagnostics for M := D2 + D3.

Definitions at r=m-2:
  D2 = p1*q1 - p0*qm
  D3 = p0*q1 - p1*q0
  M  = D2 + D3 = p1*(q1-q0) + p0*(q1-qm)

Write Q = x*G so q0=G_{r-1}, q1=G_r, qm=G_{r+1}. For any A define
  M_A(G,r) := a_{r+1}*(G_r - G_{r-1}) + a_r*(G_r - G_{r+1}).
Then M = M_P(G,r), and M_A is linear in A.

With P = sum_{S} H_S (subset expansion of prod(g_i+h_i)), we get:
  M = sum_S M_{H_S}(G,r).

This script checks whether individual subset contributions M_{H_S} are nonnegative.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from typing import Any

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import mode_index_leftmost, remove_vertices
from graph6 import parse_graph6
from indpoly import independence_poly, _polymul, _polyadd


def getc(poly: list[int], k: int) -> int:
    return poly[k] if 0 <= k < len(poly) else 0


def compute_rooted_child_factors(adj_b: list[list[int]], u_in_b: int) -> tuple[list[list[int]], list[list[int]], list[list[int]]]:
    n = len(adj_b)
    if n <= 1:
        return [], [], []

    children = [[] for _ in range(n)]
    visited = [False] * n
    visited[u_in_b] = True
    bfs_queue = [u_in_b]
    head = 0
    while head < len(bfs_queue):
        v = bfs_queue[head]
        head += 1
        for w in adj_b[v]:
            if not visited[w]:
                visited[w] = True
                children[v].append(w)
                bfs_queue.append(w)

    order: list[int] = []
    stack: list[tuple[int, bool]] = [(u_in_b, False)]
    while stack:
        v, done = stack.pop()
        if done:
            order.append(v)
            continue
        stack.append((v, True))
        for c in children[v]:
            stack.append((c, False))

    dp0: list[list[int]] = [[] for _ in range(n)]
    dp1: list[list[int]] = [[] for _ in range(n)]
    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
            continue
        a0 = [1]
        for c in children[v]:
            a0 = _polymul(a0, _polyadd(dp0[c], dp1[c]))
        dp0[v] = a0
        a1 = [1]
        for c in children[v]:
            a1 = _polymul(a1, dp0[c])
        dp1[v] = [0] + a1

    fs: list[list[int]] = []
    gs: list[list[int]] = []
    hs: list[list[int]] = []
    for c in children[u_in_b]:
        g = dp0[c]
        h = dp1[c]
        f = _polyadd(g, h)
        fs.append(f)
        gs.append(g)
        hs.append(h)
    return fs, gs, hs


def poly_prod(polys: list[list[int]]) -> list[int]:
    out = [1]
    for p in polys:
        out = _polymul(out, p)
    return out


def M_of_A(A: list[int], G: list[int], r: int) -> int:
    return getc(A, r + 1) * (getc(G, r) - getc(G, r - 1)) + getc(A, r) * (getc(G, r) - getc(G, r + 1))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--explicit-subset-max-degree", type=int, default=8)
    ap.add_argument("--out", default="results/verify_strong_c2_route_b_subset_M_2026_02_19.json")
    args = ap.parse_args()

    stats: dict[str, Any] = {
        "checked": 0,
        "M_neg": 0,
        "M_min": None,
        "M_subset_checked_trees": 0,
        "M_subset_neg_trees": 0,
        "M_subset_neg_total": 0,
        "M_subset_min": None,
        "M_identity_fail": 0,
        "wall_s": 0.0,
    }

    t0 = time.time()
    for n in range(args.min_n, args.max_n + 1):
        proc = subprocess.Popen(
            [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        assert proc.stdout is not None

        for raw in proc.stdout:
            nn, adj = parse_graph6(raw)
            if not is_dleaf_le_1(nn, adj):
                continue

            poly_t = independence_poly(nn, adj)
            m = mode_index_leftmost(poly_t)
            if m == 0 or m >= len(poly_t) or poly_t[m] == 0:
                continue

            deg = [len(nb) for nb in adj]
            leaves = [v for v in range(nn) if deg[v] == 1]
            min_parent_deg = min(deg[adj[l][0]] for l in leaves)
            leaf = min(l for l in leaves if deg[adj[l][0]] == min_parent_deg)
            support = adj[leaf][0]
            if deg[support] != 2:
                continue

            u = [x for x in adj[support] if x != leaf][0]
            b_adj = remove_vertices(adj, {leaf, support})
            if len(b_adj) == 0:
                continue
            b_poly = independence_poly(len(b_adj), b_adj)
            if m - 2 < 0 or m >= len(b_poly) or m - 1 >= len(b_poly) or b_poly[m - 1] == 0:
                continue

            keep = [v for v in range(nn) if v not in {leaf, support}]
            idx_map = {v: i for i, v in enumerate(keep)}
            u_in_b = idx_map[u]
            fs, gs, hs = compute_rooted_child_factors(b_adj, u_in_b)
            P = poly_prod(fs) if fs else [1]
            G = poly_prod(gs) if gs else [1]
            Q = [0] + G

            r = m - 2
            p0 = getc(P, r)
            p1 = getc(P, r + 1)
            q0 = getc(Q, r)
            q1 = getc(Q, r + 1)
            qm = getc(Q, r + 2)
            M = p1 * (q1 - q0) + p0 * (q1 - qm)
            stats["checked"] += 1

            if M < 0:
                stats["M_neg"] += 1
            if stats["M_min"] is None or M < stats["M_min"]:
                stats["M_min"] = M

            # Identity M == M_of_A(P,G,r)
            if M != M_of_A(P, G, r):
                stats["M_identity_fail"] += 1

            d = len(gs)
            if d <= args.explicit_subset_max_degree:
                stats["M_subset_checked_trees"] += 1
                neg_here = 0
                min_here = None
                total = 0
                for mask in range(1 << d):
                    A = [1]
                    for i in range(d):
                        A = _polymul(A, hs[i] if (mask & (1 << i)) else gs[i])
                    v = M_of_A(A, G, r)
                    total += v
                    if v < 0:
                        neg_here += 1
                    if min_here is None or v < min_here:
                        min_here = v
                if total != M:
                    stats["M_identity_fail"] += 1
                if neg_here > 0:
                    stats["M_subset_neg_trees"] += 1
                    stats["M_subset_neg_total"] += neg_here
                if stats["M_subset_min"] is None or (min_here is not None and min_here < stats["M_subset_min"]):
                    stats["M_subset_min"] = min_here

        proc.wait()
        print(
            f"n={n:2d}: checked={stats['checked']:8d} "
            f"M_neg={stats['M_neg']} subset_trees={stats['M_subset_checked_trees']} "
            f"subset_neg_trees={stats['M_subset_neg_trees']}",
            flush=True,
        )

    stats["wall_s"] = time.time() - t0
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    print("\n" + "=" * 72)
    print(f"checked: {stats['checked']}")
    print(f"M_neg: {stats['M_neg']}  M_min: {stats['M_min']}")
    print(f"subset_checked_trees: {stats['M_subset_checked_trees']}")
    print(f"subset_neg_trees: {stats['M_subset_neg_trees']}")
    print(f"subset_neg_total: {stats['M_subset_neg_total']}")
    print(f"subset_min: {stats['M_subset_min']}")
    print(f"M_identity_fail: {stats['M_identity_fail']}")
    print(f"wall_s: {stats['wall_s']:.1f}")
    print(f"wrote: {args.out}")


if __name__ == "__main__":
    main()
