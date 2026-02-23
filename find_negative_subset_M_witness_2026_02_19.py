#!/usr/bin/env python3
"""Find first tree where a subset contribution to M is negative.

M_A(G,r) = a_{r+1}(G_r-G_{r-1}) + a_r(G_r-G_{r+1})
for A = prod(h_i for i in S) prod(g_i for i not in S).
"""

from __future__ import annotations

import argparse
import subprocess

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import mode_index_leftmost, remove_vertices
from graph6 import parse_graph6
from indpoly import independence_poly, _polymul, _polyadd


def getc(poly: list[int], k: int) -> int:
    return poly[k] if 0 <= k < len(poly) else 0


def compute_rooted_child_factors(adj_b: list[list[int]], u_in_b: int):
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

    order = []
    stack = [(u_in_b, False)]
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
        a0 = [1]
        for c in children[v]:
            a0 = _polymul(a0, _polyadd(dp0[c], dp1[c]))
        dp0[v] = a0
        a1 = [1]
        for c in children[v]:
            a1 = _polymul(a1, dp0[c])
        dp1[v] = [0] + a1

    fs, gs, hs = [], [], []
    for c in children[u_in_b]:
        g = dp0[c]
        h = dp1[c]
        f = _polyadd(g, h)
        fs.append(f)
        gs.append(g)
        hs.append(h)
    return fs, gs, hs


def poly_prod(polys):
    out = [1]
    for p in polys:
        out = _polymul(out, p)
    return out


def M_of_A(A, G, r):
    return getc(A, r + 1) * (getc(G, r) - getc(G, r - 1)) + getc(A, r) * (getc(G, r) - getc(G, r + 1))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=23)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--explicit-subset-max-degree", type=int, default=8)
    args = ap.parse_args()

    checked = 0
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
            d = len(gs)
            if d > args.explicit_subset_max_degree:
                continue

            checked += 1
            r = m - 2
            G = poly_prod(gs) if gs else [1]
            P = poly_prod(fs) if fs else [1]
            M_total = M_of_A(P, G, r)

            for mask in range(1 << d):
                A = [1]
                for i in range(d):
                    A = _polymul(A, hs[i] if (mask & (1 << i)) else gs[i])
                v = M_of_A(A, G, r)
                if v < 0:
                    print("FOUND")
                    print("n", nn, "m", m, "r", r, "deg_u_B", d, "mask", mask, "value", v)
                    print("g6", raw.decode("ascii").strip())
                    print("M_total", M_total)
                    print("G_r-1,G_r,G_r+1", getc(G, r-1), getc(G, r), getc(G, r+1))
                    print("A_r,A_r+1", getc(A, r), getc(A, r+1))
                    return

        proc.wait()
        print(f"n={n}: checked={checked}", flush=True)

    print("No negative subset contributions found")


if __name__ == "__main__":
    main()
