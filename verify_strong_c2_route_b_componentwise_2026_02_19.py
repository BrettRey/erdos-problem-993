#!/usr/bin/env python3
"""Route-B diagnostics for algebraic positivity of R = cross + mismatch.

For a degree-2 support leaf decomposition:
  B = T - {l, s}
  P = dp_B[u][0] = prod_c (g_c + h_c)
  Q = dp_B[u][1] = x * prod_c g_c = x * G

At r = m-2 (where m is the leftmost mode index of I(T)), define:
  p0 = P_r, p1 = P_{r+1}, p2 = P_{r+2}
  q0 = Q_r = G_{r-1}
  q1 = Q_{r+1} = G_r
  qm = Q_{r+2} = G_{r+1}

Residual:
  R = 2*p1*q1 - p2*q0 - p0*qm + p0*q1 - p1*q0
    = D1 + D2 + D3
where
  D1 = p1*q1 - p2*q0
  D2 = p1*q1 - p0*qm
  D3 = p0*q1 - p1*q0  (mismatch wrt Q-shift form)

Key linearity: for fixed G and index r, the functional
  L(A; G, r) := 2*a_{r+1}*G_r - a_{r+2}*G_{r-1}
                - a_r*G_{r+1} + a_r*G_r - a_{r+1}*G_{r-1}
is linear in A and satisfies R = L(P; G, r).

Since P = sum_{S subset children(u)} H_S with
  H_S = (prod_{i in S} h_i) * (prod_{i notin S} g_i),
we get:
  R = sum_S L(H_S; G, r).

This script checks signs of:
1) D1, D2, D3 and useful pairwise sums.
2) Layer contributions L_j = L(E_j; G, r), where E_j = sum_{|S|=j} H_S.
3) Optional explicit per-subset contributions L(H_S; G, r) for small degree.
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
from indpoly import independence_poly, _polyadd, _polymul


def getc(poly: list[int], k: int) -> int:
    return poly[k] if 0 <= k < len(poly) else 0


def compute_rooted_child_factors(adj_b: list[list[int]], u_in_b: int) -> list[tuple[list[int], list[int], list[int]]]:
    """Return per-child factors (f, g, h) for root u_in_b in B.

    f = I(T_c) = dp0[c] + dp1[c]
    g = dp0[c]
    h = dp1[c]
    """
    n = len(adj_b)
    if n <= 1:
        return []

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
        prod0 = [1]
        for c in children[v]:
            prod0 = _polymul(prod0, _polyadd(dp0[c], dp1[c]))
        dp0[v] = prod0

        prod1 = [1]
        for c in children[v]:
            prod1 = _polymul(prod1, dp0[c])
        dp1[v] = [0] + prod1

    out = []
    for c in children[u_in_b]:
        g = dp0[c]
        h = dp1[c]
        f = _polyadd(g, h)
        out.append((f, g, h))
    return out


def poly_prod(polys: list[list[int]]) -> list[int]:
    out = [1]
    for p in polys:
        out = _polymul(out, p)
    return out


def L_of_A(A: list[int], G: list[int], r: int) -> int:
    """Linear functional in A giving R when A=P."""
    a0 = getc(A, r)
    a1 = getc(A, r + 1)
    a2 = getc(A, r + 2)
    g0 = getc(G, r - 1)
    g1 = getc(G, r)
    g2 = getc(G, r + 1)
    return 2 * a1 * g1 - a2 * g0 - a0 * g2 + a0 * g1 - a1 * g0


def explicit_subset_contribs(gs: list[list[int]], hs: list[list[int]], G: list[int], r: int) -> tuple[int, int]:
    """Return (num_negative, min_value) over all subset contributions L(H_S)."""
    d = len(gs)
    if d == 0:
        v = L_of_A([1], G, r)
        return (1 if v < 0 else 0, v)

    neg = 0
    min_v: int | None = None
    for mask in range(1 << d):
        A = [1]
        for i in range(d):
            A = _polymul(A, hs[i] if (mask & (1 << i)) else gs[i])
        v = L_of_A(A, G, r)
        if v < 0:
            neg += 1
        if min_v is None or v < min_v:
            min_v = v
    assert min_v is not None
    return neg, min_v


def layer_polys_by_h_count(gs: list[list[int]], hs: list[list[int]]) -> list[list[int]]:
    """Build E_j polynomials where E_j sums products with exactly j h-factors."""
    d = len(gs)
    # dp[j] stores polynomial E_j after processed prefix.
    dp: list[list[int]] = [[1]]
    for i in range(d):
        g = gs[i]
        h = hs[i]
        new: list[list[int]] = [[] for _ in range(len(dp) + 1)]
        for j in range(len(dp)):
            a = dp[j]
            ag = _polymul(a, g)
            ah = _polymul(a, h)
            new[j] = _polyadd(new[j], ag) if new[j] else ag
            new[j + 1] = _polyadd(new[j + 1], ah) if new[j + 1] else ah
        dp = new
    return dp


def update_min(stats: dict[str, Any], key: str, value: int, witness: dict[str, Any]) -> None:
    cur = stats.get(key)
    if cur is None or value < cur:
        stats[key] = value
        stats[key + "_witness"] = witness


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--out", default="results/verify_strong_c2_route_b_componentwise_2026_02_19.json")
    ap.add_argument(
        "--explicit-subset-max-degree",
        type=int,
        default=8,
        help="Explicitly enumerate all subsets only when deg(u in B) <= this value.",
    )
    args = ap.parse_args()

    stats: dict[str, Any] = {
        "seen": 0,
        "considered": 0,
        "checked": 0,
        "R_neg": 0,
        "D1_neg": 0,
        "D2_neg": 0,
        "D3_neg": 0,
        "D1_plus_D3_neg": 0,
        "D2_plus_D3_neg": 0,
        "cross_minus_abs_mismatch_neg": 0,
        "R_min": None,
        "D1_min": None,
        "D2_min": None,
        "D3_min": None,
        "subset_explicit_checked": 0,
        "subset_explicit_neg_contrib_trees": 0,
        "subset_explicit_total_neg_contrib": 0,
        "subset_explicit_min_contrib": None,
        "layer_neg_trees": 0,
        "layer_min_contrib": None,
        "layer_stats": {},  # j -> {checked, neg, min}
        "identity_fail_P_expand": 0,
        "identity_fail_R_linear": 0,
        "wall_s": 0.0,
    }

    t0 = time.time()
    for n in range(args.min_n, args.max_n + 1):
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert proc.stdout is not None

        for raw in proc.stdout:
            nn, adj = parse_graph6(raw)
            stats["seen"] += 1
            if not is_dleaf_le_1(nn, adj):
                continue
            stats["considered"] += 1

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
            child_factors = compute_rooted_child_factors(b_adj, u_in_b)
            fs = [f for (f, _, _) in child_factors]
            gs = [g for (_, g, _) in child_factors]
            hs = [h for (_, _, h) in child_factors]

            P = poly_prod(fs) if fs else [1]
            G = poly_prod(gs) if gs else [1]
            Q = [0] + G

            # Sanity: I(B) = P + Q
            if _polyadd(P, Q) != b_poly:
                stats["identity_fail_P_expand"] += 1

            stats["checked"] += 1
            r = m - 2

            p0 = getc(P, r)
            p1 = getc(P, r + 1)
            p2 = getc(P, r + 2)
            q0 = getc(Q, r)
            q1 = getc(Q, r + 1)
            qm = getc(Q, r + 2)

            D1 = p1 * q1 - p2 * q0
            D2 = p1 * q1 - p0 * qm
            D3 = p0 * q1 - p1 * q0
            cross = D1 + D2
            mismatch = D3
            R = cross + mismatch

            witness = {
                "n": nn,
                "m": m,
                "r": r,
                "g6": raw.decode("ascii").strip(),
                "deg_u_B": len(child_factors),
                "p0": p0,
                "p1": p1,
                "p2": p2,
                "q0": q0,
                "q1": q1,
                "qm": qm,
                "D1": D1,
                "D2": D2,
                "D3": D3,
                "R": R,
            }

            if D1 < 0:
                stats["D1_neg"] += 1
            if D2 < 0:
                stats["D2_neg"] += 1
            if D3 < 0:
                stats["D3_neg"] += 1
            if D1 + D3 < 0:
                stats["D1_plus_D3_neg"] += 1
            if D2 + D3 < 0:
                stats["D2_plus_D3_neg"] += 1
            if cross - abs(mismatch) < 0:
                stats["cross_minus_abs_mismatch_neg"] += 1
            if R < 0:
                stats["R_neg"] += 1

            update_min(stats, "R_min", R, witness)
            update_min(stats, "D1_min", D1, witness)
            update_min(stats, "D2_min", D2, witness)
            update_min(stats, "D3_min", D3, witness)

            # Linearity identity check: R == L(P; G, r).
            Lp = L_of_A(P, G, r)
            if Lp != R:
                stats["identity_fail_R_linear"] += 1

            # Layer decomposition E_j and contribution signs.
            layer_polys = layer_polys_by_h_count(gs, hs)
            layer_neg_here = False
            layer_sum = 0
            for j, Ej in enumerate(layer_polys):
                Lj = L_of_A(Ej, G, r)
                layer_sum += Lj
                key = str(j)
                if key not in stats["layer_stats"]:
                    stats["layer_stats"][key] = {"checked": 0, "neg": 0, "min": None}
                info = stats["layer_stats"][key]
                info["checked"] += 1
                if Lj < 0:
                    info["neg"] += 1
                    layer_neg_here = True
                if info["min"] is None or Lj < info["min"]:
                    info["min"] = Lj

                if stats["layer_min_contrib"] is None or Lj < stats["layer_min_contrib"]:
                    stats["layer_min_contrib"] = Lj
            if layer_sum != R:
                stats["identity_fail_R_linear"] += 1
            if layer_neg_here:
                stats["layer_neg_trees"] += 1

            # Optional explicit per-subset contributions.
            d = len(child_factors)
            if d <= args.explicit_subset_max_degree:
                neg_sub, min_sub = explicit_subset_contribs(gs, hs, G, r)
                stats["subset_explicit_checked"] += 1
                if neg_sub > 0:
                    stats["subset_explicit_neg_contrib_trees"] += 1
                    stats["subset_explicit_total_neg_contrib"] += neg_sub
                if (
                    stats["subset_explicit_min_contrib"] is None
                    or min_sub < stats["subset_explicit_min_contrib"]
                ):
                    stats["subset_explicit_min_contrib"] = min_sub

        proc.wait()
        print(
            f"n={n:2d}: checked={stats['checked']:8d} R_neg={stats['R_neg']} "
            f"D1_neg={stats['D1_neg']} D2_neg={stats['D2_neg']} D3_neg={stats['D3_neg']} "
            f"layer_neg_trees={stats['layer_neg_trees']} explicit={stats['subset_explicit_checked']}",
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
    print(f"R_neg: {stats['R_neg']}   R_min: {stats['R_min']}")
    print(
        f"D1_neg: {stats['D1_neg']} D2_neg: {stats['D2_neg']} D3_neg: {stats['D3_neg']}"
    )
    print(
        f"D1+D3 neg: {stats['D1_plus_D3_neg']}  D2+D3 neg: {stats['D2_plus_D3_neg']}"
    )
    print(f"cross-|mismatch| neg: {stats['cross_minus_abs_mismatch_neg']}")
    print(
        f"layer_neg_trees: {stats['layer_neg_trees']}  layer_min_contrib: {stats['layer_min_contrib']}"
    )
    print(
        f"explicit_subset_checked: {stats['subset_explicit_checked']}  "
        f"explicit_neg_trees: {stats['subset_explicit_neg_contrib_trees']}  "
        f"explicit_min_contrib: {stats['subset_explicit_min_contrib']}"
    )
    print(f"identity_fail_P_expand: {stats['identity_fail_P_expand']}")
    print(f"identity_fail_R_linear: {stats['identity_fail_R_linear']}")
    print(f"wall_s: {stats['wall_s']:.1f}")
    print(f"wrote: {args.out}")


if __name__ == "__main__":
    main()
