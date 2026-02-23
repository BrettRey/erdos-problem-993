#!/usr/bin/env python3
"""Verify whether Q-drop forces shift-0 in the STRONG C2 bridge setup.

For each d_leaf<=1 tree T (or all trees with --all-trees), choose canonical leaf:
  - leaf with minimum support degree, tie by smallest leaf id
and require deg(support)=2.

Let B = T-{leaf,support}, with hub u (other neighbor of support), and root B at u:
  P = dp_B[u][0]
  Q = dp_B[u][1] = x * R, where R = product(dp0[child]).

Define m = leftmost mode(I(T)). At index m-1:
  q_drop <=> q_{m-1} < q_{m-2}.

This script verifies:
1) q_drop => mode(B) = m-1 (shift-0?) on scanned frontier.
2) R and Q indexing identities at witness points.
3) Optional per-witness details for downstream analysis.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from typing import Any

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import compute_hub_polys, mode_index_leftmost, remove_vertices
from graph6 import parse_graph6
from indpoly import _polyadd, _polymul, independence_poly


def rooted_children(adj: list[list[int]], root: int) -> tuple[list[list[int]], list[int]]:
    """Return rooted children lists and parent array for a tree rooted at root."""
    n = len(adj)
    parent = [-1] * n
    children = [[] for _ in range(n)]
    seen = [False] * n
    seen[root] = True
    q = [root]
    head = 0
    while head < len(q):
        v = q[head]
        head += 1
        for w in adj[v]:
            if not seen[w]:
                seen[w] = True
                parent[w] = v
                children[v].append(w)
                q.append(w)
    return children, parent


def rooted_dp(adj: list[list[int]], root: int) -> tuple[list[list[int]], list[list[int]], list[list[int]]]:
    """Return (children, dp0, dp1) for tree rooted at root."""
    children, _ = rooted_children(adj, root)
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

        prod0 = [1]
        for c in children[v]:
            prod0 = _polymul(prod0, _polyadd(dp0[c], dp1[c]))
        dp0[v] = prod0

        prod1 = [1]
        for c in children[v]:
            prod1 = _polymul(prod1, dp0[c])
        dp1[v] = [0] + prod1

    return children, dp0, dp1


def coeff(poly: list[int], k: int) -> int:
    return poly[k] if 0 <= k < len(poly) else 0


def prod_all(polys: list[list[int]]) -> list[int]:
    out = [1]
    for p in polys:
        out = _polymul(out, p)
    return out


def make_witness_payload(
    *,
    n: int,
    g6: str,
    m: int,
    mode_b: int,
    leaf: int,
    support: int,
    u: int,
    p_poly: list[int],
    q_poly: list[int],
    b_poly: list[int],
    r_poly: list[int],
) -> dict[str, Any]:
    p0 = coeff(p_poly, m - 2)
    p1 = coeff(p_poly, m - 1)
    pm = coeff(p_poly, m)

    q0 = coeff(q_poly, m - 2)
    q1 = coeff(q_poly, m - 1)
    qm = coeff(q_poly, m)

    b0 = coeff(b_poly, m - 2)
    b1 = coeff(b_poly, m - 1)
    b2 = coeff(b_poly, m)

    db = b1 - b0
    dp = p1 - p0
    dq = q1 - q0
    shift = mode_b - (m - 1)

    r_m3 = coeff(r_poly, m - 3)
    r_m2 = coeff(r_poly, m - 2)
    r_m1 = coeff(r_poly, m - 1)

    # Q = x*R implies q_{m-2}=R_{m-3}, q_{m-1}=R_{m-2}
    q_index_ok = (q0 == r_m3) and (q1 == r_m2)
    transfer_ratio = ((-dq) / db) if (dq < 0 and db > 0) else None
    need_ratio = (p1 / b1) if b1 > 0 else None

    return {
        "n": n,
        "g6": g6,
        "mode_T": m,
        "mode_B": mode_b,
        "shift": shift,
        "leaf": leaf,
        "support": support,
        "hub_u_in_T": u,
        "p_m2": p0,
        "p_m1": p1,
        "p_m": pm,
        "q_m2": q0,
        "q_m1": q1,
        "q_m": qm,
        "b_m2": b0,
        "b_m1": b1,
        "b_m": b2,
        "delta_p": dp,
        "delta_q": dq,
        "delta_b": db,
        "q_drop": q1 < q0,
        "p1_over_b1": need_ratio,
        "transfer_ratio": transfer_ratio,
        "ratio_gap": (need_ratio - transfer_ratio) if (need_ratio is not None and transfer_ratio is not None) else None,
        "mode_R": mode_index_leftmost(r_poly) if r_poly else 0,
        "R_m3": r_m3,
        "R_m2": r_m2,
        "R_m1": r_m1,
        "delta_R_qwindow": r_m2 - r_m3,   # aligns with dq
        "delta_R_same_window": r_m1 - r_m2,  # user's "m-1 window" wording
        "delta_P_same_window": coeff(p_poly, m - 1) - coeff(p_poly, m - 2),
        "delta_S_same_window": (coeff(p_poly, m - 1) - coeff(p_poly, m - 2)) - (r_m1 - r_m2),
        "q_index_identity_ok": q_index_ok,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Verify Q-drop implies shift-0 in STRONG C2 Route A setup.")
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--mod", type=int, default=1, help="Use nauty residue classes res/mod.")
    ap.add_argument("--all-trees", action="store_true", help="Disable d_leaf<=1 filter.")
    ap.add_argument("--out", default="results/verify_strong_c2_qdrop_shift0_2026_02_19.json")
    args = ap.parse_args()

    summary: dict[str, Any] = {
        "seen": 0,
        "considered": 0,
        "checked": 0,
        "skip_mode": 0,
        "skip_non_deg2": 0,
        "skip_short_B": 0,
        "qdrop_total": 0,
        "qdrop_shift0": 0,
        "qdrop_shift1": 0,
        "qdrop_shift_other": 0,
        "qdrop_shift_nonzero": 0,
        "qdrop_q_index_identity_fail": 0,
        "qdrop_R_decreasing_same_window": 0,
        "qdrop_P_increasing_same_window": 0,
        "qdrop_R_decrease_and_P_increase_same_window": 0,
        "qdrop_witnesses": [],
        "wall_s": 0.0,
    }
    per_n: dict[str, dict[str, int]] = {}

    t0 = time.time()
    for n in range(args.min_n, args.max_n + 1):
        n_stats = {
            "seen": 0,
            "considered": 0,
            "checked": 0,
            "qdrop_total": 0,
            "qdrop_shift_nonzero": 0,
        }

        for res in range(args.mod):
            cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
            if args.mod > 1:
                cmd.append(f"{res}/{args.mod}")

            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            assert proc.stdout is not None
            for raw in proc.stdout:
                n_stats["seen"] += 1
                summary["seen"] += 1

                nn, adj = parse_graph6(raw)
                if (not args.all_trees) and (not is_dleaf_le_1(nn, adj)):
                    continue

                n_stats["considered"] += 1
                summary["considered"] += 1

                poly_t = independence_poly(nn, adj)
                m = mode_index_leftmost(poly_t)
                if m == 0 or m >= len(poly_t) or poly_t[m] == 0:
                    summary["skip_mode"] += 1
                    continue

                deg = [len(nb) for nb in adj]
                leaves = [v for v in range(nn) if deg[v] == 1]
                min_parent_deg = min(deg[adj[l][0]] for l in leaves)
                leaf = min(l for l in leaves if deg[adj[l][0]] == min_parent_deg)
                support = adj[leaf][0]
                if deg[support] != 2:
                    summary["skip_non_deg2"] += 1
                    continue

                n_stats["checked"] += 1
                summary["checked"] += 1

                u = [x for x in adj[support] if x != leaf][0]
                b_adj = remove_vertices(adj, {leaf, support})
                b_poly = independence_poly(len(b_adj), b_adj)
                if m - 2 < 0 or m >= len(b_poly) or (m - 1) >= len(b_poly) or b_poly[m - 1] == 0:
                    summary["skip_short_B"] += 1
                    continue

                keep = [v for v in range(nn) if v not in {leaf, support}]
                idx_map = {v: i for i, v in enumerate(keep)}
                u_in_b = idx_map[u]
                p_poly, q_poly = compute_hub_polys(b_adj, u_in_b)

                # Independent reconstruction of R from child dp0 factors at u.
                children, dp0, dp1 = rooted_dp(b_adj, u_in_b)
                child_F = [dp0[c] for c in children[u_in_b]]
                child_I = [_polyadd(dp0[c], dp1[c]) for c in children[u_in_b]]
                r_poly = prod_all(child_F)
                p_rebuilt = prod_all(child_I)

                if p_rebuilt != p_poly:
                    raise RuntimeError("Internal inconsistency: rebuilt P != compute_hub_polys P")

                q0 = coeff(q_poly, m - 2)
                q1 = coeff(q_poly, m - 1)
                if q1 >= q0:
                    continue

                n_stats["qdrop_total"] += 1
                summary["qdrop_total"] += 1

                g6 = raw.decode("ascii").strip()
                mode_b = mode_index_leftmost(b_poly)
                payload = make_witness_payload(
                    n=nn,
                    g6=g6,
                    m=m,
                    mode_b=mode_b,
                    leaf=leaf,
                    support=support,
                    u=u,
                    p_poly=p_poly,
                    q_poly=q_poly,
                    b_poly=b_poly,
                    r_poly=r_poly,
                )
                summary["qdrop_witnesses"].append(payload)

                if payload["shift"] == 0:
                    summary["qdrop_shift0"] += 1
                elif payload["shift"] == 1:
                    summary["qdrop_shift1"] += 1
                    n_stats["qdrop_shift_nonzero"] += 1
                    summary["qdrop_shift_nonzero"] += 1
                else:
                    summary["qdrop_shift_other"] += 1
                    n_stats["qdrop_shift_nonzero"] += 1
                    summary["qdrop_shift_nonzero"] += 1

                if not payload["q_index_identity_ok"]:
                    summary["qdrop_q_index_identity_fail"] += 1

                if payload["delta_R_same_window"] < 0:
                    summary["qdrop_R_decreasing_same_window"] += 1
                if payload["delta_P_same_window"] > 0:
                    summary["qdrop_P_increasing_same_window"] += 1
                if payload["delta_R_same_window"] < 0 and payload["delta_P_same_window"] > 0:
                    summary["qdrop_R_decrease_and_P_increase_same_window"] += 1

            proc.wait()

        per_n[str(n)] = n_stats
        print(
            f"n={n:2d}: seen={n_stats['seen']:9d} considered={n_stats['considered']:8d} "
            f"checked={n_stats['checked']:8d} qdrop={n_stats['qdrop_total']:3d} "
            f"qdrop_shift_nonzero={n_stats['qdrop_shift_nonzero']:2d}",
            flush=True,
        )

    summary["wall_s"] = time.time() - t0

    out_payload = {
        "params": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "mod": args.mod,
            "all_trees": args.all_trees,
        },
        "summary": summary,
        "per_n": per_n,
    }
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(out_payload, f, indent=2)

    print("-" * 96)
    print(
        f"TOTAL qdrop={summary['qdrop_total']} shift0={summary['qdrop_shift0']} "
        f"shift1={summary['qdrop_shift1']} shift_other={summary['qdrop_shift_other']} "
        f"shift_nonzero={summary['qdrop_shift_nonzero']} q_index_fail={summary['qdrop_q_index_identity_fail']}",
        flush=True,
    )
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
