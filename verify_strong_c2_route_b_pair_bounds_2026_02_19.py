#!/usr/bin/env python3
"""Fast scan of pairwise Route-B inequalities for R = cross + mismatch.

Tracks exact minima/witnesses for:
  R
  cross
  mismatch
  D1 = p1*q1 - p_m*q0
  D2 = p1*q1 - p0*qm
  D3 = mismatch = p0*q1 - p1*q0
  D1 + D3
  D2 + D3
  cross - |mismatch|
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
from indpoly import independence_poly


def getc(poly: list[int], k: int) -> int:
    return poly[k] if 0 <= k < len(poly) else 0


def update_min(stats: dict[str, Any], key: str, value: int, witness: dict[str, Any]) -> None:
    if stats[key] is None or value < stats[key]:
        stats[key] = value
        stats[key + "_witness"] = witness


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--out", default="results/verify_strong_c2_route_b_pair_bounds_2026_02_19.json")
    args = ap.parse_args()

    stats: dict[str, Any] = {
        "seen": 0,
        "considered": 0,
        "checked": 0,
        "R_neg": 0,
        "cross_neg": 0,
        "D1_neg": 0,
        "D2_neg": 0,
        "D3_neg": 0,
        "D1_plus_D3_neg": 0,
        "D2_plus_D3_neg": 0,
        "cross_minus_abs_mismatch_neg": 0,
        "R_min": None,
        "cross_min": None,
        "mismatch_min": None,
        "D1_min": None,
        "D2_min": None,
        "D3_min": None,
        "D1_plus_D3_min": None,
        "D2_plus_D3_min": None,
        "cross_minus_abs_mismatch_min": None,
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
            P, Q = compute_hub_polys(b_adj, u_in_b)

            r = m - 2
            p0 = getc(P, r)
            p1 = getc(P, r + 1)
            pm = getc(P, r + 2)
            q0 = getc(Q, r)
            q1 = getc(Q, r + 1)
            qm = getc(Q, r + 2)

            D1 = p1 * q1 - pm * q0
            D2 = p1 * q1 - p0 * qm
            D3 = p0 * q1 - p1 * q0
            cross = D1 + D2
            mismatch = D3
            R = cross + mismatch
            d13 = D1 + D3
            d23 = D2 + D3
            cabs = cross - abs(mismatch)

            stats["checked"] += 1
            if R < 0:
                stats["R_neg"] += 1
            if cross < 0:
                stats["cross_neg"] += 1
            if D1 < 0:
                stats["D1_neg"] += 1
            if D2 < 0:
                stats["D2_neg"] += 1
            if D3 < 0:
                stats["D3_neg"] += 1
            if d13 < 0:
                stats["D1_plus_D3_neg"] += 1
            if d23 < 0:
                stats["D2_plus_D3_neg"] += 1
            if cabs < 0:
                stats["cross_minus_abs_mismatch_neg"] += 1

            w = {
                "n": nn,
                "m": m,
                "r": r,
                "g6": raw.decode("ascii").strip(),
                "deg_u_B": len(b_adj[idx_map[u]]),
                "p0": p0,
                "p1": p1,
                "pm": pm,
                "q0": q0,
                "q1": q1,
                "qm": qm,
                "D1": D1,
                "D2": D2,
                "D3": D3,
                "cross": cross,
                "mismatch": mismatch,
                "R": R,
                "D1_plus_D3": d13,
                "D2_plus_D3": d23,
                "cross_minus_abs_mismatch": cabs,
            }

            update_min(stats, "R_min", R, w)
            update_min(stats, "cross_min", cross, w)
            update_min(stats, "mismatch_min", mismatch, w)
            update_min(stats, "D1_min", D1, w)
            update_min(stats, "D2_min", D2, w)
            update_min(stats, "D3_min", D3, w)
            update_min(stats, "D1_plus_D3_min", d13, w)
            update_min(stats, "D2_plus_D3_min", d23, w)
            update_min(stats, "cross_minus_abs_mismatch_min", cabs, w)

        proc.wait()
        print(
            f"n={n:2d}: checked={stats['checked']:8d} "
            f"R_neg={stats['R_neg']} D1_neg={stats['D1_neg']} D2_neg={stats['D2_neg']} "
            f"D3_neg={stats['D3_neg']} d13_neg={stats['D1_plus_D3_neg']} "
            f"d23_neg={stats['D2_plus_D3_neg']} cabs_neg={stats['cross_minus_abs_mismatch_neg']}",
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
    print(f"R_neg: {stats['R_neg']}  R_min: {stats['R_min']}")
    print(f"cross_neg: {stats['cross_neg']}  cross_min: {stats['cross_min']}")
    print(f"D1_neg: {stats['D1_neg']}  D1_min: {stats['D1_min']}")
    print(f"D2_neg: {stats['D2_neg']}  D2_min: {stats['D2_min']}")
    print(f"D3_neg: {stats['D3_neg']}  D3_min: {stats['D3_min']}")
    print(f"D1+D3 neg: {stats['D1_plus_D3_neg']}  min: {stats['D1_plus_D3_min']}")
    print(f"D2+D3 neg: {stats['D2_plus_D3_neg']}  min: {stats['D2_plus_D3_min']}")
    print(
        f"cross-|mismatch| neg: {stats['cross_minus_abs_mismatch_neg']}  "
        f"min: {stats['cross_minus_abs_mismatch_min']}"
    )
    print(f"wall_s: {stats['wall_s']:.1f}")
    print(f"wrote: {args.out}")


if __name__ == "__main__":
    main()
