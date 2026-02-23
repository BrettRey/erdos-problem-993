#!/usr/bin/env python3
"""Analyze the shift-1 regime in detail.

In shift-1 (b2 > b1, i.e., mode(B) = m):
  combined = (rise-neg) + b0*(b1-b2)
  = (rise-neg) - b0*(b2-b1)

So we need rise-neg >= b0*(b2-b1).

Check:
- How large is b0*(b2-b1) relative to rise-neg?
- What fraction of rise-neg is consumed by the shift-1 loss?
- Is the LC decomposition route (lc_P + lc_Q + R) better here?
"""

from __future__ import annotations

import argparse
import subprocess
import time

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import compute_hub_polys, mode_index_leftmost, remove_vertices
from graph6 import parse_graph6
from indpoly import independence_poly


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    args = ap.parse_args()

    total = 0
    shift1_count = 0
    shift1_combined_neg = 0
    max_loss_ratio = 0.0  # max b0*(b2-b1) / rise-neg among shift-1
    min_margin = None
    min_margin_witness = None

    # For shift-1: what's the min(lc_P + lc_Q) vs max(-R)?
    # (R is always >= 4, so lc_P + lc_Q + R >= lc_P + lc_Q + 4.)
    shift1_min_combined = None
    shift1_min_lc_P = None
    shift1_min_lc_Q = None
    shift1_min_R = None

    t0 = time.time()

    for n in range(args.min_n, args.max_n + 1):
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
            p_poly, q_poly = compute_hub_polys(b_adj, u_in_b)

            total += 1

            p0 = p_poly[m - 2] if m - 2 < len(p_poly) else 0
            p1 = p_poly[m - 1] if m - 1 < len(p_poly) else 0
            pm = p_poly[m] if m < len(p_poly) else 0
            q0 = q_poly[m - 2] if m - 2 >= 0 and m - 2 < len(q_poly) else 0
            q1 = q_poly[m - 1] if m - 1 < len(q_poly) else 0
            qm = q_poly[m] if m < len(q_poly) else 0

            b0 = p0 + q0
            b1 = p1 + q1
            b2 = pm + qm

            if b1 >= b2:  # shift-0
                continue

            shift1_count += 1

            rise_neg = p1 * (b1 - b0) + b1 * (q1 - q0)
            loss = b0 * (b2 - b1)
            combined = rise_neg - loss

            lc_P = p1 * p1 - pm * p0
            lc_Q = q1 * q1 - qm * q0
            R = 2 * p1 * q1 - pm * q0 - p0 * qm + p0 * q1 - p1 * q0
            combined2 = lc_P + lc_Q + R

            if combined < 0:
                shift1_combined_neg += 1

            if rise_neg > 0:
                ratio = loss / rise_neg
                if ratio > max_loss_ratio:
                    max_loss_ratio = ratio

            if min_margin is None or combined < min_margin:
                min_margin = combined
                min_margin_witness = {
                    "n": nn, "m": m,
                    "p": [p0, p1, pm], "q": [q0, q1, qm],
                    "b": [b0, b1, b2],
                    "rise_neg": rise_neg, "loss": loss,
                    "combined": combined,
                    "lc_P": lc_P, "lc_Q": lc_Q, "R": R,
                    "combined2": combined2,
                    "g6": raw.decode("ascii").strip(),
                    "deg_u": deg[u],
                }

            if shift1_min_combined is None or combined < shift1_min_combined:
                shift1_min_combined = combined
            if shift1_min_lc_P is None or lc_P < shift1_min_lc_P:
                shift1_min_lc_P = lc_P
            if shift1_min_lc_Q is None or lc_Q < shift1_min_lc_Q:
                shift1_min_lc_Q = lc_Q
            if shift1_min_R is None or R < shift1_min_R:
                shift1_min_R = R

        proc.wait()
        print(f"n={n:2d}: shift1={shift1_count:8d} max_loss_ratio={max_loss_ratio:.6f}", flush=True)

    print(f"\nTotal: {total}, Shift-1: {shift1_count}")
    print(f"Shift-1 combined neg: {shift1_combined_neg}")
    print(f"Max loss/rise-neg ratio: {max_loss_ratio:.6f}")
    print(f"Shift-1 min combined: {shift1_min_combined}")
    print(f"Shift-1 min lc_P: {shift1_min_lc_P}")
    print(f"Shift-1 min lc_Q: {shift1_min_lc_Q}")
    print(f"Shift-1 min R: {shift1_min_R}")
    print(f"\nTightest shift-1 witness:")
    if min_margin_witness:
        for k, v in min_margin_witness.items():
            print(f"  {k}: {v}")
    print(f"\nTime: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
