#!/usr/bin/env python3
"""Unified verification of STRONG C2 via the P-dominance + LC decomposition route.

Combined proof strategy:
1. lc_surplus + mismatch = lc_P + lc_Q + R  (exact identity)
2. lc_P >= 0 (P is LC, product of LC subtree IS polynomials)
3. lc_Q >= 0 (Q is LC, x times product of LC dp0 polynomials)
4. R = cross + mismatch = (2*p1*q1 - pm*q0 - p0*qm) + (p0*q1 - p1*q0)

For R >= 0, decompose into:
  R = (p1*q1 - pm*q0) + (p1*q1 - p0*qm) + (p0*q1 - p1*q0)
    = A + B' + C'  [different grouping from original A+B+C]

Or use rise-compensation: combined = (rise - neg) + b0*(b1-b2)
  where rise - neg = p1*db + b1*dq >= 0 when dq >= 0 (Regime E, 99.998% of cases).

For the 2 Q-drop cases (dq < 0):
  Need (-dq/db) <= p1/b1.
  Since p1/b1 >= 1/2 always (computationally verified: p1 >= q1 + 1),
  need (-dq/db) <= 1/2.

This script verifies the full chain on the n<=23 frontier.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import compute_hub_polys, mode_index_leftmost, remove_vertices
from graph6 import parse_graph6
from indpoly import independence_poly


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    # Counters
    total = 0
    combined_neg = 0  # main target: should be 0

    # Route 1: LC decomposition
    lc_P_neg = 0
    lc_Q_neg = 0
    R_neg = 0
    R_min = None

    # Route 2: rise-compensation
    regime_E = 0  # dq >= 0
    regime_H = 0  # dq < 0
    rise_neg_fail = 0  # rise-neg < 0
    hard_ratio_fail = 0  # (-dq/db) > p1/b1 in regime H

    # P-dominance
    p_ge_q1 = 0
    min_p_minus_q = None  # min(p1 - q1)

    # Combined
    shift0 = 0  # b1 >= b2
    shift1 = 0  # b1 < b2
    shift1_combined_neg = 0

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

            # Extract coefficients
            p0 = p_poly[m - 2] if m - 2 < len(p_poly) else 0
            p1 = p_poly[m - 1] if m - 1 < len(p_poly) else 0
            pm = p_poly[m] if m < len(p_poly) else 0

            q0 = q_poly[m - 2] if m - 2 >= 0 and m - 2 < len(q_poly) else 0
            q1 = q_poly[m - 1] if m - 1 < len(q_poly) else 0
            qm = q_poly[m] if m < len(q_poly) else 0

            b0 = p0 + q0
            b1 = p1 + q1
            b2 = pm + qm

            # --- Route 1: LC decomposition ---
            lc_P = p1 * p1 - pm * p0
            lc_Q = q1 * q1 - qm * q0
            cross = 2 * p1 * q1 - pm * q0 - p0 * qm
            mismatch = p0 * q1 - p1 * q0
            R = cross + mismatch
            combined = lc_P + lc_Q + R

            if lc_P < 0:
                lc_P_neg += 1
            if lc_Q < 0:
                lc_Q_neg += 1
            if R < 0:
                R_neg += 1
            if R_min is None or R < R_min:
                R_min = R
            if combined < 0:
                combined_neg += 1

            # --- Route 2: rise-compensation ---
            db = b1 - b0
            dq = q1 - q0
            rise_minus_neg = p1 * db + b1 * dq

            if dq >= 0:
                regime_E += 1
            else:
                regime_H += 1
                # Check ratio inequality
                if db > 0:
                    transfer = (-dq) / db
                    need = p1 / b1
                    if transfer > need:
                        hard_ratio_fail += 1

            if rise_minus_neg < 0:
                rise_neg_fail += 1

            # --- P-dominance ---
            if p1 >= q1:
                p_ge_q1 += 1
            gap = p1 - q1
            if min_p_minus_q is None or gap < min_p_minus_q:
                min_p_minus_q = gap

            # --- Shift analysis ---
            if b1 >= b2:
                shift0 += 1
            else:
                shift1 += 1
                if combined < 0:
                    shift1_combined_neg += 1

        proc.wait()
        print(f"n={n:2d}: total={total:8d} combined_neg={combined_neg} R_neg={R_neg} "
              f"rise_fail={rise_neg_fail}", flush=True)

    elapsed = time.time() - t0
    results = {
        "total": total,
        "combined_neg": combined_neg,
        "lc_P_neg": lc_P_neg,
        "lc_Q_neg": lc_Q_neg,
        "R_neg": R_neg,
        "R_min": R_min,
        "regime_E": regime_E,
        "regime_H": regime_H,
        "rise_neg_fail": rise_neg_fail,
        "hard_ratio_fail": hard_ratio_fail,
        "p_ge_q1": p_ge_q1,
        "min_p_minus_q": min_p_minus_q,
        "shift0": shift0,
        "shift1": shift1,
        "shift1_combined_neg": shift1_combined_neg,
        "wall_s": elapsed,
    }

    print(f"\n{'='*60}")
    print(f"STRONG C2 UNIFIED VERIFICATION (n <= {args.max_n})")
    print(f"{'='*60}")
    print(f"Total trees checked: {total}")
    print(f"COMBINED NEGATIVE: {combined_neg}  {'PASS' if combined_neg == 0 else 'FAIL'}")
    print()
    print(f"--- Route 1: LC decomposition ---")
    print(f"lc_P < 0: {lc_P_neg}  (P always LC: {'YES' if lc_P_neg == 0 else 'NO'})")
    print(f"lc_Q < 0: {lc_Q_neg}  (Q always LC: {'YES' if lc_Q_neg == 0 else 'NO'})")
    print(f"R < 0: {R_neg}  (R always non-neg: {'YES' if R_neg == 0 else 'NO'})")
    print(f"R min: {R_min}")
    print()
    print(f"--- Route 2: rise-compensation ---")
    print(f"Regime E (dq >= 0): {regime_E} ({100*regime_E/total:.4f}%)")
    print(f"Regime H (dq < 0):  {regime_H}")
    print(f"rise-neg < 0 (should be 0): {rise_neg_fail}")
    print(f"hard ratio fail: {hard_ratio_fail}")
    print()
    print(f"--- P-dominance at m-1 ---")
    print(f"p1 >= q1: {p_ge_q1} / {total} ({100*p_ge_q1/total:.2f}%)")
    print(f"min(p1 - q1): {min_p_minus_q}")
    print()
    print(f"--- Shift analysis ---")
    print(f"Shift-0 (b1 >= b2): {shift0}")
    print(f"Shift-1 (b1 < b2):  {shift1}")
    print(f"Shift-1 combined neg: {shift1_combined_neg}")
    print()
    print(f"Time: {elapsed:.1f}s")

    if args.out:
        with open(args.out, "w") as f:
            json.dump(results, f, indent=2)
        print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
