#!/usr/bin/env python3
"""Extract and display the 2 Q-drop witness trees and the tightest mismatch-neg cases.

For each witness, show full structural data:
- Tree structure, P, Q polynomials
- All indices' p_k, q_k, b_k values
- Mode structure, LC checks
- Rise-compensation quantities
"""

from __future__ import annotations

import argparse
import subprocess
import time

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import compute_hub_polys, mode_index_leftmost, remove_vertices
from graph6 import parse_graph6
from indpoly import independence_poly, is_log_concave


def analyze_tree(nn, adj, raw_g6):
    """Full analysis of a single tree."""
    poly_t = independence_poly(nn, adj)
    m = mode_index_leftmost(poly_t)

    deg = [len(nb) for nb in adj]
    leaves = [v for v in range(nn) if deg[v] == 1]
    min_parent_deg = min(deg[adj[l][0]] for l in leaves)
    leaf = min(l for l in leaves if deg[adj[l][0]] == min_parent_deg)
    support = adj[leaf][0]

    if deg[support] != 2:
        return None

    u = [x for x in adj[support] if x != leaf][0]
    b_adj = remove_vertices(adj, {leaf, support})
    b_poly = independence_poly(len(b_adj), b_adj)

    keep = [v for v in range(nn) if v not in {leaf, support}]
    idx_map = {v: i for i, v in enumerate(keep)}
    u_in_b = idx_map[u]
    p_poly, q_poly = compute_hub_polys(b_adj, u_in_b)

    # Coefficients at mode
    p0 = p_poly[m-2] if m-2 >= 0 and m-2 < len(p_poly) else 0
    p1 = p_poly[m-1] if m-1 < len(p_poly) else 0
    pm = p_poly[m] if m < len(p_poly) else 0

    q0 = q_poly[m-2] if m-2 >= 0 and m-2 < len(q_poly) else 0
    q1 = q_poly[m-1] if m-1 < len(q_poly) else 0
    qm = q_poly[m] if m < len(q_poly) else 0

    b0 = p0 + q0
    b1 = p1 + q1
    b2 = pm + qm

    mismatch = p0 * b1 - p1 * b0
    neg = -mismatch
    lc_surplus = b1*b1 - b2*b0
    combined = lc_surplus + mismatch
    rise = b1 * (b1 - b0)
    db = b1 - b0
    dq = q1 - q0
    rise_minus_neg = p1 * db + b1 * dq

    lc_P = p1*p1 - pm*p0
    lc_Q = q1*q1 - qm*q0
    cross = 2*p1*q1 - pm*q0 - p0*qm
    R = cross + mismatch

    mode_P = mode_index_leftmost(p_poly)
    mode_Q = mode_index_leftmost(q_poly)
    mode_B = mode_index_leftmost(b_poly)

    return {
        "g6": raw_g6,
        "n": nn,
        "m": m,
        "deg_u": deg[u],
        "I(T)": list(poly_t),
        "I(B)": list(b_poly),
        "P": list(p_poly),
        "Q": list(q_poly),
        "p0": p0, "p1": p1, "pm": pm,
        "q0": q0, "q1": q1, "qm": qm,
        "b0": b0, "b1": b1, "b2": b2,
        "mode_P": mode_P, "mode_Q": mode_Q, "mode_B": mode_B,
        "P_is_LC": is_log_concave(p_poly),
        "Q_is_LC": is_log_concave(q_poly),
        "mismatch": mismatch,
        "neg": neg,
        "lc_surplus": lc_surplus,
        "combined": combined,
        "rise_minus_neg": rise_minus_neg,
        "db": db, "dq": dq,
        "lc_P": lc_P, "lc_Q": lc_Q,
        "cross": cross, "R": R,
        "p1_minus_q1": p1 - q1,
        "p1_over_b1": p1/b1 if b1 > 0 else None,
    }


def print_analysis(info):
    """Pretty-print analysis."""
    print(f"  g6: {info['g6']}")
    print(f"  n={info['n']}, m={info['m']}, deg(u)={info['deg_u']}")
    print(f"  I(T) = {info['I(T)']}")
    print(f"  I(B) = {info['I(B)']}")
    print(f"  P    = {info['P']}")
    print(f"  Q    = {info['Q']}")
    print(f"  mode(P)={info['mode_P']}, mode(Q)={info['mode_Q']}, mode(B)={info['mode_B']}")
    print(f"  P is LC: {info['P_is_LC']}, Q is LC: {info['Q_is_LC']}")
    print(f"  At mode m={info['m']}:")
    print(f"    p0={info['p0']}, p1={info['p1']}, pm={info['pm']}")
    print(f"    q0={info['q0']}, q1={info['q1']}, qm={info['qm']}")
    print(f"    b0={info['b0']}, b1={info['b1']}, b2={info['b2']}")
    print(f"  p1 - q1 = {info['p1_minus_q1']}")
    print(f"  p1/b1 = {info['p1_over_b1']:.6f}" if info['p1_over_b1'] else "  p1/b1 = N/A")
    print(f"  mismatch = {info['mismatch']}")
    print(f"  lc_surplus = {info['lc_surplus']}")
    print(f"  combined = {info['combined']}")
    print(f"  rise-neg = {info['rise_minus_neg']}")
    print(f"  db={info['db']}, dq={info['dq']}")
    print(f"  lc_P={info['lc_P']}, lc_Q={info['lc_Q']}, cross={info['cross']}, R={info['R']}")
    print()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    args = ap.parse_args()

    q_drop_witnesses = []
    mismatch_neg_witnesses = []
    tightest_combined = []

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
            b_poly = independence_poly(len(b_adj), b_adj)
            if m - 2 < 0 or m >= len(b_poly) or m - 1 >= len(b_poly) or b_poly[m - 1] == 0:
                continue

            keep = [v for v in range(nn) if v not in {leaf, support}]
            idx_map = {v: i for i, v in enumerate(keep)}
            u_in_b = idx_map[u]
            p_poly, q_poly = compute_hub_polys(b_adj, u_in_b)

            q0 = q_poly[m-2] if m-2 >= 0 and m-2 < len(q_poly) else 0
            q1 = q_poly[m-1] if m-1 < len(q_poly) else 0
            p0 = p_poly[m-2] if m-2 >= 0 and m-2 < len(p_poly) else 0
            p1 = p_poly[m-1] if m-1 < len(p_poly) else 0

            mismatch = p0 * (p1+q1) - p1 * (p0+q0)
            combined = (p1+q1)**2 - (p_poly[m] if m < len(p_poly) else 0 + q_poly[m] if m < len(q_poly) else 0) * (p0+q0) + mismatch

            g6 = raw.decode("ascii").strip()

            if q1 < q0:
                info = analyze_tree(nn, adj, g6)
                if info:
                    q_drop_witnesses.append(info)

            if mismatch < 0:
                info = analyze_tree(nn, adj, g6)
                if info:
                    mismatch_neg_witnesses.append(info)

        proc.wait()
        print(f"n={n:2d}: q_drop={len(q_drop_witnesses)}, mismatch_neg={len(mismatch_neg_witnesses)}", flush=True)

    print(f"\n{'='*60}")
    print(f"Q-DROP WITNESSES ({len(q_drop_witnesses)}):")
    print(f"{'='*60}")
    for info in q_drop_witnesses:
        print_analysis(info)

    print(f"\n{'='*60}")
    print(f"TIGHTEST MISMATCH-NEG WITNESSES (lowest combined):")
    print(f"{'='*60}")
    mismatch_neg_witnesses.sort(key=lambda x: x["combined"])
    for info in mismatch_neg_witnesses[:10]:
        print_analysis(info)

    print(f"Time: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
