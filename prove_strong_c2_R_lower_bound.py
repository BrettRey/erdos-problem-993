#!/usr/bin/env python3
"""Investigate the lower bound on R = cross + mismatch.

R = 2*p1*q1 - pm*q0 - p0*qm + p0*q1 - p1*q0

Observations:
- R >= 4 (min at n=5)
- cross >= 3 (min at n=5)
- mismatch ranges from -87808 to large positive

Try to find what determines R. Candidates:
1. R >= some function of deg(u)
2. R >= some function of mode difference
3. R / b1 or R / (p1*q1) has a nice lower bound
4. R >= cross - |mismatch|, and cross dominates

Also check: does R = (p1-p0)*(q1-q0) + something nice?
(p1-p0)*(q1-q0) = p1*q1 - p1*q0 - p0*q1 + p0*q0

R = 2*p1*q1 - pm*q0 - p0*qm + p0*q1 - p1*q0
  = (p1*q1 - pm*q0) + (p1*q1 - p0*qm) + (p0*q1 - p1*q0)
  ... the first two are cross-mode terms

Let me check: R = (p1-p0)*(q1-q0) + p0*q0 + p1*q1 - pm*q0 - p0*qm?
(p1-p0)*(q1-q0) = p1*q1 - p1*q0 - p0*q1 + p0*q0
R - (p1-p0)*(q1-q0) = 2*p1*q1 - pm*q0 - p0*qm + p0*q1 - p1*q0
                      - p1*q1 + p1*q0 + p0*q1 - p0*q0
                     = p1*q1 - pm*q0 - p0*qm + 2*p0*q1 - p0*q0

Hmm not clean. Let me just compute statistics.
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
    min_R = None
    min_R_tree = None
    min_cross = None

    # Check various candidate lower bounds
    min_R_over_q1 = None
    min_R_over_p1 = None
    min_cross_over_q1 = None

    # Check a specific identity: R = dp*dq + ?
    # where dp = p1-p0, dq = q1-q0
    identity_check_fail = 0

    # Check: cross = p1*q1 + (p1*q1 - pm*q0 - p0*qm + p0*q0) - p0*q0 + p1*q1
    # Hmm, just check: does cross = lc_P_shifted + lc_Q_shifted + something?

    # New idea: factor R through products
    # P = prod f_c, Q = x * prod g_c where f_c = I(T_c), g_c = dp0[c]
    # For a single factor: R_1 = cross_1 + mismatch_1 where cross_1 = 2*f1*g0 - f2*g_{-1} - f0*g1
    # and g is just dp0[c] shifted... not sure this simplifies.

    # Check: R >= q0 always?
    R_ge_q0 = 0
    R_ge_p0 = 0
    R_ge_2 = 0

    # Small-R witnesses
    small_R_witnesses = []

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

            cross = 2 * p1 * q1 - pm * q0 - p0 * qm
            mismatch = p0 * q1 - p1 * q0
            R = cross + mismatch

            if min_R is None or R < min_R:
                min_R = R
                min_R_tree = {
                    "n": nn, "m": m,
                    "p": [p0, p1, pm], "q": [q0, q1, qm],
                    "R": R, "cross": cross, "mismatch": mismatch,
                    "g6": raw.decode("ascii").strip(),
                    "deg_u": deg[u],
                }

            if min_cross is None or cross < min_cross:
                min_cross = cross

            if q1 > 0:
                ratio = R / q1
                if min_R_over_q1 is None or ratio < min_R_over_q1:
                    min_R_over_q1 = ratio
            if p1 > 0:
                ratio = R / p1
                if min_R_over_p1 is None or ratio < min_R_over_p1:
                    min_R_over_p1 = ratio

            if R >= q0:
                R_ge_q0 += 1
            if R >= p0:
                R_ge_p0 += 1
            if R >= 2:
                R_ge_2 += 1

            if R <= 20:
                small_R_witnesses.append({
                    "n": nn, "m": m,
                    "p": [p0, p1, pm], "q": [q0, q1, qm],
                    "R": R, "cross": cross, "mismatch": mismatch,
                    "deg_u": deg[u],
                    "g6": raw.decode("ascii").strip(),
                })

        proc.wait()
        print(f"n={n:2d}: total={total:8d} min_R={min_R}", flush=True)

    print(f"\nTotal: {total}")
    print(f"min(R) = {min_R}")
    print(f"min(cross) = {min_cross}")
    print(f"min(R/q1) = {min_R_over_q1:.6f}" if min_R_over_q1 is not None else "")
    print(f"min(R/p1) = {min_R_over_p1:.6f}" if min_R_over_p1 is not None else "")
    print(f"R >= q0: {R_ge_q0} / {total}")
    print(f"R >= p0: {R_ge_p0} / {total}")
    print(f"R >= 2: {R_ge_2} / {total}")
    print(f"\nMin R witness: {min_R_tree}")

    print(f"\nSmall-R witnesses (R <= 20):")
    for w in sorted(small_R_witnesses, key=lambda x: x["R"])[:20]:
        print(f"  n={w['n']}, m={w['m']}, R={w['R']}, cross={w['cross']}, "
              f"mismatch={w['mismatch']}, p={w['p']}, q={w['q']}, deg_u={w['deg_u']}")

    print(f"\nTime: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
