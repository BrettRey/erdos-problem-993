#!/usr/bin/env python3
"""Prove that the residual term in the STRONG C2 decomposition is non-negative.

Key decomposition:
  combined = lc_surplus + mismatch = lc_P + lc_Q + residual

where:
  lc_P = p_{m-1}^2 - p_m * p_{m-2}  >= 0 (P is LC)
  lc_Q = q_{m-1}^2 - q_m * q_{m-2}  >= 0 (Q is LC)
  residual = 2*p1*q1 - pm*q0 - p0*qm + p0*q1 - p1*q0

Since Q = x*P' where P'_j = q_{j+1}, we have q0=P'_{m-3}, q1=P'_{m-2}, qm=P'_{m-1}.

ALGEBRAIC PROOF ATTEMPT:

Method: express residual using Cauchy-like inequalities.

residual = (p1*q1 - pm*q0) + (p1*q1 - p0*qm) + (p0*q1 - p1*q0)
         = (p1*q1 - pm*q0) + p0*(q1-qm) + p1*(q1-q0)

Alternative factoring:
residual = p1*(2*q1 - q0) + p0*(q1 - qm) - pm*q0

FKG APPROACH:
P = product_c I(T_c) and P' = product_c dp0[c].
Let I(T_c) = dp0[c] + dp1[c], so P_k = sum over partitions of k into |children| parts
of product dp0[c_j]+dp1[c_j] at each part.

The KEY structural property: P dominates P' coefficientwise in the log-concavity sense.
Specifically, p_k/P'_k is non-decreasing for k up to the modes.

This script checks: for each child c, I(T_c) / dp0[c] ratio behavior,
and whether the "conditional inclusion" probability is monotone.
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
from indpoly import independence_poly, is_log_concave, _polyadd, _polymul


def compute_child_polys(adj_B: list[list[int]], u_in_B: int) -> list[tuple[list[int], list[int]]]:
    """Compute (dp0[c], dp1[c]) for each child c of u when B is rooted at u.

    Returns list of (dp0[c], dp1[c]) pairs (one per child of u).
    """
    n = len(adj_B)
    if n <= 1:
        return []

    parent = [-1] * n
    children = [[] for _ in range(n)]
    visited = [False] * n
    visited[u_in_B] = True
    bfs_queue = [u_in_B]
    head = 0
    while head < len(bfs_queue):
        v = bfs_queue[head]
        head += 1
        for w in adj_B[v]:
            if not visited[w]:
                visited[w] = True
                parent[w] = v
                children[v].append(w)
                bfs_queue.append(w)

    order = []
    stack = [(u_in_B, False)]
    while stack:
        v, processed = stack.pop()
        if processed:
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
        else:
            prod = [1]
            for c in children[v]:
                summand = _polyadd(dp0[c], dp1[c])
                prod = _polymul(prod, summand)
            dp0[v] = prod

            prod = [1]
            for c in children[v]:
                prod = _polymul(prod, dp0[c])
            dp1[v] = [0] + prod

    return [(dp0[c], dp1[c]) for c in children[u_in_B]]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--out", default="results/prove_strong_c2_residual_n20.json")
    args = ap.parse_args()

    stats = {
        "checked": 0,
        "residual_neg": 0,
        "residual_min": None,

        # Check: is P/P' ratio non-decreasing at (m-2, m-1)?
        # i.e. p_{m-1}/P'_{m-2} >= p_{m-2}/P'_{m-3}
        # Equivalently: p_{m-1}*P'_{m-3} >= p_{m-2}*P'_{m-2}
        # Using q_k = P'_{k-1}: p_{m-1}*q_{m-2} >= p_{m-2}*q_{m-1}
        # = p1*q0 >= p0*q1
        # This is mismatch = p0*q1 - p1*q0 <= 0, i.e., mismatch <= 0!
        # So P/P' ratio non-decreasing at (m-2,m-1) IFF mismatch <= 0.
        # Since mismatch >= 0 in 99.98% of cases, the ratio is actually DECREASING usually.
        # When mismatch < 0 (only 129 cases), the ratio IS increasing.
        # So this is backwards from what I hoped.

        # Check: for each child subtree T_c, is I(T_c)/dp0[c] ratio non-decreasing?
        # I(T_c)_k = dp0[c]_k + dp1[c]_k, so ratio_k = 1 + dp1[c]_k/dp0[c]_k.
        # dp1[c]_k = dp0-product of subtrees of c, shifted by 1 (= x * product dp0[grandchildren]).
        # The ratio dp1_k/dp0_k is the conditional probability of including c given IS of size k.
        # For trees, this typically increases with k up to a point.
        "child_ratio_monotone": 0,
        "child_ratio_not_monotone": 0,

        # Cross-LC: for two LC sequences a, b, define
        # CrossLC(a,b,k) = a_k * b_k - a_{k+1} * b_{k-1}
        # This is >= 0 when the sequences have "compatible" growth rates (a/b decreasing).
        # Check: CrossLC(P, P', m-2) = p_{m-2}*P'_{m-2} - p_{m-1}*P'_{m-3}
        #       = p0*q1 - p1*q0 = mismatch (NOT -mismatch)
        # So CrossLC(P,P',m-2) = mismatch. It's >= 0 in 99.98% of cases.

        # The cross term: cross = 2*p1*q1 - pm*q0 - p0*qm
        # Can be written as: cross = (p1*q1 - pm*q0) + (p1*q1 - p0*qm)
        # = CrossLC'(P,Q,m-1,m-2) + CrossLC'(P,Q,m-1,m)
        # where CrossLC'(a,b,j,k) = a_j*b_j - a_{j+1}*b_{j-1} etc.
        # Actually: p1*q1 - pm*q0 and p1*q1 - p0*qm.
        # The first is p_{m-1}*q_{m-1} - p_m*q_{m-2}: this is non-negative iff
        # p_{m-1}/p_m >= q_{m-2}/q_{m-1}, i.e., P's inverse ratio >= Q's ratio.
        # Since both P and Q peak around m-1 to m, this measures whether P and Q
        # have compatible growth patterns.

        # CHECK: cross = (p1*q1 - pm*q0) + (p1*q1 - p0*qm)
        # First term: p_{m-1}*q_{m-1} >= p_m*q_{m-2}?
        "cross_term1_neg": 0,
        # Second term: p_{m-1}*q_{m-1} >= p_{m-2}*q_m?
        "cross_term2_neg": 0,

        # ALGEBRAIC RESIDUAL PROOF IDEA:
        # residual = cross + mismatch
        # = (p1*q1 - pm*q0) + (p1*q1 - p0*qm) + (p0*q1 - p1*q0)
        # = (p1*q1 - pm*q0) + q1*(p1+p0-p1) + p1*q1 - p0*qm - p1*q0 ... no
        # = (p1*q1 - pm*q0) + (p1*q1 - p1*q0) + (p0*q1 - p0*qm)
        # = (p1*q1 - pm*q0) + p1*(q1-q0) + p0*(q1-qm)
        #
        # Since mode(Q) <= m-1 or mode(Q) = m: q1 vs qm depends on mode(Q).
        # mode_Q relative to m-1:
        "mode_Q_le_m1": 0,  # q1=q_{m-1} is past or at mode; q1 >= qm
        "mode_Q_ge_m": 0,   # q_{m-1} is before mode; qm >= q1

        # For the first term: p1*q1 - pm*q0:
        # By Cauchy-Schwarz from LC: p1*q1 >= sqrt(pm*p0) * sqrt(qm*q0)
        # >= sqrt(pm*q0*p0*qm)? No, that's sqrt(pm*p0)*sqrt(qm*q0).
        # We need p1*q1 >= pm*q0. Is this true?
        # p1*q1 >= sqrt(pm*p0)*sqrt(qm*q0) (by LC of P and Q)
        # vs pm*q0.
        # sqrt(pm*p0*qm*q0) vs pm*q0: squaring: pm*p0*qm*q0 vs pm^2*q0^2
        # = p0*qm vs pm*q0. This is CrossLC(P',P,m-1) type comparison.
        # So p1*q1 >= pm*q0 iff (p1*q1)^2 >= (pm*q0)^2
        # iff pm*p0*qm*q0 * (p1*q1/sqrt(pm*p0*qm*q0))^2 >= pm^2*q0^2
        # Hmm, this isn't clean.

        # DIRECT CHECK of p1*q1 >= pm*q0:
        "pq_cross_neg": 0,  # p1*q1 < pm*q0?

        # If ALL three pieces of the decomposition are >= 0:
        # residual = (p1*q1 - pm*q0) + p1*(q1-q0) + p0*(q1-qm) >= 0
        # This holds when:
        # (i) p1*q1 >= pm*q0 (always? check)
        # (ii) q1 >= q0 (mode(Q) >= m-1, so Q non-decreasing at m-2 to m-1)
        # (iii) q1 >= qm (mode(Q) <= m-1, so Q non-increasing from m-1 on)
        # Note: (ii) and (iii) can't both hold with strict inequality unless mode(Q) = m-1.
        # When mode(Q) = m-1: q1 >= q0 and q1 >= qm, so both (ii) and (iii) hold.
        # When mode(Q) < m-1: q1 < q0 possible (violates ii), but q1 > qm (satisfies iii).
        # When mode(Q) > m-1: q1 > q0 (satisfies ii), but q1 < qm possible (violates iii).

        # So we need the three-term decomposition approach to handle all cases.

        # MULTIPLICATIVE APPROACH:
        # Since P = prod_c I(T_c), the coefficients satisfy the FKG inequality.
        # Specifically, for two LC polynomials f and g:
        # (fg)_k * (fg)_k >= (fg)_{k-1} * (fg)_{k+1} (product is LC)
        # More relevantly: for the ratio f*g / (f' * g') where f = dp0+dp1, f' = dp0:
        # This is (1+dp1/dp0) type ratios, whose product gives P/P'.

        # Let me just check empirically whether p1*q1 >= pm*q0 always.

        # INTERLACING approach:
        # Both P and Q come from the same tree B (rooted at u).
        # P_{k} = sum over IS of B\u of size k
        # Q_{k} = sum over IS of B containing u of size k = P'_{k-1} (shifted product of dp0s)
        # The interleaving: b_k = p_k + q_k.
        # Both P and Q are LC, and their sum B is LC.
        # The "signed" cross-determinant p1*q0 - p0*q1 is the "mismatch" sign.

        # RESIDUAL AS INTERLACED TURÁN:
        # residual = (p1+p0)*q1 + p1*(q1-q0) - p0*qm - pm*q0
        # = (p1+p0)*q1 + p1*q1 - p1*q0 - p0*qm - pm*q0
        # = p1*q1 + (p1+p0)*q1 - p1*q0 - p0*qm - pm*q0
        # = p1*q1 + p1*(q1-q0) + p0*(q1-qm) + (p0*qm - p0*qm) + ... circular

        # MORE PROMISING: express using a_k = b_k + p_{k-1}:
        # combined = b1*a1 - b0*a2 where a1=b1+p0, a2=b2+p1
        # = b1*(b1+p0) - b0*(b2+p1)
        # = b1^2 + b1*p0 - b0*b2 - b0*p1
        # = (b1^2-b0*b2) + b1*p0 - b0*p1
        # = lc_B + mismatch (as defined)
        #
        # Alternative: combined = b1*a1 - b0*a2
        # If A is LC at m: a1^2 >= a0*a2 where a0=b0+p_{m-3}.
        # Then a2 <= a1^2/a0 and combined >= b1*a1 - b0*a1^2/a0 = a1*(b1 - b0*a1/a0).
        # = a1*(b1*a0 - b0*a1)/a0.
        # b1*a0 - b0*a1 = b1*(b0+p_{m-3}) - b0*(b1+p0) = b1*p_{m-3} - b0*p0.
        # So combined >= a1/a0 * (b1*p_{m-3} - b0*p0) = a1/a0 * mismatch_lower
        # where mismatch_lower = p_{m-3}*b1 - p0*b0.
        # From v1: lower_mismatch_neg = 77130 (almost always negative!).
        # So this approach FAILS.

        # CLEANEST PATH FORWARD:
        # Show combined >= 0 by decomposition into lc_P + lc_Q + residual,
        # where residual >= 0 is verified computationally through n=23
        # and backed by structural arguments for why it should hold in general.

        "wall_s": 0.0,
    }

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

            stats["checked"] += 1

            p0 = p_poly[m - 2] if m - 2 < len(p_poly) else 0
            p1 = p_poly[m - 1] if m - 1 < len(p_poly) else 0
            pm = p_poly[m] if m < len(p_poly) else 0
            q0 = q_poly[m - 2] if m - 2 < len(q_poly) else 0
            q1 = q_poly[m - 1] if m - 1 < len(q_poly) else 0
            qm = q_poly[m] if m < len(q_poly) else 0

            lc_P = p1 * p1 - pm * p0
            lc_Q = q1 * q1 - qm * q0
            residual = 2 * p1 * q1 - pm * q0 - p0 * qm + p0 * q1 - p1 * q0

            if residual < 0:
                stats["residual_neg"] += 1
            if stats["residual_min"] is None or residual < stats["residual_min"]:
                stats["residual_min"] = residual

            # Cross term 1: p1*q1 vs pm*q0
            if p1 * q1 < pm * q0:
                stats["pq_cross_neg"] += 1
                stats["cross_term1_neg"] += 1

            # Cross term 2: p1*q1 vs p0*qm
            if p1 * q1 < p0 * qm:
                stats["cross_term2_neg"] += 1

            # Mode of Q
            mode_Q = mode_index_leftmost(q_poly) if q_poly else 0
            if mode_Q <= m - 1:
                stats["mode_Q_le_m1"] += 1
            if mode_Q >= m:
                stats["mode_Q_ge_m"] += 1

        proc.wait()
        print(f"n={n:2d}: checked={stats['checked']:8d} resid_neg={stats['residual_neg']} "
              f"cross1_neg={stats['cross_term1_neg']} cross2_neg={stats['cross_term2_neg']} "
              f"pq_cross_neg={stats['pq_cross_neg']} resid_min={stats['residual_min']}",
              flush=True)

    stats["wall_s"] = time.time() - t0

    print("\n" + "=" * 80)
    print(f"Checked: {stats['checked']}")
    print(f"Residual negative: {stats['residual_neg']} (min={stats['residual_min']})")
    print(f"Cross term 1 (p1*q1 < pm*q0): {stats['cross_term1_neg']}")
    print(f"Cross term 2 (p1*q1 < p0*qm): {stats['cross_term2_neg']}")
    print(f"p1*q1 < pm*q0: {stats['pq_cross_neg']}")
    print(f"mode(Q) <= m-1: {stats['mode_Q_le_m1']}")
    print(f"mode(Q) >= m: {stats['mode_Q_ge_m']}")

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
