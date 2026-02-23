#!/usr/bin/env python3
"""Route C structural analysis: prove p_{m-1} >= q_{m-1}.

Mathematical framework
======================
B = T - {l,s}, rooted at u with children c_1, ..., c_d.

For each child c:
  f_c = I(T_c) = dp0[c] + dp1[c]    (full subtree IS poly)
  g_c = dp0[c]                        (subtree IS poly, c excluded)
  h_c = dp1[c] = x * prod_{gc} dp0[gc]  (subtree IS poly, c included)

Key polynomials:
  P  = prod_c f_c = prod_c (g_c + h_c)
  P' = prod_c g_c
  Q  = x * P'      (so q_k = P'_{k-1})
  E  = P - P'      (excess, all coefficients >= 0)

Target: P_{m-1} >= P'_{m-2}, equivalently p_{m-1} >= q_{m-1}.

Decomposition:
  p_{m-1} - q_{m-1} = P_{m-1} - P'_{m-2}
                     = [P'_{m-1} - P'_{m-2}] + E_{m-1}

Case 1 (ascending): P'_{m-1} >= P'_{m-2}  =>  trivially p >= q.
Case 2 (descending): Need E_{m-1} >= P'_{m-2} - P'_{m-1} = descent.

This script explores several structural approaches:

Approach 1: Single-child inclusion lower bound
  When d=1: P = f_1 = g_1 + h_1, P' = g_1, E = h_1.
  Then p_{m-1} - q_{m-1} = g_1[m-1] + h_1[m-1] - g_1[m-2].
  Since h_1 = x * prod_{gc} dp0[gc], h_1[m-1] = (prod gc-excluded subtree polys)[m-2].

Approach 2: Ratio bound via vertex probability
  p_u = Z_Q(lam) / Z_B(lam) = probability u is in the IS at fugacity lam.
  If p_u <= 1/2, then Z_P >= Z_Q at lam, hence P "dominates" Q in a weighted sense.
  Key question: does p_u <= 1/2 imply p_{m-1} >= q_{m-1}?

Approach 3: Weighted mean comparison
  mu_P(lam) = sum_c mu_{f_c}(lam), mu_{P'}(lam) = sum_c mu_{g_c}(lam).
  Since f_c = g_c + h_c with h_c >= 0, mu_{f_c}(lam) >= mu_{g_c}(lam).
  So mu_P >= mu_P'. Combined with P' mode structure, this constrains
  where P' can descend relative to P.

Approach 4: First-order inclusion-exclusion
  P = prod(g_c + h_c) = P' + sum_c [h_c * prod_{j!=c} g_j] + higher-order terms.
  The first-order excess E1 = sum_c [h_c * prod_{j!=c} g_j] has coefficients:
    E1_k = sum_c (h_c * prod_{j!=c} g_j)_k.
  Each term is h_c conv (prod_{j!=c} g_j), with h_c = [0, ...] (starts at x).
  So E1_{m-1} relates to the IS structure of T_c minus c.

Approach 5: Comparison via d_leaf <= 1 constraint
  In a d_leaf <= 1 tree, each child c of u has at most one leaf neighbor.
  This constrains the structure of g_c (and hence P') near the mode.

This script verifies these approaches computationally and records diagnostics.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from fractions import Fraction

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import compute_hub_polys, mode_index_leftmost, remove_vertices
from graph6 import parse_graph6
from indpoly import independence_poly, _polyadd, _polymul


def compute_child_data(adj_B: list[list[int]], u_in_B: int):
    """Compute full DP data for each child of u in B.

    Returns:
      children_of_u: list of child indices
      dp0, dp1: dp arrays for all vertices
      f_c: list of I(T_c) for each child c
      g_c: list of dp0[c] for each child c
      h_c: list of dp1[c] for each child c
    """
    n = len(adj_B)
    if n <= 1:
        return [], {}, {}, [], [], []

    # BFS to build tree rooted at u_in_B
    children_list: list[list[int]] = [[] for _ in range(n)]
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
                children_list[v].append(w)
                bfs_queue.append(w)

    # Post-order traversal
    order = []
    stack = [(u_in_B, False)]
    while stack:
        v, processed = stack.pop()
        if processed:
            order.append(v)
            continue
        stack.append((v, True))
        for c in children_list[v]:
            stack.append((c, False))

    dp0: dict[int, list[int]] = {}
    dp1: dict[int, list[int]] = {}

    for v in order:
        if not children_list[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
        else:
            prod = [1]
            for c in children_list[v]:
                summand = _polyadd(dp0[c], dp1[c])
                prod = _polymul(prod, summand)
            dp0[v] = prod

            prod = [1]
            for c in children_list[v]:
                prod = _polymul(prod, dp0[c])
            dp1[v] = [0] + prod

    children_of_u = children_list[u_in_B]
    f_list = [_polyadd(dp0[c], dp1[c]) for c in children_of_u]
    g_list = [dp0[c] for c in children_of_u]
    h_list = [dp1[c] for c in children_of_u]

    return children_of_u, dp0, dp1, f_list, g_list, h_list


def poly_coeff(p: list[int], k: int) -> int:
    if k < 0 or k >= len(p):
        return 0
    return p[k]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    total = 0
    p_ge_q = 0
    ascending = 0
    descending = 0
    descending_compensated = 0

    # Approach 1: first-order excess E1
    min_E1_ratio = None   # min E1_{m-1} / descent (when descending)
    E1_always_compensates = 0  # E1_{m-1} >= descent

    # Approach 2: vertex probability
    max_pu_at_mode = 0.0
    min_pu_at_mode = 1.0
    pu_ge_half_count = 0
    pu_lt_half_count = 0

    # Approach 3: mean comparison
    min_mu_gap = None  # min(mu_P - mu_P') at lam_m
    max_mu_gap = None

    # Approach 4: per-child excess analysis
    min_excess_ratio = None  # min E_{m-1}/descent (full excess)
    max_descent = 0
    min_E_at_descent = None

    # Approach 5: mode comparison
    mode_P_ge_m_minus_1 = 0
    mode_Pprime_ge_m_minus_1 = 0
    mode_Pprime_ge_m_minus_2 = 0

    # Deeper structural analysis
    # For single-child (d=1) trees: is h_c[m-1] >= descent always?
    single_child_count = 0
    single_child_h_compensates = 0

    # For multi-child trees: is the first-order term sufficient?
    multi_child_count = 0
    multi_child_E1_compensates = 0

    # New: comparison of P'_{m-1} and P'_{m-2} vs weighted means
    pprime_ascending_but_mu_low = 0  # P' ascending at m-1 but mu_P' < m-1

    # Track where the excess E comes from
    min_h_contrib_ratio = None  # min sum_c h_c[m-1] / E_{m-1}

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

            keep = [v for v in range(nn) if v not in {leaf, support}]
            idx_map = {v: i for i, v in enumerate(keep)}
            u_in_b = idx_map[u]

            # Get P, Q
            p_poly, q_poly = compute_hub_polys(b_adj, u_in_b)

            # Get child data
            children_of_u, dp0_map, dp1_map, f_list, g_list, h_list = \
                compute_child_data(b_adj, u_in_b)

            d = len(children_of_u)
            total += 1

            p1 = poly_coeff(p_poly, m - 1)
            q1 = poly_coeff(q_poly, m - 1)

            # Compute P' = prod g_c
            pprime = [1]
            for g in g_list:
                pprime = _polymul(pprime, g)

            pprime_m1 = poly_coeff(pprime, m - 1)
            pprime_m2 = poly_coeff(pprime, m - 2)

            # Verify q1 = pprime_m2 (since Q = x * P')
            assert q1 == pprime_m2, f"q1={q1} != pprime_m2={pprime_m2}"

            # E = P - P' (excess)
            E_m1 = p1 - pprime_m1
            assert E_m1 >= 0, f"E_{m-1} negative"

            gap = p1 - q1
            step = pprime_m1 - pprime_m2  # ascending step of P'

            if gap >= 0:
                p_ge_q += 1

            if step >= 0:
                ascending += 1
            else:
                descent = -step
                descending += 1
                if E_m1 >= descent:
                    descending_compensated += 1

                if descent > max_descent:
                    max_descent = descent
                if min_E_at_descent is None or E_m1 < min_E_at_descent:
                    min_E_at_descent = E_m1
                if descent > 0:
                    ratio = E_m1 / descent
                    if min_excess_ratio is None or ratio < min_excess_ratio:
                        min_excess_ratio = ratio

            # --- Approach 1: First-order excess E1 ---
            # E1 = sum_c [h_c * prod_{j!=c} g_j]
            E1 = [0] * max(len(p_poly), 1)
            for idx_c in range(d):
                # prod_{j != idx_c} g_j
                prod_others = [1]
                for j in range(d):
                    if j != idx_c:
                        prod_others = _polymul(prod_others, g_list[j])
                term = _polymul(h_list[idx_c], prod_others)
                E1 = _polyadd(E1, term)

            E1_m1 = poly_coeff(E1, m - 1)

            if step < 0:
                descent = -step
                if E1_m1 >= descent:
                    E1_always_compensates += 1
                if descent > 0:
                    r = E1_m1 / descent
                    if min_E1_ratio is None or r < min_E1_ratio:
                        min_E1_ratio = r

            # --- Approach 2: vertex probability p_u ---
            lam = poly_t[m - 1] / poly_t[m] if poly_t[m] > 0 else 1.0
            z_P = sum(poly_coeff(p_poly, k) * lam**k for k in range(len(p_poly)))
            z_Q = sum(poly_coeff(q_poly, k) * lam**k for k in range(len(q_poly)))
            z_B = z_P + z_Q
            p_u = z_Q / z_B if z_B > 0 else 0.0

            if p_u >= 0.5:
                pu_ge_half_count += 1
            else:
                pu_lt_half_count += 1
            if p_u > max_pu_at_mode:
                max_pu_at_mode = p_u
            if p_u < min_pu_at_mode:
                min_pu_at_mode = p_u

            # --- Approach 3: weighted mean comparison ---
            # mu_P = sum_c mu_{f_c}, mu_P' = sum_c mu_{g_c}
            mu_P_sum = 0.0
            mu_Pprime_sum = 0.0
            for idx_c in range(d):
                f = f_list[idx_c]
                g = g_list[idx_c]
                zf = sum(poly_coeff(f, k) * lam**k for k in range(len(f)))
                zg = sum(poly_coeff(g, k) * lam**k for k in range(len(g)))
                if zf > 0:
                    mu_f = sum(k * poly_coeff(f, k) * lam**k for k in range(len(f))) / zf
                    mu_P_sum += mu_f
                if zg > 0:
                    mu_g = sum(k * poly_coeff(g, k) * lam**k for k in range(len(g))) / zg
                    mu_Pprime_sum += mu_g

            mu_gap = mu_P_sum - mu_Pprime_sum
            if min_mu_gap is None or mu_gap < min_mu_gap:
                min_mu_gap = mu_gap
            if max_mu_gap is None or mu_gap > max_mu_gap:
                max_mu_gap = mu_gap

            # --- Approach 5: mode comparison ---
            mode_P = mode_index_leftmost(p_poly)
            mode_Pprime_val = mode_index_leftmost(pprime) if pprime else 0

            if mode_P >= m - 1:
                mode_P_ge_m_minus_1 += 1
            if mode_Pprime_val >= m - 1:
                mode_Pprime_ge_m_minus_1 += 1
            if mode_Pprime_val >= m - 2:
                mode_Pprime_ge_m_minus_2 += 1

            # --- Single vs multi child ---
            if d == 1:
                single_child_count += 1
                h_m1 = poly_coeff(h_list[0], m - 1)
                if step < 0:
                    descent = -step
                    if h_m1 >= descent:
                        single_child_h_compensates += 1
            else:
                multi_child_count += 1
                if step < 0:
                    descent = -step
                    if E1_m1 >= descent:
                        multi_child_E1_compensates += 1

        proc.wait()
        elapsed = time.time() - t0
        print(f"n={n:2d}: total={total:8d} p_ge_q={p_ge_q} asc={ascending} "
              f"desc={descending} compensated={descending_compensated} "
              f"({elapsed:.1f}s)", flush=True)

    elapsed = time.time() - t0
    results = {
        "total": total,
        "p_ge_q": p_ge_q,
        "ascending": ascending,
        "descending": descending,
        "descending_compensated": descending_compensated,
        "max_descent": max_descent,
        "min_E_at_descent": min_E_at_descent,
        "min_excess_ratio_full": min_excess_ratio,

        "E1_always_compensates_descending": E1_always_compensates,
        "min_E1_ratio_descending": min_E1_ratio,

        "max_pu_at_mode_fugacity": max_pu_at_mode,
        "min_pu_at_mode_fugacity": min_pu_at_mode,
        "pu_ge_half_count": pu_ge_half_count,
        "pu_lt_half_count": pu_lt_half_count,

        "min_mu_gap": min_mu_gap,
        "max_mu_gap": max_mu_gap,

        "mode_P_ge_m_minus_1": mode_P_ge_m_minus_1,
        "mode_Pprime_ge_m_minus_1": mode_Pprime_ge_m_minus_1,
        "mode_Pprime_ge_m_minus_2": mode_Pprime_ge_m_minus_2,

        "single_child_count": single_child_count,
        "single_child_h_compensates_descending": single_child_h_compensates,
        "multi_child_count": multi_child_count,
        "multi_child_E1_compensates_descending": multi_child_E1_compensates,

        "wall_s": elapsed,
    }

    print(f"\n{'='*70}")
    print(f"ROUTE C: P-DOMINANCE STRUCTURAL ANALYSIS (n <= {args.max_n})")
    print(f"{'='*70}")
    print(f"Total trees checked: {total}")
    print(f"p_{{m-1}} >= q_{{m-1}}: {p_ge_q}/{total} "
          f"({'PASS' if p_ge_q == total else 'FAIL'})")
    print()

    print("--- Ascending/Descending at P' ---")
    print(f"P' ascending at m-1: {ascending} ({100*ascending/total:.1f}%)")
    print(f"P' descending at m-1: {descending} ({100*descending/total:.1f}%)")
    print(f"  All compensated by full E: {descending_compensated}")
    print(f"  Max descent: {max_descent}")
    print(f"  Min E at descent: {min_E_at_descent}")
    print(f"  Min excess ratio (E/descent): {min_excess_ratio}")
    print()

    print("--- Approach 1: First-order excess E1 ---")
    print(f"E1 compensates all descents: "
          f"{E1_always_compensates}/{descending}")
    print(f"Min E1 ratio at descent: {min_E1_ratio}")
    print()

    print("--- Approach 2: Vertex probability p_u ---")
    print(f"p_u < 1/2 (P dominates in Z): {pu_lt_half_count}/{total} "
          f"({100*pu_lt_half_count/total:.1f}%)")
    print(f"p_u >= 1/2: {pu_ge_half_count}/{total}")
    print(f"Max p_u: {max_pu_at_mode:.6f}")
    print(f"Min p_u: {min_pu_at_mode:.6f}")
    print()

    print("--- Approach 3: Weighted mean gap ---")
    print(f"mu_P - mu_P' min: {min_mu_gap}")
    print(f"mu_P - mu_P' max: {max_mu_gap}")
    print()

    print("--- Approach 5: Mode structure ---")
    print(f"mode(P) >= m-1: {mode_P_ge_m_minus_1}/{total} "
          f"({100*mode_P_ge_m_minus_1/total:.1f}%)")
    print(f"mode(P') >= m-1: {mode_Pprime_ge_m_minus_1}/{total} "
          f"({100*mode_Pprime_ge_m_minus_1/total:.1f}%)")
    print(f"mode(P') >= m-2: {mode_Pprime_ge_m_minus_2}/{total} "
          f"({100*mode_Pprime_ge_m_minus_2/total:.1f}%)")
    print()

    print("--- Single vs multi child ---")
    print(f"Single child (d=1): {single_child_count}")
    print(f"  h_c compensates descent: {single_child_h_compensates}")
    print(f"Multi child (d>=2): {multi_child_count}")
    print(f"  E1 compensates descent: {multi_child_E1_compensates}")
    print()

    print(f"Time: {elapsed:.1f}s")

    if args.out:
        with open(args.out, "w") as f:
            json.dump(results, f, indent=2)
        print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
