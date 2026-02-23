#!/usr/bin/env python3
"""Prove p_{m-1} >= q_{m-1} (P-dominance at mode index).

Structural argument:
  P = prod_c I(T_c), Q = x * prod_c dp0[c].
  Since I(T_c) = dp0[c] + dp1[c] with dp1[c] >= 0 coefficient-wise,
  P_k >= (prod dp0)_k for all k.

  Also q_k = (prod dp0)_{k-1}.

  So p_{m-1} >= (prod dp0)_{m-1}.
  If (prod dp0)_{m-1} >= (prod dp0)_{m-2} = q_{m-1}, done.
  This is mode(prod dp0) >= m-1, equivalently mode(Q) >= m.

  But mode(Q) >= m fails for 2 Q-drop trees.

  Alternative: We need
    P_{m-1} >= q_{m-1} = (prod dp0)_{m-2}.

  Write D_c = dp1[c] = x * prod_{grandchildren} dp0[...].
  Then I(T_c) = dp0[c] + D_c.

  Expanding the product:
    P = prod_c (dp0[c] + D_c) = prod_c dp0[c] + sum_c D_c * prod_{j!=c} dp0[j] + ...

  So P_k = (prod dp0)_k + (sum_c D_c * prod_{j!=c} dp0[j])_k + ...

  The "excess" E = P - prod dp0 has all non-negative coefficients.
  E_k = P_k - (prod dp0)_k >= 0.

  p_{m-1} - q_{m-1} = P_{m-1} - (prod dp0)_{m-2}
                     = (prod dp0)_{m-1} + E_{m-1} - (prod dp0)_{m-2}
                     = [(prod dp0)_{m-1} - (prod dp0)_{m-2}] + E_{m-1}

  The first bracket is the "ascending step" of prod dp0 at m-1.
  E_{m-1} is the excess from dp1 contributions (always >= 0).

  So p_{m-1} >= q_{m-1} iff
    E_{m-1} >= (prod dp0)_{m-2} - (prod dp0)_{m-1}   [when prod dp0 is descending at m-1]
  or trivially when prod dp0 is ascending at m-1.

This script checks:
1. Is p_{m-1} >= q_{m-1} always? (verified above: yes through n=22)
2. For cases where prod dp0 descends at m-1: how big is E_{m-1}?
3. What fraction of cases need E_{m-1} compensation?
4. Statistics on the ratio E_{m-1} / descent.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import compute_hub_polys, mode_index_leftmost, remove_vertices
from graph6 import parse_graph6
from indpoly import independence_poly, _polyadd, _polymul


def compute_child_polys(adj_B: list[list[int]], u_in_B: int):
    """Compute dp0[c] and dp1[c] for each child c of u in B.

    Returns list of (dp0_c, dp1_c) tuples, one per child.
    """
    n = len(adj_B)
    if n <= 1:
        return []

    # BFS to build tree rooted at u_in_B
    children: list[list[int]] = [[] for _ in range(n)]
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
                children[v].append(w)
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
        for c in children[v]:
            stack.append((c, False))

    dp0: list[list[int]] = [[] for _ in range(n)]
    dp1: list[list[int]] = [[] for _ in range(n)]

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

    # Return child polys for u's children
    result = []
    for c in children[u_in_B]:
        result.append((dp0[c], dp1[c]))
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    total = 0
    p_ge_q = 0
    ascending = 0  # prod dp0 ascending at m-1
    descending_compensated = 0  # descending but E compensates
    descent_count = 0
    min_gap = None  # min(p_{m-1} - q_{m-1})
    min_gap_witness = None
    min_excess_ratio = None  # min E_{m-1} / descent
    max_descent = 0
    min_E_at_descent = None

    # Also check: ratio p1/(p1+q1) when q drops
    q_drop_p_ge_q = 0  # p1>=q1 even in Q-drop cases
    q_drop_count = 0

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

            p_poly, q_poly = compute_hub_polys(b_adj, u_in_b)
            child_polys = compute_child_polys(b_adj, u_in_b)

            total += 1

            p1 = p_poly[m - 1] if m - 1 < len(p_poly) else 0
            q1 = q_poly[m - 1] if m - 1 < len(q_poly) else 0

            # Compute prod dp0
            prod_dp0 = [1]
            for dp0_c, dp1_c in child_polys:
                prod_dp0 = _polymul(prod_dp0, dp0_c)

            pd_m1 = prod_dp0[m - 1] if m - 1 < len(prod_dp0) else 0  # (prod dp0)_{m-1}
            pd_m2 = prod_dp0[m - 2] if m - 2 >= 0 and m - 2 < len(prod_dp0) else 0  # (prod dp0)_{m-2}

            # q_{m-1} should equal pd_m2 (since Q = x * prod dp0)
            assert q1 == pd_m2, f"q1={q1} != pd_m2={pd_m2}"

            # Excess: E_{m-1} = P_{m-1} - (prod dp0)_{m-1}
            E_m1 = p1 - pd_m1
            assert E_m1 >= 0, f"E_m1={E_m1} < 0"

            # p1 - q1 = (pd_m1 - pd_m2) + E_m1
            gap = p1 - q1
            step = pd_m1 - pd_m2  # ascending step of prod dp0

            if gap >= 0:
                p_ge_q += 1

            if step >= 0:
                ascending += 1
                # p1 >= q1 trivially (step + E both non-negative)
            else:
                descent = -step  # how much prod dp0 descends
                descent_count += 1
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

            if min_gap is None or gap < min_gap:
                min_gap = gap
                min_gap_witness = {
                    "n": nn, "m": m,
                    "p1": p1, "q1": q1,
                    "pd_m1": pd_m1, "pd_m2": pd_m2,
                    "E_m1": E_m1, "step": step,
                    "g6": raw.decode("ascii").strip(),
                }

            # Check Q-drop status
            q0 = q_poly[m - 2] if m - 2 >= 0 and m - 2 < len(q_poly) else 0
            if q1 < q0:
                q_drop_count += 1
                if p1 >= q1:
                    q_drop_p_ge_q += 1

        proc.wait()
        mpf_str = f"min_gap={min_gap}" if min_gap is not None else "min_gap=N/A"
        print(f"n={n:2d}: total={total:8d} {mpf_str}", flush=True)

    results = {
        "total": total,
        "p_ge_q": p_ge_q,
        "ascending": ascending,
        "descent_count": descent_count,
        "descending_compensated": descending_compensated,
        "min_gap": min_gap,
        "max_descent": max_descent,
        "min_E_at_descent": min_E_at_descent,
        "min_excess_ratio": min_excess_ratio,
        "q_drop_count": q_drop_count,
        "q_drop_p_ge_q": q_drop_p_ge_q,
        "min_gap_witness": min_gap_witness,
        "wall_s": time.time() - t0,
    }

    print(f"\n=== P-dominance analysis (n <= {args.max_n}) ===")
    print(f"Total trees: {total}")
    print(f"p1 >= q1: {p_ge_q} / {total} ({100*p_ge_q/total:.2f}%)" if total > 0 else "No trees")
    print(f"  Ascending (prod dp0 non-decreasing at m-1): {ascending}")
    print(f"  Descending: {descent_count}")
    print(f"    Compensated by E: {descending_compensated}")
    print(f"  Min gap (p1 - q1): {min_gap}")
    print(f"  Max descent: {max_descent}")
    print(f"  Min E at descent: {min_E_at_descent}")
    print(f"  Min excess ratio (E/descent): {min_excess_ratio}")
    print(f"Q-drop cases: {q_drop_count}, p1>=q1 in {q_drop_p_ge_q}")
    print(f"Witness: {min_gap_witness}")
    print(f"Time: {time.time()-t0:.1f}s")

    if args.out:
        with open(args.out, "w") as f:
            json.dump(results, f, indent=2)
        print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
