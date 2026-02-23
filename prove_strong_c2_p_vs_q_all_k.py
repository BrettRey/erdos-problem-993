#!/usr/bin/env python3
"""Check whether p_k >= q_k for ALL indices k (not just k = m-1).

If P dominates Q coefficient-wise, then STRONG C2 proof simplifies dramatically:
  p_k >= q_k for all k
  => p1/b1 = p1/(p1+q1) >= 1/2
  => hard-regime ratio inequality trivially bounded.

But even stronger: p_k >= q_k means P = Q + E where E has all non-neg coefficients.
Then b_k = p_k + q_k = 2*q_k + E_k, and
  neg = p1*q0 - p0*q1 = (q1+E1)*q0 - (q0+E0)*q1 = E1*q0 - E0*q1
     = (E1/E0 - q1/q0) * E0 * q0   [when E0, q0 > 0]

This could yield a clean proof if E is also "nicely behaved".

NOTE: Q = x * prod_c dp0[c], so q_0 = 0. And P_0 = 1 always (empty set).
So p_0 = 1, q_0 = 0: p_0 > q_0 trivially.
And p_1 = sum_c |V(T_c)| (number of vertices across subtrees),
   q_1 = #{children of u} = deg(u) - 1 in original T (or just |children| in B).
Actually Q = x * prod dp0[c], so q_1 = prod dp0[c]_0 = 1 (since dp0[c]_0 = 1 for all c).
Hmm wait: q_1 = (prod dp0[c])_0 = 1. And p_1 = sum_c (dp0[c] + dp1[c])_1 convolved.
For a single child: p_1 = dp0[c]_1 + dp1[c]_1 = |V(T_c)| + ... Actually it's more complex
for multiple children.

Let me just check computationally.
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
    p_dominates_q = 0  # p_k >= q_k for ALL k
    fail_count = 0
    first_fail = None

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

            total += 1

            # Check p_k >= q_k for all k
            max_len = max(len(p_poly), len(q_poly))
            dominated = True
            min_diff = None
            min_diff_k = None
            for k in range(max_len):
                pk = p_poly[k] if k < len(p_poly) else 0
                qk = q_poly[k] if k < len(q_poly) else 0
                diff = pk - qk
                if min_diff is None or diff < min_diff:
                    min_diff = diff
                    min_diff_k = k
                if pk < qk:
                    dominated = False
                    if first_fail is None:
                        first_fail = {
                            "n": nn, "m": m, "k": k,
                            "pk": pk, "qk": qk,
                            "g6": raw.decode("ascii").strip(),
                            "P": list(p_poly),
                            "Q": list(q_poly),
                        }
                    break

            if dominated:
                p_dominates_q += 1
            else:
                fail_count += 1

        proc.wait()
        print(f"n={n:2d}: total={total:8d} dominated={p_dominates_q:8d} fails={fail_count}", flush=True)

    print(f"\nTotal: {total}")
    print(f"P dominates Q (all k): {p_dominates_q} / {total} ({100*p_dominates_q/total:.2f}%)" if total else "")
    print(f"Failures: {fail_count}")
    if first_fail:
        print(f"First fail: {first_fail}")
    print(f"Time: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
