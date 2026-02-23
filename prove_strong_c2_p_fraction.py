#!/usr/bin/env python3
"""Check the P-fraction p_{m-1}/b_{m-1} = p_{m-1}/(p_{m-1}+q_{m-1}).

If p_{m-1}/b_{m-1} >= 1/2 always (i.e., p_{m-1} >= q_{m-1}), then
the hard-regime ratio inequality (-dq/db) <= p1/b1 reduces to
(-dq/db) <= 1/2, which needs db >= 2*(-dq).

More ambitiously, if p_{m-1}/b_{m-1} is bounded away from 0, we get
a uniform bound.

Also check: at lambda_m = b_{m-2}/b_{m-1}, what is the vertex probability P(u)?
P(u) = Z_Q(lambda) / Z_B(lambda) = Q(lambda)/B(lambda).
If P(u) <= 1/3 at lambda = lambda_m, then q1/b1 <= 1/3 and p1/b1 >= 2/3.
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

    min_p_frac = None
    min_p_frac_witness = None
    total = 0
    p_frac_ge_half = 0
    p_frac_ge_two_thirds = 0
    min_p_frac_all = None  # across ALL trees, not just mismatch-neg

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

            total += 1

            p1 = p_poly[m - 1] if m - 1 < len(p_poly) else 0
            q1 = q_poly[m - 1] if m - 1 < len(q_poly) else 0
            b1 = p1 + q1

            if b1 == 0:
                continue

            p_frac = p1 / b1

            if min_p_frac_all is None or p_frac < min_p_frac_all:
                min_p_frac_all = p_frac

            if p_frac >= 0.5:
                p_frac_ge_half += 1
            if p_frac >= 2/3:
                p_frac_ge_two_thirds += 1

            if min_p_frac is None or p_frac < min_p_frac:
                min_p_frac = p_frac
                min_p_frac_witness = {
                    "n": nn, "m": m,
                    "p1": p1, "q1": q1, "b1": b1,
                    "p_frac": p_frac,
                    "g6": raw.decode("ascii").strip(),
                    "deg_u": deg[u] if u < nn else -1,
                }

        proc.wait()
        mpf_str = f"{min_p_frac:.6f}" if min_p_frac is not None else "N/A"
        print(f"n={n:2d}: total={total:8d} min_p_frac={mpf_str}", flush=True)

    if total == 0:
        print("No trees checked.")
        return
    print(f"\nTotal trees: {total}")
    print(f"Min p-fraction: {min_p_frac:.10f}")
    print(f"p-frac >= 1/2: {p_frac_ge_half} / {total} ({100*p_frac_ge_half/total:.2f}%)")
    print(f"p-frac >= 2/3: {p_frac_ge_two_thirds} / {total} ({100*p_frac_ge_two_thirds/total:.2f}%)")
    print(f"Witness: {min_p_frac_witness}")
    print(f"Time: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
