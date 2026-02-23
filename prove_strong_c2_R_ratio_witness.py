#!/usr/bin/env python3
"""Find the tree with minimum R/p1 ratio."""

from __future__ import annotations

import subprocess

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import compute_hub_polys, mode_index_leftmost, remove_vertices
from graph6 import parse_graph6
from indpoly import independence_poly


def main() -> None:
    min_ratio = None
    min_witness = None

    for n in range(4, 21):
        cmd = ["/opt/homebrew/bin/geng", "-q", str(n), f"{n-1}:{n-1}", "-c"]
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

            p0 = p_poly[m - 2] if m - 2 < len(p_poly) else 0
            p1 = p_poly[m - 1] if m - 1 < len(p_poly) else 0
            pm = p_poly[m] if m < len(p_poly) else 0
            q0 = q_poly[m - 2] if m - 2 >= 0 and m - 2 < len(q_poly) else 0
            q1 = q_poly[m - 1] if m - 1 < len(q_poly) else 0
            qm = q_poly[m] if m < len(q_poly) else 0

            R = 2 * p1 * q1 - pm * q0 - p0 * qm + p0 * q1 - p1 * q0

            if p1 > 0:
                ratio = R / p1
                if min_ratio is None or ratio < min_ratio:
                    min_ratio = ratio
                    min_witness = {
                        "n": nn, "m": m,
                        "p": [p0, p1, pm],
                        "q": [q0, q1, qm],
                        "R": R, "R/p1": ratio,
                        "g6": raw.decode("ascii").strip(),
                        "P": list(p_poly),
                        "Q": list(q_poly),
                        "deg_u": deg[u],
                    }

        proc.wait()

    print(f"Min R/p1 = {min_ratio}")
    print(f"Witness: {min_witness}")


if __name__ == "__main__":
    main()
