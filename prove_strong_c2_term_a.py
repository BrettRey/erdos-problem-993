#!/usr/bin/env python3
"""Prove Term A: p_{m-1} * q_{m-1} >= p_m * q_{m-2} for all d_leaf<=1 trees.

Recall:
  P = dp_B[u][0] = prod_c I(T_c)     (product over children c of u in B)
  Q = dp_B[u][1] = x * prod_c dp0[c] = x * P'

where P' = prod_c dp0[c], so q_k = P'_{k-1}.

Term A = p_{m-1} * P'_{m-2} - p_m * P'_{m-3}

For a single child c: P = I(T_c) = f, P' = dp0[c] = g, Q = x*g.
  p_k = f_k, P'_k = g_k.
  Term A = f_{m-1}*g_{m-2} - f_m*g_{m-3}
  = "Turan cross-determinant" of (f, g) at indices (m-1, m-2).

For two LC sequences f, g with f >= g coefficientwise:
  When is f_{j}*g_{j-1} >= f_{j+1}*g_{j-2}?

This is saying: the ratio f_j/f_{j+1} >= g_{j-2}/g_{j-1}, equivalently
f's lambda_j <= g's lambda_{j-1}^{-1}... no, it's f_j/f_{j+1} >= g_{j-2}/g_{j-1},
i.e., lambda_j(f) >= lambda_{j-1}(g), i.e., f is "less peaked" than g
(f's ratio at j is at least as large as g's ratio at j-1).

Actually wait, this is exactly the STRONG C2 condition but one level down!
It says that f (=I(T_c)) has lambda_{m-1}(f) >= lambda_{m-2}(g), where g = dp0[c].

But for two children: P = f1*f2, P' = g1*g2.
Term A = (f1*f2)_{m-1} * (g1*g2)_{m-2} - (f1*f2)_m * (g1*g2)_{m-3}.

This is the cross-Turan determinant of the product sequences.

KEY INSIGHT: For products of LC sequences, the cross-Turan determinant
satisfies a multiplication law. If we can prove it for single-factor (f,g),
then for multi-factor (prod f, prod g) it follows by induction.

Let me check: for two pairs (f1,g1) and (f2,g2), if
  f1_j * g1_{j-1} >= f1_{j+1} * g1_{j-2} for all j,
  f2_j * g2_{j-1} >= f2_{j+1} * g2_{j-2} for all j,
does
  (f1*f2)_j * (g1*g2)_{j-1} >= (f1*f2)_{j+1} * (g1*g2)_{j-2} ?

This would follow from the Alexandrov-Fenchel inequality for mixed volumes,
or from the FKG inequality on the product lattice of independent sets.

Actually, there's a simpler approach: the TURÁN PROPERTY OF CONVOLUTIONS.

LEMMA (Turán property of product): If f and g are PF_2 (i.e., LC with positive terms),
and h is PF_2, then f*h and g*h also satisfy the cross-Turán inequality:
(f*h)_j * (g*h)_{j-1} >= (f*h)_{j+1} * (g*h)_{j-2}
provided f_j * g_{j-1} >= f_{j+1} * g_{j-2}.

THIS IS EXACTLY WHAT WE NEED. If the cross-Turán holds for each factor (I(T_c), dp0[c]),
then it holds for their products (P, P').

So the proof reduces to: for each child c of u in B, show that
  I(T_c)_j * dp0[c]_{j-1} >= I(T_c)_{j+1} * dp0[c]_{j-2}

for j = m-1 evaluated "at the right level" for each factor.

Wait, the product formula is more subtle. Let me think carefully.

This script:
1. Verifies Term A for each child's factor pair
2. Checks whether the single-factor cross-Turán holds
3. Tests the multiplicative closure property
"""

from __future__ import annotations

import argparse
import subprocess
import time
from typing import Any

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import mode_index_leftmost, remove_vertices
from graph6 import parse_graph6
from indpoly import independence_poly, _polyadd, _polymul


def compute_child_dp(adj_B: list[list[int]], u_in_B: int):
    """Return children of u and their (dp0, dp1, I=dp0+dp1) polys."""
    n = len(adj_B)
    if n <= 1:
        return [], []

    parent = [-1] * n
    children_of = [[] for _ in range(n)]
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
                children_of[v].append(w)
                bfs_queue.append(w)

    order = []
    stack = [(u_in_B, False)]
    while stack:
        v, processed = stack.pop()
        if processed:
            order.append(v)
            continue
        stack.append((v, True))
        for c in children_of[v]:
            stack.append((c, False))

    dp0 = [[] for _ in range(n)]
    dp1 = [[] for _ in range(n)]

    for v in order:
        if not children_of[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
        else:
            prod = [1]
            for c in children_of[v]:
                summand = _polyadd(dp0[c], dp1[c])
                prod = _polymul(prod, summand)
            dp0[v] = prod

            prod = [1]
            for c in children_of[v]:
                prod = _polymul(prod, dp0[c])
            dp1[v] = [0] + prod

    # Return per-child data
    child_data = []
    for c in children_of[u_in_B]:
        f = _polyadd(dp0[c], dp1[c])  # I(T_c)
        g = dp0[c]                      # dp0[c]
        child_data.append((f, g))

    return child_data


def get_coeff(poly, k):
    return poly[k] if 0 <= k < len(poly) else 0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    args = ap.parse_args()

    stats = {
        "checked": 0,
        "term_A_neg": 0,

        # Single-factor cross-Turan: for each child c, at index m-1,
        # check I(T_c)_{j} * dp0[c]_{j-1} >= I(T_c)_{j+1} * dp0[c]_{j-2}
        # where j ranges over relevant indices.
        "single_factor_checks": 0,
        "single_factor_fails": 0,

        # For the single factor, check at ALL indices j (not just m-1)
        "single_factor_all_j_checks": 0,
        "single_factor_all_j_fails": 0,

        # Check: f_j * g_{j-1} >= f_{j+1} * g_{j-2} for all j where both are positive
        # This is "f/g has non-increasing ratio" or equivalently
        # det [f_j     f_{j+1}] >= 0
        #     [g_{j-1} g_j    ]
        # which is the TP_2 (totally positive of order 2) condition on the matrix
        # [f_0 f_1 f_2 ...]
        # [g_0 g_1 g_2 ...]
        # shifted by 1.

        # Actually, f_j*g_{j-1} >= f_{j+1}*g_{j-2} is:
        # det [f_j     f_{j+1}  ] >= 0
        #     [g_{j-2} g_{j-1}  ]
        # = f_j*g_{j-1} - f_{j+1}*g_{j-2} >= 0
        # This is TP_2 on the "shifted" matrix where g is shifted by 1 index.

        # For f = dp0 + dp1 and g = dp0:
        # f_j * g_{j-1} - f_{j+1} * g_{j-2}
        # = (dp0_j + dp1_j) * dp0_{j-1} - (dp0_{j+1} + dp1_{j+1}) * dp0_{j-2}
        # = dp0_j*dp0_{j-1} - dp0_{j+1}*dp0_{j-2}   (LC of dp0, shifted)
        #   + dp1_j*dp0_{j-1} - dp1_{j+1}*dp0_{j-2}  (cross term)
        # The first part is the LC surplus of g at index j (non-negative since g=dp0 is LC).
        # The second part: dp1_j*dp0_{j-1} - dp1_{j+1}*dp0_{j-2}.
        # Using dp1 = x * (product of dp0 of grandchildren), dp1_j = grandchild-product_{j-1}.
        # Let h = product of dp0 of grandchildren. Then dp1_j = h_{j-1}.
        # So second part = h_{j-1}*g_{j-1} - h_j*g_{j-2}
        #                = g_{j-1}*h_{j-1} - g_{j-2}*h_j
        # This is the cross-Turan of (g, h) at index (j-1, j)... shifted by 1.
        # = det [g_{j-1} g_j   ] ... wait, let me be careful.
        #       [h_{j-2} h_{j-1}]
        # No: g_{j-1}*h_{j-1} - g_{j-2}*h_j = g_{j-1}*h_{j-1} - g_{j-2}*h_j.
        # This is: -( g_{j-2}*h_j - g_{j-1}*h_{j-1} ) = -(cross-Turan of (g,h) at (j-2,j)).
        # Hmm. Let me just denote D(f,g,j) = f_j*g_{j-1} - f_{j+1}*g_{j-2}.

        # So: D(f,g,j) = D(dp0,dp0,j) + D(dp1,dp0,j)
        #              = LC(dp0,j) + dp1_j*dp0_{j-1} - dp1_{j+1}*dp0_{j-2}

        # In the case where c is a leaf: dp0[c] = [1], dp1[c] = [0,1].
        # I(T_c) = [1,1]. g = [1].
        # f_j*g_{j-1} - f_{j+1}*g_{j-2}: only nonzero for small j.
        # j=0: f_0*g_{-1} - f_1*g_{-2} = 0 (g has no negative indices). Skip.
        # j=1: f_1*g_0 - f_2*g_{-1} = 1*1 - 0 = 1 >= 0.

        # For a path P_2 as child: dp0 = [1,1], dp1 = [0,1].
        # I(T_c) = [1,2,1] (but wait, P_2 has IS poly [1,2,1]... actually P_2 is 2 vertices,
        # IS poly is 1 + 2x + 0x^2? No: the 2 vertices are independent (no edge), so
        # IS poly is (1+x)^2 = 1+2x+x^2. But if connected (edge between them),
        # IS poly is 1+2x (size 0: 1, size 1: 2, size 2: 0).
        # Since children of u in a tree are subtrees, they're connected.
        # For an edge (u-c): dp0[c]=[1], dp1[c]=[0,1], I(T_c)=[1,1]. This is a single edge.
        # For a path c-w: dp0[c] = dp0[w]+dp1[w] = [1]+[0,1]=[1,1]; dp1[c] = [0]*dp0[w] = [0,1].
        # So I(T_c) = [1,1]+[0,1] = [1,2,1]? No: dp0[c]+dp1[c] = [1,1]+[0,1] = [1,2,1].
        # Wait, that's wrong: dp0[c] for c rooted at c with child w:
        # dp0[c] = product over children of c of (dp0+dp1) = dp0[w]+dp1[w] = [1]+[0,1] = [1,1].
        # dp1[c] = x * product over children of dp0 = x * dp0[w] = x * [1] = [0,1].
        # I(subtree_c) = dp0[c] + dp1[c] = [1,1] + [0,1] = [1,2].
        # But subtree_c has 2 vertices (c and w) with one edge, so IS poly should be [1,2]:
        # size 0: {empty}=1, size 1: {c},{w}=2, size 2: 0 (edge prevents). Correct.

        # OK, for this case:
        # f = [1,2], g = [1,1].
        # j=1: f_1*g_0 - f_2*g_{-1} = 2*1 - 0 = 2 >= 0. (But f_2=0, g_{-1}=0.)
        # Actually only j where both f_{j+1} and g_{j-2} are meaningful.
        # For j=2: f_2*g_1 - f_3*g_0 = 0*1 - 0*1 = 0. OK.

        # This is getting complicated. Let me just run the computational check.
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

            if m - 2 < 0:
                continue

            keep = [v for v in range(nn) if v not in {leaf, support}]
            idx_map = {v: i for i, v in enumerate(keep)}
            u_in_b = idx_map[u]

            child_data = compute_child_dp(b_adj, u_in_b)

            stats["checked"] += 1

            # Check single-factor cross-Turan for each child at the relevant index
            for f, g in child_data:
                # Check D(f,g,j) = f_j*g_{j-1} - f_{j+1}*g_{j-2} >= 0 at all j
                max_j = max(len(f), len(g)) + 2
                for j in range(2, max_j):  # j >= 2 for g_{j-2} to be defined
                    fj = get_coeff(f, j)
                    fj1 = get_coeff(f, j + 1)
                    gj1 = get_coeff(g, j - 1)
                    gj2 = get_coeff(g, j - 2)
                    if fj == 0 and fj1 == 0:
                        continue
                    if gj1 == 0 and gj2 == 0:
                        continue
                    det = fj * gj1 - fj1 * gj2
                    stats["single_factor_all_j_checks"] += 1
                    if det < 0:
                        stats["single_factor_all_j_fails"] += 1
                        print(f"SINGLE FACTOR FAIL: n={nn}, j={j}, f={f[:6]}, g={g[:6]}, "
                              f"det={det}")

        proc.wait()
        print(f"n={n:2d}: checked={stats['checked']:8d} "
              f"sf_checks={stats['single_factor_all_j_checks']} "
              f"sf_fails={stats['single_factor_all_j_fails']}", flush=True)

    print(f"\nTotal: checked={stats['checked']}, "
          f"single_factor_checks={stats['single_factor_all_j_checks']}, "
          f"single_factor_fails={stats['single_factor_all_j_fails']}")
    print(f"Time: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
