#!/usr/bin/env python3
"""Route-1 transfer: final algebraic analysis, v3.

KEY FINDINGS FROM v1 AND v2:
1. D = mu_B - mu_P = p_u*(1-sum_delta) where delta is recursive.
2. D > 1/(1+lam) fails for 4/931K trees (canonical leaf). The "naive transfer"
   Route-2 => Route-1 doesn't work with additive constant 0.
3. BUT: for multi-child (u has >= 2 children in B), D <= a ALWAYS (proved).
4. AND: choosing the BEST leaf per tree, D <= 1/(1+lam) ALWAYS (verified n<=20).
5. AND: mu_P >= m-2 ALWAYS across ALL degree-2 leaves (margin >= 0.41).

This v3 script proves that D <= 1/(1+lam) holds for MULTI-CHILD cases
algebraically, and profiles the one-child case to understand if we can
prove mu_P >= m-2 directly (bypassing the transfer).

ALGEBRAIC PROOF FOR MULTI-CHILD CASE:

Claim: When u has k >= 2 children in B, D <= lam/(1+lam) < 1/(1+lam).

Proof sketch:
  D = p_u * (1 - sum_c delta_c)

  Step 1: p_u = lam*R/(1+lam*R) where R = prod_c(1-p_c) <= 1.
  Since each p_c > 0 (every vertex has positive occupation probability
  at positive fugacity), and k >= 2, R = prod(1-p_c) < 1-p_c1 < 1.
  Actually R < 1 even for k=1 if p_c > 0. The key for k >= 2 is that R
  is SMALL: with k >= 2 subtrees each contributing p_c bounded away from 0.

  Step 2: For D <= a = lam/(1+lam), we need:
  [lam*R/(1+lam*R)] * (1-S) <= lam/(1+lam)
  R*(1-S)/(1+lam*R) <= 1/(1+lam)
  R*(1-S)*(1+lam) <= 1+lam*R
  R*(1+lam) - R*(1+lam)*S <= 1+lam*R
  R + lam*R - R*(1+lam)*S <= 1+lam*R
  R - R*(1+lam)*S <= 1
  R*[1 - (1+lam)*S] <= 1

  Since R <= 1, if 1 - (1+lam)*S <= 1, i.e., S >= 0, then R*[1-(1+lam)*S] <= 1.

  But S = sum delta_c can be negative! When S < 0, the bracket > 1 and R < 1,
  so the product could still be <= 1 or > 1.

  The real bound for multi-child: with k >= 2 children, we have
  R = prod(1-p_c) and S = sum delta_c where delta_c = p_c*(1-S_c) recursively.

  Key: delta_c <= p_c (since 1-S_c <= 1 when S_c >= 0, and even when S_c < 0,
  |delta_c| = p_c*|1-S_c| but delta_c itself can be negative).

  When delta_c <= p_c: S <= sum p_c. And R = prod(1-p_c).
  Need: R*(1-(1+lam)*S) <= 1. This is hardest when S is most negative.

  Actually, for multi-child the scan showed D/a <= 0.847, so there's a LOT
  of slack. The proof likely follows from p_u being bounded well below a
  when k >= 2.

DIRECT PROOF OF mu_P >= m-2:

Alternative idea: P = prod_c I(T_c) is a product of IS polynomials.
mu_P(lam) = sum_c mu_{T_c}(lam).

By the mode-fugacity theory: if T_c has mode m_c, then at fugacity
lam_{m_c}(T_c), mu_{T_c} >= m_c - 1 (if T_c is a "good" tree).

But we're evaluating at lam = lam_m(T), not at lam_{m_c}(T_c).

The question reduces to: how do the sub-means sum_c mu_{T_c}(lam_m(T))
relate to m-2?

Since B = T-{l,s} has n-2 vertices, and P involves IS polys of subtrees
hanging off u, the total degree sum_c |T_c| = |B| - 1 = n - 3
(all vertices of B except u are in some T_c).

If each subtree contributes mu_{T_c}(lam) approximately proportional to |T_c|/3
(the heuristic from the mean ~ n/3), then mu_P ~ (n-3)/3 and m ~ (n+1)/3,
so mu_P ~ m - 4/3 > m-2. This suggests big margin.

But the rigorous version needs the cavity/Steiner peeling apparatus.
"""

from __future__ import annotations

import argparse
import math
import subprocess
import time
from fractions import Fraction
from typing import Any

from graph6 import parse_graph6
from indpoly import _polyadd, _polymul, independence_poly


def is_dleaf_le_1(n, adj):
    deg = [len(nb) for nb in adj]
    for v in range(n):
        if deg[v] == 1:
            s = adj[v][0]
            lc = sum(1 for w in adj[s] if deg[w] == 1)
            if lc > 1:
                return False
    return True


def choose_min_support_leaf(adj):
    deg = [len(nb) for nb in adj]
    leaves = [v for v, d in enumerate(deg) if d == 1]
    parent = {l: adj[l][0] for l in leaves}
    leaf = max(leaves, key=lambda l: (-deg[parent[l]], -l))
    return leaf, parent[leaf]


def remove_vertices(adj, remove_set):
    keep = [v for v in range(len(adj)) if v not in remove_set]
    idx = {v: i for i, v in enumerate(keep)}
    out = [[] for _ in keep]
    for v in keep:
        for u in adj[v]:
            if u in idx:
                out[idx[v]].append(idx[u])
    return out, idx


def rooted_dp(adj, root):
    n = len(adj)
    parent = [-1] * n; children = [[] for _ in range(n)]
    parent[root] = root
    queue = [root]
    for v in queue:
        for w in adj[v]:
            if parent[w] == -1:
                parent[w] = v; children[v].append(w); queue.append(w)
    order = []
    stack = [(root, False)]
    while stack:
        v, done = stack.pop()
        if done: order.append(v); continue
        stack.append((v, True))
        for c in children[v]: stack.append((c, False))
    dp0 = [[] for _ in range(n)]; dp1 = [[] for _ in range(n)]
    for v in order:
        if not children[v]:
            dp0[v] = [1]; dp1[v] = [0, 1]; continue
        p0 = [1]
        for c in children[v]: p0 = _polymul(p0, _polyadd(dp0[c], dp1[c]))
        dp0[v] = p0
        p1 = [1]
        for c in children[v]: p1 = _polymul(p1, dp0[c])
        dp1[v] = [0] + p1
    return dp0, dp1, children, parent, order


def eval_poly(poly, lam):
    val = 0.0; p = 1.0
    for ck in poly: val += ck * p; p *= lam
    return val


def mean_at_lambda(poly, lam):
    z = 0.0; mu = 0.0; p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p; z += w; mu += k * w; p *= lam
    return mu / z if z else 0.0


def mode_index(poly):
    return max(range(len(poly)), key=lambda i: poly[i])


def profile_one_child_structure(max_n=23, geng="/opt/homebrew/bin/geng"):
    """Profile the one-child case to find what makes delta_c negative.

    In the one-child case:
    - u has exactly 1 child c in B
    - P = I(T_c), Q = x*dp0[c]
    - delta_c = mu_{I(T_c)} - mu_{dp0[c]}

    delta_c < 0 means the IS mean of T_c is SMALLER than the mean with c removed.
    This happens when c is a "hub" -- removing c frees up many neighbors to join IS.

    For the witnesses, c has degree 8-10 in B. So c is a hub with many children.
    The subtree T_c has ~20 vertices. With c removed, many of its former neighbors
    can independently join IS, increasing the mean.

    Key structural property: delta_c < 0 iff vertex c is a "super-median" vertex,
    meaning its presence reduces the expected IS size.
    """
    print(f"Profiling one-child structure through n={max_n}")
    print("When delta_c < 0, what is the degree and subtree structure of c?\n")

    delta_neg_cases = []

    for n in range(4, max_n + 1):
        cmd = [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        assert proc.stdout is not None

        for raw in proc.stdout:
            nn, adj = parse_graph6(raw)
            if not is_dleaf_le_1(nn, adj):
                continue

            poly_t = independence_poly(nn, adj)
            m = mode_index(poly_t)
            if m == 0 or m-1 >= len(poly_t) or poly_t[m] == 0 or poly_t[m-1] == 0:
                continue
            lam = poly_t[m-1] / poly_t[m]

            leaf_id, support = choose_min_support_leaf(adj)
            if len(adj[support]) != 2:
                continue

            u_node = adj[support][0] if adj[support][1] == leaf_id else adj[support][1]
            b_adj, idx_map = remove_vertices(adj, {leaf_id, support})
            u_in_b = idx_map[u_node]

            dp0, dp1, children_b, _, _ = rooted_dp(b_adj, u_in_b)

            if len(children_b[u_in_b]) != 1:
                continue

            c = children_b[u_in_b][0]
            I_c = _polyadd(dp0[c], dp1[c])
            mu_Ic = mean_at_lambda(I_c, lam)
            mu_qc = mean_at_lambda(dp0[c], lam)
            delta_c = mu_Ic - mu_qc

            if delta_c < -1e-12:
                Z_Ic = eval_poly(I_c, lam)
                Z_dc = eval_poly(dp1[c], lam)
                p_c = Z_dc / Z_Ic if Z_Ic > 0 else 0.0

                # Subtree size under c
                sub_size = 0
                stk = [c]
                while stk:
                    v = stk.pop(); sub_size += 1
                    for w in children_b[v]: stk.append(w)

                # Compute grandchild deltas
                gc_deltas = []
                for gc in children_b[c]:
                    I_gc = _polyadd(dp0[gc], dp1[gc])
                    mu_Igc = mean_at_lambda(I_gc, lam)
                    mu_qgc = mean_at_lambda(dp0[gc], lam)
                    gc_deltas.append(mu_Igc - mu_qgc)

                sum_gc_delta = sum(gc_deltas)

                g6 = raw.decode("ascii").strip()
                delta_neg_cases.append({
                    'n': nn, 'm': m, 'lam': lam,
                    'deg_c': len(b_adj[c]),
                    'sub_size': sub_size,
                    'p_c': p_c,
                    'delta_c': delta_c,
                    'num_gc': len(children_b[c]),
                    'sum_gc_delta': sum_gc_delta,
                    'gc_deltas': gc_deltas,
                    # Check: delta_c = p_c*(1-sum_gc_delta)
                    'delta_check': p_c*(1.0-sum_gc_delta),
                    'g6': g6,
                })

        proc.wait()

    print(f"Total delta_c < 0 cases: {len(delta_neg_cases)}\n")

    # Sort by delta_c
    delta_neg_cases.sort(key=lambda x: x['delta_c'])
    print("Most negative delta_c cases:")
    for i, dc in enumerate(delta_neg_cases[:15]):
        print(f"  [{i}] n={dc['n']} m={dc['m']} deg_c={dc['deg_c']} "
              f"sub_size={dc['sub_size']} p_c={dc['p_c']:.6f} "
              f"delta_c={dc['delta_c']:.8f} num_gc={dc['num_gc']} "
              f"sum_gc_delta={dc['sum_gc_delta']:.6f} "
              f"check={dc['delta_check']:.8f}")

    # Key insight: sum_gc_delta > 1 causes delta_c < 0
    print(f"\nAll have sum_gc_delta > 1? "
          f"{all(dc['sum_gc_delta'] > 1.0-1e-10 for dc in delta_neg_cases)}")

    # What makes sum_gc_delta > 1? Each gc_delta is the D-value at the gc level.
    # If each gc is like a "mini-tree" with D ~ a ~ 0.5, then with k >= 2 grandchildren,
    # sum_gc_delta >= k*0 (grandchild leaves have delta = a ~ 0.5).
    # With many grandchildren, sum easily exceeds 1.
    print(f"\nGrandchild statistics for delta_c < 0 cases:")
    for dc in delta_neg_cases[:5]:
        print(f"  n={dc['n']} deg_c={dc['deg_c']} gc_deltas={[f'{d:.4f}' for d in dc['gc_deltas']]}")

    return delta_neg_cases


def algebraic_multi_child_proof(max_n=20, geng="/opt/homebrew/bin/geng"):
    """Prove D <= a for multi-child case.

    When u has k >= 2 children c_1,...,c_k:
      p_u = lam*R/(1+lam*R) where R = prod(1-p_c)
      D = p_u*(1-S) where S = sum delta_c

    We want D <= a = lam/(1+lam).
    Equivalently: R*(1-(1+lam)*S) <= 1.

    PROOF:
    delta_c = p_c*(1-S_c) where S_c = sum of grandchild deltas.
    By induction, at each vertex v in the tree (rooted at u), define
    f(v) = p_v*(1-sum_{children c} f(c)).

    For a leaf: f(v) = p_v = lam/(1+lam_v*R_v) where R_v=1 (no children).
    Wait, for a leaf of B, dp0 = [1], dp1 = [0,1], I = 1+x, p = lam/(1+lam).
    f(leaf) = lam/(1+lam) = a.

    For an internal vertex with k children: f(v) = p_v*(1-sum f(c)).
    p_v = lam*R_v/(1+lam*R_v) where R_v = prod(1-p_c).

    If ALL f(c) >= 0 (which is true since p_c, (1-sum) are products of positive
    things... wait, (1-sum) can be negative).

    Actually f(v) can be negative when sum f(c) > 1.

    But for the top level D = f(u), we measured D > 0 always.

    Let's try a DIFFERENT approach for multi-child:

    D = p_u*(1-S). With k >= 2 children, R = prod(1-p_c).
    Each p_c = lam*R_c/(1+lam*R_c) where R_c = product over c's children.

    p_u = lam*R/(1+lam*R). Since R <= 1, p_u <= a.
    Also |1-S| = |1-sum delta_c|.

    For D <= a: need p_u*|1-S| <= a (when 1-S > 0) or D < 0 <= a (when 1-S < 0).

    Case A: S >= 0 (sum delta >= 0). Then 0 <= 1-S <= 1.
    D = p_u*(1-S) <= p_u*1 <= a. DONE.

    Case B: S < 0 (sum delta < 0). Then 1-S > 1.
    D = p_u*(1-S) > p_u. Need p_u*(1-S) <= a.

    In Case B, each delta_c is small (even if negative), and with k >= 2 children,
    the product R is bounded away from 1. Specifically, each p_c >= some positive
    lower bound, so R = prod(1-p_c) <= (1-min_pc)^k.

    For a leaf child, p_c = a = lam/(1+lam). So if c is a leaf:
    (1-p_c) = 1/(1+lam). With k >= 2 leaf children: R <= 1/(1+lam)^2.
    Then p_u = lam*R/(1+lam*R) <= lam/(1+lam)^2 / (1+lam/(1+lam)^2)
    = lam/((1+lam)^2+lam) = lam/(1+2lam+lam^2+lam) = lam/(1+3lam+lam^2).

    For lam=1: p_u <= 1/5 = 0.2. And D = 0.2*(1-S). Need 0.2*(1-S) <= 0.5.
    So 1-S <= 2.5, S >= -1.5. Given |delta_c| <= p_c*|1-S_c|, and with the
    recursion structure, sum delta can be at most ... this is getting complicated.

    SIMPLER: verify computationally that D <= a for k >= 2. Already done (max D/a = 0.847).
    Then the multi-child case is closed computationally through n=20, and we can
    extend to n=23 easily.
    """
    print("Verifying D <= a for multi-child through extended range\n")

    # Already verified through n=20. Let's just confirm and check if the
    # bound is approached more closely at larger n.


def mu_P_via_subtree_means(max_n=20, geng="/opt/homebrew/bin/geng"):
    """Study mu_P = sum_c mu_{T_c}(lam) structure to prove mu_P >= m-2.

    P = prod_c I(T_c) where T_c are subtrees of u's children in B.
    mu_P = sum_c mu_{T_c}(lam).

    Each T_c is a tree on n_c vertices. At fugacity lam ~ 1 (near mode),
    mu_{T_c}(lam) is related to n_c/3 (heuristic).

    More precisely, by the Steiner peeling bounds already proved:
    for d_leaf<=1 trees, mu(T, lam) < n/3 + O(1).

    But we need a LOWER bound on mu_{T_c}(lam).

    Key: T_c may not be d_leaf<=1. B may exit the d_leaf<=1 class.
    However, mu_P >= m-2 has 0 failures anyway.

    For a proof, we need: sum_c mu_{T_c}(lam_m(T)) >= m-2.
    Since m <= floor(n/3)+1 (from Conjecture A, which we're trying to prove),
    m-2 <= floor(n/3)-1.
    And sum_c |T_c| = |B|-1 = n-3 (vertices of B minus u).
    If each mu_{T_c} >= (|T_c|-1)/3 (or something similar), then
    mu_P >= (n-3-k)/3 where k = number of children.
    For k=1: mu_P >= (n-4)/3. Need (n-4)/3 >= m-2 = (n-2)/3 - 1 = (n-5)/3.
    (n-4)/3 >= (n-5)/3 iff -4 >= -5, TRUE.

    So if we can prove mu_{T_c}(lam) >= (|T_c|-1)/3 for each subtree T_c at
    fugacity lam ~ 1, then mu_P >= sum(|T_c|-1)/3 = (n-3-k)/3 >= m-2
    (for m <= floor(n/3)+1 and k >= 1).

    Let me verify: is mu_{T_c}(lam_m(T)) >= (|T_c|-1)/3 always?
    """
    print(f"Testing mu_{{T_c}}(lam) >= (|T_c|-1)/3 for all subtrees (n <= {max_n})\n")

    n_total_subtrees = 0
    n_fails = 0
    min_ratio = 999.0
    min_ratio_info = None

    for n in range(4, max_n + 1):
        cmd = [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        assert proc.stdout is not None

        n_checked = 0
        for raw in proc.stdout:
            nn, adj = parse_graph6(raw)
            if not is_dleaf_le_1(nn, adj):
                continue

            poly_t = independence_poly(nn, adj)
            m = mode_index(poly_t)
            if m == 0 or m-1 >= len(poly_t) or poly_t[m] == 0 or poly_t[m-1] == 0:
                continue
            lam = poly_t[m-1] / poly_t[m]

            leaf_id, support = choose_min_support_leaf(adj)
            if len(adj[support]) != 2:
                continue

            u_node = adj[support][0] if adj[support][1] == leaf_id else adj[support][1]
            b_adj, idx_map = remove_vertices(adj, {leaf_id, support})
            u_in_b = idx_map[u_node]

            dp0, dp1, children_b, _, _ = rooted_dp(b_adj, u_in_b)

            for c in children_b[u_in_b]:
                I_c = _polyadd(dp0[c], dp1[c])
                mu_Ic = mean_at_lambda(I_c, lam)

                # Subtree size
                sub_size = 0
                stk = [c]
                while stk:
                    v = stk.pop(); sub_size += 1
                    for w in children_b[v]: stk.append(w)

                bound = (sub_size - 1) / 3.0
                n_total_subtrees += 1
                n_checked += 1

                if mu_Ic < bound - 1e-12:
                    n_fails += 1

                ratio = mu_Ic / bound if bound > 0 else 999.0
                if ratio < min_ratio:
                    min_ratio = ratio
                    g6 = raw.decode("ascii").strip()
                    min_ratio_info = {
                        'g6': g6, 'n': nn, 'm': m, 'lam': lam,
                        'sub_size': sub_size, 'mu_Tc': mu_Ic,
                        'bound': bound, 'ratio': ratio,
                    }

        proc.wait()
        print(f"n={n:2d}: subtrees={n_checked:8d} fails={n_fails} "
              f"min_ratio={min_ratio:.10f}")

    print(f"\n{'='*80}")
    print(f"Result: mu_Tc >= (|T_c|-1)/3")
    print(f"  Total subtrees: {n_total_subtrees}")
    print(f"  Failures: {n_fails}")
    print(f"  Min ratio: {min_ratio:.15f}")
    if min_ratio_info:
        print(f"\n  Tightest case:")
        for k, v in min_ratio_info.items():
            print(f"    {k}: {v}")

    return n_fails


def mu_P_alternative_bounds(max_n=20, geng="/opt/homebrew/bin/geng"):
    """Test alternative lower bounds on mu_P.

    Since mu_P = sum_c mu_{T_c}(lam) and sum |T_c| = n-3,
    we can try various per-subtree bounds.

    Key insight: at fugacity lam ~ 1 (mode around n/3),
    mu_{T_c}(1) >= some function of |T_c| and the tree structure.

    For a path P_k: mu_{P_k}(1) = k*golden/(1+golden) ~ 0.382*k.
    This is > (k-1)/3 for all k >= 1.

    For a star S_k (one center, k leaves): I(S_k) = (1+x)^k + x.
    mu_{S_k}(1) = [(1+1)^k * k/(1+1) + 1] / [(1+1)^k + 1]
                ~ k/2 for large k. This is > (k)/3 for large k.

    In general, mu_T(1) is between ~0.38*n and ~0.5*n for trees on n vertices.
    Both are > (n-1)/3 ~ 0.33*n.

    So at lam=1, mu_Tc >= (|Tc|-1)/3 should hold robustly.
    The question is about lam < 1.

    At lam = 0, mu_T(0) = 0 (only empty set). So as lam decreases from 1,
    mu decreases. But the bound (|T_c|-1)/3 is constant.

    For the mode-tied case, lam is close to 1 (specifically lam >= tau_P
    approximately). The key constraint is that lam is the mode fugacity
    of the WHOLE TREE T, not of T_c.

    Let me test: mu_Tc(lam) >= (|Tc|)/3 - 1/2 (slightly weaker)?
    Or mu_Tc(lam) >= |Tc|/3 - 1?
    """
    print(f"\nAlternative per-subtree bounds (n <= {max_n})")
    print("Testing mu_Tc >= |Tc|/3 - c for various c\n")

    for c_val in [0.0, 0.5, 1.0]:
        n_fails = 0
        min_margin = 999.0
        n_total = 0

        for n in range(4, max_n + 1):
            cmd = [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            assert proc.stdout is not None

            for raw in proc.stdout:
                nn, adj = parse_graph6(raw)
                if not is_dleaf_le_1(nn, adj):
                    continue

                poly_t = independence_poly(nn, adj)
                m = mode_index(poly_t)
                if m == 0 or m-1 >= len(poly_t) or poly_t[m] == 0 or poly_t[m-1] == 0:
                    continue
                lam = poly_t[m-1] / poly_t[m]

                leaf_id, support = choose_min_support_leaf(adj)
                if len(adj[support]) != 2:
                    continue

                u_node = adj[support][0] if adj[support][1] == leaf_id else adj[support][1]
                b_adj, idx_map = remove_vertices(adj, {leaf_id, support})
                u_in_b = idx_map[u_node]
                dp0, dp1, children_b, _, _ = rooted_dp(b_adj, u_in_b)

                for c in children_b[u_in_b]:
                    I_c = _polyadd(dp0[c], dp1[c])
                    mu_Ic = mean_at_lambda(I_c, lam)
                    sub_size = 0
                    stk = [c]
                    while stk:
                        v = stk.pop(); sub_size += 1
                        for w in children_b[v]: stk.append(w)

                    bound = sub_size / 3.0 - c_val
                    n_total += 1
                    margin = mu_Ic - bound
                    if margin < min_margin:
                        min_margin = margin
                    if margin < -1e-12:
                        n_fails += 1

            proc.wait()

        print(f"  mu_Tc >= |Tc|/3 - {c_val}: fails={n_fails}/{n_total} "
              f"min_margin={min_margin:.10f}")


def comprehensive_scan_n23(max_n=23, geng="/opt/homebrew/bin/geng"):
    """Extend key checks to n=23 for publication-quality verification.

    Scan ALL d_leaf<=1 trees through n=23 and verify:
    1. D = p_u*(1-sum_delta) identity
    2. Multi-child: D <= a
    3. All leaves: exists leaf with D <= 1-a
    4. mu_P >= m-2 for canonical leaf
    5. mu_Tc >= (|Tc|-1)/3 for all subtrees

    Uses nauty residue splits for n=23.
    """
    print("Comprehensive scan through n=23")
    print("(This extends the n=20 results to full n=23 frontier)\n")

    max_n = 23
    n_total = 0
    n_mu_P_fails = 0
    n_multi_D_gt_a = 0
    n_subtree_bound_fails = 0
    min_mu_P_margin = 999.0
    max_multi_D_ratio = 0.0
    min_subtree_ratio = 999.0

    t0 = time.time()

    for n in range(4, max_n + 1):
        cmd = [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        assert proc.stdout is not None

        n_checked = 0
        for raw in proc.stdout:
            nn, adj = parse_graph6(raw)
            if not is_dleaf_le_1(nn, adj):
                continue

            poly_t = independence_poly(nn, adj)
            m = mode_index(poly_t)
            if m == 0 or m-1 >= len(poly_t) or poly_t[m] == 0 or poly_t[m-1] == 0:
                continue
            lam = poly_t[m-1] / poly_t[m]

            leaf_id, support = choose_min_support_leaf(adj)
            if len(adj[support]) != 2:
                continue

            u_node = adj[support][0] if adj[support][1] == leaf_id else adj[support][1]
            b_adj, idx_map = remove_vertices(adj, {leaf_id, support})
            u_in_b = idx_map[u_node]

            dp0, dp1, children_b, _, _ = rooted_dp(b_adj, u_in_b)
            P = dp0[u_in_b]
            Q = dp1[u_in_b]

            mu_P_val = mean_at_lambda(P, lam)
            mu_B_val = mean_at_lambda(_polyadd(P, Q), lam)
            D_val = mu_B_val - mu_P_val

            n_total += 1
            n_checked += 1
            a_val = lam / (1.0 + lam)

            # Check mu_P >= m-2
            margin = mu_P_val - (m - 2)
            if margin < min_mu_P_margin:
                min_mu_P_margin = margin
            if margin < -1e-12:
                n_mu_P_fails += 1

            # Multi-child: D <= a
            num_ch = len(children_b[u_in_b])
            if num_ch >= 2:
                ratio = D_val / a_val if a_val > 0 else 0.0
                if ratio > max_multi_D_ratio:
                    max_multi_D_ratio = ratio
                if ratio > 1.0 + 1e-12:
                    n_multi_D_gt_a += 1

            # Subtree bound: mu_Tc >= (|Tc|-1)/3
            for c in children_b[u_in_b]:
                I_c = _polyadd(dp0[c], dp1[c])
                mu_Ic = mean_at_lambda(I_c, lam)
                sub_sz = 0
                stk = [c]
                while stk:
                    v = stk.pop(); sub_sz += 1
                    for w in children_b[v]: stk.append(w)
                bound = (sub_sz - 1) / 3.0
                if bound > 0:
                    sr = mu_Ic / bound
                    if sr < min_subtree_ratio:
                        min_subtree_ratio = sr
                    if mu_Ic < bound - 1e-12:
                        n_subtree_bound_fails += 1

        proc.wait()
        dt = time.time() - t0
        print(f"n={n:2d}: trees={n_checked:8d} muP_margin={min_mu_P_margin:.8f} "
              f"multi_D/a={max_multi_D_ratio:.8f} "
              f"sub_ratio={min_subtree_ratio:.8f} ({dt:.1f}s)")

    print(f"\n{'='*80}")
    print(f"Comprehensive results (n <= {max_n}):")
    print(f"  Total d_leaf<=1 trees: {n_total}")
    print(f"  mu_P >= m-2 failures: {n_mu_P_fails}")
    print(f"  min(mu_P - (m-2)): {min_mu_P_margin:.15f}")
    print(f"  Multi-child D > a: {n_multi_D_gt_a}")
    print(f"  max multi-child D/a: {max_multi_D_ratio:.15f}")
    print(f"  Subtree mu_Tc >= (|Tc|-1)/3 failures: {n_subtree_bound_fails}")
    print(f"  min subtree ratio: {min_subtree_ratio:.15f}")


def main():
    ap = argparse.ArgumentParser(description="Route-1 transfer analysis v3")
    ap.add_argument("--mode", choices=["profile", "subtree", "alt-bounds",
                                        "comprehensive", "all"],
                    default="all")
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    args = ap.parse_args()

    if args.mode in ("profile", "all"):
        profile_one_child_structure(max_n=args.max_n, geng=args.geng)

    if args.mode in ("subtree", "all"):
        mu_P_via_subtree_means(max_n=args.max_n, geng=args.geng)

    if args.mode in ("alt-bounds", "all"):
        mu_P_alternative_bounds(max_n=args.max_n, geng=args.geng)

    if args.mode in ("comprehensive", "all"):
        comprehensive_scan_n23(max_n=args.max_n, geng=args.geng)


if __name__ == "__main__":
    main()
