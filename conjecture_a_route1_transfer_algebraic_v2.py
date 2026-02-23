#!/usr/bin/env python3
"""Route-1 transfer: tight algebraic analysis, v2.

This script focuses on proving the correct bound on D = mu_B - mu_P
via the identity D = p_u * (1 - sum_c delta_c), where:
  p_u = lam*R/(1+lam*R), R = prod_c(1-p_c)
  delta_c = p_c * (1 - sum_{gc} delta_{gc})  (recursive)

The bound D <= 1/(1+lam) is equivalent to:
  lam*R*[lam - (1+lam)*S] <= 1,  where S = sum_c delta_c.

Key structural finding from v1:
  - The 4 failure witnesses all have EXACTLY 1 CHILD of u in B.
  - delta_c < 0 for that child (negative delta => 1-S > 1).
  - deg(u) in B = 1 (u becomes a leaf in B).

When u has exactly 1 child c in B:
  P = I(T_c)  (the IS poly of the subtree)
  Q = x * dp0[c]  (x times the c-excluded poly)
  R = 1-p_c
  S = delta_c
  p_u = lam*(1-p_c)/(1+lam*(1-p_c))
  D = p_u*(1-delta_c)

For delta_c < 0 to happen, we need sum_{gc} delta_{gc} > 1 at the child level.

ALTERNATIVE APPROACH: Direct algebraic bound using the Φ structure.

Instead of proving D <= 1/(1+lam), we can work directly with the proven facts:
  Route 2: mu_B >= m - 1 - a  (0 failures through n=23)
  Need:    mu_P >= m - 2

From mu_B = mu_P + D:
  mu_P = mu_B - D >= (m-1-a) - D
  Need: (m-1-a) - D >= m-2, i.e. D <= 1-a.

D = p_u*(1-S) where S can be negative.
We need: p_u*(1-S) <= 1-a = 1/(1+lam).

With p_u = lam*R/(1+lam*R):
  lam*R*(1-S)/(1+lam*R) <= 1/(1+lam)
  lam*(1+lam)*R*(1-S) <= 1+lam*R

Let Psi = lam*(1+lam)*R*(1-S) - (1+lam*R).

Psi = lam*(1+lam)*R - lam*(1+lam)*R*S - 1 - lam*R
    = lam*R*[(1+lam) - (1+lam)*S - 1]  -  1
    = lam*R*[lam - (1+lam)*S]  -  1

So Psi <= 0 iff lam*R*[lam - (1+lam)*S] <= 1.

NEW IDEA: Work in terms of the Φ decomposition directly.

We know Φ_m(T) = Φ_m(A) + lam*Φ_{m-1}(B).
And Φ_m(A) = (1+lam)*Φ_m(B) + lam*Z_B + lam*Φ_{m-1}(P).

So Φ_m(T) = (1+lam)*Φ_m(B) + lam*Z_B + lam*Φ_{m-1}(P) + lam*Φ_{m-1}(B).

We want Φ_{m-1}(P) >= 0, which means mu_P >= m-2.

What if instead, we don't prove mu_P >= m-2 at all, but note that
we already have Route 2 (mu_B >= m-1-a) AND the pendant bonus gives
Φ_m(A) >= 0, AND Φ_{m-1}(B) >= 0, so Φ_m(T) >= 0 WITHOUT needing Route 1?

Actually, the whole point of Route 1 was as an alternative path. The
decomposition Φ_m(T) = Φ_m(A) + lam*Φ_{m-1}(B) already works if BOTH
terms are non-negative. And both ARE non-negative (0 failures through n=23).

So perhaps the transfer analysis is not needed at all -- the direct
decomposition already proves Conjecture A if we can prove:
  (a) Φ_m(A) >= 0
  (b) Φ_{m-1}(B) >= 0

independently. And Route 2 (mu_B >= m-1-a) implies Φ_m(A) >= 0 via the
pendant bonus, while Φ_{m-1}(B) >= 0 needs its own argument.

But for completeness, let me find the EXACT tight bound on D.

THEOREM ATTEMPT: D <= a + a*(1-a)*epsilon for some small epsilon.

Actually, let's just characterize when D > 1/(1+lam) and see if those cases
are covered by having extra slack in Route 2.
"""

from __future__ import annotations

import argparse
import subprocess
import time
from fractions import Fraction
from typing import Any

from graph6 import parse_graph6
from indpoly import _polyadd, _polymul, independence_poly


def is_dleaf_le_1(n: int, adj: list[list[int]]) -> bool:
    deg = [len(nb) for nb in adj]
    for v in range(n):
        if deg[v] == 1:
            s = adj[v][0]
            leaf_count = sum(1 for w in adj[s] if deg[w] == 1)
            if leaf_count > 1:
                return False
    return True


def choose_min_support_leaf(adj: list[list[int]]) -> tuple[int, int]:
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
        vv = idx[v]
        for u in adj[v]:
            if u in idx:
                out[vv].append(idx[u])
    return out, idx


def rooted_dp(adj, root):
    n = len(adj)
    parent = [-1] * n
    children = [[] for _ in range(n)]
    parent[root] = root
    queue = [root]
    for v in queue:
        for w in adj[v]:
            if parent[w] == -1:
                parent[w] = v
                children[v].append(w)
                queue.append(w)
    order = []
    stack = [(root, False)]
    while stack:
        v, done = stack.pop()
        if done:
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
            continue
        p0 = [1]
        for c in children[v]:
            p0 = _polymul(p0, _polyadd(dp0[c], dp1[c]))
        dp0[v] = p0
        p1 = [1]
        for c in children[v]:
            p1 = _polymul(p1, dp0[c])
        dp1[v] = [0] + p1
    return dp0, dp1, children, parent, order


def eval_poly(poly, lam):
    val = 0.0
    p = 1.0
    for ck in poly:
        val += ck * p
        p *= lam
    return val


def mean_at_lambda(poly, lam):
    z = 0.0
    mu_num = 0.0
    p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p
        z += w
        mu_num += k * w
        p *= lam
    return mu_num / z if z else 0.0


def mode_index(poly):
    return max(range(len(poly)), key=lambda i: poly[i])


def combined_analysis(max_n: int = 23, geng: str = "/opt/homebrew/bin/geng"):
    """Check the COMBINED condition: does Route-2 slack cover the D excess?

    Route-2 gives mu_B >= m-1-a with slack:
      slack_route2 = mu_B - (m-1-a)

    Route-1 needs mu_P >= m-2, which requires:
      D <= mu_B - (m-2) = (m-1-a + slack) - (m-2) = 1-a + slack

    So Route-1 fails only when D > 1-a + slack_route2, i.e.:
      D - (1-a) > slack_route2

    We check: is exact_excess := D - (1-a) always <= slack_route2?
    I.e., does mu_P actually end up >= m-2 even when D > 1-a?
    """
    print(f"Combined analysis: checking mu_P >= m-2 directly (n <= {max_n})")
    print("Also tracking D vs Route-2 slack\n")

    n_total = 0
    n_mu_P_fail = 0
    min_mu_P_margin = 999.0
    max_D = 0.0
    max_excess = -999.0

    # Track: when exact_excess > 0, how does Route-2 slack compare?
    n_D_gt_1_minus_a = 0
    n_covered_by_slack = 0

    # One-child vs multi-child statistics
    one_child_max_D = 0.0
    multi_child_max_D = 0.0

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
            if m == 0 or m - 1 >= len(poly_t) or poly_t[m] == 0 or poly_t[m - 1] == 0:
                continue
            lam = poly_t[m - 1] / poly_t[m]

            leaf, support = choose_min_support_leaf(adj)
            if len(adj[support]) != 2:
                continue

            u_node = adj[support][0] if adj[support][1] == leaf else adj[support][1]
            b_adj, idx_map = remove_vertices(adj, {leaf, support})
            u_in_b = idx_map[u_node]

            dp0, dp1, children_b, _, _ = rooted_dp(b_adj, u_in_b)
            P = dp0[u_in_b]
            Q = dp1[u_in_b]
            B_poly = _polyadd(P, Q)

            mu_P_val = mean_at_lambda(P, lam)
            mu_B_val = mean_at_lambda(B_poly, lam)
            D_val = mu_B_val - mu_P_val

            a_val = lam / (1.0 + lam)
            threshold = 1.0 - a_val
            excess = D_val - threshold

            # Route-2 slack
            route2_target = m - 1 - a_val
            route2_slack = mu_B_val - route2_target

            # mu_P margin
            mu_P_margin = mu_P_val - (m - 2)

            n_total += 1
            n_checked += 1

            if mu_P_margin < min_mu_P_margin:
                min_mu_P_margin = mu_P_margin

            if mu_P_margin < -1e-10:
                n_mu_P_fail += 1

            if D_val > max_D:
                max_D = D_val

            if excess > max_excess:
                max_excess = excess

            if excess > 1e-12:
                n_D_gt_1_minus_a += 1
                # Check if Route-2 slack covers the excess
                if route2_slack >= excess - 1e-12:
                    n_covered_by_slack += 1

            num_children = len(children_b[u_in_b])
            if num_children == 1:
                if D_val > one_child_max_D:
                    one_child_max_D = D_val
            else:
                if D_val > multi_child_max_D:
                    multi_child_max_D = D_val

        proc.wait()
        dt = time.time() - t0
        print(f"n={n:2d}: trees={n_checked:8d} "
              f"mu_P_fails={n_mu_P_fail} min_muP_margin={min_mu_P_margin:.10f} "
              f"D_excess={n_D_gt_1_minus_a} ({dt:.1f}s)")

    print(f"\n{'='*80}")
    print(f"SUMMARY (n <= {max_n}):")
    print(f"  Total d_leaf<=1 trees: {n_total}")
    print(f"  mu_P < m-2 failures: {n_mu_P_fail}")
    print(f"  min(mu_P - (m-2)): {min_mu_P_margin:.15f}")
    print(f"  max D: {max_D:.15f}")
    print(f"  max D excess over 1/(1+lam): {max_excess:.15f}")
    print(f"  D > 1/(1+lam) cases: {n_D_gt_1_minus_a}")
    print(f"  Covered by Route-2 slack: {n_covered_by_slack}")
    print(f"  One-child max D: {one_child_max_D:.15f}")
    print(f"  Multi-child max D: {multi_child_max_D:.15f}")


def one_child_detailed(max_n: int = 23, geng: str = "/opt/homebrew/bin/geng"):
    """Detailed analysis of the one-child case (u has exactly 1 child c in B).

    In this case:
      P = I(T_c), Q = x * dp0[c]
      R = 1 - p_c
      p_u = lam*(1-p_c) / (1 + lam*(1-p_c))
      D = p_u * (1 - delta_c) where delta_c = mu_{I_c} - mu_{dp0[c]}

    The child c has a potentially large subtree. delta_c depends on the
    recursive structure of that subtree.

    For the 4 failure witnesses:
      p_c is very small (~0.02-0.03): c has many neighbors, making it unlikely
      to be in the IS.
      delta_c is negative (~-0.03): removing c from its subtree INCREASES the
      mean IS size (because c is a hub that blocks many neighbors).

    When p_c is small and delta_c is negative:
      R = 1-p_c ~ 1
      p_u ~ lam/(1+lam) ~ a ~ 0.5
      1 - delta_c ~ 1 + |delta_c| > 1
      D ~ 0.5 * (1 + |delta_c|) > 0.5

    The question: can p_u*(1-delta_c) exceed 1/(1+lam)?

    Let's rewrite using exact quantities:
      p_c = lam*R_c / (1+lam*R_c) where R_c = product over c's children
      delta_c = p_c*(1-sum_{gc} delta_{gc})
      p_u = lam*(1-p_c)/(1+lam*(1-p_c))

    Since we know R_c through the child subtree, and R = 1-p_c = 1/(1+lam*R_c),
    we get:
      p_u = lam/(1+lam*R_c) / (1 + lam/(1+lam*R_c))
          = lam / (1+lam*R_c + lam)
          = lam / (1+lam*(1+R_c))

    And D = p_u*(1-delta_c) = [lam/(1+lam+lam*R_c)] * (1 - delta_c).

    For D <= 1/(1+lam), we need:
      lam*(1+lam)*(1-delta_c) <= 1+lam+lam*R_c
      lam*(1+lam) - lam*(1+lam)*delta_c <= 1+lam+lam*R_c
      lam*lam - lam*(1+lam)*delta_c <= 1+lam*R_c
      lam^2 - lam*(1+lam)*delta_c <= 1+lam*R_c

    With delta_c = mu_{I_c} - mu_{dp0[c]}:
    If delta_c >= 0, LHS <= lam^2 < 1 <= 1+lam*R_c. DONE.
    If delta_c < 0, LHS = lam^2 + lam*(1+lam)*|delta_c|.

    Need: lam^2 + lam*(1+lam)*|delta_c| <= 1+lam*R_c
    i.e. lam*(1+lam)*|delta_c| <= 1-lam^2 + lam*R_c = (1-lam)(1+lam) + lam*R_c
    i.e. lam*|delta_c| <= (1-lam) + lam*R_c/(1+lam)

    For lam close to 1: (1-lam) ~ 0, and we need:
      lam*|delta_c| <= lam*R_c/(1+lam) (approximately)
      |delta_c| <= R_c/(1+lam) ~ R_c/2

    Now delta_c = p_c*(1-S_c) where S_c = sum of grandchild deltas.
    |delta_c| = p_c*|1-S_c|.
    If S_c > 1 (sum of grandchild deltas > 1), then |delta_c| = p_c*(S_c-1).

    And R_c = prod_{gc}(1-p_{gc}).

    So the condition |delta_c| <= R_c/(1+lam) becomes:
      p_c*(S_c-1) <= R_c/(1+lam)

    But p_c = lam*R_c/(1+lam*R_c), so:
      lam*R_c*(S_c-1)/(1+lam*R_c) <= R_c/(1+lam)
      lam*(S_c-1)/(1+lam*R_c) <= 1/(1+lam)
      lam*(1+lam)*(S_c-1) <= 1+lam*R_c

    This is getting recursive. Let me just scan for the exact statistics.
    """
    print(f"One-child detailed analysis (n <= {max_n})")
    print("Focus on u-degree-1 cases and their subtree structure\n")

    max_excess = -999.0
    max_excess_info = None
    n_one_child = 0
    n_multi_child = 0
    n_delta_neg = 0
    max_abs_delta = 0.0

    # Track the quantity: lam^2 + lam*(1+lam)*|delta_c| vs 1+lam*R_c
    max_Psi = -999.0  # Psi = lam^2 + lam*(1+lam)*|delta| - 1 - lam*R
    n_Psi_positive = 0

    t0 = time.time()

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
            if m == 0 or m - 1 >= len(poly_t) or poly_t[m] == 0 or poly_t[m - 1] == 0:
                continue
            lam = poly_t[m - 1] / poly_t[m]

            leaf, support = choose_min_support_leaf(adj)
            if len(adj[support]) != 2:
                continue

            u_node = adj[support][0] if adj[support][1] == leaf else adj[support][1]
            b_adj, idx_map = remove_vertices(adj, {leaf, support})
            u_in_b = idx_map[u_node]

            dp0, dp1, children_b, _, _ = rooted_dp(b_adj, u_in_b)

            num_children = len(children_b[u_in_b])
            if num_children != 1:
                n_multi_child += 1
                continue
            n_one_child += 1

            c = children_b[u_in_b][0]
            I_c = _polyadd(dp0[c], dp1[c])
            q_c = dp0[c]

            Z_Ic = eval_poly(I_c, lam)
            Z_qc = eval_poly(q_c, lam)
            Z_dc = eval_poly(dp1[c], lam)
            if Z_Ic == 0:
                continue
            p_c = Z_dc / Z_Ic

            mu_Ic = mean_at_lambda(I_c, lam)
            mu_qc = mean_at_lambda(q_c, lam)
            delta_c = mu_Ic - mu_qc

            if delta_c < -1e-12:
                n_delta_neg += 1
            if abs(delta_c) > max_abs_delta:
                max_abs_delta = abs(delta_c)

            R_c_prod = 1.0  # prod over c's grandchildren
            for gc in children_b[c]:
                I_gc = _polyadd(dp0[gc], dp1[gc])
                Z_Igc = eval_poly(I_gc, lam)
                Z_dgc = eval_poly(dp1[gc], lam)
                if Z_Igc > 0:
                    p_gc = Z_dgc / Z_Igc
                    R_c_prod *= (1.0 - p_gc)

            # The bound quantity for one-child case
            R = 1.0 - p_c
            P_poly = dp0[u_in_b]
            Q_poly = dp1[u_in_b]
            Z_P = eval_poly(P_poly, lam)
            Z_Q = eval_poly(Q_poly, lam)
            p_u = Z_Q / (Z_P + Z_Q) if (Z_P + Z_Q) > 0 else 0.0

            D_val = p_u * (1.0 - delta_c)
            a_val = lam / (1.0 + lam)
            threshold = 1.0 - a_val
            excess = D_val - threshold

            if excess > max_excess:
                max_excess = excess
                g6 = raw.decode("ascii").strip()
                max_excess_info = {
                    'g6': g6, 'n': nn, 'm': m, 'lam': lam,
                    'p_u': p_u, 'p_c': p_c, 'delta_c': delta_c,
                    'R': R, 'R_c_prod': R_c_prod, 'D': D_val,
                    'excess': excess, 'deg_c': len(b_adj[c]),
                    'num_gc': len(children_b[c]),
                }

            # Compute Psi = lam^2 + lam*(1+lam)*|delta_c| - 1 - lam*R
            if delta_c < 0:
                Psi = lam**2 + lam*(1+lam)*abs(delta_c) - 1 - lam*R
                if Psi > max_Psi:
                    max_Psi = Psi
                if Psi > 1e-12:
                    n_Psi_positive += 1

        proc.wait()
        dt = time.time() - t0
        print(f"n={n:2d}: one_child={n_one_child} multi={n_multi_child} "
              f"delta_neg={n_delta_neg} max_excess={max_excess:.10f} "
              f"max_Psi={max_Psi:.10f} ({dt:.1f}s)")

    print(f"\n{'='*80}")
    print(f"One-child summary:")
    print(f"  One-child trees: {n_one_child}")
    print(f"  Multi-child trees: {n_multi_child}")
    print(f"  delta_c < 0: {n_delta_neg}")
    print(f"  max |delta_c|: {max_abs_delta:.15f}")
    print(f"  max excess (D - 1/(1+lam)): {max_excess:.15f}")
    print(f"  max Psi (need <= 0): {max_Psi:.15f}")
    print(f"  Psi > 0 cases: {n_Psi_positive}")

    if max_excess_info:
        print(f"\n  Worst excess case:")
        for k, v in max_excess_info.items():
            print(f"    {k}: {v}")


def multi_child_analysis(max_n: int = 23, geng: str = "/opt/homebrew/bin/geng"):
    """For multi-child cases, verify D <= a (which is stronger than D <= 1-a).

    When u has >= 2 children:
      P = prod_c I(T_c), Q = x * prod_c dp0[c]
      R = prod_c(1-p_c) < 1 (with strict inequality if any p_c > 0)
      p_u = lam*R/(1+lam*R) < lam/(1+lam) = a

    For multi-child, p_u < a strictly. And 1-S may be > 1 (if S < 0),
    but p_u is bounded away from a.

    Question: is D = p_u*(1-S) <= a for multi-child always?
    """
    print(f"\nMulti-child analysis (n <= {max_n})")
    print("Testing D <= a = lam/(1+lam) for trees where u has >= 2 children\n")

    n_multi = 0
    max_D_over_a = 0.0
    max_D = 0.0

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
            if m == 0 or m - 1 >= len(poly_t) or poly_t[m] == 0 or poly_t[m - 1] == 0:
                continue
            lam = poly_t[m - 1] / poly_t[m]

            leaf, support = choose_min_support_leaf(adj)
            if len(adj[support]) != 2:
                continue

            u_node = adj[support][0] if adj[support][1] == leaf else adj[support][1]
            b_adj, idx_map = remove_vertices(adj, {leaf, support})
            u_in_b = idx_map[u_node]

            dp0, dp1, children_b, _, _ = rooted_dp(b_adj, u_in_b)

            if len(children_b[u_in_b]) < 2:
                continue
            n_multi += 1

            P = dp0[u_in_b]
            Q = dp1[u_in_b]
            mu_P_val = mean_at_lambda(P, lam)
            mu_B_val = mean_at_lambda(_polyadd(P, Q), lam)
            D_val = mu_B_val - mu_P_val

            a_val = lam / (1.0 + lam)
            if D_val > max_D:
                max_D = D_val
            ratio = D_val / a_val if a_val > 0 else 0.0
            if ratio > max_D_over_a:
                max_D_over_a = ratio

        proc.wait()

    print(f"  Multi-child trees: {n_multi}")
    print(f"  max D: {max_D:.15f}")
    print(f"  max D/a: {max_D_over_a:.15f} (need <= 1 for D <= a)")
    if max_D_over_a <= 1.0 + 1e-12:
        print(f"  CONFIRMED: D <= a for ALL multi-child trees through n={max_n}")
    else:
        print(f"  D > a for some multi-child trees!")


def all_leaves_scan(max_n: int = 23, geng: str = "/opt/homebrew/bin/geng"):
    """Instead of canonical leaf, check ALL degree-2-support leaves.

    For the 4 witnesses, the canonical leaf has 1 child of u in B.
    Maybe a different leaf choice gives a better (smaller) D?

    Also check: for the witnesses, is there ANY leaf choice with D <= 1-a?
    """
    print(f"\nAll-leaves scan (n <= {max_n})")
    print("For each tree, check ALL degree-2-support leaves and pick the one")
    print("that minimizes D. Then check if min D <= 1-a.\n")

    n_total = 0
    n_min_D_exceeds = 0  # trees where even best leaf gives D > 1-a
    max_min_D_excess = -999.0

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
            if m == 0 or m - 1 >= len(poly_t) or poly_t[m] == 0 or poly_t[m - 1] == 0:
                continue
            lam = poly_t[m - 1] / poly_t[m]

            deg = [len(nb) for nb in adj]
            leaves = [v for v, d in enumerate(deg) if d == 1]

            best_D = 999.0
            a_val = lam / (1.0 + lam)
            threshold = 1.0 - a_val

            for leaf in leaves:
                support = adj[leaf][0]
                if deg[support] != 2:
                    continue

                u_node = adj[support][0] if adj[support][1] == leaf else adj[support][1]
                b_adj, idx_map = remove_vertices(adj, {leaf, support})
                u_in_b = idx_map[u_node]

                dp0, dp1, _, _, _ = rooted_dp(b_adj, u_in_b)
                P = dp0[u_in_b]
                Q = dp1[u_in_b]
                mu_P_val = mean_at_lambda(P, lam)
                mu_B_val = mean_at_lambda(_polyadd(P, Q), lam)
                D_val = mu_B_val - mu_P_val

                if D_val < best_D:
                    best_D = D_val

            if best_D >= 999.0:
                continue  # no deg-2 support leaves?

            n_total += 1
            n_checked += 1
            excess = best_D - threshold
            if excess > 1e-12:
                n_min_D_exceeds += 1
            if excess > max_min_D_excess:
                max_min_D_excess = excess

        proc.wait()
        dt = time.time() - t0
        print(f"n={n:2d}: trees={n_checked:8d} minD_exceeds={n_min_D_exceeds} "
              f"max_minD_excess={max_min_D_excess:.10f} ({dt:.1f}s)")

    print(f"\n{'='*80}")
    print(f"All-leaves result (n <= {max_n}):")
    print(f"  Total d_leaf<=1 trees: {n_total}")
    print(f"  Trees where EVERY leaf gives D > 1-a: {n_min_D_exceeds}")
    print(f"  Max (min-over-leaves D - 1/(1+lam)): {max_min_D_excess:.15f}")
    if n_min_D_exceeds == 0:
        print(f"  EVERY TREE has at least one leaf with D <= 1/(1+lam).")
        print(f"  This means D <= 1-a can be achieved with the right leaf choice!")


def direct_mu_P_bound(max_n: int = 23, geng: str = "/opt/homebrew/bin/geng"):
    """Directly verify mu_P >= m-2 for ALL degree-2 leaves (not just canonical).

    This is the actual Route-1 target. We know it holds (0 failures through n=23).
    But let's measure the margin per-leaf.
    """
    print(f"\nDirect mu_P >= m-2 verification with margin (n <= {max_n})")
    print("Checking ALL degree-2-support leaves\n")

    n_checks = 0
    n_fails = 0
    min_margin = 999.0
    min_margin_info = None

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
            if m == 0 or m - 1 >= len(poly_t) or poly_t[m] == 0 or poly_t[m - 1] == 0:
                continue
            lam = poly_t[m - 1] / poly_t[m]

            deg = [len(nb) for nb in adj]
            leaves = [v for v, d in enumerate(deg) if d == 1]

            for leaf in leaves:
                support = adj[leaf][0]
                if deg[support] != 2:
                    continue

                u_node = adj[support][0] if adj[support][1] == leaf else adj[support][1]
                b_adj, idx_map = remove_vertices(adj, {leaf, support})
                u_in_b = idx_map[u_node]

                dp0, dp1, children_b, _, _ = rooted_dp(b_adj, u_in_b)
                P = dp0[u_in_b]
                mu_P_val = mean_at_lambda(P, lam)
                margin = mu_P_val - (m - 2)

                n_checks += 1
                n_checked += 1

                if margin < -1e-12:
                    n_fails += 1

                if margin < min_margin:
                    min_margin = margin
                    g6 = raw.decode("ascii").strip()
                    min_margin_info = {
                        'g6': g6, 'n': nn, 'm': m, 'lam': lam,
                        'leaf': leaf, 'support': support,
                        'mu_P': mu_P_val, 'margin': margin,
                        'num_children_u': len(children_b[u_in_b]),
                    }

        proc.wait()
        dt = time.time() - t0
        print(f"n={n:2d}: checked={n_checked:8d} fails={n_fails} "
              f"min_margin={min_margin:.10f} ({dt:.1f}s)")

    print(f"\n{'='*80}")
    print(f"Direct mu_P >= m-2 result (n <= {max_n}):")
    print(f"  Total leaf-checks: {n_checks}")
    print(f"  Failures: {n_fails}")
    print(f"  Min margin: {min_margin:.15f}")
    if min_margin_info:
        print(f"\n  Tightest case:")
        for k, v in min_margin_info.items():
            print(f"    {k}: {v}")


def main():
    ap = argparse.ArgumentParser(description="Route-1 transfer analysis v2")
    ap.add_argument("--mode", choices=["combined", "one-child", "multi-child",
                                        "all-leaves", "direct", "all"],
                    default="all")
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    args = ap.parse_args()

    if args.mode in ("combined", "all"):
        combined_analysis(max_n=args.max_n, geng=args.geng)

    if args.mode in ("one-child", "all"):
        one_child_detailed(max_n=args.max_n, geng=args.geng)

    if args.mode in ("multi-child", "all"):
        multi_child_analysis(max_n=args.max_n, geng=args.geng)

    if args.mode in ("all-leaves", "all"):
        all_leaves_scan(max_n=args.max_n, geng=args.geng)

    if args.mode in ("direct", "all"):
        direct_mu_P_bound(max_n=args.max_n, geng=args.geng)


if __name__ == "__main__":
    main()
