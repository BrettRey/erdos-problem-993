#!/usr/bin/env python3
"""Algebraic exploration of the Route-2 compensation inequality.

Target: mu_B(lambda_m(T)) >= m - 3/2

This script explores several algebraic approaches:

1. Cavity-field decomposition of mu_B in terms of vertex occupation
   probabilities at fugacity lambda.

2. Relationship between lambda_m(T) and tau_B = b_{m-2}/b_{m-1}:
   - Lower bound on lambda_m - tau_B
   - Lower bound on d(mu_B)/d(lambda) = Var_B(X)/lambda

3. Product-of-means decomposition of mu_P and its implications for mu_B.

4. The pendant-shift identity: I(T) = (1+x)I(B) + xP, and what it implies
   for how mode(T) relates to mode(B) and the means.

5. Small-tree exact analysis to identify structural patterns.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from collections import Counter
from fractions import Fraction
from typing import Any

# Allow imports from project root
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from graph6 import parse_graph6
from indpoly import independence_poly


def mode_index_leftmost(poly: list[int]) -> int:
    return max(range(len(poly)), key=lambda i: poly[i])


def remove_vertices(adj: list[list[int]], remove_set: set[int]) -> list[list[int]]:
    keep = [v for v in range(len(adj)) if v not in remove_set]
    idx = {v: i for i, v in enumerate(keep)}
    out = [[] for _ in keep]
    for v in keep:
        vv = idx[v]
        for u in adj[v]:
            if u in idx:
                out[vv].append(idx[u])
    return out


def choose_canonical_deg2_leaf(adj: list[list[int]]) -> int | None:
    deg = [len(nb) for nb in adj]
    cand = []
    for v, d in enumerate(deg):
        if d != 1:
            continue
        s = adj[v][0]
        if deg[s] == 2:
            cand.append(v)
    if not cand:
        return None
    return min(cand)


def mean_at_lambda(poly: list[int], lam: float) -> float:
    z = 0.0
    mu_num = 0.0
    p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p
        z += w
        mu_num += k * w
        p *= lam
    return mu_num / z if z else 0.0


def variance_at_lambda(poly: list[int], lam: float) -> float:
    """Var_lambda(X) where X is size under Gibbs measure at fugacity lambda."""
    z = 0.0
    mu_num = 0.0
    mu2_num = 0.0
    p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p
        z += w
        mu_num += k * w
        mu2_num += k * k * w
        p *= lam
    if z == 0:
        return 0.0
    mu = mu_num / z
    mu2 = mu2_num / z
    return mu2 - mu * mu


def mean_exact(poly: list[int], lam: Fraction) -> Fraction:
    """Exact mean at rational fugacity."""
    z = Fraction(0)
    mu_num = Fraction(0)
    p = Fraction(1)
    for k, ck in enumerate(poly):
        w = Fraction(ck) * p
        z += w
        mu_num += k * w
        p *= lam
    return mu_num / z if z else Fraction(0)


def dp_tree(n: int, adj: list[list[int]], root: int = 0):
    """Compute DP tables dp[v][0], dp[v][1] for independence polynomial.

    dp[v][0] = polynomial for subtree(v) with v excluded
    dp[v][1] = polynomial for subtree(v) with v included

    Returns dict mapping vertex -> (poly_exclude, poly_include)
    """
    from indpoly import _polymul

    # BFS to get parent/children structure
    parent = [-1] * n
    children: list[list[int]] = [[] for _ in range(n)]
    visited = [False] * n
    order = []
    queue = [root]
    visited[root] = True
    idx = 0
    while idx < len(queue):
        u = queue[idx]
        idx += 1
        order.append(u)
        for v in adj[u]:
            if not visited[v]:
                visited[v] = True
                parent[v] = u
                children[u].append(v)
                queue.append(v)

    dp: dict[int, tuple[list[int], list[int]]] = {}

    # Post-order traversal
    for v in reversed(order):
        if not children[v]:
            dp[v] = ([1], [0, 1])
        else:
            excl = [1]
            incl = [1]
            for c in children[v]:
                c_excl, c_incl = dp[c]
                c_total = [0] * max(len(c_excl), len(c_incl))
                for i, x in enumerate(c_excl):
                    c_total[i] += x
                for i, x in enumerate(c_incl):
                    c_total[i] += x
                excl = _polymul(excl, c_total)
                incl = _polymul(incl, c_excl)
            # dp[v][1] = x * incl (shift by 1)
            dp[v] = (excl, [0] + incl)

    return dp, children, parent


def analyze_tree_detailed(n: int, adj: list[list[int]], leaf: int):
    """Detailed analysis of a single tree at a specific degree-2 leaf."""
    from indpoly import _polymul

    support = adj[leaf][0]
    deg = [len(nb) for nb in adj]
    assert deg[support] == 2, f"support has degree {deg[support]}, expected 2"

    # Find u (other neighbor of support)
    u = [v for v in adj[support] if v != leaf][0]

    # Compute I(T)
    poly_T = independence_poly(n, adj)
    m = mode_index_leftmost(poly_T)
    if m < 2:
        return None
    lam = poly_T[m - 1] / poly_T[m]

    # Compute B = T - {leaf, support}
    b_adj = remove_vertices(adj, {leaf, support})
    poly_B = independence_poly(len(b_adj), b_adj)

    if m - 1 >= len(poly_B) or poly_B[m - 2] <= 0 or poly_B[m - 1] <= 0:
        return None

    tau = poly_B[m - 2] / poly_B[m - 1]
    mu_B_lam = mean_at_lambda(poly_B, lam)
    mu_B_tau = mean_at_lambda(poly_B, tau)
    var_B_lam = variance_at_lambda(poly_B, lam)
    var_B_tau = variance_at_lambda(poly_B, tau)

    # Compute P and Q via DP on B rooted at u
    # Need to map u to its index in B
    keep = [v for v in range(n) if v not in {leaf, support}]
    u_in_B = keep.index(u)

    dp_B, children_B, parent_B = dp_tree(len(b_adj), b_adj, root=u_in_B)
    poly_P = dp_B[u_in_B][0]  # u excluded
    poly_Q_raw = dp_B[u_in_B][1]  # u included (has leading 0)

    # Verify: P + Q should equal I(B)
    poly_Q = poly_Q_raw
    check = [0] * max(len(poly_P), len(poly_Q))
    for i, x in enumerate(poly_P):
        check[i] += x
    for i, x in enumerate(poly_Q):
        check[i] += x
    assert check == poly_B, "P + Q != I(B)"

    mu_P_lam = mean_at_lambda(poly_P, lam)

    # Probability u is in IS at fugacity lambda (under B's Gibbs measure)
    z_P = sum(poly_P[k] * lam ** k for k in range(len(poly_P)))
    z_Q = sum(poly_Q[k] * lam ** k for k in range(len(poly_Q)))
    z_B = z_P + z_Q
    p_u = z_Q / z_B

    # mu_B = (1-p_u)*mu_P + p_u*mu_Q = mu_P + p_u (since mu_Q = mu_P_children + 1)
    # Actually mu_B = [(1-p_u)*Z_P*mu_P + p_u*Z_Q*mu_Q] / Z_B
    # Let's compute directly
    mu_Q_lam = mean_at_lambda(poly_Q, lam) if z_Q > 0 else 0.0

    # Steiner peeling: mean decomposition
    # mu_B(lam) = sum over v in B of P(v in IS | fugacity lam)
    # For the rooted DP: mu_B = mu_P + p_u

    # Compute the slope d(mu_B)/d(lambda) = Var_B(X)/lambda
    slope_at_lam = var_B_lam / lam if lam > 0 else 0
    slope_at_tau = var_B_tau / tau if tau > 0 else 0

    # Route-2 bound analysis
    target = m - 1.5
    deficit_tau = target - mu_B_tau
    gap = lam - tau
    gain = mu_B_lam - mu_B_tau
    route2_slack = mu_B_lam - target
    exact_threshold = m - 1 - lam / (1 + lam)
    exact_slack = mu_B_lam - exact_threshold

    return {
        "n": n,
        "m": m,
        "lam": lam,
        "tau": tau,
        "gap": gap,
        "poly_T": poly_T[:m + 3],
        "poly_B": poly_B[:m + 2],
        "poly_P": poly_P[:m + 1],
        "mu_B_lam": mu_B_lam,
        "mu_B_tau": mu_B_tau,
        "mu_P_lam": mu_P_lam,
        "mu_Q_lam": mu_Q_lam,
        "var_B_lam": var_B_lam,
        "var_B_tau": var_B_tau,
        "slope_at_lam": slope_at_lam,
        "slope_at_tau": slope_at_tau,
        "p_u": p_u,
        "deficit_tau": deficit_tau,
        "gain": gain,
        "route2_slack": route2_slack,
        "exact_slack": exact_slack,
        "target": target,
        "deg_sig": dict(sorted(Counter(deg).items())),
    }


def scan_structural_bounds(max_n: int = 18, geng: str = "/opt/homebrew/bin/geng"):
    """Scan for structural relationships that could yield an analytic proof."""

    print(f"Scanning structural bounds for Route-2 compensation, n=4..{max_n}")
    print("=" * 100)

    # Track worst cases for various candidate bounds
    stats = {
        "checked": 0,
        "min_route2_slack": float("inf"),
        "min_exact_slack": float("inf"),
        "min_var_ratio": float("inf"),         # Var_B(X) / (target - mu_B(tau))
        "min_slope_compensation": float("inf"),  # slope * gap - deficit
        "min_p_u": float("inf"),
        "max_p_u": float("-inf"),
        "min_mode_gap_B": float("inf"),  # mode(B) - (m-1)
        "max_deficit_tau": float("-inf"),
        # Bound: mu_B(tau) >= m - 2 + p_u(tau)
        "min_mu_B_tau_minus_m_plus_2": float("inf"),
        # Bound: mu_P >= m - 2 - epsilon for some epsilon
        "min_mu_P_minus_m_plus_2": float("inf"),
        # The mean-value theorem bound: gain >= slope_min * gap
        "min_mvt_ratio": float("inf"),  # gain / (slope_at_tau * gap) if deficit > 0
    }

    worst_tree = None

    for n in range(4, max_n + 1):
        proc = subprocess.Popen(
            [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"],
            stdout=subprocess.PIPE,
        )
        assert proc.stdout is not None

        n_count = 0
        for line in proc.stdout:
            nn, adj = parse_graph6(line)
            if not is_dleaf_le_1(nn, adj):
                continue

            leaf = choose_canonical_deg2_leaf(adj)
            if leaf is None:
                continue

            result = analyze_tree_detailed(nn, adj, leaf)
            if result is None:
                continue

            stats["checked"] += 1
            n_count += 1

            r2 = result["route2_slack"]
            if r2 < stats["min_route2_slack"]:
                stats["min_route2_slack"] = r2
                worst_tree = {**result, "g6": line.decode("ascii").strip()}

            stats["min_exact_slack"] = min(stats["min_exact_slack"], result["exact_slack"])
            stats["min_p_u"] = min(stats["min_p_u"], result["p_u"])
            stats["max_p_u"] = max(stats["max_p_u"], result["p_u"])
            stats["min_mu_P_minus_m_plus_2"] = min(
                stats["min_mu_P_minus_m_plus_2"],
                result["mu_P_lam"] - (result["m"] - 2),
            )
            stats["min_mu_B_tau_minus_m_plus_2"] = min(
                stats["min_mu_B_tau_minus_m_plus_2"],
                result["mu_B_tau"] - (result["m"] - 2),
            )

            # For deficit cases, check mean-value theorem bound
            if result["deficit_tau"] > 1e-12 and result["gap"] > 1e-12:
                stats["max_deficit_tau"] = max(stats["max_deficit_tau"], result["deficit_tau"])
                # Slope at tau gives a LOWER bound on gain (since mu_B is convex in lambda)
                mvt_gain_lb = result["slope_at_tau"] * result["gap"]
                mvt_ratio = mvt_gain_lb / result["gain"] if result["gain"] > 0 else 0
                stats["min_mvt_ratio"] = min(stats["min_mvt_ratio"], mvt_ratio)

                # Variance at tau: ratio to deficit
                if result["deficit_tau"] > 0:
                    var_ratio = result["var_B_tau"] / result["deficit_tau"]
                    stats["min_var_ratio"] = min(stats["min_var_ratio"], var_ratio)

                # Direct compensation check: slope * gap >= deficit?
                slope_comp = result["slope_at_tau"] * result["gap"] - result["deficit_tau"]
                stats["min_slope_compensation"] = min(stats["min_slope_compensation"], slope_comp)

            # mode(B) vs m-1
            mode_B = mode_index_leftmost(result["poly_B"])
            mode_gap_B = mode_B - (result["m"] - 1)
            stats["min_mode_gap_B"] = min(stats["min_mode_gap_B"], mode_gap_B)

        proc.wait()
        print(f"n={n:2d}: {n_count:8d} trees checked")

    print()
    print("=" * 100)
    print(f"Total checked: {stats['checked']}")
    print()
    print("--- Route-2 slack ---")
    print(f"  min route2_slack:    {stats['min_route2_slack']:.10f}")
    print(f"  min exact_slack:     {stats['min_exact_slack']:.10f}")
    print()
    print("--- Cavity field decomposition ---")
    print(f"  min p_u:             {stats['min_p_u']:.10f}")
    print(f"  max p_u:             {stats['max_p_u']:.10f}")
    print(f"  min mu_B(tau)-(m-2): {stats['min_mu_B_tau_minus_m_plus_2']:.10f}")
    print(f"  min mu_P(lam)-(m-2): {stats['min_mu_P_minus_m_plus_2']:.10f}")
    print(f"  min mode(B)-(m-1):   {stats['min_mode_gap_B']:.1f}")
    print()
    print("--- Mean-value theorem analysis (deficit cases) ---")
    print(f"  max deficit_tau:     {stats['max_deficit_tau']:.10f}")
    print(f"  min Var/deficit:     {stats['min_var_ratio']:.10f}")
    print(f"  min mvt_ratio:       {stats['min_mvt_ratio']:.10f}")
    print(f"  min slope*gap-def:   {stats['min_slope_compensation']:.10f}")

    if worst_tree:
        print()
        print("--- Worst tree (min route2_slack) ---")
        for k, v in worst_tree.items():
            if k != "poly_T" and k != "poly_B" and k != "poly_P":
                print(f"  {k}: {v}")
        print(f"  poly_T (truncated): {worst_tree['poly_T']}")
        print(f"  poly_B (truncated): {worst_tree['poly_B']}")


def scan_convexity_check(max_n: int = 16, geng: str = "/opt/homebrew/bin/geng"):
    """Check whether mu_B(lambda) is convex in lambda for all trees.

    If mu_B is convex, then the MVT lower bound using slope at tau is valid:
      gain >= slope_at_tau * gap

    Convexity of mu_B(lambda) = d^2(mu_B)/d(lambda)^2 >= 0.
    Since d(mu_B)/d(lambda) = Var(X)/lambda, convexity requires
    d(Var(X)/lambda)/d(lambda) >= 0.
    """
    print(f"\nChecking convexity of mu_B(lambda) for n=4..{max_n}")
    print("=" * 80)

    total = 0
    convex_fail = 0

    for n in range(4, max_n + 1):
        proc = subprocess.Popen(
            [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"],
            stdout=subprocess.PIPE,
        )
        assert proc.stdout is not None

        for line in proc.stdout:
            nn, adj = parse_graph6(line)
            # Check ALL trees for convexity, not just d_leaf <= 1
            poly = independence_poly(nn, adj)
            total += 1

            # Check convexity by numerical sampling
            # mu(lam) = sum_k k * c_k * lam^k / sum_k c_k * lam^k
            # d^2(mu)/d(lam)^2 should be >= 0 for convexity
            lam_vals = [0.1 * i for i in range(1, 20)]
            mus = [mean_at_lambda(poly, lam) for lam in lam_vals]

            # Check discrete second differences
            for i in range(1, len(mus) - 1):
                d2 = mus[i + 1] - 2 * mus[i] + mus[i - 1]
                if d2 < -1e-10:
                    convex_fail += 1
                    break

        proc.wait()

    print(f"  Total trees checked: {total}")
    print(f"  Convexity failures:  {convex_fail}")
    print()
    if convex_fail == 0:
        print("  mu_B(lambda) appears convex in lambda for all trees checked!")
        print("  This validates the mean-value theorem approach:")
        print("  gain = mu_B(lam) - mu_B(tau) >= (d mu_B/d lam)|_{tau} * (lam - tau)")
        print("       = Var_B(X;tau)/tau * (lam - tau)")
    else:
        print(f"  WARNING: {convex_fail} trees have non-convex mu_B(lambda)")


def scan_pendant_identity_analysis(max_n: int = 16, geng: str = "/opt/homebrew/bin/geng"):
    """Analyze the pendant-shift identity I(T) = (1+x)I(B) + xP
    to derive bounds on lambda_m(T) in terms of B's coefficients.

    Key identity: i_k(T) = b_k + b_{k-1} + p_{k-1}

    At the mode m of T: i_{m-1}/i_m = lambda_m
    So: (b_{m-1} + b_{m-2} + p_{m-2}) / (b_m + b_{m-1} + p_{m-1}) = lambda_m

    We want to bound lambda_m from below in terms of tau = b_{m-2}/b_{m-1}.
    """
    print(f"\nPendant identity analysis for lambda_m vs tau_B, n=4..{max_n}")
    print("=" * 80)

    min_lam_tau_ratio = float("inf")
    min_lam_minus_tau = float("inf")

    # Track the ratio (b_{m-1} + p_{m-2}) / (b_m + p_{m-1}) vs tau
    min_ratio_analysis = float("inf")

    checked = 0
    deficit_cases = 0

    for n in range(4, max_n + 1):
        proc = subprocess.Popen(
            [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"],
            stdout=subprocess.PIPE,
        )
        assert proc.stdout is not None

        for line in proc.stdout:
            nn, adj = parse_graph6(line)
            if not is_dleaf_le_1(nn, adj):
                continue

            leaf = choose_canonical_deg2_leaf(adj)
            if leaf is None:
                continue

            support = adj[leaf][0]
            u = [v for v in adj[support] if v != leaf][0]

            poly_T = independence_poly(nn, adj)
            m = mode_index_leftmost(poly_T)
            if m < 2:
                continue

            b_adj = remove_vertices(adj, {leaf, support})
            poly_B = independence_poly(len(b_adj), b_adj)

            if m - 1 >= len(poly_B) or poly_B[m - 2] <= 0 or poly_B[m - 1] <= 0:
                continue

            lam = poly_T[m - 1] / poly_T[m]
            tau = poly_B[m - 2] / poly_B[m - 1]

            checked += 1

            # From I(T) = (1+x)I(B) + xP:
            # i_k = b_k + b_{k-1} + p_{k-1}
            # lambda_m = i_{m-1}/i_m = (b_{m-1} + b_{m-2} + p_{m-2}) / (b_m + b_{m-1} + p_{m-1})

            # We know tau = b_{m-2}/b_{m-1}.
            # Let's write lambda_m in terms of tau and the P-ratios.

            b0 = poly_B[m - 2]  # b_{m-2}
            b1 = poly_B[m - 1]  # b_{m-1}
            b2 = poly_B[m] if m < len(poly_B) else 0  # b_m

            # Need P coefficients
            # Compute P via DP
            keep = [v for v in range(nn) if v not in {leaf, support}]
            u_in_B = keep.index(u)
            dp_B, _, _ = dp_tree(len(b_adj), b_adj, root=u_in_B)
            poly_P = dp_B[u_in_B][0]

            p0 = poly_P[m - 2] if m - 2 < len(poly_P) else 0  # p_{m-2}
            p1 = poly_P[m - 1] if m - 1 < len(poly_P) else 0  # p_{m-1}

            # lambda_m = (b1 + b0 + p0) / (b2 + b1 + p1)
            numer = b1 + b0 + p0
            denom = b2 + b1 + p1

            # Verify
            lam_check = numer / denom
            assert abs(lam_check - lam) < 1e-10, f"lambda mismatch: {lam_check} vs {lam}"

            # lambda_m - tau = (b1 + b0 + p0)/(b2 + b1 + p1) - b0/b1
            # = [b1(b1 + b0 + p0) - b0(b2 + b1 + p1)] / [b1(b2 + b1 + p1)]
            # = [b1^2 + b0*b1 + b1*p0 - b0*b2 - b0*b1 - b0*p1] / [b1*(b2+b1+p1)]
            # = [b1^2 - b0*b2 + b1*p0 - b0*p1] / [b1*(b2+b1+p1)]
            # = [LC_surplus(B,m-1) + mismatch(P,B)] / [b1 * i_m(T)]
            #
            # This is exactly the STRONG C2 determinant!

            lc_surplus = b1 * b1 - b0 * b2
            mismatch = b1 * p0 - b0 * p1
            combined = lc_surplus + mismatch

            lam_minus_tau = combined / (b1 * denom)
            assert abs(lam_minus_tau - (lam - tau)) < 1e-10

            min_lam_minus_tau = min(min_lam_minus_tau, lam - tau)

            mu_B_tau = mean_at_lambda(poly_B, tau)
            deficit = (m - 1.5) - mu_B_tau
            if deficit > 1e-12:
                deficit_cases += 1
                # For the MVT approach: gain >= slope_at_tau * gap
                # slope_at_tau = Var_B(X;tau)/tau
                var_tau = variance_at_lambda(poly_B, tau)
                slope_tau = var_tau / tau
                mvt_lb = slope_tau * (lam - tau)

                # Check: does MVT lower bound exceed deficit?
                compensation = mvt_lb - deficit
                min_ratio_analysis = min(min_ratio_analysis, compensation)

        proc.wait()

    print(f"  Checked: {checked}")
    print(f"  Deficit cases: {deficit_cases}")
    print(f"  min(lambda_m - tau): {min_lam_minus_tau:.10f}")
    print(f"  min(MVT_lb - deficit): {min_ratio_analysis:.10f}")
    print()
    if min_ratio_analysis > 0:
        print("  MVT lower bound ALWAYS exceeds deficit!")
        print("  This means: Var_B(tau)/tau * (lam-tau) >= (m-3/2) - mu_B(tau)")
        print("  i.e., the slope at tau is steep enough to compensate in all cases.")
    else:
        print("  MVT lower bound does NOT always suffice.")


def scan_mu_B_lower_bound_at_tau(max_n: int = 18, geng: str = "/opt/homebrew/bin/geng"):
    """Investigate whether mu_B(tau) >= m - 2 always holds.

    If mu_B(tau) >= m - 2 (which is 0.5 below the target m - 3/2),
    and the slope Var/tau is at least 0.5/gap_min, the compensation closes.

    Also check: mu_B(tau) >= m - 2 + p_u(tau) where p_u is the occupation
    probability of u.
    """
    print(f"\nLower bound on mu_B(tau) relative to m-2, n=4..{max_n}")
    print("=" * 80)

    min_mu_B_tau_minus_target = float("inf")  # target = m - 2
    min_mu_B_tau_minus_m15 = float("inf")      # target = m - 1.5
    checked = 0

    for n in range(4, max_n + 1):
        proc = subprocess.Popen(
            [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"],
            stdout=subprocess.PIPE,
        )
        assert proc.stdout is not None

        for line in proc.stdout:
            nn, adj = parse_graph6(line)
            if not is_dleaf_le_1(nn, adj):
                continue

            leaf = choose_canonical_deg2_leaf(adj)
            if leaf is None:
                continue

            poly_T = independence_poly(nn, adj)
            m = mode_index_leftmost(poly_T)
            if m < 2:
                continue

            b_adj = remove_vertices(adj, {leaf, support := adj[leaf][0]})
            poly_B = independence_poly(len(b_adj), b_adj)
            if m - 1 >= len(poly_B) or poly_B[m - 2] <= 0 or poly_B[m - 1] <= 0:
                continue

            tau = poly_B[m - 2] / poly_B[m - 1]
            mu_B_tau = mean_at_lambda(poly_B, tau)

            checked += 1
            min_mu_B_tau_minus_target = min(min_mu_B_tau_minus_target, mu_B_tau - (m - 2))
            min_mu_B_tau_minus_m15 = min(min_mu_B_tau_minus_m15, mu_B_tau - (m - 1.5))

        proc.wait()

    print(f"  Checked: {checked}")
    print(f"  min(mu_B(tau) - (m-2)):   {min_mu_B_tau_minus_target:.10f}")
    print(f"  min(mu_B(tau) - (m-1.5)): {min_mu_B_tau_minus_m15:.10f}")
    print()
    if min_mu_B_tau_minus_target > 0:
        print("  mu_B(tau) > m-2 ALWAYS holds!")
        print("  This gives a baseline: need only 0.5 more gain to reach m-1.5.")
    if min_mu_B_tau_minus_m15 > 0:
        print("  mu_B(tau) > m-1.5 ALWAYS holds (no deficit at all).")
    elif min_mu_B_tau_minus_m15 < 0:
        print(f"  max deficit at tau: {-min_mu_B_tau_minus_m15:.10f}")


def scan_phi_m_bound_direct(max_n: int = 16, geng: str = "/opt/homebrew/bin/geng"):
    """Test the direct algebraic bound approach:

    From Phi_m(T) = (1+lam)Phi_m(B) + lam*Z_B + lam*Phi_{m-1}(P) + lam*Phi_{m-1}(B)

    We need Phi_m(T) >= 0, which is equivalent to mu_T(lam_m) >= m - 1.

    The sufficient condition mu_B(lam) >= m - 1 - lam/(1+lam) comes from
    ignoring the Phi_{m-1}(P) >= 0 term. But we can use both:

    Phi_m(T) >= (1+lam)*Z_B*(mu_B - (m-1)) + lam*Z_B + lam*Z_P*(mu_P - (m-2))
             = Z_B * [(1+lam)(mu_B - m + 1) + lam] + lam*Z_P*(mu_P - m + 2)
             = Z_B * [(1+lam)*mu_B - (1+lam)(m-1) + lam] + lam*Z_P*(mu_P - m + 2)
             = Z_B * [(1+lam)*mu_B - m + 1 - lam*m + lam + lam] + ...
             = Z_B * [(1+lam)*mu_B - m*(1+lam) + (1+2lam)] + lam*Z_P*(mu_P - m + 2)

    Hmm. Let me recompute. Actually the four-term identity should be:

    Phi_m(T) = Phi_m(A) + lam*Phi_{m-1}(B)

    where Phi_m(A) = (1+lam)*Phi_m(B) + lam*Z_B + lam*Phi_{m-1}(P).

    So:
    Phi_m(T) = (1+lam)*Phi_m(B) + lam*Z_B + lam*Phi_{m-1}(P) + lam*Phi_{m-1}(B)

    Now Phi_m(B) = Z_B*(mu_B - (m-1)) and Phi_{m-1}(B) = Z_B*(mu_B - (m-2)).

    So the B-terms: (1+lam)*Z_B*(mu_B - m + 1) + lam*Z_B*(mu_B - m + 2)
    = Z_B * [(1+lam)(mu_B - m + 1) + lam*(mu_B - m + 2)]
    = Z_B * [(1+2lam)*mu_B - (1+lam)(m-1) - lam*(m-2)]
    = Z_B * [(1+2lam)*mu_B - (m-1) - lam*m + lam - lam*m + 2*lam]
    = Z_B * [(1+2lam)*mu_B - (m-1) - 2*lam*m + 3*lam]
    = Z_B * [(1+2lam)*mu_B - m + 1 - 2*lam*m + 3*lam]
    = Z_B * [(1+2lam)*mu_B - (1+2lam)*m + (1+2lam) + 2*lam]  (hmm let me be more careful)

    (1+lam)(mu_B - m + 1) + lam(mu_B - m + 2)
    = (1+lam)*mu_B - (1+lam)(m-1) + lam*mu_B - lam(m-2)
    = (1+2lam)*mu_B - (1+lam)(m-1) - lam(m-2)
    = (1+2lam)*mu_B - (m-1) - lam(m-1) - lam(m-2)
    = (1+2lam)*mu_B - (m-1) - lam(2m-3)
    = (1+2lam)*[mu_B - (m-1)] + lam(2m-2) - lam(2m-3)  -- nope, let me just expand
    = (1+2lam)*mu_B - (m-1) - 2lam*m + 3lam
    = (1+2lam)*mu_B - m + 1 - 2lam*m + 3lam
    = (1+2lam)*(mu_B - m) + (1 + 3lam)

    So: B-terms = Z_B * [(1+2lam)(mu_B - m) + 1 + 3lam]

    And the full expression:
    Phi_m(T) = Z_B * [(1+2lam)(mu_B - m) + 1 + 3lam] + lam*Z_P*(mu_P - m + 2) + lam*Z_B

    Wait, I already included the pendant bonus (lam*Z_B) in the B terms.
    Let me redo from scratch.

    Actually, I realize the clean computation is simpler. We don't need this
    for the scan, the scan just verifies numerical properties.
    """
    print(f"\nDirect Phi_m bound check, n=4..{max_n}")
    print("=" * 80)
    print("(This analysis uses the four-term identity symbolically.)")
    print("Skipping numerical scan -- see algebraic write-up.")


def main():
    import argparse

    ap = argparse.ArgumentParser(description="Route-2 compensation algebraic exploration.")
    ap.add_argument("--max-n", type=int, default=18)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument(
        "--scan",
        choices=["structural", "convexity", "pendant", "mu_B_lb", "all"],
        default="all",
    )
    args = ap.parse_args()

    if args.scan in ("structural", "all"):
        scan_structural_bounds(args.max_n, args.geng)

    if args.scan in ("convexity", "all"):
        scan_convexity_check(min(args.max_n, 16), args.geng)

    if args.scan in ("pendant", "all"):
        scan_pendant_identity_analysis(min(args.max_n, 16), args.geng)

    if args.scan in ("mu_B_lb", "all"):
        scan_mu_B_lower_bound_at_tau(args.max_n, args.geng)


if __name__ == "__main__":
    main()
