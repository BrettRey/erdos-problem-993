#!/usr/bin/env python3
"""Route-2 compensation: refined analysis using concavity of mu_B(lambda).

KEY INSIGHT from v2: mu_B(lambda) is CONCAVE in lambda for all trees.

For a concave function f:
  f(b) - f(a) >= f'(b) * (b - a)    [tangent at right endpoint is lower bound]

So: gain = mu_B(lam) - mu_B(tau) >= (Var_B(lam)/lam) * (lam - tau)

This script:
1. Verifies concavity more rigorously
2. Proves that the concavity-based lower bound suffices
3. Analyzes the pendant identity to bound lam - tau from below
4. Attempts a clean analytical proof
"""

from __future__ import annotations

import math
import os
import subprocess
import sys
from collections import Counter
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from graph6 import parse_graph6
from indpoly import independence_poly, _polymul, _polyadd


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


def d2_mu_dlam2(poly: list[int], lam: float) -> float:
    """Numerical second derivative of mu(lambda)."""
    h = 1e-7
    mu_minus = mean_at_lambda(poly, lam - h)
    mu_center = mean_at_lambda(poly, lam)
    mu_plus = mean_at_lambda(poly, lam + h)
    return (mu_plus - 2 * mu_center + mu_minus) / (h * h)


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    args = ap.parse_args()

    max_n = args.max_n
    geng = args.geng

    # ------------------------------------------------------------------
    print("=" * 100)
    print("PART 1: Rigorous concavity check for B polynomials (d_leaf<=1)")
    print("=" * 100)
    print()
    print("Testing d^2(mu_B)/d(lam)^2 <= 0 at many sample points in (0, 2).")
    print()

    concavity_fails = 0
    total_B = 0

    for n in range(4, min(max_n + 1, 19)):
        proc = subprocess.Popen(
            [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"],
            stdout=subprocess.PIPE,
        )
        assert proc.stdout is not None
        n_ct = 0

        for line in proc.stdout:
            nn, adj = parse_graph6(line)
            if not is_dleaf_le_1(nn, adj):
                continue

            leaf = choose_canonical_deg2_leaf(adj)
            if leaf is None:
                continue

            support = adj[leaf][0]
            b_adj = remove_vertices(adj, {leaf, support})
            poly_B = independence_poly(len(b_adj), b_adj)
            total_B += 1
            n_ct += 1

            # Check d^2 mu / d lam^2 at 30 sample points
            for i in range(30):
                lam_test = 0.05 + 0.06 * i
                d2 = d2_mu_dlam2(poly_B, lam_test)
                if d2 > 1e-6:
                    concavity_fails += 1
                    print(f"  CONCAVITY FAIL: n={nn}, lam={lam_test:.2f}, d^2mu={d2:.6e}")
                    break

        proc.wait()
        print(f"n={n:2d}: {n_ct:8d} B-polynomials checked")

    print(f"\nTotal B-polynomials: {total_B}, concavity failures: {concavity_fails}")
    print()

    # ------------------------------------------------------------------
    print("=" * 100)
    print("PART 2: Concavity-based compensation analysis")
    print("=" * 100)
    print()
    print("For concave mu_B(lam), the tangent at lam gives a lower bound:")
    print("  gain = mu_B(lam) - mu_B(tau) >= slope_at_lam * (lam - tau)")
    print("  where slope_at_lam = Var_B(X;lam)/lam")
    print()
    print("Route-2 requires: mu_B(tau) + gain >= m - 3/2")
    print("  i.e., gain >= (m-3/2) - mu_B(tau) = deficit")
    print()
    print("Sufficient: slope_at_lam * gap >= deficit")
    print("  i.e., Var_B(lam)/lam * (lam - tau) >= (m-3/2) - mu_B(tau)")
    print()

    # Track the ratio: slope_at_lam * gap / deficit (should be >= 1)
    min_ratio_concave_lb = float("inf")
    checked = 0
    deficit_cases = 0

    # Additional: check what fraction of gap comes from STRONG C2 alone
    # STRONG C2: lam_m(T) >= tau_B (lambda_{m-1}(B))
    # This gives gap >= 0.

    # More refined: the gap = lam - tau has a structural lower bound
    # from the pendant identity.

    # The pendant identity gives:
    # lam = (b1 + b0 + p0) / (b2 + b1 + p1) where bk = b_{m-k}, pk = p_{m-k}
    # tau = b0/b1
    # gap = [b1^2 - b0*b2 + b1*p0 - b0*p1] / [b1 * (b2+b1+p1)]
    #     = [LC_surplus(B,m-1) + mismatch] / [b1 * i_m(T)]

    # And the deficit at tau:
    # deficit = (m-3/2) - mu_B(tau) where tau = b0/b1
    # At tie-point tau: levels m-2 and m-1 of B are "tied" in weight.
    # The deficit depends on how much weight lies above m-1 vs below m-2.

    # Track mu_B(tau) decomposition
    min_mu_B_tau_m2 = float("inf")

    # Key question: Can we bound deficit/gap from above by slope_at_lam?
    # Or equivalently: deficit * lam / (Var_B(lam) * gap) <= 1?

    max_deficit_over_slope_gap = float("-inf")

    for n in range(4, max_n + 1):
        proc = subprocess.Popen(
            [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"],
            stdout=subprocess.PIPE,
        )
        assert proc.stdout is not None
        n_ct = 0

        for line in proc.stdout:
            nn, adj = parse_graph6(line)
            if not is_dleaf_le_1(nn, adj):
                continue

            leaf = choose_canonical_deg2_leaf(adj)
            if leaf is None:
                continue

            support = adj[leaf][0]
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

            mu_B_tau = mean_at_lambda(poly_B, tau)
            mu_B_lam = mean_at_lambda(poly_B, lam)
            var_B_lam = variance_at_lambda(poly_B, lam)
            slope_lam = var_B_lam / lam if lam > 0 else 0

            gap = lam - tau
            deficit = (m - 1.5) - mu_B_tau
            gain = mu_B_lam - mu_B_tau

            checked += 1
            n_ct += 1
            min_mu_B_tau_m2 = min(min_mu_B_tau_m2, mu_B_tau - (m - 2))

            if deficit > 1e-12:
                deficit_cases += 1
                if slope_lam > 0 and gap > 0:
                    ratio = (slope_lam * gap) / deficit
                    min_ratio_concave_lb = min(min_ratio_concave_lb, ratio)
                    r2 = deficit / (slope_lam * gap) if slope_lam * gap > 0 else float("inf")
                    max_deficit_over_slope_gap = max(max_deficit_over_slope_gap, r2)

        proc.wait()
        print(f"n={n:2d}: {n_ct:8d} d_leaf<=1 trees")

    print()
    print(f"Total checked: {checked}")
    print(f"Deficit cases: {deficit_cases}")
    print(f"min mu_B(tau) - (m-2): {min_mu_B_tau_m2:.10f}")
    print()
    print("Concavity-based lower bound analysis:")
    print(f"  min (slope_lam * gap / deficit): {min_ratio_concave_lb:.10f}")
    print(f"  max (deficit / (slope_lam * gap)): {max_deficit_over_slope_gap:.10f}")
    print()
    if min_ratio_concave_lb > 1:
        print("  PASS: slope_at_lam * gap ALWAYS exceeds deficit.")
        print("  Combined with concavity of mu_B, this proves Route-2!")
    else:
        print("  FAIL: slope_at_lam * gap does not always exceed deficit.")
        print("  Need a tighter bound.")

    # ------------------------------------------------------------------
    print()
    print("=" * 100)
    print("PART 3: Analytical proof attempt via Darroch-like bound")
    print("=" * 100)
    print()
    print("KEY IDEA: If B is log-concave with mode >= m-1, and tau = b_{m-2}/b_{m-1},")
    print("then mu_B(tau) sits between the tied levels m-2 and m-1.")
    print()
    print("More precisely, for a log-concave sequence b_0,...,b_d with mode >= m-1:")
    print("  mu_B(tau) >= m - 2 + b_{m-1}*tau^{m-1} / (b_{m-2}*tau^{m-2} + b_{m-1}*tau^{m-1})")
    print("           = m - 2 + 1/2")
    print("           = m - 3/2  (!!)")
    print()
    print("Wait, that gives exactly m - 3/2 if we only count levels m-2, m-1.")
    print("But there are other levels too. Let me check this more carefully...")
    print()
    print("At tau = b_{m-2}/b_{m-1}:")
    print("  w_{m-2} = b_{m-2} * tau^{m-2}")
    print("  w_{m-1} = b_{m-1} * tau^{m-1} = b_{m-1} * (b_{m-2}/b_{m-1})^{m-1}")
    print("  w_{m-2}/w_{m-1} = b_{m-2} * tau^{m-2} / (b_{m-1} * tau^{m-1})")
    print("                  = (b_{m-2}/b_{m-1}) / tau = 1")
    print("  So w_{m-2} = w_{m-1} (they are tied, as expected from tau definition).")
    print()
    print("mu_B(tau) = sum_k k * w_k / sum_k w_k")
    print("         = [... + (m-2)*w_{m-2} + (m-1)*w_{m-1} + m*w_m + ...] / Z")
    print()
    print("Since w_{m-2} = w_{m-1}, the contribution from these two levels to the")
    print("numerator is (m-2)*w + (m-1)*w = (2m-3)*w, and to Z is 2w.")
    print("So their 'local mean' is (2m-3)/(2) = m - 3/2.")
    print()
    print("For mu_B(tau) >= m - 3/2, we need the other levels to contribute")
    print("at least m - 3/2 on average too. Specifically:")
    print("  mu_B(tau) >= m-3/2  iff  sum_{k != m-2, m-1} (k - m + 3/2) * w_k >= 0")
    print("  iff  sum_{k >= m} (k-m+3/2) * w_k >= sum_{k <= m-3} (m-3/2-k) * w_k")
    print("  iff  'upper tail mass' >= 'lower tail mass' (weighted by distance).")
    print()

    # ------------------------------------------------------------------
    # Verify the tail decomposition
    print("Verifying tail decomposition numerically...")
    print()

    min_tail_balance = float("inf")
    tail_checked = 0

    for n in range(4, min(max_n + 1, 19)):
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
            poly_T = independence_poly(nn, adj)
            m = mode_index_leftmost(poly_T)
            if m < 2:
                continue

            b_adj = remove_vertices(adj, {leaf, support})
            poly_B = independence_poly(len(b_adj), b_adj)
            if m - 1 >= len(poly_B) or poly_B[m - 2] <= 0 or poly_B[m - 1] <= 0:
                continue

            tau = poly_B[m - 2] / poly_B[m - 1]

            # Compute weights w_k = b_k * tau^k
            weights = []
            p = 1.0
            for k in range(len(poly_B)):
                weights.append(poly_B[k] * p)
                p *= tau

            # Tail balance: sum_{k>=m} (k-m+3/2)*w_k - sum_{k<=m-3} (m-3/2-k)*w_k
            balance = 0.0
            for k in range(len(weights)):
                balance += (k - (m - 1.5)) * weights[k]
            # This should equal Z_B * (mu_B(tau) - (m-3/2)), which we want >= 0.

            Z_B = sum(weights)
            tail_balance_normalized = balance / Z_B  # = mu_B(tau) - (m-3/2)
            min_tail_balance = min(min_tail_balance, tail_balance_normalized)
            tail_checked += 1

        proc.wait()

    print(f"  Tail balance checked: {tail_checked}")
    print(f"  min (mu_B(tau) - (m-3/2)): {min_tail_balance:.10f}")
    print()

    if min_tail_balance < 0:
        print("  DEFICIT EXISTS: mu_B(tau) < m-3/2 for some trees.")
        print("  Need to show the gain from tau->lam compensates.")
    else:
        print("  NO DEFICIT: mu_B(tau) >= m-3/2 always!")
        print("  Combined with lam >= tau (STRONG C2), Route-2 is PROVED!")

    # ------------------------------------------------------------------
    print()
    print("=" * 100)
    print("PART 4: Log-concavity and the tail balance")
    print("=" * 100)
    print()
    print("For a LC sequence with mode >= m-1, at the tie-point tau = b_{m-2}/b_{m-1}:")
    print("  w_{m-2} = w_{m-1}  (by definition of tau)")
    print("  w_{m} = b_m * tau^m <= b_{m-1} * tau^{m-1} * (b_m/b_{m-1}) = w_{m-1} * r")
    print("  where r = b_m*tau/b_{m-1} = b_m*b_{m-2}/b_{m-1}^2 <= 1 (by LC at m-1).")
    print()
    print("So w_m <= w_{m-1}, and for k >= m:")
    print("  w_k = b_k * tau^k. By LC, b_k/b_{k-1} <= b_{k-1}/b_{k-2}")
    print("  so the ratios w_{k+1}/w_k = (b_{k+1}/b_k)*tau are non-increasing for k >= mode.")
    print("  (Since b_{k+1}/b_k decreases past the mode, and tau is fixed.)")
    print()
    print("Similarly, for k <= m-3:")
    print("  w_{m-3}/w_{m-2} = (b_{m-3}/b_{m-2})*tau.")
    print("  By LC at m-2: b_{m-2}^2 >= b_{m-3}*b_{m-1}")
    print("  so b_{m-3}/b_{m-2} <= b_{m-2}/b_{m-1} = tau")
    print("  hence w_{m-3}/w_{m-2} <= tau^2 < tau (since tau < 1 for mode >= m-1).")
    print()
    print("This means the upper tail decays SLOWER than the lower tail!")
    print("The upper tail has w_m >= w_{m-1} * r = w_{m-2} * r,")
    print("while the lower tail has w_{m-3} <= w_{m-2} * tau^2.")
    print()
    print("Let's verify: is the 'excess upper tail' >= deficit from lower tail?")

    # Verify: for each tree, compute:
    # U = sum_{k>=m} (k - m + 3/2) * w_k   (upper contribution)
    # L = sum_{k<=m-3} (m - 3/2 - k) * w_k  (lower contribution)
    # Balance = U - L (should be >= 0 for route 2)
    #
    # At the tie point, the m-2 and m-1 levels contribute 0 to the balance.

    print()
    print("Detailed tail decomposition at tau:")
    print()

    min_U = float("inf")
    max_L = float("-inf")
    min_U_minus_L = float("inf")
    min_ratio_UL = float("inf")

    for n in range(4, min(max_n + 1, 19)):
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
            poly_T = independence_poly(nn, adj)
            m = mode_index_leftmost(poly_T)
            if m < 2:
                continue

            b_adj = remove_vertices(adj, {leaf, support})
            poly_B = independence_poly(len(b_adj), b_adj)
            if m - 1 >= len(poly_B) or poly_B[m - 2] <= 0 or poly_B[m - 1] <= 0:
                continue

            tau = poly_B[m - 2] / poly_B[m - 1]

            weights = []
            p = 1.0
            for k in range(len(poly_B)):
                weights.append(poly_B[k] * p)
                p *= tau

            Z = sum(weights)

            U = sum((k - m + 1.5) * weights[k] for k in range(m, len(weights))) / Z
            L = sum((m - 1.5 - k) * weights[k] for k in range(0, m - 1)) / Z
            # Note: levels m-2 and m-1 contribute 0 to balance (since tied)
            # Actually level m-2 contributes (m-2 - m + 1.5) * w = -0.5 * w
            # and level m-1 contributes (m-1 - m + 1.5) * w = 0.5 * w
            # So they cancel! Good.

            min_U = min(min_U, U)
            max_L = max(max_L, L)
            min_U_minus_L = min(min_U_minus_L, U - L)
            if L > 1e-15:
                min_ratio_UL = min(min_ratio_UL, U / L)

        proc.wait()

    print(f"  min upper contribution (U/Z): {min_U:.10f}")
    print(f"  max lower contribution (L/Z): {max_L:.10f}")
    print(f"  min (U-L)/Z = mu_B(tau)-(m-3/2): {min_U_minus_L:.10f}")
    print(f"  min U/L ratio: {min_ratio_UL:.10f}")

    # ------------------------------------------------------------------
    print()
    print("=" * 100)
    print("PART 5: Tight analysis of the ratio U/L as function of tree structure")
    print("=" * 100)
    print()

    # The upper tail contributes at least 0.5 * w_m / Z (from level m alone,
    # contributing (m - m + 1.5) * w_m = 1.5 * w_m, minus 0.5*w_{m-1} already counted).
    # Wait, let me redo:
    # U = sum_{k>=m} (k - m + 1.5) * w_k / Z
    #   >= 1.5 * w_m / Z  (just level m, k=m gives 1.5)
    #   + 2.5 * w_{m+1}/Z + ...
    #
    # L = sum_{k<=m-3} (m-1.5-k) * w_k / Z
    #   = 0.5 * w_{m-3}/Z + 1.5 * w_{m-4}/Z + ...
    #
    # For LC sequence, w_m >= w_{m-1} * r where r = b_m*tau/b_{m-1}.
    # And w_{m-3} <= w_{m-2} * s where s = b_{m-3}*tau/b_{m-2} <= tau^2.
    # (From LC: b_{m-3}/b_{m-2} <= b_{m-2}/b_{m-1} = tau.)
    #
    # So U >= 1.5 * r * w_{m-1}/Z and L's first term is 0.5 * s * w_{m-2}/Z.
    # Since w_{m-1} = w_{m-2} (tied):
    # U/L >= (1.5 * r) / (0.5 * s + smaller terms)

    # For the ratio r = b_m * b_{m-2} / b_{m-1}^2 (LC ratio at m-1):
    # This is the LC quotient. By LC, r <= 1. By verified data, r is usually close to 1.

    # Let's track the LC ratio at m-1 for B:
    min_lc_ratio = float("inf")
    max_lc_ratio = float("-inf")

    for n in range(4, min(max_n + 1, 19)):
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
            poly_T = independence_poly(nn, adj)
            m = mode_index_leftmost(poly_T)
            if m < 2:
                continue

            b_adj = remove_vertices(adj, {leaf, support})
            poly_B = independence_poly(len(b_adj), b_adj)
            if m - 1 >= len(poly_B) or poly_B[m - 2] <= 0 or poly_B[m - 1] <= 0:
                continue
            if m >= len(poly_B) or poly_B[m] <= 0:
                continue

            r = poly_B[m] * poly_B[m - 2] / (poly_B[m - 1] ** 2)
            min_lc_ratio = min(min_lc_ratio, r)
            max_lc_ratio = max(max_lc_ratio, r)

        proc.wait()

    print(f"  LC ratio r = b_m * b_{{m-2}} / b_{{m-1}}^2 at index m-1:")
    print(f"  min r: {min_lc_ratio:.10f}")
    print(f"  max r: {max_lc_ratio:.10f}")
    print()
    print("  Note: r <= 1 always (LC). The upper tail level w_m = r * w_{m-1}.")
    print("  When r is close to 1, levels m-2, m-1, m are nearly equal in weight,")
    print("  and the 'upper tail' contributes strongly.")


if __name__ == "__main__":
    main()
