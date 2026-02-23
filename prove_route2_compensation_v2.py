#!/usr/bin/env python3
"""Route-2 compensation: deeper analysis of mu_B convexity/concavity,
MVT bounds, and analytical proof attempts.

Key question: Why does mu_B(lambda_m(T)) >= m - 3/2?

Strategy: Decompose into
  (1) mu_B(tau) >= m - 2 + f(tau)   [baseline at tie-point]
  (2) mu_B(lam) - mu_B(tau) >= g(lam - tau)  [gain from fugacity increase]

where (1) + (2) >= m - 3/2.
"""

from __future__ import annotations

import math
import os
import subprocess
import sys
from collections import Counter

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
    """Compute d^2(mu)/d(lam)^2 numerically.

    We know d(mu)/d(lam) = Var(X)/lam.
    d^2(mu)/d(lam)^2 = d(Var/lam)/d(lam).

    Var(X) = E[X^2] - (E[X])^2 where expectations are under Gibbs(lam).
    """
    h = 1e-7
    mu_minus = mean_at_lambda(poly, lam - h)
    mu_center = mean_at_lambda(poly, lam)
    mu_plus = mean_at_lambda(poly, lam + h)
    return (mu_plus - 2 * mu_center + mu_minus) / (h * h)


def mu_is_concave_check(poly: list[int], lam_low: float, lam_high: float, npts: int = 50) -> bool:
    """Check if mu(lam) is concave on [lam_low, lam_high] by sampling d^2 mu / d lam^2."""
    for i in range(npts):
        lam = lam_low + (lam_high - lam_low) * (i + 0.5) / npts
        d2 = d2_mu_dlam2(poly, lam)
        if d2 > 1e-8:  # positive second derivative = convex, not concave
            return False
    return True


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=18)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    args = ap.parse_args()

    max_n = args.max_n
    geng = args.geng

    print("=" * 100)
    print("PART 1: Concavity of mu_B(lambda)")
    print("=" * 100)
    print()
    print("d(mu)/d(lam) = Var(X)/lam")
    print("If mu is CONCAVE, d^2(mu)/d(lam)^2 <= 0,")
    print("meaning the slope decreases, and the MVT bound using slope at")
    print("the LEFT endpoint (tau < lam) gives an UPPER bound on gain.")
    print()
    print("If mu is CONVEX, d^2(mu)/d(lam)^2 >= 0,")
    print("meaning slope increases, and MVT at left gives a LOWER bound on gain.")
    print()

    # Check concavity/convexity of mu for all trees
    n_convex = 0
    n_concave = 0
    n_mixed = 0
    n_total = 0
    n_concave_in_range = 0  # concave specifically in [tau, lam] interval

    for n in range(4, min(max_n + 1, 17)):
        proc = subprocess.Popen(
            [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"],
            stdout=subprocess.PIPE,
        )
        assert proc.stdout is not None

        for line in proc.stdout:
            nn, adj = parse_graph6(line)
            poly = independence_poly(nn, adj)
            n_total += 1

            # Sample d^2 mu / d lam^2 at many points
            lams = [0.1 + 0.1 * i for i in range(15)]
            d2s = [d2_mu_dlam2(poly, lam) for lam in lams]

            all_neg = all(d2 < 1e-8 for d2 in d2s)
            all_pos = all(d2 > -1e-8 for d2 in d2s)

            if all_neg and not all_pos:
                n_concave += 1
            elif all_pos and not all_neg:
                n_convex += 1
            else:
                n_mixed += 1

        proc.wait()

    print(f"Trees checked: {n_total}")
    print(f"  Concave: {n_concave}")
    print(f"  Convex:  {n_convex}")
    print(f"  Mixed:   {n_mixed}")
    print()

    # ------------------------------------------------------------------
    print("=" * 100)
    print("PART 2: mu_B(tau) baseline and slope analysis (d_leaf<=1)")
    print("=" * 100)
    print()

    # For d_leaf<=1 trees, track:
    # (a) mu_B(tau) - (m - 2)
    # (b) Var_B(tau) / tau = slope of mu_B at tau
    # (c) lambda_m - tau (the gap)
    # (d) Whether slope_at_tau * gap >= deficit

    stats = {
        "checked": 0,
        "deficit_cases": 0,
        "min_mu_B_tau_m2": float("inf"),   # mu_B(tau) - (m-2)
        "min_slope_tau": float("inf"),      # Var_B(tau)/tau
        "max_slope_tau": float("-inf"),
        "min_gap": float("inf"),
        "max_deficit": float("-inf"),
        "min_gain": float("inf"),
        "min_r2_slack": float("inf"),

        # For the concavity-based UPPER bound on gain:
        # If mu is concave, gain <= slope_at_tau * gap
        # So the MVT upper bound should exceed actual gain
        "max_gain_over_mvt_ub": float("-inf"),  # gain / (slope_at_tau * gap)

        # Combined bound: mu_B(tau) + slope_at_tau * gap - (m - 3/2)
        # = [mu_B(tau) - (m-2)] + slope_at_tau * gap - 1/2
        # If this is > 0 with concavity, we're done.
        "min_combined_lb": float("inf"),

        # Alternative: mu_B(tau) + min(slope_at_lam, slope_at_tau) * gap
        # For concave functions, slope decreases, so min is at lam (right).
        # For the lower bound we'd use slope at right endpoint for concave fn.
        "min_combined_lb_concave": float("inf"),  # using slope at lam
    }

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

            mu_B_lam = mean_at_lambda(poly_B, lam)
            mu_B_tau = mean_at_lambda(poly_B, tau)
            var_B_tau = variance_at_lambda(poly_B, tau)
            var_B_lam = variance_at_lambda(poly_B, lam)

            slope_tau = var_B_tau / tau if tau > 0 else 0
            slope_lam = var_B_lam / lam if lam > 0 else 0

            gap = lam - tau
            gain = mu_B_lam - mu_B_tau
            deficit = (m - 1.5) - mu_B_tau
            r2_slack = mu_B_lam - (m - 1.5)

            stats["checked"] += 1
            n_ct += 1
            stats["min_mu_B_tau_m2"] = min(stats["min_mu_B_tau_m2"], mu_B_tau - (m - 2))
            stats["min_slope_tau"] = min(stats["min_slope_tau"], slope_tau)
            stats["max_slope_tau"] = max(stats["max_slope_tau"], slope_tau)
            stats["min_gap"] = min(stats["min_gap"], gap)
            stats["min_r2_slack"] = min(stats["min_r2_slack"], r2_slack)

            if deficit > 1e-12:
                stats["deficit_cases"] += 1
                stats["max_deficit"] = max(stats["max_deficit"], deficit)
                stats["min_gain"] = min(stats["min_gain"], gain)

                # MVT analysis with concavity:
                # For concave mu, gain >= slope_at_lam * gap (slope at RIGHT endpoint)
                mvt_ub = slope_tau * gap  # upper bound if concave
                mvt_lb = slope_lam * gap   # lower bound if concave
                if mvt_ub > 0:
                    stats["max_gain_over_mvt_ub"] = max(
                        stats["max_gain_over_mvt_ub"], gain / mvt_ub
                    )

                # Combined: using slope at tau (upper bound if concave)
                combined = mu_B_tau + mvt_ub - (m - 1.5)
                stats["min_combined_lb"] = min(stats["min_combined_lb"], combined)

                # Using slope at lam (lower bound if concave)
                combined_c = mu_B_tau + mvt_lb - (m - 1.5)
                stats["min_combined_lb_concave"] = min(stats["min_combined_lb_concave"], combined_c)

        proc.wait()
        print(f"n={n:2d}: {n_ct:8d} d_leaf<=1 trees")

    print()
    print(f"Total checked: {stats['checked']}")
    print(f"Deficit cases: {stats['deficit_cases']}")
    print()
    print("Baseline:")
    print(f"  min mu_B(tau) - (m-2):      {stats['min_mu_B_tau_m2']:.10f}")
    print(f"  max deficit (m-3/2 - mu_B(tau)): {stats['max_deficit']:.10f}")
    print()
    print("Slope of mu_B at tau:")
    print(f"  min slope_at_tau:            {stats['min_slope_tau']:.10f}")
    print(f"  max slope_at_tau:            {stats['max_slope_tau']:.10f}")
    print()
    print("Gap and gain:")
    print(f"  min gap (lam - tau):         {stats['min_gap']:.10f}")
    print(f"  min gain in deficit cases:   {stats['min_gain']:.10f}")
    print()
    print("Route-2 outcome:")
    print(f"  min r2_slack:                {stats['min_r2_slack']:.10f}")
    print()
    print("MVT analysis:")
    print(f"  max gain/(slope_tau*gap):    {stats['max_gain_over_mvt_ub']:.10f}")
    print(f"    (if < 1 everywhere, mu is CONVEX and slope_tau*gap is UPPER bound)")
    print(f"    (if > 1 somewhere, slope_tau*gap underestimates)")
    print(f"  min combined (tau+slope_tau*gap): {stats['min_combined_lb']:.10f}")
    print(f"  min combined concave (tau+slope_lam*gap): {stats['min_combined_lb_concave']:.10f}")

    # ------------------------------------------------------------------
    print()
    print("=" * 100)
    print("PART 3: Mode-tie consequence -- mu_B(tau) at the B-tie-point")
    print("=" * 100)
    print()
    print("At tau = b_{m-2}/b_{m-1}, mu_B(tau) is the weighted mean of B")
    print("at the fugacity where b_{m-2}*tau^{m-2} = b_{m-1}*tau^{m-1} in weight.")
    print("This is the tie-point between levels m-2 and m-1 of B.")
    print()
    print("Fact: If mode(B) >= m-1 (verified 0 failures), then b_{m-1} >= b_{m-2},")
    print("so tau = b_{m-2}/b_{m-1} <= 1.")
    print()
    print("Claim: mu_B(tau) >= m - 2 + 1/2 always (i.e., deficit <= 0.5).")
    print(f"  Observed: min mu_B(tau) - (m-2) = {stats['min_mu_B_tau_m2']:.10f}")
    print(f"  This means mu_B(tau) - (m-3/2) >= {stats['min_mu_B_tau_m2'] - 0.5:.10f}")

    # ------------------------------------------------------------------
    print()
    print("=" * 100)
    print("PART 4: Direct algebraic analysis via four-term identity")
    print("=" * 100)
    print()

    # The key algebraic identity (proved in bridge notes):
    # Phi_m(T; lam) = (1+lam)*Phi_m(B; lam) + lam*Z_B(lam)
    #               + lam*Phi_{m-1}(P; lam) + lam*Phi_{m-1}(B; lam)
    #
    # = Z_B * [(1+lam)(mu_B - m+1) + lam + lam*(mu_B - m+2)]
    #   + lam*Z_P*(mu_P - m+2)
    #
    # Let me compute the B-part coefficient of Z_B:
    # (1+lam)(mu_B - m + 1) + lam + lam(mu_B - m + 2)
    # = (1+lam)(mu_B - m + 1) + lam(mu_B - m + 2) + lam
    # = (1+2lam)(mu_B) - (1+lam)(m-1) - lam(m-2) + lam
    # = (1+2lam)*mu_B - (m-1) - lam(m-1) - lam(m-2) + lam
    # = (1+2lam)*mu_B - (m-1) - lam(2m-3) + lam
    # = (1+2lam)*mu_B - (m-1) - 2lam*m + 3lam + lam
    # = (1+2lam)*mu_B - (m-1) - 2lam*m + 4lam
    #
    # Hmm, let me just verify numerically.

    print("For the four-term identity, the Route-2 condition mu_B >= m - 3/2")
    print("is exactly what makes the Phi_m(A) term nonneg (ignoring Phi_{m-1}(P)).")
    print()
    print("The sufficient condition from Phi_m(A) >= 0 alone:")
    print("  (1+lam)*(mu_B - (m-1))*Z_B + lam*Z_B + lam*Z_P*(mu_P-(m-2)) >= 0")
    print("  = Z_B*[(1+lam)(mu_B - m + 1) + lam] + lam*Z_P*(mu_P - m + 2)")
    print()
    print("Since Phi_{m-1}(P) >= 0 (verified, 0 failures), the condition reduces to:")
    print("  (1+lam)(mu_B - m + 1) + lam >= 0")
    print("  mu_B >= m - 1 - lam/(1+lam)")
    print("  mu_B >= m - 3/2 + 1/2 - lam/(1+lam)")
    print()
    print("Since 1/2 - lam/(1+lam) = 1/(2(1+lam)) > 0 for all lam > 0,")
    print("this is STRONGER than mu_B >= m - 3/2.")
    print()
    print("But we also have the FOURTH term: lam*Phi_{m-1}(B) = lam*Z_B*(mu_B - m + 2).")
    print("If mu_B >= m - 2 (verified, min slack 0.416), this term is nonnegative!")
    print()
    print("So the full Phi_m(T) >= 0 condition with ALL four terms is:")
    print("  Z_B * [(1+lam)(mu_B - m+1) + lam + lam(mu_B - m+2)] + lam*Z_P*(mu_P-m+2) >= 0")
    print("  Z_B * [(1+2lam)*mu_B - (1+lam)(m-1) - lam(m-2) + lam] + lam*Z_P*(mu_P-m+2) >= 0")
    print()
    print("The B coefficient simplifies to:")
    print("  (1+2lam)*mu_B - (m-1)(1+lam) - lam(m-2) + lam")
    print("  = (1+2lam)*mu_B - m + 1 - lam*m + lam - lam*m + 2*lam + lam")
    print("  = (1+2lam)*mu_B - m + 1 - 2*lam*m + 4*lam")
    print("  = (1+2lam)*(mu_B - m) + (1 + 4lam)")
    print()
    print("So Phi_m(T) = Z_B * [(1+2lam)(mu_B - m) + 1 + 4lam] + lam*Z_P*(mu_P - m + 2)")
    print()
    print("Setting this >= 0 and ignoring the P term (which is >= 0):")
    print("  (1+2lam)(mu_B - m) + 1 + 4lam >= 0")
    print("  mu_B >= m - (1+4lam)/(1+2lam)")

    for lam_test in [0.5, 0.7, 0.9, 1.0]:
        bound = lam_test + (1 + 4 * lam_test) / (1 + 2 * lam_test)
        # mu_B >= m - (1+4lam)/(1+2lam), which is m - bound from m
        req = - (1 + 4 * lam_test) / (1 + 2 * lam_test)
        print(f"  lam={lam_test}: mu_B >= m + {req:.4f} = m - {-req:.4f}")

    print()
    print("Compare with Route-2 target: mu_B >= m - 3/2 = m - 1.5")
    print("The four-term bound gives mu_B >= m - (1+4lam)/(1+2lam).")
    print("  At lam=1: m - 5/3 = m - 1.6667 (WEAKER than m - 1.5)")
    print("  At lam=0.5: m - 3/2 = m - 1.5 (SAME)")
    print("  At lam=0: m - 1 (STRONGER)")
    print()
    print("So the four-term identity gives a WEAKER sufficient condition when lam > 0.5!")
    print("The pendant bonus (lam*Z_B) makes the biggest difference in the Phi_m(A) term.")
    print("Route-2 uses only Phi_m(A) >= 0, not Phi_m(A) + lam*Phi_{m-1}(B) >= 0.")


if __name__ == "__main__":
    main()
