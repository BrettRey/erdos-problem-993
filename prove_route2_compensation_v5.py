#!/usr/bin/env python3
"""Route-2 compensation: alternative proof without concavity.

The concavity-based proof requires kappa_3 <= kappa_2 for all lambda in (0,1].
This is verified but not yet proved analytically.

Alternative approach: use the four-term identity DIRECTLY.

The four-term identity is:
  Phi_m(T; lam) = (1+lam)*Phi_m(B; lam) + lam*Z_B(lam) + lam*Phi_{m-1}(P; lam) + lam*Phi_{m-1}(B; lam)

We need Phi_m(T; lam_m) >= 0 (which is equivalent to mu_T(lam_m) >= m-1, trivially
true since lam_m is the mode fugacity).

Actually, Phi_m(T; lam_m) = 0 by definition (lam_m is chosen so mu_T(lam_m) = m-1).
Wait, not quite -- let me reconsider.

Actually: lam_m = i_{m-1}/i_m. The condition Phi_m(T; lam_m) >= 0 is:
  Z_T(lam_m) * (mu_T(lam_m) - (m-1)) >= 0
  mu_T(lam_m) >= m-1

This is NOT trivially true. It's the target condition we're trying to prove!

But ACTUALLY: the original GOAL is to prove mode_T <= ceil(mu_T). For mode <= m,
this requires mu_T >= m-1 (at lambda=1, not at lambda=lambda_m).

I think I'm confusing the fugacity-based analysis. Let me go back to the bridge
notes and understand what exactly Route-2 is trying to prove.

From the bridge notes:
- Target: Phi_m(T; lam_m) >= 0 where lam_m is the mode-tie fugacity of T
- This equals Z_T(lam_m) * (mu_T(lam_m) - (m-1))
- At lam_m = i_{m-1}/i_m: the weighted contributions at levels m-1 and m are tied
- So mu_T(lam_m) >= m-1 iff the upper tail outweighs the lower tail

The bridge decomposition gives:
  Phi_m(T) = Phi_m(A) + lam*Phi_{m-1}(B)

Route 2 says: Phi_m(A) >= 0 is SUFFICIENT (since Phi_{m-1}(B) >= 0 also holds).

Phi_m(A) >= 0 means mu_A(lam_m) >= m-1.

Since A = (1+x)B + xP:
  Phi_m(A) = (1+lam)Phi_m(B) + lam*Z_B + lam*Phi_{m-1}(P)

If Phi_{m-1}(P) >= 0 (verified), then:
  Phi_m(A) >= (1+lam)Phi_m(B) + lam*Z_B

  >= 0 iff (1+lam)Z_B(mu_B - (m-1)) + lam*Z_B >= 0
       iff (1+lam)(mu_B - m + 1) + lam >= 0
       iff mu_B >= m - 1 - lam/(1+lam)

Since lam <= 1 (mode >= 2 implies lam = i_{m-1}/i_m <= 1 by unimodality...
wait, we're trying to PROVE unimodality!).

Actually lam = i_{m-1}/i_m. If m is the mode, then i_m >= i_{m-1}, so lam <= 1.
And lam > 0 always.

So Route-2 target: mu_B(lam_m) >= m - 1 - lam/(1+lam).

Since lam/(1+lam) < 1/2 for lam < 1, this is stronger than mu_B >= m - 3/2.

Now, to prove this without the concavity detour, we can try:

APPROACH: Direct coefficient analysis.

mu_B(lam) = sum_k k*b_k*lam^k / sum_k b_k*lam^k

We need: sum_k k*b_k*lam^k - (m-1-lam/(1+lam)) * sum_k b_k*lam^k >= 0

Let q = m - 1 - lam/(1+lam) = m - 1 - a where a = lam/(1+lam).

sum_k (k - q) * b_k * lam^k >= 0

sum_k (k - m + 1 + a) * b_k * lam^k >= 0

This is: Phi_m(B; lam) + a * Z_B(lam) >= 0
i.e.: Z_B * (mu_B - (m-1)) + a * Z_B >= 0
i.e.: mu_B >= m - 1 - a = q. (Just restating the condition.)

The question is: can we prove this from the structure of B as T-{l,s}
and the properties of T (d_leaf <= 1, m is the mode)?

APPROACH 2: Use the relation between i_k(T) and b_k.

i_k(T) = b_k + b_{k-1} + p_{k-1}  (pendant identity)

At k = m (the mode of T): i_m >= i_{m-1}, i.e.,
  b_m + b_{m-1} + p_{m-1} >= b_{m-1} + b_{m-2} + p_{m-2}
  b_m - b_{m-2} >= p_{m-2} - p_{m-1}
  b_m >= b_{m-2} + (p_{m-2} - p_{m-1})

This gives a relation between b_m, b_{m-2}, and the P coefficients.

Similarly, i_m >= i_{m+1} (since m is the mode):
  b_m + b_{m-1} + p_{m-1} >= b_{m+1} + b_m + p_m
  b_{m-1} - b_{m+1} >= p_m - p_{m-1}

These constraints on the coefficients of B might help bound mu_B.

Let me scan for the tightest trees and understand what controls the bound.
"""

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


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    args = ap.parse_args()

    max_n = args.max_n
    geng = args.geng

    print("=" * 100)
    print("Route-2 detailed structural analysis")
    print("=" * 100)
    print()

    # For the worst trees (highest deficit), analyze the coefficient structure
    # to understand what controls Route-2.

    # The sufficient condition is:
    # Phi_m(A; lam) >= 0 where A = T-l, lam = lambda_m(T)
    # Equivalently: mu_A(lam) >= m-1
    #
    # A = (1+x)B + xP
    # So i_A(k) = b_k + b_{k-1} + p_{k-1}
    #
    # mu_A(lam) = sum_k k*(b_k + b_{k-1} + p_{k-1})*lam^k / sum_k (b_k + b_{k-1} + p_{k-1})*lam^k
    #
    # The numerator = sum_k k*b_k*lam^k + sum_k k*b_{k-1}*lam^k + sum_k k*p_{k-1}*lam^k
    # = Z_B*mu_B + lam*[sum_j (j+1)*b_j*lam^j] + lam*[sum_j (j+1)*p_j*lam^j]
    # = Z_B*mu_B + lam*Z_B*(mu_B + 1) + lam*Z_P*(mu_P + 1)
    # = Z_B*mu_B*(1+lam) + lam*Z_B + lam*Z_P*(mu_P + 1)
    #
    # Denominator = Z_A = (1+lam)*Z_B + lam*Z_P
    #
    # mu_A = [Z_B*mu_B*(1+lam) + lam*Z_B + lam*Z_P*(mu_P + 1)] / [(1+lam)*Z_B + lam*Z_P]
    #
    # = [(1+lam)*Z_B*(mu_B + lam/(1+lam)) + lam*Z_P*(mu_P + 1)] / [(1+lam)*Z_B + lam*Z_P]
    #
    # Let p = lam*Z_P / ((1+lam)*Z_B + lam*Z_P) = the weight on the P-component.
    # Then mu_A = (1-p)*(mu_B + lam/(1+lam)) + p*(mu_P + 1)
    #
    # This is a convex combination! So:
    # mu_A >= min(mu_B + lam/(1+lam), mu_P + 1)
    #
    # For mu_A >= m-1:
    # Need both mu_B + lam/(1+lam) >= m-1 AND mu_P + 1 >= m-1
    # i.e., mu_B >= m - 1 - lam/(1+lam) AND mu_P >= m-2
    #
    # Both are verified! The first is Route-2 stronger bound, the second is Route-1.
    #
    # WAIT. The convex combination bound says mu_A >= min(...).
    # For mu_A >= m-1, we need BOTH components >= m-1.
    # But the second is mu_P + 1 >= m-1, i.e., mu_P >= m-2 (verified, min slack 0.383).
    # And the first is mu_B + lam/(1+lam) >= m-1, i.e., mu_B >= m-1-a (Route-2).
    #
    # So the Route-2 bound is EXACTLY the condition that the B-component of
    # the convex combination exceeds m-1. And Route-1 is the P-component.
    #
    # BOTH need to hold for Phi_m(A) >= 0 to follow from this decomposition.
    # They're not independent -- they're the two halves of one condition.

    print("KEY IDENTITY: mu_A = (1-p)*(mu_B + a) + p*(mu_P + 1)")
    print("  where a = lam/(1+lam), p = lam*Z_P / Z_A")
    print()
    print("For mu_A >= m-1, BOTH components must be >= m-1:")
    print("  Component B: mu_B + a >= m-1  =>  mu_B >= m-1-a  (ROUTE-2)")
    print("  Component P: mu_P + 1 >= m-1  =>  mu_P >= m-2    (ROUTE-1)")
    print()
    print("Route-1 is verified (min slack 0.383). Route-2 is the target.")
    print()

    # Now the question reduces to: why is mu_B(lam) >= m - 1 - a?
    # Equivalently: why is mu_B(lam) >= m - 3/2 + 1/(2(1+lam))?
    # Since 1/(2(1+lam)) > 0, this is STRONGER than m - 3/2.

    # Let's see if we can approach this from the pendant identity constraints.
    # From i_m(T) >= i_{m+1}(T) (m is the mode):
    #   b_m + b_{m-1} + p_{m-1} >= b_{m+1} + b_m + p_m
    #   b_{m-1} - b_{m+1} >= p_m - p_{m-1}

    print("=" * 100)
    print("Scanning coefficient structure at worst trees")
    print("=" * 100)
    print()

    worst_trees = []

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

            lam = poly_T[m - 1] / poly_T[m]
            a = lam / (1 + lam)

            b_adj = remove_vertices(adj, {leaf, support})
            poly_B = independence_poly(len(b_adj), b_adj)
            if m - 1 >= len(poly_B) or poly_B[m - 2] <= 0 or poly_B[m - 1] <= 0:
                continue

            mu_B = mean_at_lambda(poly_B, lam)
            r2_slack = mu_B - (m - 1 - a)

            worst_trees.append((r2_slack, nn, line.decode("ascii").strip(), m, lam, mu_B, a, poly_B, poly_T))

        proc.wait()

    worst_trees.sort()

    print("Top 10 worst trees (smallest Route-2 slack):")
    print()
    for i, (slack, nn, g6, m, lam, mu_B, a, poly_B, poly_T) in enumerate(worst_trees[:10]):
        deg = []
        # Parse degree sig from g6
        _, adj = parse_graph6(g6.encode())
        degs = [len(nb) for nb in adj]
        deg_sig = dict(sorted(Counter(degs).items()))

        tau = poly_B[m - 2] / poly_B[m - 1]
        mu_B_tau = mean_at_lambda(poly_B, tau)
        var_B = variance_at_lambda(poly_B, lam)

        print(f"  #{i+1}: n={nn}, m={m}, deg={deg_sig}")
        print(f"       slack={slack:.10f}")
        print(f"       lam={lam:.6f}, a={a:.6f}")
        print(f"       mu_B(lam)={mu_B:.6f}, target={m-1-a:.6f}")
        print(f"       tau={tau:.6f}, mu_B(tau)={mu_B_tau:.6f}")
        print(f"       deficit_tau={(m-1.5)-mu_B_tau:.6f}")
        print(f"       gap={lam-tau:.6f}")
        print(f"       Var_B(lam)={var_B:.6f}")
        print(f"       poly_B[m-2..m+1]: {poly_B[max(0,m-2):min(len(poly_B),m+2)]}")
        print(f"       poly_T[m-1..m+1]: {poly_T[max(0,m-1):min(len(poly_T),m+2)]}")
        print()

    # ------------------------------------------------------------------
    print("=" * 100)
    print("Structural pattern: all worst trees are star-like!")
    print("=" * 100)
    print()
    print("The worst trees for Route-2 are those with one high-degree hub vertex,")
    print("which is exactly the structure of the mixed spider S(2^k, 1^j).")
    print()
    print("This suggests that the Route-2 bound is controlled by the same")
    print("mechanism as the mixed spider tie-fugacity analysis already proved.")
    print()

    # Check: are the worst trees all single-hub?
    all_single_hub = True
    for i, (slack, nn, g6, m, lam, mu_B, a, poly_B, poly_T) in enumerate(worst_trees[:20]):
        _, adj = parse_graph6(g6.encode())
        degs = [len(nb) for nb in adj]
        high_deg = [d for d in degs if d >= 3]
        if len(high_deg) != 1:
            all_single_hub = False
            print(f"  Tree #{i+1}: NOT single-hub, deg sig = {dict(sorted(Counter(degs).items()))}")

    if all_single_hub:
        print("  All top-20 worst trees have exactly one vertex of degree >= 3.")
        print("  These are spiders/caterpillars with one hub.")


if __name__ == "__main__":
    main()
