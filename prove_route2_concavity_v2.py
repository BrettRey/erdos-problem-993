#!/usr/bin/env python3
"""Concavity: does mixing (I(T) = dp[0] + dp[1]) preserve kappa_3 <= kappa_2?

Key finding from v1: kappa_3 <= kappa_2 does NOT hold for general LC sequences.
It DOES hold for all IS polynomials of trees tested (n <= 16, 32,505 trees).
The max ratio kappa_3/kappa_2 is < 0.961 for all tree families checked.

The product (convolution) of polynomials PRESERVES kappa_3 <= kappa_2
because cumulants add for independent random variables.

The crucial question: does the mixture I(T) = dp[v][0] + dp[v][1] preserve it?

This is a mixing of two distributions, not a convolution. Cumulants of mixtures
are NOT simply the sums. So this needs a separate analysis.
"""

from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def moments_at_lambda(poly: list[int], lam: float):
    """Compute mu, kappa_2, kappa_3 at fugacity lambda."""
    z = 0.0
    s1 = 0.0
    s2 = 0.0
    s3 = 0.0
    p = 1.0
    for k in range(len(poly)):
        w = poly[k] * p
        z += w
        s1 += k * w
        s2 += k * k * w
        s3 += k * k * k * w
        p *= lam
    if z == 0:
        return 0.0, 0.0, 0.0
    mu = s1 / z
    mu2 = s2 / z
    mu3 = s3 / z
    kappa_2 = mu2 - mu * mu
    kappa_3 = mu3 - 3 * mu * mu2 + 2 * mu ** 3
    return mu, kappa_2, kappa_3


def polyadd(a, b):
    out = [0] * max(len(a), len(b))
    for i in range(len(a)):
        out[i] += a[i]
    for i in range(len(b)):
        out[i] += b[i]
    return out


def main():
    print("=" * 80)
    print("Does mixing (P + Q) preserve kappa_3 <= kappa_2?")
    print("=" * 80)
    print()

    # The IS polynomial DP gives:
    #   dp[v][0] = product of (dp[c][0] + dp[c][1]) over children c
    #   dp[v][1] = x * product of dp[c][0] over children c
    #   I(subtree_v) = dp[v][0] + dp[v][1]
    #
    # At each vertex, I(T_v) is the sum of dp[v][0] and dp[v][1].
    # dp[v][0] is a product (convolution) -- preserves kappa_3 <= kappa_2.
    # dp[v][1] = x * product -- shift by 1 preserves it, product preserves it.
    #
    # But I(T_v) = dp[0] + dp[1] is a MIXTURE.
    # For a mixture f = g + h (with nonneg coefficients):
    #   Z_f = Z_g + Z_h
    #   mu_f = (Z_g/Z_f)*mu_g + (Z_h/Z_f)*mu_h = p*mu_g + q*mu_h
    #   where p = Z_g/Z_f, q = Z_h/Z_f
    #
    # The variance of the mixture:
    #   kappa_2(f) = p*kappa_2(g) + q*kappa_2(h) + p*q*(mu_g - mu_h)^2
    #
    # The third cumulant:
    #   E_f[X^3] = p*E_g[X^3] + q*E_h[X^3]
    #   mu_f = p*mu_g + q*mu_h
    #
    #   kappa_3(f) = E_f[(X-mu_f)^3]
    #             = p*E_g[(X-mu_f)^3] + q*E_h[(X-mu_f)^3]
    #
    #   Let d_g = mu_g - mu_f = q*(mu_g - mu_h)
    #       d_h = mu_h - mu_f = -p*(mu_g - mu_h)
    #   Let D = mu_g - mu_h
    #
    #   E_g[(X-mu_f)^3] = E_g[(X-mu_g+d_g)^3]
    #                   = kappa_3(g) + 3*kappa_2(g)*d_g + d_g^3
    #
    #   Similarly for h.
    #
    #   kappa_3(f) = p*[kappa_3(g) + 3*kappa_2(g)*qD + (qD)^3]
    #             + q*[kappa_3(h) + 3*kappa_2(h)*(-pD) + (-pD)^3]
    #
    #   = p*kappa_3(g) + q*kappa_3(h)
    #     + 3pqD*[kappa_2(g) - kappa_2(h)]... hmm let me be more careful.
    #
    # Actually let me just compute directly:
    #   kappa_3(f) = p*kappa_3(g) + q*kappa_3(h)
    #              + 3*p*kappa_2(g)*qD + 3*q*kappa_2(h)*(-pD)
    #              + p*(qD)^3 + q*(-pD)^3
    #
    #   = p*kappa_3(g) + q*kappa_3(h)
    #     + 3pqD*(kappa_2(g) - kappa_2(h))
    #     + pq^3*D^3 - qp^3*D^3
    #
    #   = p*kappa_3(g) + q*kappa_3(h)
    #     + 3pqD*(kappa_2(g) - kappa_2(h))
    #     + pqD^3*(q^2 - p^2)
    #
    #   = p*kappa_3(g) + q*kappa_3(h)
    #     + 3pqD*(kappa_2(g) - kappa_2(h))
    #     + pqD^3*(q - p)(q + p)
    #
    #   Since q + p = 1:
    #   = p*kappa_3(g) + q*kappa_3(h)
    #     + 3pqD*(kappa_2(g) - kappa_2(h))
    #     + pqD^3*(q - p)

    # Let me just verify numerically for tree-like IS polynomials.

    # Test with small trees built from scratch
    from indpoly import independence_poly, _polymul

    # Build all small trees and check the ratio kappa_3/kappa_2 at the
    # dp[v][0], dp[v][1], and I(T) = dp[v][0]+dp[v][1] level.

    import subprocess
    from graph6 import parse_graph6

    geng = "/opt/homebrew/bin/geng"

    max_ratio_dp0 = 0.0
    max_ratio_dp1 = 0.0
    max_ratio_sum = 0.0
    total = 0

    for n in range(2, 17):
        proc = subprocess.Popen(
            [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"],
            stdout=subprocess.PIPE,
        )
        assert proc.stdout is not None

        for line in proc.stdout:
            nn, adj = parse_graph6(line)
            total += 1

            # Root at 0 and compute DP
            parent = [-1] * nn
            children = [[] for _ in range(nn)]
            visited = [False] * nn
            visited[0] = True
            bfs = [0]
            head = 0
            while head < len(bfs):
                v = bfs[head]
                head += 1
                for u in adj[v]:
                    if not visited[u]:
                        visited[u] = True
                        parent[u] = v
                        children[v].append(u)
                        bfs.append(u)

            # Post-order
            order = []
            stack = [(0, False)]
            while stack:
                v, done = stack.pop()
                if done:
                    order.append(v)
                    continue
                stack.append((v, True))
                for c in children[v]:
                    stack.append((c, False))

            dp0 = [None] * nn
            dp1 = [None] * nn

            for v in order:
                if not children[v]:
                    dp0[v] = [1]
                    dp1[v] = [0, 1]
                else:
                    prod_total = [1]
                    prod_excl = [1]
                    for c in children[v]:
                        c_sum = polyadd(dp0[c], dp1[c])
                        prod_total = _polymul(prod_total, c_sum)
                        prod_excl = _polymul(prod_excl, dp0[c])
                    dp0[v] = prod_total
                    dp1[v] = [0] + prod_excl

            # Check kappa_3/kappa_2 at root
            root_dp0 = dp0[0]
            root_dp1 = dp1[0]
            root_sum = polyadd(root_dp0, root_dp1)

            for i in range(100):
                lam = 0.01 + 0.01 * i
                _, k2_0, k3_0 = moments_at_lambda(root_dp0, lam)
                _, k2_1, k3_1 = moments_at_lambda(root_dp1, lam)
                _, k2_s, k3_s = moments_at_lambda(root_sum, lam)

                if k2_0 > 1e-15:
                    max_ratio_dp0 = max(max_ratio_dp0, k3_0 / k2_0)
                if k2_1 > 1e-15:
                    max_ratio_dp1 = max(max_ratio_dp1, k3_1 / k2_1)
                if k2_s > 1e-15:
                    max_ratio_sum = max(max_ratio_sum, k3_s / k2_s)

        proc.wait()
        print(f"n={n:2d}: max k3/k2 for dp0={max_ratio_dp0:.6f}, dp1={max_ratio_dp1:.6f}, sum={max_ratio_sum:.6f}")

    print()
    print(f"Total trees: {total}")
    print(f"max k3/k2 for dp[v][0]: {max_ratio_dp0:.6f}")
    print(f"max k3/k2 for dp[v][1]: {max_ratio_dp1:.6f}")
    print(f"max k3/k2 for I(T):     {max_ratio_sum:.6f}")
    print()

    if max_ratio_sum < 1:
        print("CONFIRMED: kappa_3 < kappa_2 for ALL IS polynomials of trees (n<=16).")
        print("The mixing (dp0 + dp1) PRESERVES the property for IS polynomials.")
        print()
        print("Structural explanation:")
        print("  1. For a leaf v: dp0 = [1], dp1 = [0,1].")
        print("     I = [1,1]. kappa_3 = 0 = kappa_2. (Bernoulli)")
        print()
        print("  2. Convolution: dp0[v] = prod (dp0[c] + dp1[c]) = prod I(T_c).")
        print("     By cumulant additivity, kappa_3(dp0) = sum kappa_3(I(T_c)).")
        print("     Since each kappa_3(I(T_c)) <= kappa_2(I(T_c)):")
        print("       kappa_3(dp0) <= sum kappa_2(I(T_c)) = kappa_2(dp0).")
        print("     So dp0 inherits kappa_3 <= kappa_2 from children. QED for dp0.")
        print()
        print("  3. dp1[v] = x * prod dp0[c]. Shift by 1 preserves cumulants.")
        print("     kappa_3(dp1) = sum kappa_3(dp0[c]) <= sum kappa_2(dp0[c]) = kappa_2(dp1).")
        print("     QED for dp1.")
        print()
        print("  4. I(T_v) = dp0[v] + dp1[v] is a MIXTURE.")
        print("     The mixture formula gives additional positive terms in kappa_2")
        print("     (from the mean difference), which helps maintain kappa_3 <= kappa_2.")
        print("     This needs a case-by-case analysis...")


if __name__ == "__main__":
    main()
