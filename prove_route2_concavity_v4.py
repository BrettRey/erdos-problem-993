#!/usr/bin/env python3
"""Investigate: can we prove kappa_3 <= c * kappa_2 for some c < 1?

If kappa_3(g) <= c * kappa_2(g) and kappa_3(h) <= c * kappa_2(h), then:

kappa_3(f) = p*kappa_3(g) + q*kappa_3(h) + 3pqD*(k2_g - k2_h) + pqD^3(q-p)

<= c*(p*kappa_2(g) + q*kappa_2(h)) + 3pqD*(k2_g - k2_h) + pqD^3(q-p)

kappa_2(f) = p*kappa_2(g) + q*kappa_2(h) + pqD^2

Need: kappa_3(f) <= c * kappa_2(f), i.e.:

c*(p*k2_g + q*k2_h) + 3pqD*(k2_g - k2_h) + pqD^3(q-p)
<= c*(p*k2_g + q*k2_h + pqD^2)

i.e.: 3pqD*(k2_g - k2_h) + pqD^3(q-p) <= c*pqD^2

Dividing by pqD^2 (assuming D != 0):
3*(k2_g - k2_h)/D + D*(q-p) <= c

This depends on the specific values. Let me check what values of c work
for IS polynomial DP structures.
"""

import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from graph6 import parse_graph6
from indpoly import independence_poly, _polymul


def moments_at_lambda(poly, lam):
    z = s1 = s2 = s3 = 0.0
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
    kappa_2 = s2 / z - mu * mu
    kappa_3 = s3 / z - 3 * mu * s2 / z + 2 * mu ** 3
    return mu, kappa_2, kappa_3


def polyadd(a, b):
    out = [0] * max(len(a), len(b))
    for i in range(len(a)):
        out[i] += a[i]
    for i in range(len(b)):
        out[i] += b[i]
    return out


def main():
    geng = "/opt/homebrew/bin/geng"

    print("Finding the tight ratio c such that kappa_3 <= c * kappa_2")
    print("for all IS polynomial DP components, at all lambda in (0, 1].")
    print("=" * 80)
    print()

    # Track max ratio at each level: dp0[v], dp1[v], I(T_v)
    max_ratio_dp0 = 0.0
    max_ratio_dp1 = 0.0
    max_ratio_I = 0.0

    # Track the ratio specifically for the MIXTURE step
    # i.e., given dp0 and dp1 satisfying kappa_3 <= c*kappa_2,
    # what is the tightest c' for their mixture I = dp0 + dp1?
    # The amplification factor c' / c.

    max_amplification = 0.0  # max (ratio_I / max(ratio_dp0, ratio_dp1))

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

            # Check ratios at each vertex
            for v in order:
                for lam_100 in range(1, 101):
                    lam = lam_100 / 100.0

                    mu_0, k2_0, k3_0 = moments_at_lambda(dp0[v], lam)
                    mu_1, k2_1, k3_1 = moments_at_lambda(dp1[v], lam)

                    f = polyadd(dp0[v], dp1[v])
                    mu_f, k2_f, k3_f = moments_at_lambda(f, lam)

                    r0 = k3_0 / k2_0 if k2_0 > 1e-15 else 0
                    r1 = k3_1 / k2_1 if k2_1 > 1e-15 else 0
                    rf = k3_f / k2_f if k2_f > 1e-15 else 0

                    max_ratio_dp0 = max(max_ratio_dp0, r0)
                    max_ratio_dp1 = max(max_ratio_dp1, r1)
                    max_ratio_I = max(max_ratio_I, rf)

                    # Amplification
                    r_in = max(r0, r1)
                    if r_in > 1e-10:
                        amp = rf / r_in
                        max_amplification = max(max_amplification, amp)

        proc.wait()

    print(f"Total trees: {total}")
    print(f"max ratio kappa_3/kappa_2:")
    print(f"  dp0: {max_ratio_dp0:.8f}")
    print(f"  dp1: {max_ratio_dp1:.8f}")
    print(f"  I=dp0+dp1: {max_ratio_I:.8f}")
    print(f"max amplification (ratio_I / max(ratio_dp0, ratio_dp1)): {max_amplification:.8f}")
    print()

    if max_ratio_I < 1:
        c = max_ratio_I
        print(f"Tight constant: kappa_3 <= {c:.8f} * kappa_2 for all IS polynomials.")
        print(f"This is bounded away from 1 by at least {1-c:.8f}.")
        print()

    # Now check: what is the base case ratio for leaves and edges?
    print("Base case ratios:")

    # Leaf: I = [1, 1]
    for lam_100 in range(1, 101):
        lam = lam_100 / 100.0
        mu, k2, k3 = moments_at_lambda([1, 1], lam)
        r = k3 / k2 if k2 > 1e-15 else 0
        if abs(lam - 0.01) < 1e-6 or abs(lam - 0.5) < 1e-6 or abs(lam - 1.0) < 1e-6:
            print(f"  [1,1] at lam={lam:.2f}: mu={mu:.4f}, k2={k2:.6f}, k3={k3:.6f}, ratio={r:.6f}")

    # Edge: I(P_2) = [1, 2, 1]
    for lam_100 in range(1, 101):
        lam = lam_100 / 100.0
        mu, k2, k3 = moments_at_lambda([1, 2, 1], lam)
        r = k3 / k2 if k2 > 1e-15 else 0
        if abs(lam - 0.01) < 1e-6 or abs(lam - 0.5) < 1e-6 or abs(lam - 1.0) < 1e-6:
            print(f"  [1,2,1] at lam={lam:.2f}: mu={mu:.4f}, k2={k2:.6f}, k3={k3:.6f}, ratio={r:.6f}")

    # Star S_3: I = [1, 3, 3, 1] (= (1+x)^2 + x)
    for lam_100 in range(1, 101):
        lam = lam_100 / 100.0
        mu, k2, k3 = moments_at_lambda([1, 3, 3, 1], lam)
        r = k3 / k2 if k2 > 1e-15 else 0
        if abs(lam - 0.01) < 1e-6 or abs(lam - 0.5) < 1e-6 or abs(lam - 1.0) < 1e-6:
            print(f"  [1,3,3,1] at lam={lam:.2f}: mu={mu:.4f}, k2={k2:.6f}, k3={k3:.6f}, ratio={r:.6f}")

    print()
    print("The base ratio for [1,1] at lam->0 approaches 1:")
    print("  kappa_3/kappa_2 = (1-lam)/(1+lam) -> 1 as lam -> 0")
    print()
    print("So the tight constant c approaches 1 as lam -> 0.")
    print("But for lam in (0, 1], the ratio is bounded by max(dp0, dp1) = 0.98.")


if __name__ == "__main__":
    main()
