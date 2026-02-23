#!/usr/bin/env python3
"""Investigate concavity of mu(lambda) for IS polynomials.

The key question: for the Gibbs measure on {0,...,d} with weights
  w_k(lambda) = c_k * lambda^k  (c_k = IS polynomial coefficients)

is kappa_3(lambda) <= kappa_2(lambda) = Var(X; lambda) for all lambda in (0,1]?

Here kappa_3 = E[(X-mu)^3] is the third central moment (= third cumulant).

This is equivalent to d^2(mu)/d(lambda)^2 <= 0 (concavity of the mean).

For log-concave (LC) distributions on integers, there are known bounds on
skewness. We explore:

1. Whether kappa_3 <= kappa_2 holds for ALL LC distributions (not just IS polys)
2. The structural reason it holds for IS polynomials of trees
3. Connection to the product structure (IS poly of tree = product over edges)
"""

from __future__ import annotations

import math
import os
import subprocess
import sys
import itertools

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from indpoly import independence_poly


def moments_at_lambda(poly: list[int], lam: float):
    """Compute mu, kappa_2, kappa_3 at fugacity lambda."""
    z = 0.0
    s1 = 0.0
    s2 = 0.0
    s3 = 0.0
    p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p
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


def test_lc_counterexample():
    """Search for LC sequences where kappa_3 > kappa_2 at some lambda in (0,1].

    If we find one, concavity is NOT a consequence of LC alone.
    """
    print("=" * 80)
    print("TEST 1: Does kappa_3 <= kappa_2 hold for ALL LC sequences?")
    print("=" * 80)
    print()

    # Generate random LC sequences and check
    import random
    random.seed(42)

    n_tested = 0
    n_violated = 0

    # Test 1a: Binomial-like sequences (always LC)
    print("Binomial coefficients C(n,k):")
    for n in range(2, 20):
        poly = [math.comb(n, k) for k in range(n + 1)]
        for i in range(100):
            lam = 0.01 + 0.01 * i
            mu, k2, k3 = moments_at_lambda(poly, lam)
            n_tested += 1
            if k3 > k2 + 1e-10:
                n_violated += 1
                print(f"  VIOLATION: C({n},k), lam={lam:.2f}, k2={k2:.6f}, k3={k3:.6f}")
    print(f"  Tested: {n_tested}, violations: {n_violated}")
    print()

    # Test 1b: Geometric-like LC sequences: a_k = r^k for various r
    n_tested = 0
    n_violated = 0
    print("Geometric sequences a_k = r^k (always LC):")
    for d in range(3, 15):
        for r_num in range(1, 20):
            r = r_num / 10.0
            poly = [int(round(r ** k * 1000)) for k in range(d + 1)]
            # Check LC
            is_lc = True
            for k in range(1, len(poly) - 1):
                if poly[k] * poly[k] < poly[k - 1] * poly[k + 1]:
                    is_lc = False
                    break
            if not is_lc:
                continue

            for i in range(100):
                lam = 0.01 + 0.01 * i
                mu, k2, k3 = moments_at_lambda(poly, lam)
                n_tested += 1
                if k3 > k2 + 1e-10:
                    n_violated += 1
                    if n_violated <= 5:
                        print(f"  VIOLATION: r={r}, d={d}, lam={lam:.2f}, k2={k2:.6f}, k3={k3:.6f}")
    print(f"  Tested: {n_tested}, violations: {n_violated}")
    print()

    # Test 1c: Ultra-log-concave sequences: a_k / C(n,k) is LC
    n_tested = 0
    n_violated = 0
    print("Ultra-log-concave sequences (a_k/C(n,k) decreasing):")
    for n in range(4, 12):
        for trial in range(100):
            # Generate ULC sequence
            poly = [1]
            ratio = random.uniform(0.5, 1.0)
            for k in range(1, n + 1):
                next_val = int(round(poly[-1] * ratio * math.comb(n, k) / math.comb(n, k - 1)))
                if next_val < 1:
                    next_val = 1
                poly.append(next_val)
                ratio *= random.uniform(0.7, 1.0)

            for i in range(100):
                lam = 0.01 + 0.01 * i
                mu, k2, k3 = moments_at_lambda(poly, lam)
                n_tested += 1
                if k3 > k2 + 1e-10:
                    n_violated += 1
                    if n_violated <= 5:
                        print(f"  VIOLATION: n={n}, poly={poly[:6]}..., lam={lam:.2f}")
                        print(f"    k2={k2:.6f}, k3={k3:.6f}, ratio={k3/k2:.6f}")
    print(f"  Tested: {n_tested}, violations: {n_violated}")
    print()

    # Test 1d: Arbitrary LC sequences with small support
    n_tested = 0
    n_violated = 0
    print("Exhaustive LC sequences with support {0,1,2,3}:")
    for a0 in range(1, 6):
        for a1 in range(1, 20):
            for a2 in range(1, 20):
                if a1 * a1 < a0 * a2:
                    continue
                for a3 in range(1, 20):
                    if a2 * a2 < a1 * a3:
                        continue
                    poly = [a0, a1, a2, a3]
                    for i in range(100):
                        lam = 0.01 + 0.01 * i
                        mu, k2, k3 = moments_at_lambda(poly, lam)
                        n_tested += 1
                        if k3 > k2 + 1e-10:
                            n_violated += 1
                            if n_violated <= 10:
                                print(f"  VIOLATION: {poly}, lam={lam:.2f}")
                                print(f"    mu={mu:.4f}, k2={k2:.6f}, k3={k3:.6f}")
    print(f"  Tested: {n_tested}, violations: {n_violated}")
    print()

    # Test 1e: Extend to support {0,...,4}
    n_tested = 0
    n_violated = 0
    print("Exhaustive LC sequences with support {0,...,4} (sampling):")
    import random
    random.seed(123)
    for _ in range(100000):
        d = 4
        poly = [1]
        for k in range(1, d + 1):
            max_val = min(20, poly[-1] * poly[-1] // (poly[-2] if k >= 2 else poly[-1]))
            if max_val < 1:
                max_val = 1
            poly.append(random.randint(1, max_val))
        # Verify LC
        is_lc = True
        for k in range(1, len(poly) - 1):
            if poly[k] * poly[k] < poly[k - 1] * poly[k + 1]:
                is_lc = False
                break
        if not is_lc:
            continue

        for i in range(100):
            lam = 0.01 + 0.01 * i
            mu, k2, k3 = moments_at_lambda(poly, lam)
            n_tested += 1
            if k3 > k2 + 1e-10:
                n_violated += 1
                if n_violated <= 5:
                    print(f"  VIOLATION: {poly}, lam={lam:.2f}")
                    print(f"    mu={mu:.4f}, k2={k2:.6f}, k3={k3:.6f}")
    print(f"  Tested: {n_tested}, violations: {n_violated}")
    print()

    return n_violated == 0


def test_product_structure():
    """Check if the product structure of IS polynomials helps.

    For a tree, the IS polynomial has the structure:
      I(T) = dp[root][0] + dp[root][1]

    where dp[v][0] = prod_c (dp[c][0] + dp[c][1]) and dp[v][1] = x * prod_c dp[c][0].

    The IS polynomial of a star S_n is:
      I(S_n) = (1+x)^{n-1} + x

    For a path P_n:
      I(P_n) satisfies the recurrence I(P_n) = I(P_{n-1}) + x * I(P_{n-2})

    For a spider S(k1, k2, ...):
      I(spider) = prod_i I(P_{ki}) + x * prod_i dp[endpoint_i][0]

    For these specific forms, can we prove kappa_3 <= kappa_2?
    """
    print("=" * 80)
    print("TEST 2: Concavity for specific tree families")
    print("=" * 80)
    print()

    # Stars: I(S_n) = (1+x)^{n-1} + x
    print("Stars S_n:")
    for n in range(2, 30):
        poly = [0] * (n + 1)
        for k in range(n):
            poly[k] += math.comb(n - 1, k)
        poly[1] += 1  # add x term
        max_ratio = 0.0
        for i in range(100):
            lam = 0.01 + 0.01 * i
            mu, k2, k3 = moments_at_lambda(poly, lam)
            if k2 > 1e-15:
                max_ratio = max(max_ratio, k3 / k2)
        print(f"  S_{n}: max(k3/k2) = {max_ratio:.6f}" + (" OK" if max_ratio < 1 else " FAIL"))
    print()

    # Paths: I(P_n) via recurrence
    print("Paths P_n:")
    polys = {1: [1, 1], 2: [1, 2]}
    for n in range(3, 30):
        p1 = polys[n - 1]
        p2 = polys[n - 2]
        # I(P_n) = I(P_{n-1}) + x * I(P_{n-2})
        poly = list(p1) + [0] * max(0, len(p2) + 1 - len(p1))
        for k in range(len(p2)):
            poly[k + 1] += p2[k]
        polys[n] = poly

    for n in range(2, 30):
        poly = polys[n]
        max_ratio = 0.0
        for i in range(100):
            lam = 0.01 + 0.01 * i
            mu, k2, k3 = moments_at_lambda(poly, lam)
            if k2 > 1e-15:
                max_ratio = max(max_ratio, k3 / k2)
        print(f"  P_{n}: max(k3/k2) = {max_ratio:.6f}" + (" OK" if max_ratio < 1 else " FAIL"))
    print()

    # Caterpillars: spine of length L with pendant edges
    print("Caterpillars with long spine:")
    for spine_len in [5, 10, 15, 20]:
        for n_pendants in [0, 1, 2, spine_len // 2]:
            # Build caterpillar
            n = spine_len + n_pendants
            adj = [[] for _ in range(n)]
            # Spine: 0-1-2-...(spine_len-1)
            for i in range(spine_len - 1):
                adj[i].append(i + 1)
                adj[i + 1].append(i)
            # Pendants at middle vertices
            next_v = spine_len
            for i in range(n_pendants):
                target = spine_len // 2 + i % (spine_len // 2)
                if target >= spine_len:
                    target = spine_len - 1
                adj[target].append(next_v)
                adj[next_v].append(target)
                next_v += 1

            poly = independence_poly(n, adj)
            max_ratio = 0.0
            for i in range(100):
                lam = 0.01 + 0.01 * i
                mu, k2, k3 = moments_at_lambda(poly, lam)
                if k2 > 1e-15:
                    max_ratio = max(max_ratio, k3 / k2)
            print(f"  Cat(spine={spine_len}, pend={n_pendants}): max(k3/k2) = {max_ratio:.6f}"
                  + (" OK" if max_ratio < 1 else " FAIL"))
    print()


def test_convolution_preservation():
    """Does kappa_3 <= kappa_2 survive polynomial multiplication?

    If P and Q both satisfy kappa_3 <= kappa_2, does P*Q?

    For independent random variables X ~ P, Y ~ Q:
      kappa_2(X+Y) = kappa_2(X) + kappa_2(Y)
      kappa_3(X+Y) = kappa_3(X) + kappa_3(Y)

    So if kappa_3(X) <= kappa_2(X) and kappa_3(Y) <= kappa_2(Y), then
    kappa_3(X+Y) <= kappa_2(X+Y).

    This means the product of polynomials preserving kappa_3 <= kappa_2
    also preserves it! This is the KEY structural fact for IS polynomials.
    """
    print("=" * 80)
    print("TEST 3: Convolution preservation of kappa_3 <= kappa_2")
    print("=" * 80)
    print()
    print("If X ~ P and Y ~ Q are independent, then:")
    print("  kappa_2(X+Y) = kappa_2(X) + kappa_2(Y)")
    print("  kappa_3(X+Y) = kappa_3(X) + kappa_3(Y)")
    print()
    print("So kappa_3 <= kappa_2 for P and Q implies kappa_3 <= kappa_2 for P*Q.")
    print("This means the property is PRESERVED under convolution (polynomial multiply).")
    print()
    print("For IS polynomials of trees:")
    print("  P = prod_c I(subtree_c)  (product of subtree IS polynomials)")
    print("  If each I(T_c) satisfies kappa_3 <= kappa_2, then so does P.")
    print()
    print("And I(B) = P + Q where Q = x * prod_c dp[c][0].")
    print("Mixing preserves kappa_3 <= kappa_2 if both P and Q have it.")
    print("But Q = x * (product), which is a shifted product.")
    print("The shift by 1 doesn't affect the cumulants (translation invariance).")
    print("Wait, it does affect them -- cumulants of X+1:")
    print("  kappa_2(X+1) = kappa_2(X)")
    print("  kappa_3(X+1) = kappa_3(X)")
    print("So shifting by 1 preserves the inequality!")
    print()
    print("Therefore: if all subtree IS polynomials satisfy kappa_3 <= kappa_2,")
    print("then both dp[v][0] and dp[v][1] do, and therefore I(B) = P + Q")
    print("also does (as a mixture -- need to check mixture separately).")
    print()

    # Check: does MIXING preserve kappa_3 <= kappa_2?
    # I(B) = P + Q is NOT a product; it's a mixture.
    # For a mixture f = alpha*g + beta*h:
    #   kappa_2(f) = alpha*Var(g) + beta*Var(h) + alpha*beta*(mu_g - mu_h)^2
    #   kappa_3(f) = ... (more complex)
    #
    # The answer is: mixture does NOT necessarily preserve kappa_3 <= kappa_2.
    # But the specific structure I(B) = P + Q might.

    print("Testing: does mixing preserve kappa_3 <= kappa_2?")
    print("Take P = [1,2,1] and Q = [0,1,2,1] (both satisfy kappa_3 <= kappa_2).")

    P = [1, 2, 1]
    Q = [0, 1, 2, 1]
    PQ = [P[k] + (Q[k] if k < len(Q) else 0) for k in range(max(len(P), len(Q)))]

    for lam_100 in range(1, 101):
        lam = lam_100 / 100
        mu_p, k2_p, k3_p = moments_at_lambda(P, lam)
        mu_q, k2_q, k3_q = moments_at_lambda(Q, lam)
        mu_pq, k2_pq, k3_pq = moments_at_lambda(PQ, lam)

        if k3_pq > k2_pq + 1e-10:
            print(f"  Mixture violation at lam={lam:.2f}!")
            break
    else:
        print("  No mixture violations for this example.")
    print()

    # Check many random mixtures
    import random
    random.seed(42)
    mix_violations = 0
    mix_tested = 0

    for _ in range(10000):
        d1 = random.randint(2, 6)
        d2 = random.randint(2, 6)
        # Generate LC sequences
        P = [random.randint(1, 10)]
        for k in range(1, d1):
            max_next = P[-1] * P[-1] // (P[-2] if k >= 2 else P[-1])
            P.append(random.randint(1, max(1, max_next)))
        Q = [0] + [random.randint(1, 10)]
        for k in range(2, d2 + 1):
            max_next = Q[-1] * Q[-1] // (Q[-2] if Q[-2] > 0 else Q[-1])
            Q.append(random.randint(1, max(1, max_next)))

        # Check both are LC
        ok = True
        for seq in [P, Q]:
            for k in range(1, len(seq) - 1):
                if seq[k] > 0 and seq[k] * seq[k] < seq[k - 1] * seq[k + 1]:
                    ok = False
                    break
        if not ok:
            continue

        # Check both satisfy kappa_3 <= kappa_2
        both_ok = True
        for seq in [P, Q]:
            for lam_100 in range(1, 101):
                lam = lam_100 / 100
                _, k2, k3 = moments_at_lambda(seq, lam)
                if k3 > k2 + 1e-10:
                    both_ok = False
                    break
            if not both_ok:
                break
        if not both_ok:
            continue

        # Check mixture
        PQ = [0] * max(len(P), len(Q))
        for k in range(len(P)):
            PQ[k] += P[k]
        for k in range(len(Q)):
            PQ[k] += Q[k]

        mix_tested += 1
        for lam_100 in range(1, 101):
            lam = lam_100 / 100
            _, k2_pq, k3_pq = moments_at_lambda(PQ, lam)
            if k3_pq > k2_pq + 1e-10:
                mix_violations += 1
                print(f"  Mixture VIOLATION: P={P}, Q={Q}, lam={lam:.2f}")
                print(f"    k2={k2_pq:.6f}, k3={k3_pq:.6f}")
                break

    print(f"  Random mixtures tested: {mix_tested}, violations: {mix_violations}")


def main():
    all_lc_ok = test_lc_counterexample()
    if all_lc_ok:
        print("\n*** kappa_3 <= kappa_2 holds for ALL tested LC sequences! ***")
        print("This suggests it may be a general property of LC distributions.")
    else:
        print("\n*** kappa_3 <= kappa_2 does NOT hold for all LC sequences. ***")
        print("Need a more specific property of IS polynomials.")

    print()
    test_product_structure()
    test_convolution_preservation()


if __name__ == "__main__":
    main()
