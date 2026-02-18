#!/usr/bin/env python3
"""Explore how μ/n changes under tree surgery operations.

Goal: understand which structural changes increase μ/n, and whether
S(2^k) is the extremal d_leaf ≤ 1 tree (maximizing μ/n).

Operations to test:
1. Edge contraction: contract an interior edge (merge two vertices)
2. Arm shortening: shorten a long arm by 1 (remove leaf, merge parent)
3. Arm splitting: split one long arm into two shorter arms at the hub

If μ/n always increases when we make the tree "more spider-like",
then S(2^k) is extremal and μ < n/3 follows from the algebraic proof.
"""

import subprocess
import sys
import time

sys.path.insert(0, ".")
from indpoly import independence_poly


def parse_graph6(s):
    s = s.strip()
    if not s:
        return 0, []
    idx = 0
    if ord(s[0]) - 63 < 63:
        n = ord(s[0]) - 63; idx = 1
    else:
        idx = 1; n = 0
        for _ in range(3):
            n = n * 64 + (ord(s[idx]) - 63); idx += 1
    adj = [[] for _ in range(n)]
    bits = []
    for c in s[idx:]:
        val = ord(c) - 63
        for b in range(5, -1, -1):
            bits.append((val >> b) & 1)
    bit_idx = 0
    for j in range(n):
        for i in range(j):
            if bit_idx < len(bits) and bits[bit_idx]:
                adj[i].append(j); adj[j].append(i)
            bit_idx += 1
    return n, adj


def poly_mean(poly):
    total = sum(poly)
    weighted = sum(k * poly[k] for k in range(len(poly)))
    return weighted / total if total else 0.0


def make_spider(arms):
    """Build a spider with given arm lengths. Returns (n, adj)."""
    n = 1 + sum(arms)  # hub + arms
    adj = [[] for _ in range(n)]
    v = 1
    for a in arms:
        # Connect arm to hub (vertex 0)
        adj[0].append(v)
        adj[v].append(0)
        for i in range(a - 1):
            adj[v].append(v + 1)
            adj[v + 1].append(v)
            v += 1
        v += 1
    return n, adj


def check_dleaf(n, adj):
    leaves = set()
    for v in range(n):
        if len(adj[v]) == 1:
            leaves.add(v)
    for v in range(n):
        if v in leaves:
            continue
        d_leaf = sum(1 for u in adj[v] if u in leaves)
        if d_leaf >= 2:
            return False
    return True


def main():
    print("MEAN MONOTONICITY UNDER TREE SURGERY", flush=True)
    print("=" * 70, flush=True)
    print(flush=True)

    # Test 1: Among spiders with fixed n, which maximizes μ/n?
    print("TEST 1: Among d_leaf<=1 spiders with fixed n, which max μ/n?",
          flush=True)
    print("-" * 70, flush=True)

    for n in [13, 17, 21, 25, 31]:
        best_mu = 0
        best_config = None

        # Generate all spider configs with this n and d_leaf <= 1
        # Arms must be >= 2 (to avoid d_leaf >= 2 at hub), except at most 1
        # arm of length 1
        total_arm = n - 1
        configs_tested = 0

        # Try different numbers of arms
        for d in range(2, total_arm + 1):
            if total_arm < d:
                break

            # Enumerate: j arms of length 1 (j <= 1 for d_leaf <= 1)
            for j in [0, 1]:
                remaining = total_arm - j
                long_arms = d - j
                if long_arms <= 0:
                    continue
                if remaining < 2 * long_arms:
                    continue  # can't make all arms length >= 2

                # Distribute remaining among long_arms, each >= 2
                # Try uniform + slight variation
                base = remaining // long_arms
                if base < 2:
                    continue
                extra = remaining - base * long_arms

                # Config: extra arms of length (base+1), rest of length base
                arms = [base + 1] * extra + [base] * (long_arms - extra)
                arms += [1] * j

                tn, tadj = make_spider(arms)
                if tn != n:
                    continue
                if not check_dleaf(tn, tadj):
                    continue

                poly = independence_poly(tn, tadj)
                mu = poly_mean(poly)
                configs_tested += 1

                if mu > best_mu:
                    best_mu = mu
                    arms_sorted = sorted(arms, reverse=True)
                    best_config = (d, arms_sorted[:5], mu)

        if best_config:
            d, arms_show, mu = best_config
            ratio = mu / (n / 3)
            print(f"  n={n}: best μ={mu:.6f}, μ/(n/3)={ratio:.6f}, "
                  f"d={d}, arms={arms_show} ({configs_tested} configs)",
                  flush=True)

    print(flush=True)

    # Test 2: For each n, compare the exhaustive-best d_leaf<=1 tree
    # against the best spider. Is the best always a spider?
    print("TEST 2: Best d_leaf<=1 tree vs best spider at each n",
          flush=True)
    print("-" * 70, flush=True)

    for n in range(8, 19):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        best_mu = 0
        best_info = None
        best_is_spider = False

        for line in lines:
            tn, adj_data = parse_graph6(line)
            if not check_dleaf(tn, adj_data):
                continue

            poly = independence_poly(tn, adj_data)
            mu = poly_mean(poly)

            if mu > best_mu:
                best_mu = mu
                deg_seq = sorted([len(adj_data[v]) for v in range(tn)],
                                 reverse=True)
                # Check if it's a spider (one vertex with degree >= 3,
                # all others degree <= 2)
                high_deg = sum(1 for d in deg_seq if d >= 3)
                best_is_spider = (high_deg <= 1)
                best_info = (deg_seq[:6], best_is_spider)

        ratio = best_mu / (n / 3)
        ds, is_sp = best_info
        sp_str = "SPIDER" if is_sp else "NOT SPIDER"
        print(f"  n={n:2d}: best_μ={best_mu:.6f}, ratio={ratio:.6f}, "
              f"deg={ds}, {sp_str}", flush=True)

    print(flush=True)

    # Test 3: For spiders, how does arm length affect μ?
    # Fix n, vary arm distribution. All arms length 2 vs mixed.
    print("TEST 3: Arm length distribution effect on μ (fixed n)", flush=True)
    print("-" * 70, flush=True)

    for n in [13, 21, 31]:
        total = n - 1
        configs = []

        # All arms length 2 (if possible)
        if total % 2 == 0:
            k = total // 2
            arms = [2] * k
            tn, tadj = make_spider(arms)
            if check_dleaf(tn, tadj):
                poly = independence_poly(tn, tadj)
                configs.append(("S(2^%d)" % k, poly_mean(poly)))

        # All length 3
        if total % 3 == 0:
            k = total // 3
            arms = [3] * k
            tn, tadj = make_spider(arms)
            if check_dleaf(tn, tadj):
                poly = independence_poly(tn, tadj)
                configs.append(("S(3^%d)" % k, poly_mean(poly)))

        # One long arm, rest length 2
        for long_count in [1, 2]:
            rest = (total - 2 * long_count)
            if rest >= 2 and rest % 2 == 0:
                k = rest // 2
                L = total - 2 * k
                if long_count == 1 and L >= 2:
                    arms = [L] + [2] * k
                    tn, tadj = make_spider(arms)
                    if check_dleaf(tn, tadj):
                        poly = independence_poly(tn, tadj)
                        configs.append(("S(%d,2^%d)" % (L, k),
                                        poly_mean(poly)))

        # Mixed: some 2s, some 3s
        for k2 in range(0, total // 2 + 1):
            rem = total - 2 * k2
            if rem > 0 and rem % 3 == 0:
                k3 = rem // 3
                arms = [3] * k3 + [2] * k2
                tn, tadj = make_spider(arms)
                if tn == n and check_dleaf(tn, tadj):
                    poly = independence_poly(tn, tadj)
                    configs.append(("S(3^%d,2^%d)" % (k3, k2),
                                    poly_mean(poly)))

        configs.sort(key=lambda x: -x[1])
        print(f"  n={n} (total_arm={total}):", flush=True)
        for name, mu in configs[:8]:
            print(f"    {name:25s}: μ={mu:.6f}, μ/(n/3)={mu/(n/3):.6f}",
                  flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
