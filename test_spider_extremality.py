#!/usr/bin/env python3
"""Test spider extremality: is μ(spider) >= μ(T) for all d_leaf <= 1
trees T on the same number of vertices?

For each n, find:
1. The d_leaf ≤ 1 tree maximizing μ
2. The d_leaf ≤ 1 spider maximizing μ
3. Check they're the same (or spider is at least as good)

Also test: does collapsing the core (contracting internal edges toward a
single hub) always increase μ?
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


def is_spider(n, adj):
    """Check if tree is a spider (at most one vertex of degree >= 3)."""
    high_deg = sum(1 for v in range(n) if len(adj[v]) >= 3)
    return high_deg <= 1


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


def get_arm_lengths(n, adj):
    """For a spider, return the arm lengths."""
    # Find hub (highest degree vertex)
    hub = max(range(n), key=lambda v: len(adj[v]))
    if len(adj[hub]) < 2:
        return [n - 1]  # path

    arms = []
    for start in adj[hub]:
        length = 1
        prev, curr = hub, start
        while len(adj[curr]) == 2:
            nxt = [u for u in adj[curr] if u != prev][0]
            prev, curr = curr, nxt
            length += 1
        arms.append(length)
    return sorted(arms, reverse=True)


def main():
    print("SPIDER EXTREMALITY TEST", flush=True)
    print("=" * 70, flush=True)
    print(flush=True)
    print("Question: Among d_leaf<=1 trees on n vertices,", flush=True)
    print("do spiders always maximize μ?", flush=True)
    print(flush=True)

    for n in range(5, 21):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        best_overall_mu = 0
        best_overall_info = None
        best_spider_mu = 0
        best_spider_info = None
        n_dleaf1 = 0

        for line in lines:
            tn, adj_data = parse_graph6(line)
            if not check_dleaf(tn, adj_data):
                continue

            n_dleaf1 += 1
            poly = independence_poly(tn, adj_data)
            mu = poly_mean(poly)
            spider = is_spider(tn, adj_data)

            if mu > best_overall_mu:
                best_overall_mu = mu
                ds = sorted([len(adj_data[v]) for v in range(tn)],
                            reverse=True)
                best_overall_info = (ds[:6], spider)

            if spider and mu > best_spider_mu:
                best_spider_mu = mu
                arms = get_arm_lengths(tn, adj_data)
                best_spider_info = (arms, mu)

        elapsed = time.time() - t0
        ratio = best_overall_mu / (n / 3)
        spider_ratio = best_spider_mu / (n / 3)
        _, is_sp = best_overall_info
        gap = best_spider_mu - best_overall_mu

        status = "SPIDER=BEST" if abs(gap) < 1e-12 else (
            "SPIDER BETTER" if gap > 0 else "NON-SPIDER BETTER!")

        arms_str = ""
        if best_spider_info:
            arms_str = f"arms={best_spider_info[0]}"

        print(f"n={n:2d}: {n_dleaf1:5d} trees | "
              f"best_μ={best_overall_mu:.6f} (ratio={ratio:.4f}) "
              f"spider?={is_sp} | "
              f"best_spider_μ={best_spider_mu:.6f} | "
              f"{status} | {arms_str} ({elapsed:.1f}s)", flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
