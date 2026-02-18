#!/usr/bin/env python3
"""Investigate why high-mode trees don't admit 1-Private IS below mode.

For trees with mode > floor(n/3)+1, check ALL maximal IS (not just below mode)
to see if 1-Private IS exist at all, and if so, at what sizes.

Key question: Is the issue that:
(a) 1-Private IS don't exist at ALL in high-mode trees? Or
(b) 1-Private IS exist but only at sizes >= mode?
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


def find_all_maximal_is(n, adj):
    nbr = [set(adj[v]) for v in range(n)]
    result = []
    def backtrack(v, current, forbidden):
        if v == n:
            s = frozenset(current)
            for w in range(n):
                if w not in s and not (nbr[w] & s):
                    return
            result.append(s)
            return
        backtrack(v + 1, current, forbidden)
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])
    backtrack(0, [], set())
    return result


def main():
    print("HIGH-MODE TREES: 1-PRIVATE IS ANALYSIS", flush=True)
    print("=" * 60, flush=True)
    print(flush=True)

    for n in range(5, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        threshold = n // 3 + 1  # floor(n/3) + 1

        n_highmode_trees = 0
        n_1priv_any = 0       # 1-Private IS at any size in high-mode trees
        n_1priv_below = 0     # 1-Private IS below mode in high-mode trees
        n_1priv_at_mode = 0   # 1-Private IS at mode
        n_1priv_above = 0     # 1-Private IS above mode
        min_1priv_k = None
        max_mode_highmode = 0
        # Track: what sizes do 1-Private IS appear at?
        size_dist = {}

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])

            if mode <= threshold:
                continue
            n_highmode_trees += 1
            max_mode_highmode = max(max_mode_highmode, mode)

            nbr = [set(adj_data[v]) for v in range(tn)]
            all_mis = find_all_maximal_is(tn, adj_data)

            for s in all_mis:
                k = len(s)
                # Check all priv <= 1
                all_leq1 = True
                for u in s:
                    priv_u = 0
                    for v in adj_data[u]:
                        if v not in s and len(nbr[v] & s) == 1:
                            priv_u += 1
                    if priv_u > 1:
                        all_leq1 = False
                        break
                if not all_leq1:
                    continue

                # This IS has all priv <= 1
                n_1priv_any += 1
                size_dist[k] = size_dist.get(k, 0) + 1

                if k < mode:
                    n_1priv_below += 1
                elif k == mode:
                    n_1priv_at_mode += 1
                else:
                    n_1priv_above += 1

                if min_1priv_k is None or k < min_1priv_k:
                    min_1priv_k = k

        elapsed = time.time() - t0

        if n_highmode_trees > 0:
            print(f"n={n:2d}: {n_highmode_trees} high-mode trees (mode>{threshold}), "
                  f"max_mode={max_mode_highmode}", flush=True)
            if n_1priv_any > 0:
                print(f"       1-Priv IS: {n_1priv_any} total "
                      f"(below={n_1priv_below}, at={n_1priv_at_mode}, above={n_1priv_above}), "
                      f"min_k={min_1priv_k}", flush=True)
                sorted_sizes = sorted(size_dist.items())
                print(f"       Size dist: {dict(sorted_sizes)}", flush=True)
            else:
                print(f"       NO 1-Private IS at any size!", flush=True)
        else:
            print(f"n={n:2d}: 0 high-mode trees ({elapsed:.1f}s)", flush=True)

        print(f"       ({elapsed:.1f}s)", flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
