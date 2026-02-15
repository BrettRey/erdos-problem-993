#!/usr/bin/env python3
"""Verify: the 1-Private Mode Conjecture gap is FULLY closed by Part 1.

Gemini's Part 1 gives: any 1-Private maximal IS has k >= (n+1)/3.
Since k is an integer, k >= ceil((n+1)/3) = floor(n/3) + 1.

Part 2 handles mode <= floor(n/3) + 1 immediately (k >= mode).
Part 3 (empirical) handles mode >= floor(n/3) + 2.

This script checks: are there ANY trees where mode >= floor(n/3) + 2
AND a 1-Private maximal IS exists? If not, Part 3 is vacuously true
and the proof is COMPLETE.

Also checks: can we close the gap purely analytically?
"""

import subprocess
import sys
import time
from collections import defaultdict

sys.path.insert(0, ".")
from indpoly import independence_poly


def parse_graph6(s):
    s = s.strip()
    if not s:
        return 0, []
    idx = 0
    if ord(s[0]) - 63 < 63:
        n = ord(s[0]) - 63
        idx = 1
    else:
        idx = 1
        n = 0
        for _ in range(3):
            n = n * 64 + (ord(s[idx]) - 63)
            idx += 1
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
                adj[i].append(j)
                adj[j].append(i)
            bit_idx += 1
    return n, adj


def find_all_maximal_is(n, adj):
    """Find all maximal independent sets."""
    nbr = [set(adj[v]) for v in range(n)]
    result = []

    def backtrack(v, current, forbidden):
        if v == n:
            # Check maximality
            s = frozenset(current)
            for w in range(n):
                if w not in s and not (nbr[w] & s):
                    return  # not maximal
            result.append(s)
            return
        backtrack(v + 1, current, forbidden)
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])

    backtrack(0, [], set())
    return result


def is_1_private(s, n, nbr):
    """Check if IS s is 1-Private (every vertex has exactly 1 private neighbor)."""
    for u in s:
        priv_count = 0
        for v in nbr[u]:
            if v not in s:
                # Check if v's only S-neighbor is u
                s_nbrs = nbr[v] & s
                if len(s_nbrs) == 1:
                    priv_count += 1
        if priv_count != 1:
            return False
    return True


def main():
    print("=" * 70)
    print("VERIFYING 1-PRIVATE MODE CONJECTURE GAP CLOSURE")
    print("=" * 70)
    print()
    print("Part 1 gives: 1-Private maximal IS has k >= floor(n/3) + 1")
    print("This directly covers mode <= floor(n/3) + 1.")
    print("Gap: need mode >= floor(n/3) + 2 trees to have no 1-Private IS.")
    print()

    for n in range(5, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        threshold = n // 3 + 2  # mode >= floor(n/3) + 2

        n_trees = len(lines)
        n_high_mode = 0
        n_high_mode_with_1priv = 0
        n_1priv_total = 0
        max_mode_seen = 0

        for line in lines:
            tn, adj_data = parse_graph6(line)
            nbr = [set(adj_data[v]) for v in range(tn)]
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])
            max_mode_seen = max(max_mode_seen, mode)

            if mode >= threshold:
                n_high_mode += 1
                # Check for 1-Private maximal IS
                all_mis = find_all_maximal_is(tn, adj_data)
                for s in all_mis:
                    if is_1_private(s, tn, nbr):
                        n_1priv_total += 1
                        n_high_mode_with_1priv += 1
                        print(f"  FOUND: n={tn} mode={mode} threshold={threshold}"
                              f" |S|={len(s)} S={sorted(s)}")

        elapsed = time.time() - t0
        print(f"n={n:2d}: {n_trees:6d} trees, "
              f"max_mode={max_mode_seen}, threshold={threshold}, "
              f"high_mode={n_high_mode:5d}, "
              f"with_1priv={n_high_mode_with_1priv} "
              f"({elapsed:.1f}s)")

    print()
    print("=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print()
    print("If high_mode_with_1priv = 0 for all n, then Part 3 is vacuously true")
    print("and the 1-Private Mode Conjecture is FULLY PROVED (no empirical gap).")


if __name__ == "__main__":
    main()
