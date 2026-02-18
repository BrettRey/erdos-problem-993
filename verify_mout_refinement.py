#!/usr/bin/env python3
"""Check if m_out refinement closes the 1-Private gap analytically.

Part 1 gives k >= (n+1)/3 using m_out >= 0.
If m_out >= 1 (internal edges in V\S exist), we get k >= (n+2)/3.
If m_out = 0 (S is a color class), then k >= floor(n/2).

Question: does this fully close the gap?
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


def is_1_private(s, n, nbr):
    for u in s:
        priv_count = 0
        for v in nbr[u]:
            if v not in s:
                s_nbrs = nbr[v] & s
                if len(s_nbrs) == 1:
                    priv_count += 1
        if priv_count != 1:
            return False
    return True


def compute_mout(s, n, adj, nbr):
    """Count edges within V\S."""
    vs = set(range(n)) - s
    count = 0
    for u in range(n):
        for v in adj[u]:
            if u < v and u in vs and v in vs:
                count += 1
    return count


def main():
    print("m_out REFINEMENT CHECK")
    print("=" * 70)

    for n in range(5, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_trees = len(lines)
        n_1priv_total = 0
        mout_zero = 0
        mout_positive = 0
        max_mode = 0

        for line in lines:
            tn, adj_data = parse_graph6(line)
            nbr = [set(adj_data[v]) for v in range(tn)]
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])
            max_mode = max(max_mode, mode)

            all_mis = find_all_maximal_is(tn, adj_data)
            for s in all_mis:
                if is_1_private(s, tn, nbr):
                    n_1priv_total += 1
                    k = len(s)
                    mout = compute_mout(s, tn, adj_data, nbr)
                    if mout == 0:
                        mout_zero += 1
                    else:
                        mout_positive += 1

        elapsed = time.time() - t0
        bound_mout0 = n // 2  # color class: k >= floor(n/2)
        bound_mout1 = -(-( n + 2) // 3)  # ceil((n+2)/3)
        print(f"n={n:2d}: {n_trees:6d} trees, {n_1priv_total:5d} 1-Private IS, "
              f"mout=0: {mout_zero:4d}, mout>=1: {mout_positive:4d}, "
              f"max_mode={max_mode}, "
              f"bound(m=0)={bound_mout0}, bound(m>=1)={bound_mout1} "
              f"({elapsed:.1f}s)")

    print()
    print("If max_mode <= bound for ALL 1-Private IS cases, gap is closed.")
    print("bound(m=0): k >= floor(n/2) (color class)")
    print("bound(m>=1): k >= ceil((n+2)/3)")


if __name__ == "__main__":
    main()
