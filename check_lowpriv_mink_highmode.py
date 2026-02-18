#!/usr/bin/env python3
"""For high-mode trees: what is the minimum k among maximal IS with all priv <= 1?

Key question: is min_k >= mode for high-mode trees?
If so, then low-priv IS exist but are never below mode.
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
        for b in range(5, -1, -1):\
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


def compute_priv(u, s, nbr):
    count = 0
    for v in nbr[u]:
        if v not in s:
            if len(nbr[v] & s) == 1:
                count += 1
    return count


def main():
    print("HIGH-MODE TREES: MIN k AMONG LOW-PRIV IS")
    print("=" * 70)
    print()

    for n in range(5, 17):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        threshold = n // 3 + 2

        entries = []  # (tree_idx, mode, min_lp_k, max_lp_k, n_lp, n_mis)

        for idx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])

            if mode < threshold:
                continue

            nbr = [set(adj_data[v]) for v in range(tn)]
            all_mis = find_all_maximal_is(tn, adj_data)

            lp_sizes = []
            for s in all_mis:
                k = len(s)
                max_priv = max(compute_priv(u, s, nbr) for u in s)
                if max_priv <= 1:
                    lp_sizes.append(k)

            if lp_sizes:
                entries.append((idx, mode, min(lp_sizes), max(lp_sizes),
                                len(lp_sizes), len(all_mis), line.strip()))

        elapsed = time.time() - t0
        n_high = sum(1 for _ in [1 for l in lines
                      if (lambda p=independence_poly(*parse_graph6(l)):
                          max(range(len(p)), key=lambda i: p[i]) >= threshold)()])

        # Simpler count
        n_with_lp = len(entries)
        if entries:
            all_ok = all(e[2] >= e[1] for e in entries)
            min_gap = min(e[2] - e[1] for e in entries)
            print(f"n={n:2d}: {n_with_lp} high-mode trees with low-priv IS, "
                  f"min(min_k - mode) = {min_gap}, "
                  f"all min_k >= mode: {all_ok} ({elapsed:.1f}s)")
            if min_gap <= 2:
                for e in entries[:3]:
                    print(f"  tree {e[0]}: mode={e[1]}, lp_k=[{e[2]},{e[3]}], "
                          f"n_lp={e[4]}/{e[5]} MIS")
        else:
            print(f"n={n:2d}: no high-mode trees with low-priv IS ({elapsed:.1f}s)")

    print()
    print("=" * 70)
    print("If 'all min_k >= mode' is True for all n:")
    print("  Low-priv IS in high-mode trees are always at or above mode.")
    print("  Combined with low-mode analytic bound: PNP proved.")
    print("=" * 70)


if __name__ == "__main__":
    main()
