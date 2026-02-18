#!/usr/bin/env python3
"""Do high-mode trees have ANY maximal IS with all priv <= 1?

If not, that's a clean structural statement: high-mode trees categorically
exclude low-priv maximal IS. Combined with the low-mode analytic bound,
this gives the full 1-Private Mode Conjecture.

Also reports: for each high-mode tree, what is the minimum max_priv
across all maximal IS? This shows how far high-mode trees are from
having a low-priv IS.
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


def compute_priv(u, s, nbr):
    count = 0
    for v in nbr[u]:
        if v not in s:
            if len(nbr[v] & s) == 1:
                count += 1
    return count


def main():
    print("HIGH-MODE TREES: DO ANY HAVE LOW-PRIV MAXIMAL IS?")
    print("=" * 70)
    print()
    print(f"{'n':>3} {'trees':>7} {'high':>6} {'hm_MIS':>8} {'hm_lp':>6} "
          f"{'min_mp':>7} {'min_mp_bm':>10} {'time':>6}")
    print("  (hm_lp = low-priv MIS in high-mode trees)")
    print("  (min_mp = min max_priv across ALL MIS in high-mode trees)")
    print("  (min_mp_bm = min max_priv among MIS BELOW mode)")
    print("-" * 60)

    for n in range(5, 17):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        threshold = n // 3 + 2

        n_trees = len(lines)
        n_high = 0
        hm_mis_total = 0
        hm_lowpriv = 0
        global_min_maxpriv = float('inf')
        global_min_maxpriv_bm = float('inf')

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])

            if mode < threshold:
                continue

            n_high += 1
            nbr = [set(adj_data[v]) for v in range(tn)]
            all_mis = find_all_maximal_is(tn, adj_data)
            hm_mis_total += len(all_mis)

            for s in all_mis:
                k = len(s)
                max_priv = max(compute_priv(u, s, nbr) for u in s)
                global_min_maxpriv = min(global_min_maxpriv, max_priv)
                if k < mode:
                    global_min_maxpriv_bm = min(global_min_maxpriv_bm, max_priv)
                if max_priv <= 1:
                    hm_lowpriv += 1

        elapsed = time.time() - t0
        mmp = global_min_maxpriv if global_min_maxpriv < float('inf') else '-'
        mmpb = global_min_maxpriv_bm if global_min_maxpriv_bm < float('inf') else '-'
        print(f"{n:3d} {n_trees:7d} {n_high:6d} {hm_mis_total:8d} {hm_lowpriv:6d} "
              f"{str(mmp):>7} {str(mmpb):>10} {elapsed:6.1f}s")

    print()
    print("=" * 70)
    print("If hm_lp = 0 for all n:")
    print("  High-mode trees categorically exclude low-priv maximal IS.")
    print("  Combined with low-mode analytic bound: full PNP proof.")
    print("=" * 70)


if __name__ == "__main__":
    main()
