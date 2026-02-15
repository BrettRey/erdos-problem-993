#!/usr/bin/env python3
"""Check: if S is a maximal IS in a tree with all priv(u) <= 1, is |S| >= floor(n/2)?

This is the key claim needed to close the high-mode gap.
If true for ALL maximal IS (not just below mode), then:
  all priv <= 1 => k >= floor(n/2) >= mode(I(T))
which proves PNP (the Private Neighbor Property) analytically.

Also checks:
- all priv <= 1 => k >= ceil(n/2)?  (even stronger)
- minimum k among all-priv-<=1 IS
- m_out distribution for these IS
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
    """Count external private neighbors of u in S."""
    count = 0
    for v in nbr[u]:
        if v not in s:
            if len(nbr[v] & s) == 1:
                count += 1
    return count


def compute_m_out(n, adj, s):
    """Count edges within V\S."""
    count = 0
    for v in range(n):
        if v not in s:
            for w in adj[v]:
                if w > v and w not in s:
                    count += 1
    return count


def main():
    print("LOW-PRIV MAXIMAL IS: SIZE vs FLOOR(n/2)")
    print("=" * 70)
    print()
    print("For each tree, find ALL maximal IS with max_priv <= 1.")
    print("Check if |S| >= floor(n/2).")
    print()
    print(f"{'n':>3} {'trees':>7} {'MIS_tot':>8} {'low_priv':>9} {'min_k':>6} "
          f"{'n//2':>5} {'violate':>8} {'min_mout':>9} {'max_mout':>9} {'time':>6}")
    print("-" * 76)

    for n in range(4, 17):  # n<=16: MIS enumeration tractable
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        half = n // 2
        n_trees = len(lines)
        total_mis = 0
        n_low_priv = 0
        min_k_low_priv = float('inf')
        n_violate = 0  # low-priv IS with k < floor(n/2)
        min_mout = float('inf')
        max_mout = 0

        # For detailed reporting of violations
        violations = []

        for line in lines:
            tn, adj_data = parse_graph6(line)
            nbr = [set(adj_data[v]) for v in range(tn)]
            all_mis = find_all_maximal_is(tn, adj_data)
            total_mis += len(all_mis)

            for s in all_mis:
                k = len(s)
                max_priv = max(compute_priv(u, s, nbr) for u in s)
                if max_priv <= 1:
                    n_low_priv += 1
                    min_k_low_priv = min(min_k_low_priv, k)
                    mo = compute_m_out(tn, adj_data, s)
                    min_mout = min(min_mout, mo)
                    max_mout = max(max_mout, mo)
                    if k < half:
                        n_violate += 1
                        if len(violations) < 5:
                            poly = independence_poly(tn, adj_data)
                            mode = max(range(len(poly)), key=lambda i: poly[i])
                            violations.append((n, k, half, mode, mo, line.strip()))

        elapsed = time.time() - t0
        mk = min_k_low_priv if min_k_low_priv < float('inf') else '-'
        mmo = min_mout if min_mout < float('inf') else '-'
        xmo = max_mout if max_mout > 0 else '-'
        print(f"{n:3d} {n_trees:7d} {total_mis:8d} {n_low_priv:9d} {str(mk):>6} "
              f"{half:5d} {n_violate:8d} {str(mmo):>9} {str(xmo):>9} {elapsed:6.1f}s")

        for v in violations:
            print(f"  VIOLATION: n={v[0]}, k={v[1]}, n//2={v[2]}, mode={v[3]}, m_out={v[4]}, g6={v[5]}")

    print()
    print("=" * 70)
    print("KEY CLAIMS:")
    print("1. 'low_priv' = count of maximal IS where ALL vertices have priv <= 1")
    print("2. 'violate' = count where such IS has |S| < floor(n/2)")
    print("3. If violate = 0 for all n: 'all priv <= 1 => k >= floor(n/2)' holds")
    print("4. Combined with 'mode <= floor(n/2)': gives k >= mode (PNP proved!)")
    print("=" * 70)


if __name__ == "__main__":
    main()
