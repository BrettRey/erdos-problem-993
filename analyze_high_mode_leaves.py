#!/usr/bin/env python3
"""Deeper analysis: for high-mode trees, how many leaf-neighbors does the
max-degree vertex have? This determines whether including it in any IS
automatically gives priv >= 2.

Key claim to prove: In high-mode trees, every maximal IS below mode
contains a vertex with >= 2 leaf-neighbors (hence priv >= 2).

Approach:
1. For each high-mode tree, find the max-degree vertex v.
2. Count leaf-neighbors of v: leaf_deg(v).
3. Check: if v is in any maximal IS S of size k < mode, priv(v) >= leaf_deg(v) >= 2.
4. If v is NOT in S: all leaf-neighbors must be in S -> |S| >= leaf_deg(v).
   Is leaf_deg(v) >= mode? If so, excluding v makes S too large (>= mode).

Also check the "support vertex" property: a support vertex has >= 1 leaf-neighbor.
If every vertex in every maximal IS below mode is either a support (priv >= 1 from leaf)
or has other sources of private neighbors, we can argue PNP.
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


def leaf_degree(v, n, adj):
    """Count leaf-neighbors of v (neighbors with degree 1)."""
    return sum(1 for w in adj[v] if len(adj[w]) == 1)


def main():
    print("HIGH-MODE TREE: LEAF STRUCTURE ANALYSIS")
    print("=" * 70)
    print()

    for n in range(8, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        threshold = n // 3 + 2  # high mode

        n_high = 0
        min_hub_leaf_deg = float('inf')
        max_hub_leaf_deg = 0

        # For each high-mode tree, check if excluding the hub forces IS above mode
        hub_exclusion_forces_above = 0
        hub_inclusion_gives_priv2 = 0
        both_hold = 0

        # Check actual maximal IS below mode
        n_mis_below = 0
        n_mis_with_hub = 0
        n_mis_without_hub = 0
        min_hub_priv = float('inf')

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])

            if mode < threshold:
                continue

            n_high += 1
            nbr = [set(adj_data[v]) for v in range(tn)]

            # Find max-degree vertex (the "hub")
            hub = max(range(tn), key=lambda v: len(adj_data[v]))
            hub_deg = len(adj_data[hub])
            ld = leaf_degree(hub, tn, adj_data)
            min_hub_leaf_deg = min(min_hub_leaf_deg, ld)
            max_hub_leaf_deg = max(max_hub_leaf_deg, ld)

            # Check: does excluding hub force IS above mode?
            excl_forces = (ld >= mode)  # all leaves must be in IS
            if excl_forces:
                hub_exclusion_forces_above += 1

            # Check: does including hub give priv >= 2?
            incl_priv2 = (ld >= 2)
            if incl_priv2:
                hub_inclusion_gives_priv2 += 1

            if excl_forces and incl_priv2:
                both_hold += 1

            # Enumerate maximal IS below mode (only for small n to save time)
            if n <= 16:
                all_mis = find_all_maximal_is(tn, adj_data)
                for s in all_mis:
                    k = len(s)
                    if k < mode:
                        n_mis_below += 1
                        if hub in s:
                            n_mis_with_hub += 1
                            # Compute priv(hub) in this IS
                            priv_hub = sum(1 for v in nbr[hub]
                                          if v not in s and len(nbr[v] & s) == 1)
                            min_hub_priv = min(min_hub_priv, priv_hub)
                        else:
                            n_mis_without_hub += 1

        elapsed = time.time() - t0
        if n_high == 0:
            print(f"n={n:2d}: {n_high:3d} high-mode trees ({elapsed:.1f}s)")
            continue

        print(f"n={n:2d}: {n_high:3d} high-mode trees, "
              f"hub_leaf_deg=[{min_hub_leaf_deg},{max_hub_leaf_deg}], "
              f"excl_forces_above={hub_exclusion_forces_above}/{n_high}, "
              f"incl_gives_priv2={hub_inclusion_gives_priv2}/{n_high}, "
              f"both={both_hold}/{n_high} "
              f"({elapsed:.1f}s)")
        if n <= 16 and n_mis_below > 0:
            print(f"  MIS below mode: {n_mis_below} total, "
                  f"with_hub={n_mis_with_hub}, without_hub={n_mis_without_hub}, "
                  f"min_hub_priv={min_hub_priv}")

    print()
    print("=" * 70)
    print("PROOF STRATEGY:")
    print("If 'both' = n_high for all n, then for every high-mode tree:")
    print("  - Including the hub in IS gives priv(hub) >= 2 (from leaf-neighbors)")
    print("  - Excluding the hub forces IS size >= mode (all leaves must be in IS)")
    print("Either way, no IS below mode can have all priv <= 1.")
    print("=" * 70)


if __name__ == "__main__":
    main()
