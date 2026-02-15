#!/usr/bin/env python3
"""
Analyze gap trees: trees where sum L_v < mode but avoiding MLSVs
still gives IS >= mode.

Focus: understand the "extra cost" beyond forced leaves.

Key insight to test: when S avoids MLSV v and includes v's leaf-children,
those leaf-children dominate v. But they DON'T dominate v's non-leaf
neighbors. So v's non-leaf neighbors need OTHER S-vertices to dominate
them. This creates a cascade.

Specifically: if v is an MLSV with degree d and r leaf-children (r >= 2),
and v is not in S, then:
- v's r leaf-children are in S (forced by maximality)
- v is dominated (by its leaf-children)
- v's (d - r) non-leaf neighbors need domination from S

Each non-leaf neighbor w of v has deg(w) >= 2. If w is not in S,
w needs some S-neighbor other than v (since v not in S). So the non-leaf
neighborhood of v creates additional demands on S.

Let me trace this cascade for gap trees.
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


def compute_mode(poly):
    if not poly:
        return 0
    return max(range(len(poly)), key=lambda i: poly[i])


def find_all_maximal_is(n, adj):
    nbr = [set(adj[v]) for v in range(n)]
    results = []

    def backtrack(v, current, forbidden):
        if v == n:
            s = frozenset(current)
            can_extend = any(u not in s and not (nbr[u] & s) for u in range(n))
            if not can_extend:
                results.append(s)
            return
        backtrack(v + 1, current, forbidden)
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])

    backtrack(0, [], set())
    return results


def main():
    max_n = 16

    print("=" * 75)
    print("GAP TREE ANALYSIS: why avoiding MLSVs costs extra")
    print("=" * 75)

    for n in range(5, max_n + 1):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        gap_count = 0

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            deg = [len(adj_data[v]) for v in range(n)]
            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)

            mlsvs = {}
            for v in range(n):
                lc = sum(1 for w in adj_data[v] if deg[w] == 1)
                if lc >= 2:
                    mlsvs[v] = lc

            if not mlsvs:
                continue

            total_forced = sum(mlsvs.values())
            if total_forced >= mode:
                continue

            # Check for below-mode IS
            all_mis = find_all_maximal_is(n, adj_data)
            below = [s for s in all_mis if len(s) < mode]
            if not below:
                continue

            gap_count += 1

            # Analyze the smallest avoiding IS
            avoiding = [s for s in all_mis if not (s & set(mlsvs))]
            if not avoiding:
                continue

            min_avoid_is = min(avoiding, key=len)
            min_avoid_sz = len(min_avoid_is)

            # How much is the "extra cost" beyond forced leaves?
            leaves_in_avoid = set()
            for v in mlsvs:
                for w in adj_data[v]:
                    if deg[w] == 1:
                        leaves_in_avoid.add(w)

            forced_in = min_avoid_is & leaves_in_avoid
            extra_in = min_avoid_is - leaves_in_avoid
            extra_cost = len(extra_in)

            if gap_count <= 3:  # Print a few examples per n
                print(f"\nn={n} tree={tidx} mode={mode} forced_leaves={total_forced} "
                      f"min_avoid_sz={min_avoid_sz}")
                print(f"  MLSVs: {mlsvs}")
                print(f"  Forced leaves in IS: {len(forced_in)}/{total_forced}")
                print(f"  Extra non-leaf vertices: {extra_cost}")
                print(f"  Total: {min_avoid_sz} = {len(forced_in)} forced + "
                      f"{extra_cost} extra")
                print(f"  Deg sequence: {sorted(deg, reverse=True)}")

                # Trace: for each non-leaf neighbor of each MLSV, how is
                # it dominated?
                nbr = [set(adj_data[v]) for v in range(n)]
                for v, lc in mlsvs.items():
                    non_leaf_nbrs = [w for w in adj_data[v] if deg[w] >= 2]
                    print(f"  MLSV {v} (deg={deg[v]}, {lc} leaf-children, "
                          f"{len(non_leaf_nbrs)} non-leaf nbrs):")
                    for w in non_leaf_nbrs:
                        dom = nbr[w] & min_avoid_is
                        print(f"    non-leaf nbr {w} (deg={deg[w]}): "
                              f"dominated by {sorted(dom)}")

                # What's the smallest below-mode IS? (containing MLSV)
                below_with_mlsv = [s for s in below if s & set(mlsvs)]
                if below_with_mlsv:
                    smallest_below = min(below_with_mlsv, key=len)
                    print(f"  Smallest below-mode IS with MLSV: size={len(smallest_below)} "
                          f"S={sorted(smallest_below)}")

        elapsed = time.time() - t0
        print(f"n={n:2d}: {gap_count} gap trees ({elapsed:.1f}s)", flush=True)

    # KEY THEORETICAL OBSERVATION:
    print(f"\n{'='*75}")
    print("THEORETICAL ANALYSIS")
    print(f"{'='*75}")
    print("""
For a tree T with MLSV v (r >= 2 leaf-children), consider two cases:

INCLUDE v: v contributes 1 to |S|, excludes all deg(v) neighbors.
  Net contribution: 1 vertex, blocks deg(v) positions.

EXCLUDE v: v's r leaf-children are forced into S (maximality).
  These r vertices block only v (each leaf has degree 1).
  v's (deg(v) - r) non-leaf neighbors are NOT blocked by the leaves.
  They still need domination, requiring additional S-vertices.

  Net: r vertices (from leaves) + cascade from non-leaf neighbors.

The "swap cost" of excluding v instead of including v:
  Extra vertices: at least r - 1 (r leaves vs 1 hub)
  Fewer blocked: deg(v) - 1 blocked by v, vs ~1 blocked by leaves

For the proof, we don't need to track the full cascade.
We just need: |S_avoid| >= |S_include| + (r - 1) for each avoided MLSV.

If the minimum maximal IS containing at least one MLSV has size < mode,
and the minimum maximal IS avoiding all MLSVs has size >= mode, then
every maximal IS of size < mode must contain an MLSV. QED.

This is exactly what we've verified empirically.

The difficulty is proving |S_avoid| >= mode in general. The cascade
makes this hard to bound precisely.
""")

    # Alternative approach: SINGLE MLSV
    print(f"\n{'='*75}")
    print("SINGLE-MLSV ANALYSIS")
    print(f"{'='*75}")
    print("Perhaps we don't need to avoid ALL MLSVs.")
    print("Question: for each tree, does avoiding any SINGLE MLSV push |S| >= mode?")
    print("If so, the proof is: S below mode must contain every MLSV it can,")
    print("and since it can contain at least one...")

    for n in range(5, max_n + 1):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            deg = [len(adj_data[v]) for v in range(n)]
            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)

            mlsvs = {}
            for v in range(n):
                lc = sum(1 for w in adj_data[v] if deg[w] == 1)
                if lc >= 2:
                    mlsvs[v] = lc

            if not mlsvs:
                continue

            # Check: does avoiding a single MLSV force |S| >= mode?
            all_mis = find_all_maximal_is(n, adj_data)
            below = [s for s in all_mis if len(s) < mode]
            if not below:
                continue

            for v in mlsvs:
                # Find min maximal IS that avoids just v
                avoiding_v = [s for s in all_mis if v not in s]
                if avoiding_v:
                    min_av = min(len(s) for s in avoiding_v)
                    if min_av < mode:
                        # Can avoid this single MLSV and stay below mode
                        # But these IS must contain ANOTHER MLSV
                        av_below = [s for s in avoiding_v if len(s) < mode]
                        all_contain_other = all(
                            any(w in s for w in mlsvs if w != v)
                            for s in av_below
                        )
                        if not all_contain_other:
                            print(f"  *** PROBLEM: n={n} tree={tidx} "
                                  f"avoids MLSV {v}, below mode, "
                                  f"no other MLSV! ***")

    print("  (No problems reported = avoiding one MLSV below mode always")
    print("   includes another MLSV.)")


if __name__ == "__main__":
    main()
