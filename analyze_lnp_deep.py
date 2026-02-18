#!/usr/bin/env python3
"""
Deep analysis of the LNP mechanism for the 5002 "gap" trees.

These trees have:
  - sum L_v < mode (forced leaf count doesn't reach mode)
  - Maximal IS below mode that DO contain an MLSV

Questions:
1. For these gap trees, what do the maximal IS below mode look like?
2. What's the structure? (Multiple MLSVs? Leaves of different MLSVs
   in S, plus the MLSV itself?)
3. Can we find a combinatorial identity that explains why small maximal
   IS must include an MLSV even when forced leaves < mode?

Also: let's check a simpler structural property.

DOMINATION ARGUMENT:
If S is maximal IS (hence dominating), every vertex not in S has a
neighbor in S. Consider an MLSV v with leaf-children l_1, ..., l_r
(r >= 2). If v not in S, then each l_i (degree 1, sole neighbor v)
has no S-neighbor, so l_i must be in S.

Now the question is: if we put all leaf-children of all MLSVs into S,
and then extend to a maximal IS, what size do we get?

Alternative: count the "cost" of avoiding an MLSV.

For each MLSV v with r leaf-children:
  - If v in S: v is in S, contributes 1 to |S|, and excludes all
    neighbors of v from S.
  - If v not in S: all r leaf-children are in S, contribute r to |S|,
    and the leaf-children exclude only v from S (each has degree 1).

So including v "costs" 1 vertex but "blocks" deg(v) neighbors.
Excluding v "costs" r >= 2 vertices but only "blocks" v (via the leaves).

When r >= 2, excluding v is MORE EXPENSIVE (r vs 1) and LESS BLOCKING.
So excluding v tends to make S LARGER.

Hypothesis: for each MLSV v avoided, |S| increases by at least r - 1.
Total increase: sum (r_v - 1) for avoided MLSVs.

If sum (r_v - 1) >= mode - k_min where k_min is the minimum maximal IS
containing all MLSVs, then avoiding any MLSV pushes S above mode.

Let me check this more carefully.
"""

import subprocess
import sys
import time
from collections import Counter

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
    """Find ALL maximal independent sets."""
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
    max_n = 14  # Keep it fast for the deep analysis

    print("=" * 75)
    print("DEEP LNP MECHANISM: swap cost analysis")
    print("=" * 75)

    # For each tree with a gap case, analyze:
    # 1. The "swap cost" of avoiding each MLSV
    # 2. Whether this cost pushes |S| >= mode
    # 3. The actual maximal IS that avoid MLSVs: their sizes

    total_gap_trees = 0
    cost_explains = 0  # cost argument works
    cost_fails = 0

    # Also: for trees where forced < mode, find the min-size maximal IS
    # that avoids ALL MLSVs
    min_avoid_sizes = Counter()

    for n in range(5, max_n + 1):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            deg = [len(adj_data[v]) for v in range(n)]
            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)

            # Find MLSVs
            mlsvs = set()
            mlsv_lc = {}
            for v in range(n):
                lc = sum(1 for w in adj_data[v] if deg[w] == 1)
                if lc >= 2:
                    mlsvs.add(v)
                    mlsv_lc[v] = lc

            if not mlsvs:
                continue

            total_forced = sum(mlsv_lc.values())
            if total_forced >= mode:
                continue  # Not a gap tree

            # This is a gap tree: forced < mode
            # Find all maximal IS
            all_mis = find_all_maximal_is(n, adj_data)

            # Separate into: those below mode, those at/above mode
            below_mode = [s for s in all_mis if len(s) < mode]
            if not below_mode:
                continue  # No below-mode IS

            total_gap_trees += 1

            # Check: do all below-mode IS contain an MLSV?
            all_contain_mlsv = all(bool(s & mlsvs) for s in below_mode)
            if not all_contain_mlsv:
                print(f"  *** LNP FAILURE: n={n} tree={tidx} ***")
                continue

            # Find min size of maximal IS avoiding all MLSVs
            avoiding = [s for s in all_mis if not (s & mlsvs)]
            if avoiding:
                min_avoid = min(len(s) for s in avoiding)
                min_avoid_sizes[min_avoid - mode] += 1
            else:
                # No maximal IS avoids all MLSVs!
                min_avoid_sizes["no_avoid"] += 1

            # Cost analysis: for each MLSV v, the "swap cost" is
            # L_v - 1 (avoiding v adds L_v vertices instead of 1)
            swap_cost = sum(lc - 1 for lc in mlsv_lc.values())

            # If the smallest maximal IS including all MLSVs has
            # size s_min, then avoiding MLSVs gives IS of size
            # >= s_min + swap_cost. If s_min + swap_cost >= mode, done.

            # Find min size of maximal IS containing ALL MLSVs
            containing_all = [s for s in all_mis if mlsvs <= s]
            if containing_all:
                min_contain = min(len(s) for s in containing_all)
            else:
                # Can't contain all MLSVs (they might be adjacent)
                # Some MLSVs are adjacent, so we can't have all in IS
                min_contain = None

        elapsed = time.time() - t0
        print(f"n={n:2d} ({len(lines):5d} trees, {elapsed:5.1f}s)", flush=True)

    print(f"\n{'='*75}")
    print("GAP TREE ANALYSIS")
    print(f"{'='*75}")
    print(f"Gap trees (forced < mode, has below-mode IS): {total_gap_trees}")

    print(f"\nMin avoid IS size - mode distribution:")
    for k in sorted(min_avoid_sizes.keys(), key=lambda x: (isinstance(x, str), x)):
        print(f"  {k}: {min_avoid_sizes[k]}")

    # Now let's focus on the KEY observation:
    # For gap trees, the avoiding IS exist but they're all >= mode.
    # This is the structural fact we need to prove.

    print(f"\n{'='*75}")
    print("FOCUSED CHECK: avoiding IS always >= mode?")
    print(f"{'='*75}")

    # Re-check: for each gap tree, is min(|S| : S maximal IS, S ∩ MLSV = ∅) >= mode?
    gap_avoid_ge_mode = 0
    gap_avoid_lt_mode = 0
    gap_no_avoid = 0
    gap_examples = []

    for n in range(5, max_n + 1):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            deg = [len(adj_data[v]) for v in range(n)]
            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)

            mlsvs = set()
            mlsv_lc = {}
            for v in range(n):
                lc = sum(1 for w in adj_data[v] if deg[w] == 1)
                if lc >= 2:
                    mlsvs.add(v)
                    mlsv_lc[v] = lc

            if not mlsvs:
                continue

            total_forced = sum(mlsv_lc.values())
            if total_forced >= mode:
                continue

            # Check for below-mode IS
            nbr = [set(adj_data[v]) for v in range(n)]
            has_below = False

            def check_bt(v, current, forbidden):
                nonlocal has_below
                if has_below:
                    return
                if v == n:
                    k = len(current)
                    if k >= mode:
                        return
                    s = frozenset(current)
                    can_ext = any(u not in s and not (nbr[u] & s)
                                  for u in range(n))
                    if not can_ext:
                        has_below = True
                    return
                check_bt(v + 1, current, forbidden)
                if v not in forbidden and not has_below:
                    check_bt(v + 1, current + [v], forbidden | nbr[v])

            check_bt(0, [], set())
            if not has_below:
                continue

            # Find min avoiding IS
            all_mis = find_all_maximal_is(n, adj_data)
            avoiding = [s for s in all_mis if not (s & mlsvs)]

            if not avoiding:
                gap_no_avoid += 1
            else:
                min_avoid = min(len(s) for s in avoiding)
                if min_avoid >= mode:
                    gap_avoid_ge_mode += 1
                else:
                    gap_avoid_lt_mode += 1
                    if len(gap_examples) < 5:
                        gap_examples.append({
                            "n": n, "tree": tidx, "mode": mode,
                            "min_avoid": min_avoid,
                            "mlsvs": mlsv_lc,
                            "g6": line.strip(),
                            "avoiding_below": [sorted(s) for s in avoiding
                                                if len(s) < mode],
                        })

    print(f"No avoiding IS exists: {gap_no_avoid}")
    print(f"Min avoiding IS >= mode: {gap_avoid_ge_mode}")
    print(f"Min avoiding IS < mode: {gap_avoid_lt_mode}")

    if gap_avoid_lt_mode == 0:
        print("\n*** CONFIRMED: for gap trees, all MLSV-avoiding maximal IS")
        print("    have size >= mode. This is equivalent to the LNP. ***")
    else:
        print(f"\n*** GAP: {gap_avoid_lt_mode} trees where avoiding IS < mode! ***")
        for ge in gap_examples:
            print(f"  n={ge['n']} mode={ge['mode']} min_avoid={ge['min_avoid']}")
            print(f"    MLSVs: {ge['mlsvs']}")
            print(f"    Avoiding below mode: {ge['avoiding_below']}")


if __name__ == "__main__":
    main()
