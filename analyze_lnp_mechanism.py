#!/usr/bin/env python3
"""
Analyze the mechanism behind the Leaf-Neighbor Property.

Key insight from prove_leaf_neighbor.py: Route C' covers 100%.
But the "avoid all MLSVs" hypothesis fails — IS that avoid MLSVs
can be large. So the structural reason is subtler.

New approach: for each maximal IS S below mode, look at what happens
structurally. Specifically:

1. Every MLSV v has >= 2 leaf-children. If v not in S (because S is
   maximal), then all leaf-children of v must be in S (domination).
   This "absorbs" many vertices into S as leaves.

2. If S avoids ALL MLSVs, then for each MLSV v, S contains all
   leaf-children of v. This makes |S| large (many leaves forced in).

3. Hypothesis: the number of leaves forced into S by avoiding MLSVs
   pushes |S| above mode. That is:
   sum over MLSVs v of (leaf-children of v) > mode - alpha_avoiding

Key test: for each tree, compute:
  - L_forced(MLSV) = total leaf-children of all MLSVs
  - alpha_non_leaf = max IS among non-leaf vertices only
  - Is L_forced + alpha_non_leaf >= mode? (Would force S >= mode)

Actually, simpler: if S avoids all MLSVs, then S must contain all
leaf-children of every MLSV (by maximality/domination). How many
vertices is that?

Wait, that's wrong. S is independent. The leaf-children of an MLSV v
are all adjacent to v. If v is not in S, the leaf-children are NOT
forbidden by v (since v is not in S). But the leaf-children are only
forced into S by maximality if they have no OTHER neighbor in S.
A leaf has degree 1, so its only neighbor is v. If v not in S,
then the leaf has no S-neighbor, so by maximality the leaf MUST be in S.

So: if S avoids all MLSVs and S is maximal, then for each MLSV v,
ALL leaf-children of v are in S.

Now: these leaf-children are an independent set (they're all leaves
adjacent to v, and leaves of degree 1 can't be adjacent to each other
since that would give v two edges to the same pair = not a tree).
Wait, actually leaf-children of v are all adjacent to v, and leaves
have degree 1, so they can't be adjacent to each other. So they form
an independent set, and they're all in S.

Count: let L_v = number of leaf-children of v. For MLSV v, L_v >= 2.
If S avoids all MLSVs, S contains all leaves that are children of MLSVs.

These leaves may overlap if two MLSVs share a leaf-child — but that's
impossible since a leaf has degree 1, so it has exactly one neighbor.

Total forced leaves = sum_{v in MLSV} L_v.

Now, S also contains other vertices (non-leaves, non-MLSVs). But S
is independent, so non-MLSV vertices in S must not be adjacent to any
other S-vertex.

Key question: is sum_{v in MLSV} L_v >= mode?
If so, then |S| >= sum L_v >= mode, contradicting k < mode.

Let's check.
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


def main():
    max_n = 18  # push a bit further

    print("=" * 75)
    print("MECHANISM ANALYSIS: forced leaves vs mode")
    print("=" * 75)

    total_trees = 0
    forced_ge_mode = 0  # sum L_v >= mode (would prove the lemma)
    forced_lt_mode = 0
    counterexamples = []  # trees where forced < mode but LNP still holds

    # Also check: for trees with maximal IS below mode
    trees_with_below = 0
    forced_ge_mode_among_below = 0

    for n in range(5, max_n + 1):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_forced_ge = 0
        n_forced_lt = 0
        n_below = 0
        n_below_forced_ge = 0

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            deg = [len(adj_data[v]) for v in range(n)]
            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)
            total_trees += 1

            # MLSVs and their leaf-children
            mlsv_leaf_count = {}
            total_forced = 0
            for v in range(n):
                lc = sum(1 for w in adj_data[v] if deg[w] == 1)
                if lc >= 2:
                    mlsv_leaf_count[v] = lc
                    total_forced += lc

            if total_forced >= mode:
                forced_ge_mode += 1
                n_forced_ge += 1
            else:
                forced_lt_mode += 1
                n_forced_lt += 1

            # Does this tree have maximal IS below mode?
            nbr = [set(adj_data[v]) for v in range(n)]
            has_below = False

            def check_backtrack(v, current, forbidden):
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
                check_backtrack(v + 1, current, forbidden)
                if v not in forbidden and not has_below:
                    check_backtrack(v + 1, current + [v],
                                    forbidden | nbr[v])

            check_backtrack(0, [], set())

            if has_below:
                trees_with_below += 1
                n_below += 1
                if total_forced >= mode:
                    forced_ge_mode_among_below += 1
                    n_below_forced_ge += 1
                else:
                    if len(counterexamples) < 20:
                        counterexamples.append({
                            "n": n, "tree": tidx, "g6": line.strip(),
                            "mode": mode, "total_forced": total_forced,
                            "mlsv_leaf_count": mlsv_leaf_count,
                            "deg": deg,
                        })

        elapsed = time.time() - t0
        print(f"n={n:2d} ({len(lines):5d} trees, {elapsed:5.1f}s): "
              f"forced>=mode: {n_forced_ge:5d}  "
              f"forced<mode: {n_forced_lt:4d}  "
              f"has_below: {n_below:4d}  "
              f"below+forced>=mode: {n_below_forced_ge:4d}",
              flush=True)

    print(f"\n{'='*75}")
    print("SUMMARY")
    print(f"{'='*75}")
    print(f"Total trees: {total_trees}")
    print(f"sum L_v >= mode (all trees): {forced_ge_mode}/{total_trees} "
          f"({100*forced_ge_mode/total_trees:.1f}%)")
    print(f"Trees with maximal IS below mode: {trees_with_below}")
    print(f"sum L_v >= mode (among trees with below-mode IS): "
          f"{forced_ge_mode_among_below}/{trees_with_below} "
          f"({100*forced_ge_mode_among_below/trees_with_below:.1f}%)")

    gap = trees_with_below - forced_ge_mode_among_below
    print(f"\nGap (below-mode IS but forced < mode): {gap}")

    if counterexamples:
        print(f"\nCounterexamples (forced < mode but has below-mode IS):")
        for ce in counterexamples:
            print(f"  n={ce['n']} mode={ce['mode']} forced={ce['total_forced']} "
                  f"mlsv_lc={ce['mlsv_leaf_count']} g6={ce['g6']}")
    else:
        print("\n*** sum L_v >= mode for ALL trees with below-mode IS! ***")
        print("This means: if S avoids all MLSVs, S must contain all their")
        print("leaf-children, giving |S| >= sum L_v >= mode.")
        print("Contradiction with k < mode. QED.")


if __name__ == "__main__":
    main()
