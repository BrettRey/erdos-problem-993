#!/usr/bin/env python3
"""
Test a direct relationship: alpha(T \ MLSV_neighborhood) < mode.

If we remove all MLSVs and their neighborhoods from the tree, the
remaining graph has a maximum IS. If that max IS < mode, then any
IS that avoids all MLSVs (restricted to T minus neighborhoods) can't
reach mode.

Wait, this isn't quite right. Let me think again.

The key structural fact: in a tree, the independence number alpha
and the mode of I(T) are related. For most trees, mode is much
less than alpha. The question is about SMALL IS (below mode) and
whether they must include MLSVs.

NEW APPROACH: Consider the contrapositive.

Claim: if T has an MLSV v, and S is a maximal IS with v not in S,
then |S| >= L_v + f(T,v) where f captures the cascade cost.

For the simplest case: v has r leaf-children l_1,...,l_r and
d = deg(v) - r non-leaf neighbors w_1,...,w_d.

If v not in S:
  - All l_i in S (maximality, since each has only neighbor v)
  - v is dominated by its leaf-children
  - Each w_j has v as a neighbor, but v not in S. So w_j needs
    another S-neighbor or w_j itself is in S.

For each path P leaving v through a non-leaf neighbor w:
  P goes through vertices of degree 2 until it hits a branch.
  Along this path, every other vertex must be in S (domination
  in a path forces every other vertex).

So: each non-leaf arm of v contributes roughly arm_length/2
additional vertices to S.

Let A_1, ..., A_d be the non-leaf arm lengths from v.
Then |S| >= r + sum(ceil(A_i/2)) roughly.

And if we include v: |S| can be as small as 1 + sum(floor(A_i/2)).
Difference: (r - 1) + d * (1/2) roughly.

This is getting complicated. Let me just verify a clean combinatorial
relationship.

CLEAN TEST: for each tree T with MLSV set M, compute:
  alpha_minus = max IS size using only V \ (M ∪ N(M))
  alpha_leaves = |leaf-children of all MLSVs|
  Their sum alpha_minus + alpha_leaves vs mode

If alpha_minus + alpha_leaves < mode for all trees, we're done:
any IS avoiding all MLSVs can have at most alpha_minus + alpha_leaves
vertices... wait, that's not right either. Vertices in N(M) \ M
can be in S if they avoid M.

Let me just verify the cleanest possible statement.

CLEAN TEST 2: for each tree, compute:
  max_IS_avoiding_all_MLSVs = max |S| where S maximal IS, S ∩ MLSVs = ∅
  Check: this >= mode? (If always >= mode, the LNP holds.)
  This is what analyze_lnp_deep already checked. Confirmed: always >= mode.

So the statement "max IS avoiding MLSVs >= mode" is EQUIVALENT to the LNP.
We need to prove this.

Let me try a LP-based bound. In a tree:
  - IS is a set with no two adjacent vertices
  - Maximal: every non-IS vertex has a neighbor in IS
  - Size: minimize |S| subject to IS + domination

For a tree, the min dominating IS = min maximal IS = i(T) (independent
domination number). This is well-studied!

KEY INSIGHT: i(T) is the minimum maximal IS size.
If i(T) >= mode, then every maximal IS has size >= mode, and there
are NO maximal IS below mode. The LNP is vacuously true.

If i(T) < mode, there exist maximal IS below mode. The LNP says
they all contain an MLSV.

What's the independent domination number of T restricted to V \ MLSVs?
If i(T[V \ MLSVs]) >= mode (forcing all leaf-children in too),
then avoiding MLSVs gives IS >= mode.

Hmm, but T[V \ MLSVs] is a forest, not a tree. And it includes the
leaf-children which are forced.

Let me try yet another approach: just measure i(T \ MLSVs) directly
where the IS is in V and must be maximal in T (not just in the subgraph).
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


def find_min_maximal_is_avoiding(n, adj, forbidden):
    """Find the MINIMUM size of a maximal IS in T that avoids `forbidden`."""
    nbr = [set(adj[v]) for v in range(n)]
    best = [n + 1]  # upper bound

    def backtrack(v, current, blocked):
        if v == n:
            s = frozenset(current)
            can_extend = any(u not in s and u not in forbidden
                             and not (nbr[u] & s) for u in range(n))
            # Wait: maximality is about the FULL graph. S must be dominating.
            # Every vertex not in S must have a neighbor in S.
            is_dominating = all(
                any(w in s for w in nbr[u]) for u in range(n) if u not in s
            )
            if not can_extend and is_dominating:
                # Actually, for maximal IS in the full graph:
                # no vertex can be added (independence), and every vertex
                # is either in S or has a neighbor in S (domination).
                # The second condition (domination) is automatic for maximal IS.
                # Wait: if we restrict to avoiding forbidden, can_extend only
                # checks non-forbidden vertices. But maximality in the full graph
                # means: S is independent AND maximal (no IS strictly contains S).
                # If u is not in S, u in forbidden, and u has no S-neighbor,
                # then S ∪ {u} might be independent (if u not in forbidden would
                # extend it). But since we're avoiding forbidden, we don't add u.
                # The resulting S might not be maximal in the full graph.
                #
                # Let me re-check: a maximal IS that avoids forbidden means:
                # S is independent, S ∩ forbidden = ∅, and for all v ∉ S:
                #   either v ∈ forbidden, or v has a neighbor in S.
                # This makes S maximal among IS that avoid forbidden.
                # But is S maximal in the FULL graph? Not necessarily:
                # a forbidden vertex with no S-neighbor could extend S.
                #
                # For the LNP, we need S to be maximal in the full graph.
                # So let me enforce: no vertex (including forbidden) can extend.
                pass

            # Maximal in full graph: no vertex outside S can be added
            can_extend_full = any(
                u not in s and not (nbr[u] & s) for u in range(n)
            )
            if not can_extend_full:
                best[0] = min(best[0], len(current))
            return

        # Skip v
        backtrack(v + 1, current, blocked)
        # Include v if not forbidden and not blocked
        if v not in forbidden and v not in blocked:
            backtrack(v + 1, current + [v], blocked | nbr[v])

    backtrack(0, [], set())
    return best[0] if best[0] <= n else None


def main():
    max_n = 16

    print("=" * 75)
    print("MIN MAXIMAL IS AVOIDING MLSVs vs MODE")
    print("=" * 75)
    print("If min_maximal_IS_avoiding_all_MLSVs >= mode always holds,")
    print("then the LNP follows immediately.")

    total_checked = 0
    always_ge_mode = 0
    fails = 0
    fail_examples = []

    for n in range(5, max_n + 1):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_checked = 0
        n_ge = 0
        n_fail = 0

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            deg = [len(adj_data[v]) for v in range(n)]
            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)

            mlsvs = set()
            for v in range(n):
                lc = sum(1 for w in adj_data[v] if deg[w] == 1)
                if lc >= 2:
                    mlsvs.add(v)

            if not mlsvs:
                continue

            total_checked += 1
            n_checked += 1

            min_avoid = find_min_maximal_is_avoiding(n, adj_data, mlsvs)

            if min_avoid is None:
                # No maximal IS avoids all MLSVs (every maximal IS includes one)
                always_ge_mode += 1
                n_ge += 1
            elif min_avoid >= mode:
                always_ge_mode += 1
                n_ge += 1
            else:
                fails += 1
                n_fail += 1
                if len(fail_examples) < 10:
                    fail_examples.append({
                        "n": n, "tree": tidx, "mode": mode,
                        "min_avoid": min_avoid,
                        "g6": line.strip(),
                    })

        elapsed = time.time() - t0
        print(f"n={n:2d} ({len(lines):5d} trees, {elapsed:5.1f}s): "
              f"checked={n_checked:5d}  >=mode: {n_ge:5d}  <mode: {n_fail:3d}",
              flush=True)

    print(f"\n{'='*75}")
    print("RESULT")
    print(f"{'='*75}")
    print(f"Trees with MLSVs: {total_checked}")
    print(f"min IS avoiding MLSVs >= mode: {always_ge_mode}/{total_checked}")
    print(f"min IS avoiding MLSVs < mode: {fails}/{total_checked}")

    if fails == 0:
        print("\n*** CONFIRMED: for every tree T with MLSVs,")
        print("    the minimum maximal IS avoiding all MLSVs has size >= mode(I(T)). ***")
        print()
        print("PROOF STRUCTURE:")
        print("  1. If T has no MLSVs, every internal vertex has <= 1 leaf-child.")
        print("     Show: such trees have no maximal IS below mode. (Separate lemma.)")
        print("  2. If T has MLSVs, every maximal IS avoiding all MLSVs has")
        print("     size >= mode. Therefore any maximal IS below mode must")
        print("     contain at least one MLSV. Since MLSV has >= 2 leaf-children,")
        print("     and S is independent, those children are outside S = leaf-neighbors.")
        print("     This gives the Leaf-Neighbor Property. QED (modulo step 2 proof).")
    else:
        print(f"\nFailed for {fails} trees:")
        for fe in fail_examples:
            print(f"  n={fe['n']} mode={fe['mode']} min_avoid={fe['min_avoid']} "
                  f"g6={fe['g6']}")


if __name__ == "__main__":
    main()
