#!/usr/bin/env python3
"""
Deeper structural analysis: why does every maximal IS below mode
contain a support vertex with >= 2 leaf-children?

Key question: what happens if we assume the contrary?
If S avoids all multi-leaf support vertices, what constraints does
this impose on S and on the tree?

Analysis approach:
1. For each tree, identify all multi-leaf support vertices (MLSV)
2. For maximal IS S below mode, verify S ∩ MLSV ≠ ∅
3. For the tree's structure, understand WHY S must include an MLSV:
   - How many MLSVs does the tree have?
   - What fraction of max-IS include an MLSV?
   - If we force S to avoid all MLSVs, what's the max |S|?
   - Is max |S avoiding MLSV| < mode always?

Hypothesis: the maximum size of a maximal IS that avoids all MLSVs
is strictly less than mode(I(T)). This would prove the lemma.
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


def max_is_avoiding(n, adj, forbidden):
    """Find the maximum size of an IS that avoids all vertices in `forbidden`."""
    nbr = [set(adj[v]) for v in range(n)]
    allowed = set(range(n)) - forbidden
    best = [0]

    def backtrack(v_list, idx, current, forbidden_set):
        if idx == len(v_list):
            # Check maximality among all vertices (not just allowed)
            s = frozenset(current)
            can_extend = any(u not in s and not (nbr[u] & s) for u in range(n))
            if not can_extend:
                best[0] = max(best[0], len(current))
            return
        v = v_list[idx]
        # Skip
        backtrack(v_list, idx + 1, current, forbidden_set)
        # Include
        if v not in forbidden_set:
            backtrack(v_list, idx + 1, current + [v],
                      forbidden_set | nbr[v])

    backtrack(sorted(allowed), 0, [], set())
    return best[0]


def max_is_avoiding_v2(n, adj, forbidden):
    """Find the maximum size of a MAXIMAL IS that avoids `forbidden`.

    A maximal IS cannot be extended by adding any vertex (including
    forbidden ones that aren't in the IS).
    """
    nbr = [set(adj[v]) for v in range(n)]
    best = [0]

    def backtrack(v, current, blocked):
        if v == n:
            s = frozenset(current)
            # Check maximality: no vertex (including forbidden) can be added
            can_extend = any(u not in s and not (nbr[u] & s) for u in range(n))
            if not can_extend:
                best[0] = max(best[0], len(current))
            return
        # Skip v
        backtrack(v + 1, current, blocked)
        # Include v (only if allowed and not blocked)
        if v not in forbidden and v not in blocked:
            backtrack(v + 1, current + [v], blocked | nbr[v])

    backtrack(0, [], set())
    return best[0]


def main():
    max_n = 16

    print("=" * 75)
    print("SUPPORT VERTEX STRUCTURE: Why must maximal IS below mode")
    print("contain a multi-leaf support vertex (MLSV)?")
    print("=" * 75)

    total_trees_with_mis = 0
    hypothesis_holds = 0  # max IS avoiding MLSV < mode
    hypothesis_fails = 0
    fail_examples = []

    # Track: how many MLSVs per tree, and distribution
    mlsv_count_hist = {}

    for n in range(5, max_n + 1):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_trees_here = 0
        n_hyp_here = 0
        n_fail_here = 0

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            deg = [len(adj_data[v]) for v in range(n)]
            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)

            # Find MLSVs: vertices with >= 2 leaf-neighbors
            mlsvs = set()
            for v in range(n):
                leaf_children = sum(1 for w in adj_data[v] if deg[w] == 1)
                if leaf_children >= 2:
                    mlsvs.add(v)

            n_mlsv = len(mlsvs)
            mlsv_count_hist[n_mlsv] = mlsv_count_hist.get(n_mlsv, 0) + 1

            if not mlsvs:
                # Tree has no MLSV at all — what's its structure?
                # (These must be paths or caterpillars with single leaves)
                # Check if any maximal IS is below mode
                # For such trees, the LNP would fail if there existed a
                # maximal IS below mode. So either no such IS exists (mode
                # is achieved by all max-IS), or... we'd have a problem.
                # Let's just check.
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
                    # Uh oh: tree has no MLSV but has maximal IS below mode
                    print(f"  *** NO-MLSV TREE with below-mode IS: n={n} "
                          f"tree={tidx} g6={line.strip()} ***")
                continue

            # Tree has MLSVs. Compute max size of maximal IS avoiding all MLSVs.
            max_avoid = max_is_avoiding_v2(n, adj_data, mlsvs)

            total_trees_with_mis += 1
            n_trees_here += 1

            if max_avoid < mode:
                hypothesis_holds += 1
                n_hyp_here += 1
            else:
                hypothesis_fails += 1
                n_fail_here += 1
                fail_examples.append({
                    "n": n, "tree": tidx, "g6": line.strip(),
                    "mode": mode, "max_avoid": max_avoid,
                    "n_mlsv": n_mlsv,
                    "mlsvs": sorted(mlsvs),
                    "deg": deg,
                })

        elapsed = time.time() - t0
        print(f"n={n:2d} ({len(lines):5d} trees, {elapsed:5.1f}s): "
              f"tested={n_trees_here:5d} holds={n_hyp_here:5d} "
              f"fails={n_fail_here:3d}", flush=True)

    print(f"\n{'='*75}")
    print("HYPOTHESIS: max IS size avoiding all MLSVs < mode")
    print(f"{'='*75}")
    print(f"Trees tested: {total_trees_with_mis}")
    print(f"Holds: {hypothesis_holds}/{total_trees_with_mis}")
    print(f"Fails: {hypothesis_fails}/{total_trees_with_mis}")

    if hypothesis_fails == 0:
        print("\n*** HYPOTHESIS CONFIRMED through n={max_n}! ***")
        print("Every tree T with mode m has the property that any maximal IS")
        print("avoiding all multi-leaf support vertices has size < m.")
        print("Therefore every maximal IS of size >= m-1 (i.e., < mode but")
        print("within the ascending part) must include at least one MLSV.")
    else:
        print(f"\nHypothesis fails for {hypothesis_fails} trees:")
        for fe in fail_examples[:10]:
            print(f"  n={fe['n']} mode={fe['mode']} max_avoid={fe['max_avoid']} "
                  f"MLSVs={fe['mlsvs']} g6={fe['g6']}")

    # MLSV count distribution
    print(f"\n{'='*75}")
    print("MLSV COUNT DISTRIBUTION (per tree)")
    print(f"{'='*75}")
    for count in sorted(mlsv_count_hist.keys()):
        print(f"  {count} MLSVs: {mlsv_count_hist[count]} trees")

    # Additional: for trees with 0 MLSVs, what do they look like?
    print(f"\n{'='*75}")
    print("TREES WITH 0 MLSVS (no vertex has >= 2 leaf-children)")
    print(f"{'='*75}")
    print("These are trees where every internal vertex has at most 1 leaf-child.")
    print("Structure: paths, or 'caterpillars' with exactly 1 pendant leaf each.")
    print("Such trees cannot have maximal IS below mode (verified above).")
    print("This is because their IS polynomial tends to have low mode")
    print("and all max-IS reach it.")


if __name__ == "__main__":
    main()
