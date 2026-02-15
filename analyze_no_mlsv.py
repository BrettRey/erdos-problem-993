#!/usr/bin/env python3
"""
Analyze trees with no MLSV (no vertex has >= 2 leaf-children).

Claim: such trees have no maximal IS below mode.

Structure: every internal vertex has at most 1 leaf-child.
This means the tree is either:
  - A path (all vertices have degree <= 2)
  - A "bare caterpillar" where each spine vertex has at most 1 pendant
  - More complex: every branching vertex has at most 1 leaf

Let's verify the claim and characterize the trees.
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


def find_min_maximal_is(n, adj):
    """Find the minimum size maximal IS (independent domination number)."""
    nbr = [set(adj[v]) for v in range(n)]
    best = [n + 1]

    def backtrack(v, current, forbidden):
        if v == n:
            s = frozenset(current)
            can_extend = any(u not in s and not (nbr[u] & s) for u in range(n))
            if not can_extend:
                best[0] = min(best[0], len(current))
            return
        if len(current) >= best[0]:  # prune
            return
        backtrack(v + 1, current, forbidden)
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])

    backtrack(0, [], set())
    return best[0]


def classify_no_mlsv(n, adj):
    """Classify a no-MLSV tree."""
    deg = [len(adj[v]) for v in range(n)]
    max_deg = max(deg)

    if max_deg <= 2:
        return "path"

    # Check: is it a caterpillar with exactly 1 pendant per spine vertex?
    leaves = {v for v in range(n) if deg[v] == 1}
    internal = [v for v in range(n) if deg[v] >= 2]
    internal_set = set(internal)

    # Spine: internal vertices form a path?
    internal_degs = {}
    for v in internal:
        internal_degs[v] = sum(1 for w in adj[v] if w in internal_set)

    if all(d <= 2 for d in internal_degs.values()):
        # Internal vertices form a path
        leaf_per_internal = {}
        for v in internal:
            leaf_per_internal[v] = sum(1 for w in adj[v] if w in leaves)
        # No MLSV means each has <= 1 leaf
        assert all(lc <= 1 for lc in leaf_per_internal.values())
        n_with_leaf = sum(1 for lc in leaf_per_internal.values() if lc == 1)
        return f"caterpillar({n_with_leaf}_pendants)"

    return "branched"


def main():
    max_n = 20  # Can go higher since these are simpler trees

    print("=" * 75)
    print("NO-MLSV TREES: do they have maximal IS below mode?")
    print("=" * 75)

    total_no_mlsv = 0
    has_below_mode = 0
    no_below_mode = 0
    i_eq_mode = 0  # min maximal IS = mode
    i_gt_mode = 0  # min maximal IS > mode

    type_counts = Counter()

    for n in range(3, max_n + 1):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_no_mlsv = 0
        n_below = 0

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            deg = [len(adj_data[v]) for v in range(n)]

            # Check no MLSV
            is_no_mlsv = True
            for v in range(n):
                lc = sum(1 for w in adj_data[v] if deg[w] == 1)
                if lc >= 2:
                    is_no_mlsv = False
                    break

            if not is_no_mlsv:
                continue

            total_no_mlsv += 1
            n_no_mlsv += 1

            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)
            i_T = find_min_maximal_is(n, adj_data)

            tree_type = classify_no_mlsv(n, adj_data)
            type_counts[tree_type] += 1

            if i_T < mode:
                has_below_mode += 1
                n_below += 1
                print(f"  *** BELOW-MODE in no-MLSV tree: "
                      f"n={n} tree={tidx} i(T)={i_T} mode={mode} "
                      f"type={tree_type} g6={line.strip()} ***")
            elif i_T == mode:
                i_eq_mode += 1
            else:
                i_gt_mode += 1

        elapsed = time.time() - t0
        if n_no_mlsv > 0:
            print(f"n={n:2d}: {n_no_mlsv:5d} no-MLSV trees, "
                  f"{n_below} with below-mode IS ({elapsed:.1f}s)", flush=True)

    print(f"\n{'='*75}")
    print("RESULT")
    print(f"{'='*75}")
    print(f"Total no-MLSV trees: {total_no_mlsv}")
    print(f"With maximal IS below mode: {has_below_mode}")
    print(f"i(T) = mode: {i_eq_mode}")
    print(f"i(T) > mode: {i_gt_mode}")

    if has_below_mode == 0:
        print("\n*** CONFIRMED: no-MLSV trees never have maximal IS below mode. ***")
        print("This means: in a no-MLSV tree, every maximal IS has size >= mode.")
        print("Equivalently: i(T) >= mode for all no-MLSV trees.")

    print(f"\nTree type distribution:")
    for tt, cnt in type_counts.most_common():
        print(f"  {tt}: {cnt}")

    # Stronger check: is i(T) >= mode for all no-MLSV trees?
    print(f"\ni(T) >= mode for all no-MLSV trees: "
          f"{has_below_mode == 0}")

    # Even stronger: is i(T) = mode always? (i.e., min maximal IS = mode)
    print(f"i(T) = mode for all no-MLSV trees: "
          f"{i_gt_mode == 0}")


if __name__ == "__main__":
    main()
