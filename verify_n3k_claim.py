#!/usr/bin/env python3
"""
Detailed investigation of the key claim:
  For any tree T on n vertices and maximal IS S of size k < mode(I(T)),
  we have n >= 3k.

This is equivalent to pigeonhole giving priv(u) >= 2 for some u in S.

Track:
1. The ratio n/(3k) for all maximal IS below mode (how much slack?)
2. The tightest cases: what trees and IS achieve n close to 3k?
3. Whether the claim extends beyond n=16
"""

import subprocess
import sys
import time
from collections import defaultdict

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


def find_maximal_is_below_mode(n, adj, mode):
    """Find all maximal IS of size < mode. Returns list of (frozenset, size)."""
    nbr = [set(adj[v]) for v in range(n)]
    results = []

    def backtrack(v, current, forbidden):
        if v == n:
            k = len(current)
            if k >= mode:
                return
            s = frozenset(current)
            can_extend = any(u not in s and not (nbr[u] & s) for u in range(n))
            if not can_extend:
                results.append((s, k))
            return
        backtrack(v + 1, current, forbidden)
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])

    backtrack(0, [], set())
    return results


def main():
    max_n = 18  # Push as far as practical

    print("=" * 70)
    print("KEY CLAIM: maximal IS of size k < mode => n >= 3k")
    print("=" * 70)

    # Track tightest cases globally
    tightest_ratio = float("inf")  # n / (3k), want this > 1
    tightest_excess = float("inf")  # n - 3k, want this > 0
    tightest_example = None
    total_checked = 0
    violations = 0

    # Distribution of n - 3k values
    excess_hist = defaultdict(int)
    # Distribution of k values
    k_hist = defaultdict(int)

    for n in range(3, max_n + 1):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]
        n_trees = len(lines)

        n_checked = 0
        n_violations = 0
        n_tightest_excess = float("inf")
        n_tightest_example = None

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)

            # For n <= 16, enumerate all maximal IS
            if n <= 16:
                max_is_list = find_maximal_is_below_mode(n, adj_data, mode)
                for s, k in max_is_list:
                    n_checked += 1
                    total_checked += 1
                    excess = n - 3 * k
                    excess_hist[excess] += 1
                    k_hist[k] += 1

                    if excess < 0:
                        n_violations += 1
                        violations += 1

                    ratio = n / (3 * k) if k > 0 else float("inf")
                    if excess < n_tightest_excess:
                        n_tightest_excess = excess
                        n_tightest_example = {
                            "n": n, "tree": tidx, "g6": line.strip(),
                            "k": k, "mode": mode, "excess": excess,
                            "ratio": ratio, "s": sorted(s),
                            "poly": poly,
                        }
                    if excess < tightest_excess:
                        tightest_excess = excess
                        tightest_example = n_tightest_example
                    if ratio < tightest_ratio:
                        tightest_ratio = ratio
            else:
                # For n > 16, we can't enumerate all IS efficiently
                # But we can check mode vs n/3 as a weaker proxy
                # and check specific families
                pass

        elapsed = time.time() - t0
        if n <= 16:
            print(f"n={n:2d} ({n_trees:>7d} trees, {elapsed:5.1f}s): "
                  f"maxIS_below_mode={n_checked:>5d}  "
                  f"violations={n_violations}  "
                  f"tightest_excess={n_tightest_excess}", flush=True)
            if n_tightest_example and n_tightest_excess <= 4:
                ex = n_tightest_example
                print(f"       k={ex['k']} mode={ex['mode']} poly={ex['poly'][:8]}... "
                      f"IS={ex['s']}", flush=True)
        else:
            print(f"n={n:2d} ({n_trees:>7d} trees, {elapsed:5.1f}s): (skipped IS enum)", flush=True)

    print(f"\n\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")

    print(f"\nTotal maximal IS below mode checked: {total_checked}")
    print(f"Violations (n < 3k): {violations}")
    print(f"\nTightest ratio n/(3k): {tightest_ratio:.4f}")
    print(f"Tightest excess n - 3k: {tightest_excess}")
    if tightest_example:
        ex = tightest_example
        print(f"  At: n={ex['n']} tree={ex['tree']} k={ex['k']} mode={ex['mode']}")
        print(f"  IS: {ex['s']}")
        print(f"  Poly: {ex['poly']}")
        print(f"  g6: {ex['g6']}")

    print(f"\nExcess (n - 3k) distribution:")
    for exc in sorted(excess_hist.keys()):
        cnt = excess_hist[exc]
        pct = 100 * cnt / total_checked
        bar = "#" * min(50, int(pct))
        print(f"  n-3k={exc:3d}: {cnt:6d} ({pct:5.1f}%) {bar}")

    print(f"\nk distribution:")
    for k in sorted(k_hist.keys()):
        cnt = k_hist[k]
        pct = 100 * cnt / total_checked
        print(f"  k={k}: {cnt:6d} ({pct:5.1f}%)")

    # Key theoretical question: WHY does n >= 3k hold?
    print(f"\n{'='*70}")
    print("THEORETICAL ANALYSIS")
    print(f"{'='*70}")
    print("""
The claim "maximal IS of size k < mode(I(T)) => n >= 3k" is equivalent to:
  "if k >= ceil(n/3), then every maximal IS of size k has k >= mode"

This says: large maximal ISes (size >= n/3) are always at or above the mode.

Known bounds:
  - independence number alpha(T) >= ceil(n/2) for trees (bipartite)
  - mode(I(T)) <= alpha(T)
  - independent domination number i(T) = min size of maximal IS
  - For trees: gamma(T) <= n/3 (Ore's bound on domination number)
  - i(T) >= gamma(T)

The claim connects the mode to the independent domination number:
  if i(T) < mode(I(T)), then i(T) < n/3.
Equivalently: i(T) >= n/3 => i(T) >= mode(I(T)).
""")


if __name__ == "__main__":
    main()
