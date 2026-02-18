#!/usr/bin/env python3
"""
Verify two key properties needed for the Private Neighbor Bound proof:

1. mode(I(T)) > (n-1)/3 for all trees T on n vertices
   (This connects "k < mode" to "n >= 3k", enabling pigeonhole on private neighbors)

2. Tightness of P >= n - 2k + 1 bound
   (Check how tight the bound is across all maximal IS below mode)
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


def enumerate_maximal_is_below_mode(n, adj, mode):
    """Enumerate maximal independent sets of size < mode."""
    nbr = [set(adj[v]) for v in range(n)]
    results = []

    def backtrack(v, current, forbidden):
        if v == n:
            k = len(current)
            if k >= mode:
                return
            s = frozenset(current)
            # Check maximality
            can_extend = any(u not in s and not (nbr[u] & s) for u in range(n))
            if not can_extend:
                results.append(s)
            return
        backtrack(v + 1, current, forbidden)
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])

    backtrack(0, [], set())
    return results


def compute_private_neighbors(n, adj, s):
    """For a dominating IS S, compute total private neighbors P
    and per-vertex private neighbor counts."""
    nbr = [set(adj[v]) for v in range(n)]
    priv_count = {}  # u -> number of private neighbors
    total_P = 0

    for u in s:
        priv = 0
        for v in nbr[u]:
            if v not in s:
                # v is a private neighbor of u if u is v's only S-neighbor
                s_neighbors_of_v = nbr[v] & s
                if len(s_neighbors_of_v) == 1:
                    priv += 1
        priv_count[u] = priv
        total_P += priv

    return total_P, priv_count


def main():
    max_n = 20

    print("=" * 70)
    print("VERIFICATION: mode > (n-1)/3 and Private Neighbor Bound tightness")
    print("=" * 70)

    # Part 1: Check mode > (n-1)/3
    print("\n--- Part 1: mode(I(T)) vs (n-1)/3 ---\n")

    mode_violations = []
    min_mode_ratio = float("inf")
    min_mode_example = None

    # Part 2: Tightness of P >= n - 2k + 1
    print("--- Part 2: P >= n - 2k + 1 tightness ---\n")

    min_slack = float("inf")
    min_slack_example = None
    total_maximal_checked = 0
    # Track cases where some u in S has priv(u) >= 2
    swap_guaranteed = 0  # cases where pigeonhole guarantees >= 2 private nbrs for some u

    for n in range(3, max_n + 1):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]
        n_trees = len(lines)

        n_mode_violations = 0
        n_min_ratio = float("inf")
        n_min_example = None

        n_min_slack = float("inf")
        n_min_slack_example = None
        n_maximal = 0
        n_swap_ok = 0

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)

            # Part 1: check mode > (n-1)/3
            threshold = (n - 1) / 3
            ratio = mode / threshold if threshold > 0 else float("inf")
            if mode <= threshold:
                n_mode_violations += 1
                mode_violations.append({
                    "n": n, "tree": tidx, "g6": line.strip(),
                    "mode": mode, "threshold": threshold, "poly": poly,
                })
            if ratio < n_min_ratio:
                n_min_ratio = ratio
                n_min_example = {
                    "n": n, "tree": tidx, "g6": line.strip(),
                    "mode": mode, "threshold": threshold,
                }

            # Part 2: check P >= n - 2k + 1 tightness (only for small n, expensive)
            if n <= 16:
                maximal_sets = enumerate_maximal_is_below_mode(n, adj_data, mode)
                for s in maximal_sets:
                    k = len(s)
                    n_maximal += 1
                    total_maximal_checked += 1

                    total_P, priv_count = compute_private_neighbors(n, adj_data, s)
                    bound = n - 2 * k + 1
                    slack = total_P - bound

                    # Check if pigeonhole gives swap: some u has priv(u) >= 2
                    max_priv = max(priv_count.values()) if priv_count else 0
                    if max_priv >= 2:
                        n_swap_ok += 1
                        swap_guaranteed += 1

                    if slack < n_min_slack:
                        n_min_slack = slack
                        n_min_slack_example = {
                            "n": n, "tree": tidx, "g6": line.strip(),
                            "k": k, "P": total_P, "bound": bound,
                            "slack": slack, "priv_dist": dict(priv_count),
                            "max_priv": max_priv,
                        }
                    if slack < min_slack:
                        min_slack = slack
                        min_slack_example = n_min_slack_example

        if n_min_ratio < min_mode_ratio:
            min_mode_ratio = n_min_ratio
            min_mode_example = n_min_example

        elapsed = time.time() - t0
        print(f"n={n:2d} ({n_trees:>8d} trees, {elapsed:6.1f}s): "
              f"mode_violations={n_mode_violations}  "
              f"min_ratio={n_min_ratio:.3f}  ", end="", flush=True)
        if n <= 16:
            pct = 100 * n_swap_ok / max(n_maximal, 1)
            print(f"maximal_IS={n_maximal}  swap_ok={n_swap_ok} ({pct:.0f}%)  "
                  f"min_slack={n_min_slack}", flush=True)
        else:
            print("(P-bound check skipped for n>16)", flush=True)

    # Summary
    print(f"\n\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")

    print(f"\n--- mode > (n-1)/3 ---")
    if mode_violations:
        print(f"VIOLATIONS FOUND: {len(mode_violations)}")
        for v in mode_violations[:10]:
            print(f"  n={v['n']} tree={v['tree']} mode={v['mode']} "
                  f"threshold={v['threshold']:.2f} poly={v['poly']}")
    else:
        print(f"NO VIOLATIONS through n={max_n}")
    print(f"Tightest case: n={min_mode_example['n']} tree={min_mode_example['tree']} "
          f"mode={min_mode_example['mode']} threshold={min_mode_example['threshold']:.2f} "
          f"ratio={min_mode_ratio:.4f}")

    print(f"\n--- P >= n - 2k + 1 tightness (n<=16) ---")
    print(f"Total maximal IS checked: {total_maximal_checked}")
    print(f"Pigeonhole guarantees swap (some priv(u)>=2): {swap_guaranteed}/{total_maximal_checked} "
          f"({100*swap_guaranteed/max(total_maximal_checked,1):.1f}%)")
    print(f"Global min slack: {min_slack}")
    if min_slack_example:
        ex = min_slack_example
        print(f"  At: n={ex['n']} tree={ex['tree']} k={ex['k']} "
              f"P={ex['P']} bound={ex['bound']} slack={ex['slack']}")
        print(f"  Private neighbor distribution: {ex['priv_dist']}")
        print(f"  Max priv for any u: {ex['max_priv']}")

    # Key question: does the bound + pigeonhole always give a swap?
    if swap_guaranteed == total_maximal_checked:
        print(f"\n*** PIGEONHOLE ALWAYS WORKS: every maximal IS below mode "
              f"has some u with priv(u) >= 2 ***")
    else:
        n_fail = total_maximal_checked - swap_guaranteed
        print(f"\n*** PIGEONHOLE FAILS for {n_fail}/{total_maximal_checked} cases ***")
        print("  (swap still works via other mechanism — see augmented matching results)")


if __name__ == "__main__":
    main()
