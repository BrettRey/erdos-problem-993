#!/usr/bin/env python3
"""
Deep analysis of the leaf-neighbor property.

Questions:
1. For the tightest cases (where only ONE vertex has ≥2 leaf-neighbors),
   what do the trees look like?
2. How does |L_out| (leaves outside S) compare to k?
   - If |L_out| >> k, pigeonhole suffices.
   - If |L_out| ≤ k, need structural concentration.
3. What's the relationship between tree leaf-count and mode?
   Trees with few leaves (paths) — do they ever have maximal IS below mode?

Key structural insight to test:
  For a tree T with ℓ leaves and a maximal IS S of size k < mode(I(T)):
  |L_out| = ℓ - |L ∩ S|.
  If the tree has enough leaves relative to k, pigeonhole works.
  Hypothesis: ℓ ≥ k + 2 always holds when k < mode.
"""

import subprocess
import sys
import time
from collections import defaultdict, Counter

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
    max_n = 16

    print("=" * 70)
    print("LEAF-NEIGHBOR PROPERTY: structural analysis")
    print("=" * 70)

    total = 0
    # Cases where exactly one vertex has ≥2 leaf-nbrs
    tight_count = 0
    tight_cases = []
    # Distribution of |L_out| - k
    lout_minus_k_hist = Counter()
    # Distribution of max_leaf_nbrs
    max_leaf_hist = Counter()
    # Does |L_out| >= k+1 always? (pigeonhole route)
    pigeonhole_works = 0
    pigeonhole_fails = 0
    # Leaf count vs k
    leaf_minus_k_hist = Counter()

    for n in range(5, max_n + 1):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            deg = [len(adj_data[v]) for v in range(n)]
            leaves = {v for v in range(n) if deg[v] == 1}
            n_leaves = len(leaves)

            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)

            max_is_list = find_maximal_is_below_mode(n, adj_data, mode)
            for s, k in max_is_list:
                total += 1

                L_in = s & leaves
                L_out = leaves - s
                n_lout = len(L_out)

                leaf_minus_k_hist[n_leaves - k] += 1
                lout_minus_k_hist[n_lout - k] += 1

                # Compute per-vertex leaf-neighbor counts
                nbr = [set(adj_data[v]) for v in range(n)]
                leaf_nbrs = {}
                for u in s:
                    leaf_nbrs[u] = sum(1 for v in nbr[u] if v in L_out)

                max_ln = max(leaf_nbrs.values()) if leaf_nbrs else 0
                max_leaf_hist[max_ln] += 1

                # How many vertices have ≥2?
                n_with_2plus = sum(1 for u in s if leaf_nbrs[u] >= 2)
                if n_with_2plus == 1:
                    tight_count += 1
                    if len(tight_cases) < 10:
                        tight_cases.append({
                            "n": n, "tree": tidx, "k": k, "mode": mode,
                            "n_leaves": n_leaves, "n_lout": n_lout,
                            "n_lin": len(L_in),
                            "leaf_nbrs": {u: leaf_nbrs[u] for u in sorted(s)},
                            "s": sorted(s),
                        })

                if n_lout >= k + 1:
                    pigeonhole_works += 1
                else:
                    pigeonhole_fails += 1

        elapsed = time.time() - t0
        print(f"n={n:2d} ({len(lines)} trees, {elapsed:.1f}s): "
              f"total={total} tight={tight_count} pigeonhole_fails={pigeonhole_fails}",
              flush=True)

    print(f"\n\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")

    print(f"\nTotal cases: {total}")
    print(f"|L_out| >= k+1 (pigeonhole works): {pigeonhole_works}/{total} "
          f"({100*pigeonhole_works/total:.1f}%)")
    print(f"|L_out| <= k (pigeonhole insufficient): {pigeonhole_fails}/{total} "
          f"({100*pigeonhole_fails/total:.1f}%)")

    print(f"\nMax leaf-neighbor distribution:")
    for m in sorted(max_leaf_hist.keys()):
        cnt = max_leaf_hist[m]
        print(f"  max_leaf_nbrs={m}: {cnt:6d} ({100*cnt/total:.1f}%)")

    print(f"\n|L_out| - k distribution (need ≥1 for pigeonhole):")
    for d in sorted(lout_minus_k_hist.keys()):
        cnt = lout_minus_k_hist[d]
        print(f"  |L_out|-k={d:3d}: {cnt:6d} ({100*cnt/total:.1f}%)")

    print(f"\nTight cases (exactly 1 vertex with ≥2 leaf-nbrs): "
          f"{tight_count}/{total} ({100*tight_count/total:.1f}%)")
    if tight_cases:
        for tc in tight_cases:
            print(f"  n={tc['n']} k={tc['k']} mode={tc['mode']} "
                  f"leaves={tc['n_leaves']} L_out={tc['n_lout']} L_in={tc['n_lin']}")
            print(f"    S={tc['s']}  leaf_nbrs={tc['leaf_nbrs']}")


if __name__ == "__main__":
    main()
