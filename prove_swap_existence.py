#!/usr/bin/env python3
"""
Attempt to prove: every maximal IS of size k < mode(I(T)) in a tree has
some vertex u with priv(u) >= 2.

Strategy: find an alternative to pigeonhole that works even for n < 3k.

Key insight from violation analysis: when n < 3k, the IS always contains
a high-degree vertex whose leaf-neighbors are all private.

Hypothesis: In a tree, for any maximal IS S, the vertex u ∈ S with
maximum degree in the tree has at least 2 private neighbors.

Sub-hypothesis: For u ∈ S with deg(u) = d, priv(u) >= d - (k-1),
because at most k-1 of u's neighbors can be shared (each shared neighbor
needs a distinct other S-vertex at distance 2).

If d - (k-1) >= 2, i.e., d >= k+1, then we have the swap.

Question: does every maximal IS in a tree always have a vertex with
deg >= k+1?

Check: sum of degrees in S >= n-k (domination). Average deg >= (n-k)/k.
For avg >= k+1, need n >= k^2 + 2k. That's too strong for large k.

Alternative: does max deg in S >= k+1 when k < mode?

Let's check empirically.
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


def compute_priv_detail(n, adj, s):
    """For each u in S, compute priv(u) and shared(u)."""
    nbr = [set(adj[v]) for v in range(n)]
    deg = [len(adj[v]) for v in range(n)]
    k = len(s)

    result = {}
    for u in s:
        priv = 0
        shared = 0
        leaf_nbrs = 0
        for v in nbr[u]:
            if v not in s:
                s_nbrs = nbr[v] & s
                if len(s_nbrs) == 1:
                    priv += 1
                else:
                    shared += 1
                if deg[v] == 1:
                    leaf_nbrs += 1
        result[u] = {
            "deg": deg[u],
            "priv": priv,
            "shared": shared,
            "leaf_nbrs": leaf_nbrs,
        }

    return result


def main():
    max_n = 16

    print("=" * 70)
    print("SWAP EXISTENCE: structural analysis")
    print("=" * 70)

    # Test multiple hypotheses
    total = 0
    max_deg_ge_k1 = 0  # max deg(u) >= k+1
    has_leaf_nbr_2 = 0  # some u has >= 2 leaf-neighbors
    priv_bound_works = 0  # d - (k-1) >= 2 for some u (i.e., d >= k+1)
    max_priv_ge_2 = 0  # actual max priv >= 2

    # Track the cases where max_deg < k+1
    tight_cases = []

    for n in range(5, max_n + 1):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)
            deg = [len(adj_data[v]) for v in range(n)]

            max_is_list = find_maximal_is_below_mode(n, adj_data, mode)
            for s, k in max_is_list:
                total += 1
                detail = compute_priv_detail(n, adj_data, s)

                max_d = max(detail[u]["deg"] for u in s)
                max_p = max(detail[u]["priv"] for u in s)
                max_leaf = max(detail[u]["leaf_nbrs"] for u in s)

                if max_d >= k + 1:
                    max_deg_ge_k1 += 1
                if max_leaf >= 2:
                    has_leaf_nbr_2 += 1
                if max_d >= k + 1:  # same as priv_bound
                    priv_bound_works += 1
                if max_p >= 2:
                    max_priv_ge_2 += 1

                if max_d < k + 1:
                    tight_cases.append({
                        "n": n, "tree": tidx, "k": k, "mode": mode,
                        "max_deg_in_S": max_d,
                        "detail": {u: detail[u] for u in sorted(s)},
                        "s": sorted(s),
                        "g6": line.strip(),
                    })

        print(f"n={n:2d}: cumulative total={total}, max_deg>=k+1: {max_deg_ge_k1} "
              f"leaf_nbrs>=2: {has_leaf_nbr_2}  max_priv>=2: {max_priv_ge_2}",
              flush=True)

    print(f"\n\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"Total maximal IS below mode: {total}")
    print(f"Hypothesis: max deg(u) >= k+1:          {max_deg_ge_k1}/{total} ({100*max_deg_ge_k1/total:.2f}%)")
    print(f"Hypothesis: some u has >= 2 leaf-nbrs:   {has_leaf_nbr_2}/{total} ({100*has_leaf_nbr_2/total:.2f}%)")
    print(f"Actual: max priv(u) >= 2:                {max_priv_ge_2}/{total} ({100*max_priv_ge_2/total:.2f}%)")

    if tight_cases:
        print(f"\nCases where max deg in S < k+1 ({len(tight_cases)} cases):")
        for tc in tight_cases[:20]:
            print(f"  n={tc['n']} k={tc['k']} mode={tc['mode']} "
                  f"max_deg={tc['max_deg_in_S']} S={tc['s']}")
            for u in sorted(tc['detail']):
                d = tc['detail'][u]
                print(f"    u={u}: deg={d['deg']} priv={d['priv']} "
                      f"shared={d['shared']} leaf_nbrs={d['leaf_nbrs']}")
    else:
        print("\n*** max deg in S >= k+1 ALWAYS holds! ***")

    # Additional analysis: what's the relationship between deg(u) and shared(u)?
    # The bound priv(u) >= deg(u) - shared(u) is exact.
    # shared(u) <= k-1 (each shared neighbor needs a distinct other S-vertex)
    # But can shared(u) < k-1?
    print(f"\n{'='*70}")
    print("SHARED NEIGHBOR ANALYSIS")
    print(f"{'='*70}")

    # For the tightest case (highest-deg u with most shared), compute
    # the actual shared structure
    print("\nFor each u ∈ S with max priv: priv = deg - shared")
    print("shared(u) = # of u's neighbors with ≥2 S-neighbors")
    print("By tree structure: shared(u) ≤ # of other S-vertices at distance 2 from u")
    print("(because each shared neighbor of u must be adjacent to some w ∈ S, w ≠ u,")
    print(" and the path u-v-w has length 2 in the tree)")
    print()
    print("Implication: priv(u) ≥ deg(u) - |{w ∈ S : dist(u,w) = 2}|")
    print("            = deg(u) - (# S-vertices at distance 2 from u)")
    print()
    print("For the swap: need priv(u) ≥ 2, so need deg(u) > |{w : dist(u,w)=2, w∈S}| + 1")


if __name__ == "__main__":
    main()
