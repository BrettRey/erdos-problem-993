#!/usr/bin/env python3
"""Check: for maximal IS below mode, does some vertex in S always have >= 2 leaf-children?

A "leaf-child" of u in S means: v is adjacent to u, deg(v) = 1, v not in S.
Since v has degree 1 and its only neighbor u is in S, v is automatically a private neighbor.

So "2 leaf-children" => priv >= 2.

More generally, check: for u in S with deg(u) = d, how many children of u are
leaves (in the rooted-at-u sense, i.e., neighbors of u that are leaves)?

This is stronger than the degree bound because leaf-neighbors are ALWAYS private
(they can't be shared - they have only one neighbor).
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


def find_all_maximal_is(n, adj):
    nbr = [set(adj[v]) for v in range(n)]
    result = []
    def backtrack(v, current, forbidden):
        if v == n:
            s = frozenset(current)
            for w in range(n):
                if w not in s and not (nbr[w] & s):
                    return
            result.append(s)
            return
        backtrack(v + 1, current, forbidden)
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])
    backtrack(0, [], set())
    return result


def main():
    print("LEAF-CHILDREN BOUND FOR PNP")
    print("=" * 70)
    print()
    print("For each maximal IS S below mode: max over u in S of")
    print("  leaf_children(u) = |{v : v adj u, deg(v)=1, v not in S}|")
    print("If this is always >= 2, leaf-children give priv >= 2.")
    print()
    print(f"{'n':>3} {'trees':>7} {'MIS_bm':>8} {'ph_ok':>7} {'deg_ok':>7} "
          f"{'leaf2':>7} {'ph|dg|lf':>8} {'uncov':>6} {'time':>6}")
    print("-" * 68)

    for n in range(5, 21):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_trees = len(lines)
        n_mis_bm = 0
        n_ph = 0
        n_deg = 0
        n_leaf2 = 0      # max leaf-children in S >= 2
        n_any = 0         # any of the three
        n_uncov = 0       # none of the three

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])
            nbr = [set(adj_data[v]) for v in range(tn)]

            # Identify leaves
            is_leaf = [len(adj_data[v]) == 1 for v in range(tn)]

            all_mis = find_all_maximal_is(tn, adj_data)
            for s in all_mis:
                k = len(s)
                if k >= mode:
                    continue

                n_mis_bm += 1

                ph = (tn >= 3 * k)
                max_deg = max(len(adj_data[u]) for u in s)
                deg = (max_deg >= k + 1)

                # Leaf-children: for each u in S, count leaf neighbors NOT in S
                max_lc = 0
                for u in s:
                    lc = sum(1 for v in adj_data[u] if is_leaf[v] and v not in s)
                    if lc > max_lc:
                        max_lc = lc
                leaf2 = (max_lc >= 2)

                if ph: n_ph += 1
                if deg: n_deg += 1
                if leaf2: n_leaf2 += 1
                if ph or deg or leaf2: n_any += 1
                if not ph and not deg and not leaf2:
                    n_uncov += 1

        elapsed = time.time() - t0
        print(f"{n:3d} {n_trees:7d} {n_mis_bm:8d} {n_ph:7d} {n_deg:7d} "
              f"{n_leaf2:7d} {n_any:8d} {n_uncov:6d} {elapsed:6.1f}s")

        if n >= 17 and n_uncov > 0:
            print(f"  *** {n_uncov} UNCOVERED at n={n}!")

    print()
    print("=" * 70)
    print("'leaf2' = max leaf-children in S >= 2 (leaf-children are always private)")
    print("If uncov = 0: disjunction covers all cases, PNP proved!")
    print("=" * 70)


if __name__ == "__main__":
    main()
