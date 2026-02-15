#!/usr/bin/env python3
"""Check the degree bound for PNP across ALL trees.

Degree Bound Lemma: In a tree, for any u in S (maximal IS of size k):
  priv(u) >= deg(u) - (k-1)

So if max_deg_in_S >= k+1, some vertex has priv >= 2.

Question: For ALL maximal IS below mode (all trees, not just high-mode),
does max_deg_in_S >= k+1?

Also check three separate explanations:
1. Global pigeonhole: n >= 3k
2. Degree bound: max_deg_in_S >= k+1
3. Either of the above
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
    print("DEGREE BOUND FOR PNP: ALL TREES")
    print("=" * 70)
    print()
    print("Degree Bound: priv(u) >= deg(u) - (k-1). If max_deg >= k+1: priv >= 2.")
    print()
    print(f"{'n':>3} {'trees':>7} {'MIS_bm':>8} {'ph_ok':>7} {'deg_ok':>7} "
          f"{'either':>7} {'neither':>8} {'time':>6}")
    print("-" * 58)

    for n in range(5, 17):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_trees = len(lines)
        n_mis_bm = 0
        n_ph = 0      # global pigeonhole works (n >= 3k)
        n_deg = 0      # degree bound works (max_deg_in_S >= k+1)
        n_either = 0   # at least one works
        n_neither = 0  # neither works (PNP unexplained)

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])
            nbr = [set(adj_data[v]) for v in range(tn)]

            all_mis = find_all_maximal_is(tn, adj_data)
            for s in all_mis:
                k = len(s)
                if k >= mode:
                    continue

                n_mis_bm += 1

                ph = (tn >= 3 * k)
                max_deg = max(len(adj_data[u]) for u in s)
                deg = (max_deg >= k + 1)

                if ph:
                    n_ph += 1
                if deg:
                    n_deg += 1
                if ph or deg:
                    n_either += 1
                if not ph and not deg:
                    n_neither += 1

        elapsed = time.time() - t0
        print(f"{n:3d} {n_trees:7d} {n_mis_bm:8d} {n_ph:7d} {n_deg:7d} "
              f"{n_either:7d} {n_neither:8d} {elapsed:6.1f}s")

    print()
    print("=" * 70)
    print("If 'neither' = 0 for all n:")
    print("  The disjunction 'n >= 3k OR max_deg_in_S >= k+1' covers ALL cases.")
    print("  Both are analytically provable -> PNP proved!")
    print("=" * 70)


if __name__ == "__main__":
    main()
