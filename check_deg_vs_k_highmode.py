#!/usr/bin/env python3
"""For high-mode trees: does every maximal IS below mode contain a vertex with deg > k?

If YES: then priv(u) >= deg(u) - (k-1) >= 2 (degree bound lemma), proving PNP.
This would be the analytic closure of Part 3.

Equivalently: for IS below mode in high-mode trees, is max_deg_in_S > k always?
We already know this fails at n=17 (43 gap cases with max_deg = k).
But those gap cases have k=6, mode=7, and mode > ceil((n+1)/3) = 6.
So mode = 7 > 6 = ceil(18/3). These ARE high-mode trees.

So the answer is NO: high-mode trees CAN have IS below mode with max_deg = k.
But PNP still holds in all 43 cases (priv 4-6).

Let me instead check: for ALL maximal IS below mode (not just high-mode),
what's the relationship between max_deg_in_S and k?
And more importantly: when max_deg_in_S = k (the gap), what's the actual
max_priv? Is it always >= 2?

The key insight from the 43 gap cases: even when deg(u) = k (giving priv >= 1
from degree bound), the ACTUAL shared count is much less than k-1, because
the distance-2 S-vertices are sparse. In fact, shared(u) = |{w in S\{u}: dist(u,w)=2}|.
In a tree, the number of S-vertices at distance exactly 2 from u is at most
the number of non-leaf neighbors of u (each non-leaf nbr can have at most one
S-vertex behind it).

So: shared(u) <= deg(u) - leaf_nbrs_outside_S(u)
And: priv(u) >= deg(u) - shared(u) >= leaf_nbrs_outside_S(u)

If u has >= 2 leaf-neighbors outside S: priv >= 2 directly!

But we know leaf support fails (spider at n=18). So we need something else.

Let's compute: for each IS below mode, what's the actual max over u in S of
(deg(u) - |{w in S\{u}: dist(u,w)=2}|)?
This is the "distance-2 private neighbor count".
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
    print("DISTANCE-2 REFINEMENT OF DEGREE BOUND")
    print("=" * 70)
    print()
    print("For u in S: shared(u) = |{w in S\\{u}: dist(u,w)=2}|")
    print("priv(u) = deg(u) - shared(u)")
    print("If max(priv) >= 2 for all IS below mode: PNP proved!")
    print()
    print("Key question: in the gap cases (deg_in_S = k, n < 3k),")
    print("is dist2_S always < deg - 1?")
    print()

    for n in range(5, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_mis_bm = 0
        max_shared_over_all = 0
        n_gap = 0  # cases where degree bound gives only priv >= 1
        n_dist2_closes = 0  # gap cases where dist-2 refinement gives priv >= 2
        n_dist2_fails = 0  # gap cases where even dist-2 refinement gives only priv >= 1
        worst_min_priv = float('inf')

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

                # For each u in S: compute exact shared count
                # shared(u) = |{v in N(u) \ S : |N(v) & S| >= 2}|
                # priv(u) = |{v in N(u) \ S : |N(v) & S| == 1}|
                best_priv = 0
                for u in s:
                    priv_u = 0
                    shared_u = 0
                    for v in adj_data[u]:
                        if v not in s:
                            s_nbrs = len(nbr[v] & s)
                            if s_nbrs == 1:
                                priv_u += 1
                            else:
                                shared_u += 1
                    best_priv = max(best_priv, priv_u)
                    max_shared_over_all = max(max_shared_over_all, shared_u)

                worst_min_priv = min(worst_min_priv, best_priv)

                # Is this a gap case?
                max_deg = max(len(adj_data[u]) for u in s)
                ph = (tn >= 3 * k)
                deg_ok = (max_deg >= k + 1)
                if not ph and not deg_ok:
                    n_gap += 1
                    if best_priv >= 2:
                        n_dist2_closes += 1
                    else:
                        n_dist2_fails += 1

        elapsed = time.time() - t0
        gap_str = ""
        if n_gap > 0:
            gap_str = f" gap={n_gap}, closed={n_dist2_closes}, open={n_dist2_fails}"
        print(f"n={n:2d}: {n_mis_bm:6d} IS_bm, min(max_priv)={worst_min_priv},{gap_str} ({elapsed:.1f}s)")

    print()
    print("=" * 70)
    print("If min(max_priv) >= 2 for all n: PNP holds everywhere.")
    print("The 'gap' cases are where neither pigeonhole nor degree bound suffices.")
    print("If open=0: the actual priv is always >= 2 even in gap cases.")
    print("=" * 70)


if __name__ == "__main__":
    main()
