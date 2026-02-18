#!/usr/bin/env python3
"""
Verify the Private Neighbor Property (PNP) for n up to 18:
For any tree T and any maximal IS S of size k < mode(I(T)),
some vertex u in S has priv(u) >= 2.

Also check: does the pigeonhole argument (P >= n-2k+1, some u has
priv >= ceil(P/k) >= 2 when n >= 3k) cover all cases? What about
the n < 3k cases?

For n < 3k cases: check if max_priv >= 2 via degree concentration.
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
    max_n = 18

    print("=" * 75)
    print("PRIVATE NEIGHBOR PROPERTY (PNP) VERIFICATION")
    print("=" * 75)

    total = 0
    pnp_holds = 0
    pnp_fails = 0
    pigeonhole_covers = 0  # n >= 3k
    n_lt_3k = 0
    n_lt_3k_max_priv_ge2 = 0

    for n in range(5, max_n + 1):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_total = 0
        n_pnp = 0
        n_fail = 0
        n_pigeonhole = 0
        n_sub3k = 0
        n_sub3k_ok = 0

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            nbr = [set(adj_data[v]) for v in range(tn)]
            poly = independence_poly(tn, adj_data)
            mode = compute_mode(poly)

            mis_below = find_maximal_is_below_mode(tn, adj_data, mode)

            for s, k in mis_below:
                total += 1
                n_total += 1

                # Compute max_priv
                max_p = 0
                for u in s:
                    p = sum(1 for v in nbr[u]
                            if v not in s and len(nbr[v] & s) == 1)
                    max_p = max(max_p, p)

                if max_p >= 2:
                    pnp_holds += 1
                    n_pnp += 1
                else:
                    pnp_fails += 1
                    n_fail += 1
                    print(f"  *** PNP FAILURE: n={tn} tree={tidx} "
                          f"k={k} mode={mode} max_priv={max_p} ***")

                # Pigeonhole check
                if tn >= 3 * k:
                    pigeonhole_covers += 1
                    n_pigeonhole += 1
                else:
                    n_lt_3k += 1
                    n_sub3k += 1
                    if max_p >= 2:
                        n_lt_3k_max_priv_ge2 += 1
                        n_sub3k_ok += 1

        elapsed = time.time() - t0
        print(f"n={n:2d} ({len(lines):6d} trees, {elapsed:6.1f}s): "
              f"cases={n_total:6d}  "
              f"PNP_ok={n_pnp:6d}  "
              f"PNP_fail={n_fail:3d}  "
              f"pigeonhole={n_pigeonhole:6d}  "
              f"n<3k={n_sub3k:4d}(ok={n_sub3k_ok})",
              flush=True)

    print(f"\n{'='*75}")
    print("SUMMARY")
    print(f"{'='*75}")
    print(f"Total maximal IS below mode: {total}")
    print(f"PNP holds (max_priv >= 2): {pnp_holds}/{total} "
          f"({100*pnp_holds/total:.4f}%)")
    print(f"PNP fails: {pnp_fails}/{total}")
    print(f"\nPigeonhole covers (n >= 3k): {pigeonhole_covers}/{total} "
          f"({100*pigeonhole_covers/total:.2f}%)")
    print(f"n < 3k cases: {n_lt_3k} (max_priv >= 2: {n_lt_3k_max_priv_ge2})")

    if pnp_fails == 0:
        print("\n*** PNP CONFIRMED through n=18! ***")
        print("Every maximal IS of size k < mode has some u with priv(u) >= 2.")
        print("This directly gives a valid swap (remove u, add 2 private nbrs)")
        print("since in a tree, private neighbors of the same vertex are")
        print("never adjacent (acyclicity).")


if __name__ == "__main__":
    main()
