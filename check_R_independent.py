#!/usr/bin/env python3
"""Check: in pendant trees (maximal IS with all priv <= 1), is R always independent?

If yes: the mode bound proof via the expansion property is complete.
If no: need to handle R-non-independent cases separately.
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
    print("R-INDEPENDENCE CHECK FOR PENDANT TREES")
    print("=" * 60)

    for n in range(5, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_pendant_is = 0
        n_R_indep = 0
        n_R_not_indep = 0
        max_eR = 0

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

                # Check all priv <= 1
                all_leq1 = True
                privates = set()
                for u in s:
                    priv_u = 0
                    for v in adj_data[u]:
                        if v not in s:
                            s_nbrs = len(nbr[v] & s)
                            if s_nbrs == 1:
                                priv_u += 1
                                privates.add(v)
                    if priv_u > 1:
                        all_leq1 = False
                        break

                if not all_leq1:
                    continue

                n_pendant_is += 1

                # Identify R
                R = set(range(tn)) - s - privates
                # Check if R is independent
                eR = 0
                for r in R:
                    for r2 in adj_data[r]:
                        if r2 in R and r2 > r:
                            eR += 1
                max_eR = max(max_eR, eR)
                if eR == 0:
                    n_R_indep += 1
                else:
                    n_R_not_indep += 1

        elapsed = time.time() - t0
        if n_pendant_is > 0:
            print(f"n={n:2d}: {n_pendant_is:5d} pendant IS, "
                  f"R_indep={n_R_indep}, R_not_indep={n_R_not_indep}, "
                  f"max_eR={max_eR} ({elapsed:.1f}s)")
        else:
            print(f"n={n:2d}: no pendant IS below mode ({elapsed:.1f}s)")


if __name__ == "__main__":
    main()
