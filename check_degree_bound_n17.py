#!/usr/bin/env python3
"""Check degree bound for PNP at n=17-18, high-mode trees only.

At n=17, hub inclusion fails (hub exclusion bound < mode for some trees).
Does the degree bound (max_deg_in_S >= k+1) still cover all IS below mode?
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


def compute_priv(u, s, nbr):
    count = 0
    for v in nbr[u]:
        if v not in s:
            if len(nbr[v] & s) == 1:
                count += 1
    return count


def main():
    print("DEGREE BOUND CHECK FOR n=17,18 (HIGH-MODE TREES ONLY)")
    print("=" * 70)
    print()

    for n in [17, 18]:
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        threshold = n // 3 + 2

        n_trees = len(lines)
        n_high = 0
        n_mis_bm = 0
        n_ph = 0
        n_deg = 0
        n_either = 0
        n_neither = 0

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])

            if mode < threshold:
                continue

            n_high += 1
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
                    # Report details
                    max_priv = max(compute_priv(u, s, nbr) for u in s)
                    print(f"  UNCOVERED: n={tn}, k={k}, mode={mode}, "
                          f"max_deg_in_S={max_deg}, max_priv={max_priv}")

        elapsed = time.time() - t0
        print(f"n={n}: {n_high} high-mode trees, {n_mis_bm} IS below mode, "
              f"ph={n_ph}, deg={n_deg}, either={n_either}, neither={n_neither} "
              f"({elapsed:.1f}s)")

    print()
    print("=" * 70)


if __name__ == "__main__":
    main()
