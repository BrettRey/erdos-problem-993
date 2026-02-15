#!/usr/bin/env python3
"""Analyze the ACTUAL shared neighbor count vs the worst-case bound.

The degree bound says: priv(u) >= deg(u) - shared(u), where shared(u) <= k-1.
But in trees, shared(u) = |{w in S \ {u} : dist(u,w) = 2}| (each shared neighbor
maps to a distinct S-vertex at distance 2).

What's the actual max_shared for the gap cases? If we can prove
shared(u) <= deg(u) - 2 when u is the high-degree vertex in the gap cases,
we'd get priv >= 2.

More generally: for the vertex u in S with max degree, what is shared(u)?
And can we relate it to k more tightly?
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
    print("SHARED NEIGHBOR ANALYSIS FOR GAP CASES")
    print("=" * 70)
    print()
    print("For the max-deg vertex u in S: shared(u) = |{v adj u : v not in S, |N(v) & S| >= 2}|")
    print("priv(u) = deg(u) - shared(u)")
    print("Degree bound says shared(u) <= k-1. What's the actual shared(u)?")
    print()

    for n in range(5, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        max_shared_ratio = 0.0  # max(shared / (k-1)) over all gap cases
        max_shared_val = 0
        n_gap = 0

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

                # Find max-deg vertex in S
                max_deg_v = max(s, key=lambda u: len(adj_data[u]))
                u = max_deg_v
                d = len(adj_data[u])

                # Compute shared(u)
                shared = 0
                priv = 0
                for v in adj_data[u]:
                    if v not in s:
                        if len(nbr[v] & s) >= 2:
                            shared += 1
                        elif len(nbr[v] & s) == 1:
                            priv += 1

                if d <= k:  # gap case: degree bound gives only priv >= 1
                    n_gap += 1
                    if k > 1:
                        ratio = shared / (k - 1)
                        max_shared_ratio = max(max_shared_ratio, ratio)
                    max_shared_val = max(max_shared_val, shared)

        elapsed = time.time() - t0
        if n_gap > 0:
            print(f"n={n:2d}: {n_gap:5d} gap cases, max_shared={max_shared_val}, "
                  f"max_shared/(k-1)={max_shared_ratio:.3f} ({elapsed:.1f}s)")
        else:
            print(f"n={n:2d}:     0 gap cases ({elapsed:.1f}s)")

    print()
    print("=" * 70)
    print("If max_shared/(k-1) < 1 always: the degree bound is conservative.")
    print("The key question: can we prove shared(u) <= deg(u) - 2 for the")
    print("max-deg vertex in S, when k < mode and deg(u) <= k?")
    print("Equivalently: can we prove dist2_S(u) <= deg(u) - 2?")
    print("This would give priv(u) >= 2.")
    print("=" * 70)


if __name__ == "__main__":
    main()
