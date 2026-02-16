#!/usr/bin/env python3
"""Check: do ALL 1-Private IS in ALL trees have k > n/2?

If yes, combined with mode <= floor((n-1)/2), this gives k > mode always.
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
        n = ord(s[0]) - 63; idx = 1
    else:
        idx = 1; n = 0
        for _ in range(3):
            n = n * 64 + (ord(s[idx]) - 63); idx += 1
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
                adj[i].append(j); adj[j].append(i)
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
    print("1-PRIVATE IS vs n/2 CHECK", flush=True)
    print("=" * 60, flush=True)
    print(flush=True)

    total_1priv = 0
    total_above_half = 0
    total_at_half = 0
    total_below_half = 0
    min_ratio = None  # min k/(n/2)

    for n in range(5, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_1priv = 0
        n_above = 0
        n_at = 0
        n_below = 0
        local_min_k = None
        violations = []

        for line in lines:
            tn, adj_data = parse_graph6(line)
            nbr = [set(adj_data[v]) for v in range(tn)]
            all_mis = find_all_maximal_is(tn, adj_data)

            for s in all_mis:
                k = len(s)
                # Check all priv <= 1
                all_leq1 = True
                for u in s:
                    priv_u = 0
                    for v in adj_data[u]:
                        if v not in s and len(nbr[v] & s) == 1:
                            priv_u += 1
                    if priv_u > 1:
                        all_leq1 = False
                        break
                if not all_leq1:
                    continue

                n_1priv += 1
                half = tn / 2

                if k > half:
                    n_above += 1
                elif k == half:
                    n_at += 1
                else:
                    n_below += 1
                    poly = independence_poly(tn, adj_data)
                    mode = max(range(len(poly)), key=lambda i: poly[i])
                    violations.append((tn, k, half, mode, line.strip()))

                if local_min_k is None or k < local_min_k:
                    local_min_k = k

        elapsed = time.time() - t0
        total_1priv += n_1priv
        total_above_half += n_above
        total_at_half += n_at
        total_below_half += n_below

        if local_min_k is not None and min_ratio is None:
            min_ratio = local_min_k / (n / 2)
        elif local_min_k is not None:
            min_ratio = min(min_ratio, local_min_k / (n / 2))

        detail = f"min_k={local_min_k}" if local_min_k is not None else "none"
        print(f"n={n:2d}: {n_1priv:5d} 1-Priv IS, "
              f"above={n_above}, at={n_at}, below={n_below}, "
              f"{detail}, n/2={n/2:.1f} ({elapsed:.1f}s)", flush=True)

        for v in violations[:3]:
            print(f"  VIOLATION: n={v[0]} k={v[1]} n/2={v[2]} mode={v[3]}", flush=True)

    print(flush=True)
    print(f"TOTAL: {total_1priv} 1-Priv IS", flush=True)
    print(f"  above n/2: {total_above_half}", flush=True)
    print(f"  at n/2:    {total_at_half}", flush=True)
    print(f"  below n/2: {total_below_half}", flush=True)
    print(f"  min ratio k/(n/2): {min_ratio}", flush=True)


if __name__ == "__main__":
    main()
