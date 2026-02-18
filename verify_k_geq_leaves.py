#!/usr/bin/env python3
"""Verify: k >= ℓ (leaves) for all 1-Private IS.

Structural argument: every leaf of T is either in S (leaf whose support ∉ S)
or has its support vertex in S with d_leaf = 1 (since priv ≤ 1 forces this).
So k = |S_leaves| + |S_non-leaves| >= ℓ (since every leaf is accounted for).

Actually the argument is more subtle. Let me verify computationally first.
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
    print("VERIFY: k >= ℓ FOR 1-PRIVATE IS", flush=True)
    print("=" * 60, flush=True)
    print(flush=True)

    total_1priv = 0
    total_violations = 0
    min_diff = None  # min (k - ℓ)

    for n in range(5, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_1priv = 0
        n_viol = 0
        local_min_diff = None

        for line in lines:
            tn, adj_data = parse_graph6(line)
            nbr = [set(adj_data[v]) for v in range(tn)]
            leaves = sum(1 for v in range(tn) if len(adj_data[v]) == 1)
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
                diff = k - leaves
                if local_min_diff is None or diff < local_min_diff:
                    local_min_diff = diff
                if diff < 0:
                    n_viol += 1
                    poly = independence_poly(tn, adj_data)
                    mode = max(range(len(poly)), key=lambda i: poly[i])
                    print(f"  VIOLATION: n={tn} k={k} ℓ={leaves} mode={mode}", flush=True)

        elapsed = time.time() - t0
        total_1priv += n_1priv
        total_violations += n_viol
        if local_min_diff is not None:
            if min_diff is None or local_min_diff < min_diff:
                min_diff = local_min_diff

        md_str = f"min(k-ℓ)={local_min_diff}" if local_min_diff is not None else "N/A"
        print(f"n={n:2d}: {n_1priv:6d} 1-Priv IS, viol={n_viol}, {md_str} ({elapsed:.1f}s)", flush=True)

    print(flush=True)
    print(f"TOTAL: {total_1priv} 1-Priv IS, violations: {total_violations}", flush=True)
    print(f"Min (k - ℓ): {min_diff}", flush=True)


if __name__ == "__main__":
    main()
