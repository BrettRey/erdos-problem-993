#!/usr/bin/env python3
"""Find high-mode trees with smallest ℓ - mode gap.

Goal: understand what trees make the Leaf-Mode Inequality tightest,
to inform a proof strategy.
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


def degree_seq(n, adj):
    return sorted([len(adj[v]) for v in range(n)], reverse=True)


def main():
    print("TIGHT LEAF-MODE CASES IN HIGH-MODE TREES", flush=True)
    print("=" * 60, flush=True)
    print(flush=True)

    for n in range(5, 21):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        threshold = n // 3 + 1

        tight_cases = []  # (gap, n, ℓ, mode, deg_seq, g6)

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])

            if mode <= threshold:
                continue

            leaves = sum(1 for v in range(tn) if len(adj_data[v]) == 1)
            gap = leaves - mode

            if len(tight_cases) < 5 or gap < tight_cases[-1][0]:
                ds = degree_seq(tn, adj_data)
                tight_cases.append((gap, tn, leaves, mode, ds, line.strip()))
                tight_cases.sort(key=lambda x: x[0])
                tight_cases = tight_cases[:5]

        elapsed = time.time() - t0
        if tight_cases:
            print(f"n={n:2d}: tightest cases (ℓ-mode): ({elapsed:.1f}s)", flush=True)
            for g, tn, lv, m, ds, g6 in tight_cases:
                # Classify tree type
                max_deg = ds[0]
                n_leaves = lv
                n_internal = tn - lv
                print(f"  gap={g}: ℓ={lv}, mode={m}, deg_seq={ds[:6]}... "
                      f"max_deg={max_deg}, n_internal={n_internal}", flush=True)
        else:
            print(f"n={n:2d}: no high-mode trees ({elapsed:.1f}s)", flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
