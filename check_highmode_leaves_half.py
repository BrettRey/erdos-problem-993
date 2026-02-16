#!/usr/bin/env python3
"""Check: do high-mode trees (mode > floor(n/3)+1) always have ℓ > n/2?

If yes, then mode ≤ (n-1)/2 < n/2 < ℓ, giving mode < ℓ directly.
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


def main():
    print("HIGH-MODE TREES: ℓ > n/2 CHECK", flush=True)
    print("=" * 60, flush=True)
    print(flush=True)

    for n in range(5, 23):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        threshold = n // 3 + 1
        half = n / 2

        n_high = 0
        n_viol = 0
        min_leaves = None
        worst_mode = None

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])

            if mode <= threshold:
                continue

            n_high += 1
            leaves = sum(1 for v in range(tn) if len(adj_data[v]) == 1)

            if min_leaves is None or leaves < min_leaves:
                min_leaves = leaves
                worst_mode = mode

            if leaves <= half:
                n_viol += 1
                ds = sorted([len(adj_data[v]) for v in range(tn)], reverse=True)
                print(f"  VIOLATION n={tn}: ℓ={leaves}, n/2={half}, mode={mode}, "
                      f"deg={ds[:6]}", flush=True)

        elapsed = time.time() - t0
        ml_str = f"min_ℓ={min_leaves}" if min_leaves is not None else "N/A"
        print(f"n={n:2d}: {n_high:6d} high-mode, viol={n_viol}, {ml_str}, "
              f"n/2={half}, threshold={threshold} ({elapsed:.1f}s)", flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
