#!/usr/bin/env python3
"""Check: does ℓ (number of leaves) >= mode for all high-mode trees?

High-mode means mode > floor(n/3) + 1.
If ℓ >= mode for high-mode trees, and k >= ℓ for 1-Private IS, then k >= mode.
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
    print("LEAVES vs MODE CHECK FOR HIGH-MODE TREES", flush=True)
    print("=" * 60, flush=True)
    print(flush=True)

    # Also check ALL trees, not just high-mode
    total_trees = 0
    total_highmode = 0
    violations_high = 0
    violations_all = 0
    min_gap_high = None  # min (ℓ - mode) for high-mode
    min_gap_all = None   # min (ℓ - mode) for all

    for n in range(5, 21):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        threshold = n // 3 + 1
        n_trees = 0
        n_high = 0
        n_viol_high = 0
        n_viol_all = 0
        local_min_gap_high = None
        local_min_gap_all = None
        worst_tree = None

        for line in lines:
            tn, adj_data = parse_graph6(line)
            n_trees += 1

            # Count leaves
            leaves = sum(1 for v in range(tn) if len(adj_data[v]) == 1)

            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])

            gap_all = leaves - mode
            if local_min_gap_all is None or gap_all < local_min_gap_all:
                local_min_gap_all = gap_all

            if gap_all < 0:
                n_viol_all += 1

            if mode > threshold:
                n_high += 1
                gap = leaves - mode
                if local_min_gap_high is None or gap < local_min_gap_high:
                    local_min_gap_high = gap
                    worst_tree = (tn, leaves, mode, line.strip())
                if gap < 0:
                    n_viol_high += 1

        elapsed = time.time() - t0
        total_trees += n_trees
        total_highmode += n_high
        violations_high += n_viol_high
        violations_all += n_viol_all
        if local_min_gap_high is not None:
            if min_gap_high is None or local_min_gap_high < min_gap_high:
                min_gap_high = local_min_gap_high
        if local_min_gap_all is not None:
            if min_gap_all is None or local_min_gap_all < min_gap_all:
                min_gap_all = local_min_gap_all

        mg_str = f"min_gap={local_min_gap_high}" if local_min_gap_high is not None else "N/A"
        print(f"n={n:2d}: {n_trees:7d} trees, {n_high:5d} high-mode, "
              f"viol_high={n_viol_high}, viol_all={n_viol_all}, "
              f"{mg_str} ({elapsed:.1f}s)", flush=True)

        if n_viol_high > 0 and worst_tree:
            t = worst_tree
            print(f"  WORST: n={t[0]} ℓ={t[1]} mode={t[2]}", flush=True)

    print(flush=True)
    print(f"SUMMARY:", flush=True)
    print(f"  Total trees: {total_trees}", flush=True)
    print(f"  High-mode trees: {total_highmode}", flush=True)
    print(f"  Violations (ℓ < mode) in high-mode: {violations_high}", flush=True)
    print(f"  Violations (ℓ < mode) in ALL trees: {violations_all}", flush=True)
    print(f"  Min gap (high-mode): {min_gap_high}", flush=True)
    print(f"  Min gap (all): {min_gap_all}", flush=True)


if __name__ == "__main__":
    main()
