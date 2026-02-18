#!/usr/bin/env python3
"""Check whether mode(I(T)) <= floor(n/2) for all trees.

If true, this combined with the color-class argument (m_out = 0 => k >= floor(n/2))
could close the high-mode gap in the 1-Private proof.

Also computes:
- mode vs floor(n/2) for all trees
- The maximum mode/floor(n/2) ratio
- Whether any tree has mode > floor(n/2)
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


def main():
    print("MODE vs FLOOR(n/2) CHECK FOR ALL TREES")
    print("=" * 70)
    print()
    print(f"{'n':>3} {'trees':>8} {'max_mode':>8} {'n//2':>5} {'exceed':>7} "
          f"{'mode=n//2':>9} {'high_mode':>10} {'hm_exceed':>10} {'time':>6}")
    print("-" * 70)

    for n in range(3, 21):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        half = n // 2
        threshold = n // 3 + 2  # high mode: mode > floor(n/3) + 1

        n_trees = len(lines)
        max_mode = 0
        n_exceed = 0       # trees with mode > floor(n/2)
        n_equal = 0         # trees with mode = floor(n/2)
        n_high = 0          # high-mode trees
        n_hm_exceed = 0     # high-mode trees with mode > floor(n/2)

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])
            max_mode = max(max_mode, mode)

            if mode > half:
                n_exceed += 1
            if mode == half:
                n_equal += 1
            if mode >= threshold:
                n_high += 1
                if mode > half:
                    n_hm_exceed += 1

        elapsed = time.time() - t0
        print(f"{n:3d} {n_trees:8d} {max_mode:8d} {half:5d} {n_exceed:7d} "
              f"{n_equal:9d} {n_high:10d} {n_hm_exceed:10d} {elapsed:6.1f}s")

    print()
    print("=" * 70)
    print("INTERPRETATION:")
    print("If 'exceed' = 0 for all n: mode <= floor(n/2) universally.")
    print("If 'hm_exceed' = 0: high-mode trees never have mode > floor(n/2).")
    print()
    print("PROOF STRATEGY:")
    print("If all priv <= 1 with k = ceil((n+1)/3):")
    print("  n ≡ 2 mod 3: m_out = 0, V\\S independent, S is color class, k >= floor(n/2)")
    print("  n ≡ 1 mod 3: m_out <= 1")
    print("  n ≡ 0 mod 3: m_out <= 2")
    print("Color class bound k >= floor(n/2) >= mode would close the high-mode case.")
    print("=" * 70)


if __name__ == "__main__":
    main()
