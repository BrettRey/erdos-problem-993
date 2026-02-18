#!/usr/bin/env python3
"""At n=17, if all priv <= 1 and k=6: S is a color class of the tree.

Check: for trees on 17 vertices with a color class of size 6,
what is the mode? If mode <= 6 always, then k=6 >= mode and there's
no PNP issue.

More generally: for any tree T on n vertices with bipartition (A, B)
where |A| = a, |B| = n-a: what's the maximum possible mode?
If mode <= min(a, n-a) always: then color class IS always have k >= mode.
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


def get_bipartition(n, adj):
    """Return the sizes of the two color classes."""
    color = [-1] * n
    color[0] = 0
    queue = [0]
    for v in queue:
        for w in adj[v]:
            if color[w] == -1:
                color[w] = 1 - color[v]
                queue.append(w)
    a = sum(1 for c in color if c == 0)
    return min(a, n - a), max(a, n - a)


def main():
    print("COLOR CLASS SIZE vs MODE")
    print("=" * 70)
    print()
    print("For each tree: compute mode and smaller color class size (min_class).")
    print("Check: is mode > min_class ever? If so, color class bound fails.")
    print()

    for n in [17, 18, 19, 20]:
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_trees = len(lines)
        n_exceed = 0
        min_gap = float('inf')  # min(min_class - mode)
        max_mode_at_exceed = 0

        # Count by min_class size
        class_mode_exceed = {}

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])
            min_class, max_class = get_bipartition(tn, adj_data)
            gap = min_class - mode
            min_gap = min(min_gap, gap)

            if mode > min_class:
                n_exceed += 1
                max_mode_at_exceed = max(max_mode_at_exceed, mode)
                key = min_class
                if key not in class_mode_exceed:
                    class_mode_exceed[key] = 0
                class_mode_exceed[key] += 1

        elapsed = time.time() - t0
        print(f"n={n:2d}: {n_trees:7d} trees, min(min_class - mode) = {min_gap}, "
              f"exceed = {n_exceed} ({elapsed:.1f}s)")
        if class_mode_exceed:
            for mc in sorted(class_mode_exceed):
                print(f"  min_class={mc}: {class_mode_exceed[mc]} trees with mode > min_class")

    print()
    print("=" * 70)
    print("If min(min_class - mode) >= 0 for all n:")
    print("  mode <= min_class for every tree.")
    print("  Since all-priv-<=1 forces S = color class, k >= min_class >= mode.")
    print("  This would prove the 1-Private Mode Conjecture!")
    print("=" * 70)


if __name__ == "__main__":
    main()
