#!/usr/bin/env python3
"""Check Leaf-Mode Inequality restricted to trees with all d_leaf <= 1.

d_leaf(v) = number of leaf-children of v.
"All d_leaf <= 1" means every internal vertex has at most 1 leaf neighbor.

If the Leaf-Mode Inequality holds for this subclass, then combined with
the "forced inclusion" bound for trees with some d_leaf >= 2, we get
a complete proof of PNP.
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
    print("LEAF-MODE INEQUALITY FOR d_leaf <= 1 TREES", flush=True)
    print("=" * 60, flush=True)
    print(flush=True)
    print("Checking: mode <= max(ℓ, floor(n/3)+1) restricted to trees",
          flush=True)
    print("where every internal vertex has at most 1 leaf-child.", flush=True)
    print(flush=True)

    for n in range(5, 23):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        threshold = n // 3 + 1
        n_dleaf1 = 0  # trees with all d_leaf <= 1
        n_highmode_dleaf1 = 0  # high-mode trees with all d_leaf <= 1
        n_viol = 0  # violations
        min_gap = None

        for line in lines:
            tn, adj_data = parse_graph6(line)

            # Check d_leaf condition
            all_dleaf_le1 = True
            leaves = set()
            for v in range(tn):
                if len(adj_data[v]) == 1:
                    leaves.add(v)

            for v in range(tn):
                if v in leaves:
                    continue
                # Count leaf-children of v
                d_leaf_v = sum(1 for u in adj_data[v] if u in leaves)
                if d_leaf_v >= 2:
                    all_dleaf_le1 = False
                    break

            if not all_dleaf_le1:
                continue

            n_dleaf1 += 1
            ell = len(leaves)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])

            if mode <= threshold:
                continue

            n_highmode_dleaf1 += 1
            gap = ell - mode

            if min_gap is None or gap < min_gap:
                min_gap = gap

            if gap < 0:
                n_viol += 1
                ds = sorted([len(adj_data[v]) for v in range(tn)],
                            reverse=True)
                print(f"  VIOLATION n={tn}: ℓ={ell}, mode={mode}, gap={gap}, "
                      f"deg={ds[:6]}", flush=True)

        elapsed = time.time() - t0
        mg_str = f"min_gap={min_gap}" if min_gap is not None else "N/A"
        print(f"n={n:2d}: {n_dleaf1:6d} d_leaf≤1 trees, "
              f"{n_highmode_dleaf1:4d} high-mode, viol={n_viol}, "
              f"{mg_str} ({elapsed:.1f}s)", flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
