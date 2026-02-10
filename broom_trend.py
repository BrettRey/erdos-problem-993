#!/usr/bin/env python3
"""Study the near-miss ratio trend for broom trees as n grows.

Broom(p, s) = path of p vertices + star of s leaves at one end.
The targeted search found brooms are the closest to violating unimodality
(nm=0.9917 at n=500). Does nm → 1 as n → ∞?
"""

import sys
import time

from indpoly import independence_poly, is_unimodal, near_miss_ratio

from targeted import make_broom


def main():
    # Test broom(p, s) for fixed p and growing s
    # From targeted.py results, broom(33, 467) had the best nm ratio
    # Also try broom(13, s) and broom(42, s) which appeared in top 20

    for path_len in [13, 22, 33, 42, 50]:
        print(f"\n=== Broom({path_len}, s) trend ===", flush=True)
        print(f"{'s':>8} {'n':>8} {'alpha':>8} {'nm_ratio':>14} {'unimodal':>10} {'time':>8}", flush=True)

        prev_nm = 0
        for s in list(range(50, 1001, 50)) + list(range(1000, 5001, 500)):
            n = path_len + s
            t0 = time.time()
            n_tree, adj = make_broom(path_len, s)
            poly = independence_poly(n_tree, adj)
            uni = is_unimodal(poly)
            nm, nm_pos = near_miss_ratio(poly)
            elapsed = time.time() - t0

            if not uni:
                print(f"{'*'*60}", flush=True)
                print(f"*** COUNTEREXAMPLE: broom({path_len},{s}), n={n} ***", flush=True)
                print(f"*** poly = {poly} ***", flush=True)
                print(f"{'*'*60}", flush=True)
                sys.exit(0)

            trend = "↑" if nm > prev_nm else "↓" if nm < prev_nm else "="
            print(f"{s:>8} {n:>8} {len(poly)-1:>8} {nm:>14.10f} {str(uni):>10} {elapsed:>7.1f}s {trend}", flush=True)
            prev_nm = nm

            # If it's getting too slow, bail out
            if elapsed > 120:
                print(f"  (skipping larger s for path_len={path_len}, too slow)", flush=True)
                break


if __name__ == "__main__":
    main()
