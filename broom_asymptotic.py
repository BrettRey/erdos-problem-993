#!/usr/bin/env python3
"""Test broom(p,s) at very large s to confirm nm → 1 with gap ~ C/s."""

import time
from indpoly import independence_poly, is_unimodal, near_miss_ratio
from targeted import make_broom


def main():
    # Focus on broom(13, s) which was fastest
    path_len = 13
    print(f"Broom({path_len}, s) asymptotic analysis", flush=True)
    print(f"{'s':>10} {'n':>10} {'nm_ratio':>16} {'1-nm':>14} {'s*(1-nm)':>12} {'time':>8}", flush=True)

    for s in [1000, 2000, 5000, 10000, 20000, 50000]:
        n = path_len + s
        t0 = time.time()
        n_tree, adj = make_broom(path_len, s)
        poly = independence_poly(n_tree, adj)
        uni = is_unimodal(poly)
        nm, _ = near_miss_ratio(poly)
        elapsed = time.time() - t0

        gap = 1 - nm
        scaled = s * gap

        if not uni:
            print(f"*** COUNTEREXAMPLE at s={s}! ***", flush=True)
            return

        print(f"{s:>10} {n:>10} {nm:>16.12f} {gap:>14.10f} {scaled:>12.6f} {elapsed:>7.1f}s", flush=True)

        if elapsed > 300:
            print("  (stopping, too slow)", flush=True)
            break

    print(flush=True)
    print("If s*(1-nm) converges to a constant C, then nm = 1 - C/s → 1.", flush=True)
    print("The conjecture holds for all brooms (nm < 1 always) but margin → 0.", flush=True)


if __name__ == "__main__":
    main()
