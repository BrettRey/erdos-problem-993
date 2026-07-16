#!/usr/bin/env python3
"""Geometric scaling of the depth-3 phase-stacking champions.

Same protocol as scripts/valley_scaling_probe.py: scale copy counts
(A, B) proportionally at fixed gadget shapes (s, u, t), track the exact
valley margin V, the deficit law n*(1-V), dip position b/alpha, and the
rebound-window flag. Exact arithmetic (Kronecker-packed).
"""

import sys
import time

sys.path.insert(0, ".")

from scratch_depth3_valley_20260715 import tree_poly, size  # noqa: E402
from scripts.valley_search import valley_score  # noqa: E402

CHAMPS = [
    # (A, s, u, B, t), scale factors
    ((13, 2, 4, 89, 12), (1, 2, 4, 8)),
    ((34, 14, 2, 13, 2), (1, 2, 4, 8)),
    ((34, 9, 2, 34, 12), (1, 2, 4)),
    ((5, 9, 3, 89, 12), (1, 2, 4)),
]


def main():
    print(f"{'(A,s,u,B,t)':>24} {'n':>6} {'alpha':>6} {'V':>14} "
          f"{'1-V':>11} {'n(1-V)':>9} {'b/alpha':>8} {'win':>5} {'sec':>7}")
    for (A, s, u, B, t), ks in CHAMPS:
        for k in ks:
            Ak, Bk = A * k, B * k
            n = size(Ak, s, u, Bk, t)
            t0 = time.time()
            poly = tree_poly(Ak, s, u, Bk, t)
            vs = valley_score(poly)
            dt = time.time() - t0
            alpha = len(poly) - 1
            d = 1.0 - vs["ratio"]
            print(f"{str((Ak, s, u, Bk, t)):>24} {n:>6} {alpha:>6} "
                  f"{vs['ratio']:>14.10f} {d:>11.3e} {n*d:>9.4f} "
                  f"{vs['pos']/alpha:>8.4f} {str(vs['window']):>5s} "
                  f"{dt:>7.1f}", flush=True)
            if vs["witness"]:
                print(f"*** WITNESS *** {(Ak, s, u, Bk, t)}")
                print("poly=", poly, flush=True)
                return
        print(flush=True)


if __name__ == "__main__":
    main()
