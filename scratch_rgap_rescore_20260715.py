#!/usr/bin/env python3
"""Corrected re-scoring: descent-thresholded rebound metric R_gap.

R_gap(T; theta) = max over b with prefixmax(b) >= (1+theta)*i_b of
suffixmax(b)/i_b, with the rise distance c-b of the achieving pair.
R_gap > 1 at any theta > 0 implies a non-unimodality witness. Unlike the
raw valley margin V, R_gap is immune to the flat-mode lattice artifact
(V rewards trees whose mode lands near an r=1 lattice crossing).

Re-scores (a) the depth-3 grid, (b) all 2026-07-15 campaign champions.
"""

import sys
import time
from fractions import Fraction

sys.path.insert(0, ".")

from scratch_depth3_valley_20260715 import tree_poly, size  # noqa: E402
from scripts.valley_scaling_probe import bouquet_total, dumbbell_total  # noqa: E402
from scripts.valley_search import bouquet_poly  # noqa: E402


def fdiv(a, b):
    return float(Fraction(a, b))


def rgap(poly, thetas=(0.001, 0.01)):
    m = len(poly)
    prefs = [0] * m
    pref = poly[0]
    for b in range(1, m):
        prefs[b] = pref
        if poly[b] > pref:
            pref = poly[b]
    # suffix max value and its nearest achieving index
    suffs = [0] * m
    sidx = [0] * m
    sv, si = poly[-1], m - 1
    for b in range(m - 2, -1, -1):
        suffs[b], sidx[b] = sv, si
        if poly[b] > sv:
            sv, si = poly[b], b
    out = {}
    for th in thetas:
        best = None
        for b in range(1, m - 1):
            # exact integer comparison: prefs[b] >= (1+th) * poly[b]
            if prefs[b] * 1000 >= round((1 + th) * 1000) * poly[b]:
                key = (suffs[b], poly[b])
                if best is None or suffs[b] * best[1] > best[0] * poly[b]:
                    best = (suffs[b], poly[b], b, sidx[b])
        if best:
            R = fdiv(best[0], best[1])
            out[th] = {"R": R, "b": best[2], "c": best[3],
                       "witness": best[0] > best[1]}
        else:
            out[th] = None
    return out


def report(label, n, poly):
    r = rgap(poly)
    parts = [f"{label:52s} n={n:>5}"]
    for th, rec in r.items():
        if rec is None:
            parts.append(f" th={th}: no qualifying descent")
        else:
            wit = " *** WITNESS ***" if rec["witness"] else ""
            parts.append(f" th={th}: R={rec['R']:.8f} n(1-R)={n*(1-rec['R']):8.2f} "
                         f"b={rec['b']} c-b={rec['c']-rec['b']}{wit}")
    print(" | ".join(parts), flush=True)
    return r


def main():
    print("=== campaign champions, corrected metric ===")
    g = ((2,) * 3, (2,) * 77, tuple([2] * 75 + [3] * 3))
    report("bouquet S(2^3)+S(2^77)+S(2^75+3^3)", 323, bouquet_poly(list(g)))
    report("hybrid 34xS(2^2)+7xS(2^10)", 318,
           bouquet_total([{2: 2}] * 34 + [{2: 10}] * 7))
    report("hybrid 272xS(2^2)+56xS(2^10)", 2537,
           bouquet_total([{2: 2}] * 272 + [{2: 10}] * 56))
    report("dumbbell 48xS(3^6)--2--48xS(3^6)", 1827,
           dumbbell_total([{3: 6}] * 48, [{3: 6}] * 48, 2))
    report("dumbbell 224xS(2^7)--2--128xS(3^8)", 6563,
           dumbbell_total([{2: 7}] * 224, [{3: 8}] * 128, 2))

    print("\n=== depth-3 grid re-scored by R_gap(0.001) ===", flush=True)
    t0 = time.time()
    best = []
    for u in (2, 3, 4, 6):
        for s in (2, 3, 4, 6, 9, 14):
            for A in (1, 2, 3, 5, 8, 13, 21, 34):
                for t in (2, 4, 8, 12):
                    for B in (0, 1, 2, 5, 13, 34, 89):
                        n = size(A, s, u, B, t)
                        if not (120 <= n <= 1400):
                            continue
                        poly = tree_poly(A, s, u, B, t)
                        r = rgap(poly, thetas=(0.001,))[0.001]
                        if r is None:
                            continue
                        if r["witness"]:
                            print(f"*** WITNESS *** {(A,s,u,B,t)}")
                            print("poly=", poly, flush=True)
                            return
                        best.append((r["R"], (A, s, u, B, t), n, r))
    best.sort(reverse=True)
    for R, args, n, r in best[:15]:
        print(f"{str(args):>22} n={n:>5} R={R:.8f} n(1-R)={n*(1-R):8.2f} "
              f"b={r['b']} c-b={r['c']-r['b']}")
    print(f"grid rescore: {len(best)} configs in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
