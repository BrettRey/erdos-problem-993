#!/usr/bin/env python3
"""Find counterexamples to H_CW => H_R among unconstrained tree-polynomial pairs.

This search uses pairs (f,g) where:
  - f is a tree independence polynomial on n vertices,
  - g is a tree independence polynomial on n-1 vertices,
  - g_k <= f_k coefficientwise on common support.

It does NOT enforce leaf-realizable linkage (f = g + x q with q from the same T).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from indpoly import independence_poly
from trees import trees


def first_descent(seq: list[int]) -> int:
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            return i
    return len(seq)


def hr_failure_index(f: list[int], g: list[int]) -> int:
    d_g = first_descent(g)
    kmax = min(d_g - 2, min(len(f), len(g)) - 2)
    for k in range(kmax + 1):
        if g[k + 1] * f[k] > g[k] * f[k + 1]:
            return k
    return -1


def hcw_holds_on_grid(f: list[int], g: list[int], grid_count: int = 49, tol: float = 1e-12) -> bool:
    d_g = first_descent(g)
    for i in range(grid_count):
        lam = 10 ** (-6 + 12 * i / (grid_count - 1))
        fk = [c * (lam**k) for k, c in enumerate(f)]
        z = sum(fk)
        p = [v / z for v in fk]
        mu = sum(k * p[k] for k in range(len(p)))
        if mu > d_g - 1 + tol:
            continue
        r = [g[k] / f[k] for k in range(len(g))]
        er = sum(r[k] * p[k] for k in range(len(r)))
        ekr = sum(k * r[k] * p[k] for k in range(len(r)))
        cov = ekr - mu * er
        if cov > tol:
            return False
    return True


def unique_tree_polys(n: int) -> list[list[int]]:
    seen: set[tuple[int, ...]] = set()
    out: list[list[int]] = []
    for _, adj in trees(n, backend="geng"):
        p = tuple(independence_poly(n, adj))
        if p not in seen:
            seen.add(p)
            out.append(list(p))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=6)
    ap.add_argument("--max-n", type=int, default=12)
    ap.add_argument("--grid-count", type=int, default=49)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    poly_by_n: dict[int, list[list[int]]] = {}
    for n in range(2, args.max_n + 1):
        poly_by_n[n] = unique_tree_polys(n)
        print(f"n={n}: unique tree polynomials={len(poly_by_n[n])}")

    result: dict[str, Any] = {
        "claim": "H_CW does not imply H_R even when f and g are individual tree polynomials",
        "min_n": args.min_n,
        "max_n": args.max_n,
        "grid_count": args.grid_count,
        "found": False,
        "counterexample": None,
        "stats_by_n": [],
    }

    for n in range(args.min_n, args.max_n + 1):
        Fs = poly_by_n[n]
        Gs = poly_by_n[n - 1]
        tested = 0
        hr_fail = 0
        for f in Fs:
            for g in Gs:
                if len(g) > len(f):
                    continue
                if any(g[k] > f[k] for k in range(len(g))):
                    continue
                tested += 1
                bad_k = hr_failure_index(f, g)
                if bad_k < 0:
                    continue
                hr_fail += 1
                if hcw_holds_on_grid(f, g, grid_count=args.grid_count, tol=args.tol):
                    result["found"] = True
                    result["counterexample"] = {
                        "n": n,
                        "hr_bad_k": bad_k,
                        "f": f,
                        "g": g,
                        "d_g": first_descent(g),
                        "note": "Pair is not guaranteed leaf-realizable from a common tree.",
                    }
                    result["stats_by_n"].append({"n": n, "tested_pairs": tested, "hr_fail_pairs": hr_fail})
                    if args.out:
                        with open(args.out, "w", encoding="utf-8") as fobj:
                            json.dump(result, fobj, indent=2)
                    print(json.dumps(result, indent=2))
                    return
        result["stats_by_n"].append({"n": n, "tested_pairs": tested, "hr_fail_pairs": hr_fail})
        print(f"n={n}: tested_pairs={tested}, hr_fail_pairs={hr_fail}")

    if args.out:
        with open(args.out, "w", encoding="utf-8") as fobj:
            json.dump(result, fobj, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
