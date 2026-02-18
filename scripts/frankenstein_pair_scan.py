#!/usr/bin/env python3
"""Scan synthetic recurrences h = A + xB over tree-polynomial libraries.

This is a relaxed inverse search:
  - A is a tree independence polynomial on n vertices.
  - B is a tree independence polynomial on m vertices (m < n).
  - h := A + xB.

We look for valleys in h. Any hit is a "Frankenstein" target profile to later
test for actual realizability constraints B = I(T-N[v]) from a common tree T.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from indpoly import independence_poly, is_unimodal
from trees import trees


def add_shifted(a: list[int], b: list[int]) -> list[int]:
    size = max(len(a), len(b) + 1)
    out = [0] * size
    for i, v in enumerate(a):
        out[i] += v
    for i, v in enumerate(b):
        out[i + 1] += v
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def first_valley(h: list[int]) -> int:
    for k in range(1, len(h) - 1):
        if h[k - 1] > h[k] < h[k + 1]:
            return k
    return -1


def unique_tree_polys(n: int, backend: str) -> list[list[int]]:
    seen: set[tuple[int, ...]] = set()
    out: list[list[int]] = []
    for _, adj in trees(n, backend=backend):
        p = tuple(independence_poly(n, adj))
        if p in seen:
            continue
        seen.add(p)
        out.append(list(p))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-a", type=int, default=8, help="Minimum |A|-tree size")
    ap.add_argument("--max-a", type=int, default=16, help="Maximum |A|-tree size")
    ap.add_argument("--min-b", type=int, default=2, help="Minimum |B|-tree size")
    ap.add_argument("--max-b", type=int, default=15, help="Maximum |B|-tree size")
    ap.add_argument("--backend", default="geng", choices=["geng", "networkx", "auto"])
    ap.add_argument("--require-b-le-a", action="store_true", help="Filter pairs with B_k <= A_k on shared support")
    ap.add_argument("--stop-on-first", action="store_true")
    ap.add_argument("--out", default="results/frankenstein_pair_scan.json")
    args = ap.parse_args()

    t0 = time.time()
    poly_by_n: dict[int, list[list[int]]] = {}
    for n in range(args.min_b, args.max_a + 1):
        poly_by_n[n] = unique_tree_polys(n, args.backend)
        print(f"n={n}: unique tree polynomials={len(poly_by_n[n])}")

    result: dict[str, Any] = {
        "min_a": args.min_a,
        "max_a": args.max_a,
        "min_b": args.min_b,
        "max_b": args.max_b,
        "backend": args.backend,
        "require_b_le_a": args.require_b_le_a,
        "tested_pairs": 0,
        "frankenstein_valleys": 0,
        "first_examples": [],
        "stats_by_a": [],
        "elapsed_s": None,
    }

    for a_n in range(args.min_a, args.max_a + 1):
        As = poly_by_n[a_n]
        row = {"a_n": a_n, "tested_pairs": 0, "valley_pairs": 0}
        for b_n in range(args.min_b, min(args.max_b, a_n - 1) + 1):
            Bs = poly_by_n[b_n]
            for a in As:
                for b in Bs:
                    if args.require_b_le_a:
                        if any(k >= len(a) or b[k] > a[k] for k in range(len(b))):
                            continue
                    row["tested_pairs"] += 1
                    result["tested_pairs"] += 1
                    h = add_shifted(a, b)
                    if is_unimodal(h):
                        continue
                    v = first_valley(h)
                    row["valley_pairs"] += 1
                    result["frankenstein_valleys"] += 1
                    if len(result["first_examples"]) < 20:
                        result["first_examples"].append(
                            {
                                "a_n": a_n,
                                "b_n": b_n,
                                "A": a,
                                "B": b,
                                "H": h,
                                "valley_index": v,
                            }
                        )
                    if args.stop_on_first:
                        result["stats_by_a"].append(row)
                        result["elapsed_s"] = time.time() - t0
                        with open(args.out, "w", encoding="utf-8") as f:
                            json.dump(result, f, indent=2)
                        print(f"wrote {args.out}")
                        return
        result["stats_by_a"].append(row)
        print(f"a_n={a_n}: tested={row['tested_pairs']} valleys={row['valley_pairs']}")

    result["elapsed_s"] = time.time() - t0
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
