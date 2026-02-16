#!/usr/bin/env python3
"""Partitioned mode-alignment scan.

Checks mode(I(T-w)) <= d(I(T)) for one geng partition (res/mod).
This is intended for distributed runs across multiple agents/processes.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from indpoly import _polymul, independence_poly
from trees import trees_geng_raw


def mode_index(seq: list[int]) -> int:
    if not seq:
        return -1
    maxv = max(seq)
    idx = -1
    for i, v in enumerate(seq):
        if v == maxv:
            idx = i
    return idx


def first_descent(seq: list[int]) -> int:
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            return i
    return len(seq)


def forest_poly_removed(adj: list[list[int]], w: int) -> list[int]:
    n = len(adj)
    remaining = [i for i in range(n) if i != w]
    if not remaining:
        return [1]
    rem_set = set(remaining)
    seen: set[int] = set()
    out = [1]
    for start in remaining:
        if start in seen:
            continue
        comp: list[int] = []
        stack = [start]
        seen.add(start)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v in rem_set and v not in seen:
                    seen.add(v)
                    stack.append(v)
        mapping = {old: i for i, old in enumerate(comp)}
        cadj: list[list[int]] = [[] for _ in range(len(comp))]
        for old in comp:
            i = mapping[old]
            for v in adj[old]:
                j = mapping.get(v)
                if j is not None:
                    cadj[i].append(j)
        out = _polymul(out, independence_poly(len(comp), cadj))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--res", type=int, required=True)
    ap.add_argument("--mod", type=int, required=True)
    ap.add_argument("--progress-every", type=int, default=50000)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    if args.mod <= 0:
        raise ValueError("--mod must be > 0")
    if not (0 <= args.res < args.mod):
        raise ValueError("--res must satisfy 0 <= res < mod")

    if not args.out:
        args.out = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "results",
            f"mode_alignment_n{args.n}_res{args.res}_mod{args.mod}.json",
        )

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    t0 = time.time()

    trees_checked = 0
    vertex_cases = 0
    failures = 0
    equality_count = 0
    max_mode_gap = -(10**9)
    max_abs_mode_diff = 0
    sample_failures: list[dict[str, Any]] = []
    sample_equalities: list[dict[str, Any]] = []

    print(
        f"partition scan n={args.n} res={args.res}/{args.mod} "
        f"(progress every {args.progress_every:,} trees)"
    )

    for tree_n, adj, raw in trees_geng_raw(args.n, res=args.res, mod=args.mod):
        trees_checked += 1
        g6 = raw.decode("ascii", errors="replace")
        f = independence_poly(tree_n, adj)
        mode_f = mode_index(f)
        d_f = first_descent(f)

        for w in range(tree_n):
            vertex_cases += 1
            g = forest_poly_removed(adj, w)
            mode_g = mode_index(g)
            gap = mode_g - d_f
            if gap > max_mode_gap:
                max_mode_gap = gap
            abs_diff = abs(mode_g - mode_f)
            if abs_diff > max_abs_mode_diff:
                max_abs_mode_diff = abs_diff

            if gap > 0:
                failures += 1
                if len(sample_failures) < 5:
                    sample_failures.append(
                        {
                            "graph6": g6,
                            "vertex": w,
                            "mode_Tw": mode_g,
                            "d_IT": d_f,
                            "mode_IT": mode_f,
                        }
                    )
            elif gap == 0:
                equality_count += 1
                if len(sample_equalities) < 5:
                    sample_equalities.append(
                        {
                            "graph6": g6,
                            "vertex": w,
                            "mode_Tw": mode_g,
                            "d_IT": d_f,
                            "mode_IT": mode_f,
                        }
                    )

        if args.progress_every > 0 and trees_checked % args.progress_every == 0:
            elapsed = time.time() - t0
            print(
                f"  trees={trees_checked:,} vertex_cases={vertex_cases:,} "
                f"failures={failures:,} max_gap={max_mode_gap} "
                f"max|diff|={max_abs_mode_diff} elapsed={elapsed:.1f}s"
            )

    elapsed = time.time() - t0
    result = {
        "claim": "mode(I(T-w)) <= d(I(T))",
        "n": args.n,
        "res": args.res,
        "mod": args.mod,
        "trees": trees_checked,
        "vertex_cases": vertex_cases,
        "failures": failures,
        "equality_count": equality_count,
        "max_mode_gap": max_mode_gap,
        "max_abs_mode_diff": max_abs_mode_diff,
        "sample_failures": sample_failures,
        "sample_equalities": sample_equalities,
        "wall_time_seconds": round(elapsed, 1),
    }
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    status = "PASS" if failures == 0 else "FAIL"
    print(
        f"done [{status}] trees={trees_checked:,} vertex_cases={vertex_cases:,} "
        f"failures={failures:,} max_gap={max_mode_gap} "
        f"max|diff|={max_abs_mode_diff} elapsed={elapsed:.1f}s"
    )
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()

