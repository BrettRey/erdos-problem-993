#!/usr/bin/env python3
"""Exhaustive equivalence check for tight mode characterization.

Checks, for every tree T and vertex w:

  E(T,w): mode(I(T-w)) = d(I(T))

against the compact condition set:

  C(T,w): alpha(T) = alpha(T-w)
          and d(I(T)) = mode(I(T)) + 1
          and deg_T(w) >= 2.

Outputs a JSON report containing either:
  - counterexamples to equivalence (both directions, if present), or
  - a confirmed-equivalence summary.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from multiprocessing import Pool
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from indpoly import _polymul, independence_poly
from trees import trees_geng_raw


# OEIS A000055: number of unlabeled trees on n vertices.
A000055: dict[int, int] = {
    1: 1,
    2: 1,
    3: 1,
    4: 2,
    5: 3,
    6: 6,
    7: 11,
    8: 23,
    9: 47,
    10: 106,
    11: 235,
    12: 551,
    13: 1301,
    14: 3159,
    15: 7741,
    16: 19320,
    17: 48629,
    18: 123867,
    19: 317955,
    20: 823065,
}


def mode_index(seq: list[int]) -> int:
    """Return the last index achieving the maximum coefficient."""
    if not seq:
        return -1
    maxv = max(seq)
    idx = -1
    for i, val in enumerate(seq):
        if val == maxv:
            idx = i
    return idx


def first_descent(seq: list[int]) -> int:
    """First index i with seq[i] < seq[i-1]; len(seq) if none."""
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            return i
    return len(seq)


def forest_poly_removed(adj: list[list[int]], w: int) -> list[int]:
    """Compute I(T-w) by removing vertex w and multiplying component polys."""
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


def mismatch_record(
    *,
    n: int,
    g6: str,
    w: int,
    lhs: bool,
    rhs: bool,
    mode_tw: int,
    d_it: int,
    mode_it: int,
    alpha_t: int,
    alpha_tw: int,
    deg_w: int,
) -> dict[str, Any]:
    """Create a compact mismatch witness record."""
    cond_alpha = alpha_t == alpha_tw
    cond_d_mode = d_it == mode_it + 1
    cond_deg = deg_w >= 2
    return {
        "n": n,
        "graph6": g6,
        "vertex": w,
        "lhs_mode_Tw_eq_d_IT": lhs,
        "rhs_compact_conditions": rhs,
        "mode_ITw": mode_tw,
        "d_IT": d_it,
        "mode_IT": mode_it,
        "alpha_T": alpha_t,
        "alpha_Tw": alpha_tw,
        "deg_w": deg_w,
        "cond_alpha_preserved": cond_alpha,
        "cond_d_equals_mode_plus_1": cond_d_mode,
        "cond_deg_w_ge_2": cond_deg,
    }


def worker_partition(args: tuple[int, int, int, int]) -> dict[str, Any]:
    """Check one geng partition n,res/mod."""
    n, res, mod, sample_cap = args

    tree_count = 0
    vertex_cases = 0

    lhs_true_rhs_false = 0
    lhs_false_rhs_true = 0
    lhs_true = 0
    rhs_true = 0
    both_true = 0
    both_false = 0

    samples_lhs_true_rhs_false: list[dict[str, Any]] = []
    samples_lhs_false_rhs_true: list[dict[str, Any]] = []

    for tree_n, adj, raw in trees_geng_raw(n, res=res, mod=mod):
        tree_count += 1
        g6 = raw.decode("ascii", errors="replace")

        f = independence_poly(tree_n, adj)
        d_f = first_descent(f)
        mode_f = mode_index(f)
        alpha_t = len(f) - 1
        degs = [len(nei) for nei in adj]

        for w in range(tree_n):
            vertex_cases += 1
            g = forest_poly_removed(adj, w)
            mode_g = mode_index(g)
            alpha_tw = len(g) - 1

            lhs = (mode_g == d_f)
            rhs = (
                alpha_t == alpha_tw
                and d_f == mode_f + 1
                and degs[w] >= 2
            )

            if lhs:
                lhs_true += 1
            if rhs:
                rhs_true += 1
            if lhs and rhs:
                both_true += 1
            elif (not lhs) and (not rhs):
                both_false += 1
            elif lhs and (not rhs):
                lhs_true_rhs_false += 1
                if len(samples_lhs_true_rhs_false) < sample_cap:
                    samples_lhs_true_rhs_false.append(
                        mismatch_record(
                            n=tree_n,
                            g6=g6,
                            w=w,
                            lhs=lhs,
                            rhs=rhs,
                            mode_tw=mode_g,
                            d_it=d_f,
                            mode_it=mode_f,
                            alpha_t=alpha_t,
                            alpha_tw=alpha_tw,
                            deg_w=degs[w],
                        )
                    )
            else:
                lhs_false_rhs_true += 1
                if len(samples_lhs_false_rhs_true) < sample_cap:
                    samples_lhs_false_rhs_true.append(
                        mismatch_record(
                            n=tree_n,
                            g6=g6,
                            w=w,
                            lhs=lhs,
                            rhs=rhs,
                            mode_tw=mode_g,
                            d_it=d_f,
                            mode_it=mode_f,
                            alpha_t=alpha_t,
                            alpha_tw=alpha_tw,
                            deg_w=degs[w],
                        )
                    )

    return {
        "tree_count": tree_count,
        "vertex_cases": vertex_cases,
        "lhs_true_rhs_false": lhs_true_rhs_false,
        "lhs_false_rhs_true": lhs_false_rhs_true,
        "lhs_true": lhs_true,
        "rhs_true": rhs_true,
        "both_true": both_true,
        "both_false": both_false,
        "samples_lhs_true_rhs_false": samples_lhs_true_rhs_false,
        "samples_lhs_false_rhs_true": samples_lhs_false_rhs_true,
    }


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Exhaustively test tight-mode equivalence characterization."
    )
    ap.add_argument("--min-n", type=int, default=1)
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--samples-per-direction", type=int, default=20)
    ap.add_argument(
        "--out",
        default="results/tight_mode_equivalence_n20.json",
        help="Output JSON artifact path.",
    )
    args = ap.parse_args()

    if args.max_n < args.min_n:
        raise ValueError("--max-n must be >= --min-n")
    if args.workers <= 0:
        raise ValueError("--workers must be > 0")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    started = time.time()

    total_trees = 0
    total_vertex_cases = 0
    lhs_true_rhs_false = 0
    lhs_false_rhs_true = 0
    lhs_true = 0
    rhs_true = 0
    both_true = 0
    both_false = 0

    samples_lhs_true_rhs_false: list[dict[str, Any]] = []
    samples_lhs_false_rhs_true: list[dict[str, Any]] = []
    per_n: dict[str, Any] = {}

    for n in range(args.min_n, args.max_n + 1):
        if n == 1:
            # Single tree with one vertex.
            tree_count = 1
            vertex_cases = 1

            # T has I=[1,1], so mode=1, d=2 (no descent => len=2).
            # T-w is empty forest with I=[1], mode=0, alpha=0, deg(w)=0.
            local_lhs = False
            local_rhs = False
            local_lhs_true_rhs_false = 0
            local_lhs_false_rhs_true = 0
            local_lhs_true = 0
            local_rhs_true = 0
            local_both_true = 0
            local_both_false = 1
        else:
            tasks = [
                (n, res, args.workers, args.samples_per_direction)
                for res in range(args.workers)
            ]

            tree_count = 0
            vertex_cases = 0
            local_lhs_true_rhs_false = 0
            local_lhs_false_rhs_true = 0
            local_lhs_true = 0
            local_rhs_true = 0
            local_both_true = 0
            local_both_false = 0

            with Pool(args.workers) as pool:
                for rec in pool.imap_unordered(worker_partition, tasks):
                    tree_count += rec["tree_count"]
                    vertex_cases += rec["vertex_cases"]
                    local_lhs_true_rhs_false += rec["lhs_true_rhs_false"]
                    local_lhs_false_rhs_true += rec["lhs_false_rhs_true"]
                    local_lhs_true += rec["lhs_true"]
                    local_rhs_true += rec["rhs_true"]
                    local_both_true += rec["both_true"]
                    local_both_false += rec["both_false"]

                    if len(samples_lhs_true_rhs_false) < args.samples_per_direction:
                        need = args.samples_per_direction - len(samples_lhs_true_rhs_false)
                        samples_lhs_true_rhs_false.extend(
                            rec["samples_lhs_true_rhs_false"][:need]
                        )
                    if len(samples_lhs_false_rhs_true) < args.samples_per_direction:
                        need = args.samples_per_direction - len(samples_lhs_false_rhs_true)
                        samples_lhs_false_rhs_true.extend(
                            rec["samples_lhs_false_rhs_true"][:need]
                        )

        expected = A000055.get(n)
        if expected is not None and tree_count != expected:
            raise RuntimeError(
                f"Tree count mismatch at n={n}: got {tree_count:,}, expected {expected:,}"
            )

        per_n[str(n)] = {
            "trees": tree_count,
            "vertex_cases": vertex_cases,
            "lhs_true_rhs_false": local_lhs_true_rhs_false,
            "lhs_false_rhs_true": local_lhs_false_rhs_true,
            "lhs_true": local_lhs_true,
            "rhs_true": local_rhs_true,
            "both_true": local_both_true,
            "both_false": local_both_false,
            "equivalence_holds": (
                local_lhs_true_rhs_false == 0 and local_lhs_false_rhs_true == 0
            ),
        }

        total_trees += tree_count
        total_vertex_cases += vertex_cases
        lhs_true_rhs_false += local_lhs_true_rhs_false
        lhs_false_rhs_true += local_lhs_false_rhs_true
        lhs_true += local_lhs_true
        rhs_true += local_rhs_true
        both_true += local_both_true
        both_false += local_both_false

        print(
            f"n={n}: trees={tree_count:,} vertex_cases={vertex_cases:,} "
            f"lhs_true_rhs_false={local_lhs_true_rhs_false:,} "
            f"lhs_false_rhs_true={local_lhs_false_rhs_true:,}"
        )

    mismatches = lhs_true_rhs_false + lhs_false_rhs_true
    confirmed = mismatches == 0

    report = {
        "claim": (
            "mode(I(T-w)) = d(I(T)) iff "
            "(alpha(T)=alpha(T-w) and d(I(T))=mode(I(T))+1 and deg(w)>=2)"
        ),
        "min_n": args.min_n,
        "max_n": args.max_n,
        "workers": args.workers,
        "total_trees": total_trees,
        "total_vertex_cases": total_vertex_cases,
        "contingency": {
            "lhs_true_rhs_true": both_true,
            "lhs_false_rhs_false": both_false,
            "lhs_true_rhs_false": lhs_true_rhs_false,
            "lhs_false_rhs_true": lhs_false_rhs_true,
            "lhs_true_total": lhs_true,
            "rhs_true_total": rhs_true,
        },
        "equivalence_confirmed": confirmed,
        "counterexamples": {
            "lhs_true_rhs_false_examples": samples_lhs_true_rhs_false,
            "lhs_false_rhs_true_examples": samples_lhs_false_rhs_true,
        },
        "per_n_counts": per_n,
        "wall_time_seconds": round(time.time() - started, 2),
    }

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    status = "CONFIRMED" if confirmed else "COUNTEREXAMPLES FOUND"
    print(f"\n[{status}] wrote {args.out}")


if __name__ == "__main__":
    main()
