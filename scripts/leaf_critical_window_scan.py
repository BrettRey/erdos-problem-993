#!/usr/bin/env python3
"""Scan leaf-step critical-window inequalities.

For each leaf w in T with neighbor u, define:
  f = I(T)
  g = I(T-w)
  q = I(T-N[w]) = I((T-w)-u)

This script mines candidate invariants around d(g), the first descent index of g.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from indpoly import _polymul, independence_poly
from trees import trees


def coeff(seq: list[int], i: int) -> int:
    if 0 <= i < len(seq):
        return seq[i]
    return 0


def delta_at(seq: list[int], i: int) -> int:
    return coeff(seq, i + 1) - coeff(seq, i)


def first_descent(seq: list[int]) -> int:
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            return i
    return len(seq)


def encode_graph6_small(adj: list[list[int]]) -> str:
    n = len(adj)
    if n >= 63:
        return "<n>=63"
    aset = [set(nei) for nei in adj]
    bits: list[int] = []
    for j in range(1, n):
        sj = aset[j]
        for i in range(j):
            bits.append(1 if i in sj else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [chr(n + 63)]
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k : k + 6]:
            v = (v << 1) | b
        out.append(chr(v + 63))
    return "".join(out)


def forest_poly_removed_set(adj: list[list[int]], removed: set[int]) -> list[int]:
    n = len(adj)
    remaining = [i for i in range(n) if i not in removed]
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
            x = stack.pop()
            comp.append(x)
            for y in adj[x]:
                if y in rem_set and y not in seen:
                    seen.add(y)
                    stack.append(y)
        mapping = {old: i for i, old in enumerate(comp)}
        cadj: list[list[int]] = [[] for _ in range(len(comp))]
        for old in comp:
            i = mapping[old]
            for y in adj[old]:
                j = mapping.get(y)
                if j is not None:
                    cadj[i].append(j)
        out = _polymul(out, independence_poly(len(comp), cadj))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=2)
    ap.add_argument("--max-n", type=int, default=18)
    ap.add_argument("--backend", default="auto", choices=["auto", "geng", "networkx"])
    ap.add_argument("--progress-every", type=int, default=10000)
    ap.add_argument("--max-witnesses", type=int, default=5)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    stats: dict[str, Any] = {
        "min_n": args.min_n,
        "max_n": args.max_n,
        "backend": args.backend,
        "trees": 0,
        "leaf_cases": 0,
        "gap_values": {},
        "h_delta_q_dminus1_nonpos_failures": 0,
        "h_boundary_transfer_failures": 0,
        "h_prewindow_nonneg_failures": 0,
        "h_argmin_at_boundary_failures": 0,
        "h_pref_window_nonincreasing_failures": 0,
        "gap_vs_dq_dminus1_sign": {},
        "gap_vs_boundary_margin_sign": {},
        "witnesses": {
            "delta_q_dminus1_nonpos": [],
            "boundary_transfer": [],
            "prewindow_nonneg": [],
            "argmin_at_boundary": [],
            "pref_window_nonincreasing": [],
        },
        "by_n": [],
    }

    for n in range(args.min_n, args.max_n + 1):
        n_trees = 0
        n_leaf_cases = 0
        n_gap_values: dict[str, int] = {}
        n_h_dq = 0
        n_h_boundary = 0
        n_h_prewindow = 0
        n_h_argmin = 0
        n_h_noninc = 0
        n_gap_vs_dq: dict[str, int] = {}
        n_gap_vs_boundary: dict[str, int] = {}

        for _, adj in trees(n, backend=args.backend):
            n_trees += 1
            f = independence_poly(n, adj)
            leaves = [w for w in range(n) if len(adj[w]) == 1]
            if not leaves:
                continue

            for w in leaves:
                n_leaf_cases += 1
                u = adj[w][0]
                g = forest_poly_removed_set(adj, {w})
                q = forest_poly_removed_set(adj, {w, u})

                d_f = first_descent(f)
                d_g = first_descent(g)

                gap = d_f - d_g
                n_gap_values[str(gap)] = n_gap_values.get(str(gap), 0) + 1
                stats["gap_values"][str(gap)] = stats["gap_values"].get(str(gap), 0) + 1

                dq_dminus1 = delta_at(q, d_g - 1)
                dq_sign = "pos" if dq_dminus1 > 0 else ("zero" if dq_dminus1 == 0 else "neg")
                key_gap_dq = f"{gap}:{dq_sign}"
                n_gap_vs_dq[key_gap_dq] = n_gap_vs_dq.get(key_gap_dq, 0) + 1
                stats["gap_vs_dq_dminus1_sign"][key_gap_dq] = (
                    stats["gap_vs_dq_dminus1_sign"].get(key_gap_dq, 0) + 1
                )

                # H1: Delta q_{d_g-1} <= 0
                if dq_dminus1 > 0:
                    n_h_dq += 1
                    if len(stats["witnesses"]["delta_q_dminus1_nonpos"]) < args.max_witnesses:
                        stats["witnesses"]["delta_q_dminus1_nonpos"].append(
                            {
                                "n": n,
                                "graph6": encode_graph6_small(adj),
                                "leaf": w,
                                "u": u,
                                "d_f": d_f,
                                "d_g": d_g,
                                "delta_q_dminus1": dq_dminus1,
                            }
                        )

                # H2: boundary transfer at prewindow edge k=d_g-2, i.e.
                # Delta f_{d_g-2} >= 0. This is the final index that can force
                # d(f) < d(g) in first-descent indexing.
                # Delta f_k = Delta g_k + Delta q_{k-1}.
                boundary_margin = None
                if d_g - 2 >= 0:
                    boundary_margin = delta_at(g, d_g - 2) + delta_at(q, d_g - 3)
                    b_sign = (
                        "pos" if boundary_margin > 0 else ("zero" if boundary_margin == 0 else "neg")
                    )
                    key_gap_b = f"{gap}:{b_sign}"
                    n_gap_vs_boundary[key_gap_b] = n_gap_vs_boundary.get(key_gap_b, 0) + 1
                    stats["gap_vs_boundary_margin_sign"][key_gap_b] = (
                        stats["gap_vs_boundary_margin_sign"].get(key_gap_b, 0) + 1
                    )

                if boundary_margin is not None and boundary_margin < 0:
                    n_h_boundary += 1
                    if len(stats["witnesses"]["boundary_transfer"]) < args.max_witnesses:
                        stats["witnesses"]["boundary_transfer"].append(
                            {
                                "n": n,
                                "graph6": encode_graph6_small(adj),
                                "leaf": w,
                                "u": u,
                                "d_f": d_f,
                                "d_g": d_g,
                                "delta_g_dminus2": delta_at(g, d_g - 2),
                                "delta_q_dminus3": delta_at(q, d_g - 3),
                            }
                        )

                # H3: all prewindow Delta f_k >= 0 for k <= d_g-2
                min_val = None
                argmin_k = None
                pref_noninc = True
                prev_val = None
                for k in range(max(0, d_g - 1)):
                    val = delta_at(g, k) + delta_at(q, k - 1)
                    if min_val is None or val < min_val:
                        min_val = val
                        argmin_k = k
                    if prev_val is not None and val > prev_val:
                        pref_noninc = False
                    prev_val = val

                if min_val is not None and min_val < 0:
                    n_h_prewindow += 1
                    if len(stats["witnesses"]["prewindow_nonneg"]) < args.max_witnesses:
                        stats["witnesses"]["prewindow_nonneg"].append(
                            {
                                "n": n,
                                "graph6": encode_graph6_small(adj),
                                "leaf": w,
                                "u": u,
                                "d_f": d_f,
                                "d_g": d_g,
                                "min_pref_delta_f": min_val,
                                "argmin_k": argmin_k,
                            }
                        )

                # H4: argmin over prewindow k<=d_g-2 occurs at boundary k=d_g-2
                if d_g > 1 and argmin_k is not None and argmin_k != d_g - 2:
                    n_h_argmin += 1
                    if len(stats["witnesses"]["argmin_at_boundary"]) < args.max_witnesses:
                        stats["witnesses"]["argmin_at_boundary"].append(
                            {
                                "n": n,
                                "graph6": encode_graph6_small(adj),
                                "leaf": w,
                                "u": u,
                                "d_f": d_f,
                                "d_g": d_g,
                                "argmin_k": argmin_k,
                                "boundary_k": d_g - 2,
                            }
                        )

                # H5: prewindow Delta f_k nonincreasing in k<d_g
                if not pref_noninc:
                    n_h_noninc += 1
                    if len(stats["witnesses"]["pref_window_nonincreasing"]) < args.max_witnesses:
                        stats["witnesses"]["pref_window_nonincreasing"].append(
                            {
                                "n": n,
                                "graph6": encode_graph6_small(adj),
                                "leaf": w,
                                "u": u,
                                "d_f": d_f,
                                "d_g": d_g,
                            }
                        )

            if args.progress_every > 0 and n_trees % args.progress_every == 0:
                print(
                    f"n={n} trees={n_trees:,} leaf_cases={n_leaf_cases:,} "
                    f"h_prewindow_fail={n_h_prewindow:,}"
                )

        stats["trees"] += n_trees
        stats["leaf_cases"] += n_leaf_cases
        stats["h_delta_q_dminus1_nonpos_failures"] += n_h_dq
        stats["h_boundary_transfer_failures"] += n_h_boundary
        stats["h_prewindow_nonneg_failures"] += n_h_prewindow
        stats["h_argmin_at_boundary_failures"] += n_h_argmin
        stats["h_pref_window_nonincreasing_failures"] += n_h_noninc

        stats["by_n"].append(
            {
                "n": n,
                "trees": n_trees,
                "leaf_cases": n_leaf_cases,
                "gap_values": n_gap_values,
                "h_delta_q_dminus1_nonpos_failures": n_h_dq,
                "h_boundary_transfer_failures": n_h_boundary,
                "h_prewindow_nonneg_failures": n_h_prewindow,
                "h_argmin_at_boundary_failures": n_h_argmin,
                "h_pref_window_nonincreasing_failures": n_h_noninc,
                "gap_vs_dq_dminus1_sign": n_gap_vs_dq,
                "gap_vs_boundary_margin_sign": n_gap_vs_boundary,
            }
        )

        print(
            f"done n={n}: trees={n_trees:,} leaf_cases={n_leaf_cases:,} "
            f"h_prewindow_fail={n_h_prewindow:,}"
        )

    status = "PASS" if stats["h_prewindow_nonneg_failures"] == 0 else "FAIL"
    print(
        f"{status}: trees={stats['trees']:,} leaf_cases={stats['leaf_cases']:,} "
        f"h_prewindow_fail={stats['h_prewindow_nonneg_failures']:,}"
    )

    if args.out:
        with open(args.out, "w", encoding="utf-8") as fobj:
            json.dump(stats, fobj, indent=2)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
