#!/usr/bin/env python3
"""Scan global vs windowed covariance sign in leaf cases.

Global claim (false empirically): mu_g(lambda) <= mu_f(lambda) for all lambda.
Windowed claim (candidate): if mu_f(lambda) <= d(g)-1, then mu_g(lambda) <= mu_f(lambda).

Here:
  f = I(T), g = I(T-w) for leaf w,
  d(g) = first descent index of g.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from indpoly import _polymul, independence_poly
from trees import trees


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
    for t in range(0, len(bits), 6):
        v = 0
        for b in bits[t : t + 6]:
            v = (v << 1) | b
        out.append(chr(v + 63))
    return "".join(out)


def first_descent(seq: list[int]) -> int:
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            return i
    return len(seq)


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


def poly_profile(seq: list[int]) -> list[tuple[int, float]]:
    return [(k, math.log(c)) for k, c in enumerate(seq) if c > 0]


def mean_size(profile: list[tuple[int, float]], loglam: float) -> float:
    vals = [logc + k * loglam for k, logc in profile]
    m = max(vals)
    z = 0.0
    num = 0.0
    for (k, _), v in zip(profile, vals):
        w = math.exp(v - m)
        z += w
        num += k * w
    return num / z


def build_lambda_grid(min_exp: float, max_exp: float, count: int) -> list[float]:
    if count <= 1:
        return [1.0]
    out = []
    for i in range(count):
        t = i / (count - 1)
        e = min_exp + (max_exp - min_exp) * t
        out.append(10.0**e)
    if 1.0 not in out:
        out.append(1.0)
        out.sort()
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=2)
    ap.add_argument("--max-n", type=int, default=18)
    ap.add_argument("--backend", default="auto", choices=["auto", "geng", "networkx"])
    ap.add_argument("--lambda-min-exp", type=float, default=-6.0)
    ap.add_argument("--lambda-max-exp", type=float, default=6.0)
    ap.add_argument("--lambda-count", type=int, default=49)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--progress-every", type=int, default=10000)
    ap.add_argument("--max-witnesses", type=int, default=8)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    lambdas = build_lambda_grid(args.lambda_min_exp, args.lambda_max_exp, args.lambda_count)
    loglams = [math.log(x) for x in lambdas]

    stats: dict[str, Any] = {
        "global_claim": "mu_g(lambda) <= mu_f(lambda) for all tested lambda",
        "window_claim": "if mu_f(lambda) <= d(g)-1 then mu_g(lambda) <= mu_f(lambda)",
        "min_n": args.min_n,
        "max_n": args.max_n,
        "backend": args.backend,
        "lambda_grid": lambdas,
        "tol": args.tol,
        "trees": 0,
        "leaf_cases": 0,
        "all_checks": 0,
        "all_violations": 0,
        "window_checks": 0,
        "window_violations": 0,
        "max_mu_gap_all": -1e300,
        "max_mu_gap_window": -1e300,
        "all_witnesses": [],
        "window_witnesses": [],
        "by_n": [],
    }

    for n in range(args.min_n, args.max_n + 1):
        n_trees = 0
        n_leaf = 0
        n_all_checks = 0
        n_all_viol = 0
        n_win_checks = 0
        n_win_viol = 0
        n_max_all = -1e300
        n_max_win = -1e300

        for _, adj in trees(n, backend=args.backend):
            n_trees += 1
            f = independence_poly(n, adj)
            pf = poly_profile(f)

            leaves = [w for w in range(n) if len(adj[w]) == 1]
            for w in leaves:
                n_leaf += 1
                g = forest_poly_removed_set(adj, {w})
                pg = poly_profile(g)
                d_g = first_descent(g)

                for lam, loglam in zip(lambdas, loglams):
                    mu_f = mean_size(pf, loglam)
                    mu_g = mean_size(pg, loglam)
                    gap = mu_g - mu_f

                    n_all_checks += 1
                    if gap > n_max_all:
                        n_max_all = gap
                    if gap > stats["max_mu_gap_all"]:
                        stats["max_mu_gap_all"] = gap
                    if gap > args.tol:
                        n_all_viol += 1
                        if len(stats["all_witnesses"]) < args.max_witnesses:
                            stats["all_witnesses"].append(
                                {
                                    "n": n,
                                    "graph6": encode_graph6_small(adj),
                                    "leaf": w,
                                    "lambda": lam,
                                    "mu_f": mu_f,
                                    "mu_g": mu_g,
                                    "mu_gap": gap,
                                }
                            )

                    if mu_f <= d_g - 1 + args.tol:
                        n_win_checks += 1
                        if gap > n_max_win:
                            n_max_win = gap
                        if gap > stats["max_mu_gap_window"]:
                            stats["max_mu_gap_window"] = gap
                        if gap > args.tol:
                            n_win_viol += 1
                            if len(stats["window_witnesses"]) < args.max_witnesses:
                                stats["window_witnesses"].append(
                                    {
                                        "n": n,
                                        "graph6": encode_graph6_small(adj),
                                        "leaf": w,
                                        "d_g": d_g,
                                        "lambda": lam,
                                        "mu_f": mu_f,
                                        "mu_g": mu_g,
                                        "mu_gap": gap,
                                    }
                                )

            if args.progress_every > 0 and n_trees % args.progress_every == 0:
                print(
                    f"n={n} trees={n_trees:,} leaf_cases={n_leaf:,} "
                    f"all={n_all_checks:,}/{n_all_viol:,} "
                    f"window={n_win_checks:,}/{n_win_viol:,}"
                )

        stats["trees"] += n_trees
        stats["leaf_cases"] += n_leaf
        stats["all_checks"] += n_all_checks
        stats["all_violations"] += n_all_viol
        stats["window_checks"] += n_win_checks
        stats["window_violations"] += n_win_viol
        stats["by_n"].append(
            {
                "n": n,
                "trees": n_trees,
                "leaf_cases": n_leaf,
                "all_checks": n_all_checks,
                "all_violations": n_all_viol,
                "window_checks": n_win_checks,
                "window_violations": n_win_viol,
                "max_mu_gap_all": n_max_all,
                "max_mu_gap_window": n_max_win,
            }
        )

        print(
            f"done n={n}: trees={n_trees:,} leaf_cases={n_leaf:,} "
            f"all={n_all_checks:,}/{n_all_viol:,} "
            f"window={n_win_checks:,}/{n_win_viol:,}"
        )

    print(
        "summary: "
        f"all={stats['all_checks']:,}/{stats['all_violations']:,}, "
        f"window={stats['window_checks']:,}/{stats['window_violations']:,}, "
        f"max_all={stats['max_mu_gap_all']:.3e}, "
        f"max_window={stats['max_mu_gap_window']:.3e}"
    )

    if args.out:
        with open(args.out, "w", encoding="utf-8") as fobj:
            json.dump(stats, fobj, indent=2)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
