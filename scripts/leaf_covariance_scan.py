#!/usr/bin/env python3
"""Scan hard-core covariance sign for leaf vertices.

For a tree T and leaf w:
  f(x) = I(T;x), g(x) = I(T-w;x), Theta = g/f.

In the hard-core model at fugacity lambda:
  E[|S|]          = lambda f'(lambda)/f(lambda) =: mu_f(lambda)
  E[|S| | w notin S] = lambda g'(lambda)/g(lambda) =: mu_g(lambda)
  Cov(1_{w notin S}, |S|) = Theta(lambda) * (mu_g(lambda) - mu_f(lambda))
                          = lambda Theta'(lambda).

So covariance sign is determined by mu_g - mu_f.
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
    prof: list[tuple[int, float]] = []
    for k, c in enumerate(seq):
        if c > 0:
            prof.append((k, math.log(c)))
    return prof


def mean_size_from_profile(profile: list[tuple[int, float]], loglam: float) -> float:
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
    # guarantee exact 1.0 present
    if 1.0 not in out:
        out.append(1.0)
        out.sort()
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=2)
    ap.add_argument("--max-n", type=int, default=16)
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
        "claim": "For leaf w and all tested lambda, Cov(1_{w notin S}, |S|) <= 0",
        "equivalent": "mu_g(lambda) <= mu_f(lambda), where f=I(T), g=I(T-w)",
        "min_n": args.min_n,
        "max_n": args.max_n,
        "backend": args.backend,
        "lambda_grid": lambdas,
        "tol": args.tol,
        "trees": 0,
        "leaf_cases": 0,
        "lambda_checks": 0,
        "violations": 0,
        "max_mu_gap": -1e300,
        "witnesses": [],
        "by_n": [],
    }

    for n in range(args.min_n, args.max_n + 1):
        n_trees = 0
        n_leaf_cases = 0
        n_lambda_checks = 0
        n_viol = 0
        n_max_gap = -1e300

        for _, adj in trees(n, backend=args.backend):
            n_trees += 1
            f = independence_poly(n, adj)
            pf = poly_profile(f)
            leaves = [w for w in range(n) if len(adj[w]) == 1]
            if not leaves:
                continue

            for w in leaves:
                n_leaf_cases += 1
                g = forest_poly_removed_set(adj, {w})
                pg = poly_profile(g)

                for lam, loglam in zip(lambdas, loglams):
                    n_lambda_checks += 1
                    mu_f = mean_size_from_profile(pf, loglam)
                    mu_g = mean_size_from_profile(pg, loglam)
                    gap = mu_g - mu_f
                    if gap > n_max_gap:
                        n_max_gap = gap
                    if gap > stats["max_mu_gap"]:
                        stats["max_mu_gap"] = gap
                    if gap > args.tol:
                        n_viol += 1
                        if len(stats["witnesses"]) < args.max_witnesses:
                            stats["witnesses"].append(
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
                        break

            if args.progress_every > 0 and n_trees % args.progress_every == 0:
                print(
                    f"n={n} trees={n_trees:,} leaf_cases={n_leaf_cases:,} "
                    f"lambda_checks={n_lambda_checks:,} violations={n_viol:,} "
                    f"max_mu_gap={n_max_gap:.3e}"
                )

        stats["trees"] += n_trees
        stats["leaf_cases"] += n_leaf_cases
        stats["lambda_checks"] += n_lambda_checks
        stats["violations"] += n_viol
        stats["by_n"].append(
            {
                "n": n,
                "trees": n_trees,
                "leaf_cases": n_leaf_cases,
                "lambda_checks": n_lambda_checks,
                "violations": n_viol,
                "max_mu_gap": n_max_gap,
            }
        )

        print(
            f"done n={n}: trees={n_trees:,} leaf_cases={n_leaf_cases:,} "
            f"lambda_checks={n_lambda_checks:,} violations={n_viol:,} "
            f"max_mu_gap={n_max_gap:.3e}"
        )

    status = "PASS" if stats["violations"] == 0 else "FAIL"
    print(
        f"{status}: trees={stats['trees']:,} leaf_cases={stats['leaf_cases']:,} "
        f"lambda_checks={stats['lambda_checks']:,} violations={stats['violations']:,} "
        f"max_mu_gap={stats['max_mu_gap']:.3e}"
    )

    if args.out:
        with open(args.out, "w", encoding="utf-8") as fobj:
            json.dump(stats, fobj, indent=2)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
