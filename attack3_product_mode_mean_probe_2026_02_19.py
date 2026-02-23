#!/usr/bin/env python3
"""Probe mode-vs-mean for products of tree independence polynomials.

This is a broader class study than the bridge decomposition:
  P(x) = prod_j I(T_j; x), with factors taken from non-isomorphic trees
  up to a chosen size bound.

Checks:
  mode(P) >= floor(mu_P(lambda))
for:
  1) a user-provided lambda grid, and
  2) optional tie fugacities lambda_k = p_{k-1}/p_k (with cutoff).

The probe helps separate what is true specifically at mode-tie lambdas
(typically <= 1 in this project) from behavior at larger fugacities.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import subprocess
import time
from typing import Any

from graph6 import parse_graph6
from indpoly import _polymul, independence_poly


def mode_index_leftmost(poly: list[int]) -> int:
    return max(range(len(poly)), key=lambda i: poly[i])


def mean_at_lambda(poly: list[int], lam: float) -> float:
    z = 0.0
    mu_num = 0.0
    p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p
        z += w
        mu_num += k * w
        p *= lam
    return mu_num / z if z else float("nan")


def parse_lambda_grid(s: str) -> list[float]:
    out: list[float] = []
    for tok in s.split(","):
        t = tok.strip()
        if not t:
            continue
        out.append(float(t))
    if not out:
        raise ValueError("lambda grid is empty")
    return out


def collect_unique_tree_polys(max_tree_n: int, geng: str) -> list[dict[str, Any]]:
    """Collect unique tree IS polynomials up to max_tree_n with one representative each."""
    reps: dict[tuple[int, ...], dict[str, Any]] = {}
    # n=1 manually
    p1 = (1, 1)
    reps[p1] = {"n": 1, "poly": [1, 1], "rep_g6": "@"}

    for n in range(2, max_tree_n + 1):
        cmd = [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        assert proc.stdout is not None
        for raw in proc.stdout:
            nn, adj = parse_graph6(raw)
            poly = tuple(independence_poly(nn, adj))
            if poly not in reps:
                reps[poly] = {
                    "n": nn,
                    "poly": list(poly),
                    "rep_g6": raw.decode("ascii").strip(),
                }
        proc.wait()
    out = list(reps.values())
    out.sort(key=lambda e: (e["n"], len(e["poly"]), e["poly"]))
    return out


def iter_combos(pool: list[dict[str, Any]], max_factors: int, max_total_vertices: int):
    """Yield non-empty nondecreasing index tuples with bounded size budget."""
    nvals = [e["n"] for e in pool]

    def rec(start: int, depth: int, nsum: int, cur: list[int]):
        if depth > 0:
            yield tuple(cur)
        if depth == max_factors:
            return
        for i in range(start, len(pool)):
            nn = nvals[i]
            if nsum + nn > max_total_vertices:
                continue
            cur.append(i)
            yield from rec(i, depth + 1, nsum + nn, cur)
            cur.pop()

    yield from rec(0, 0, 0, [])


def main() -> None:
    ap = argparse.ArgumentParser(description="Probe mode>=floor(mu_lambda) for products of tree IS polynomials.")
    ap.add_argument("--max-tree-n", type=int, default=8)
    ap.add_argument("--max-factors", type=int, default=3)
    ap.add_argument("--max-total-vertices", type=int, default=12)
    ap.add_argument("--lambda-grid", default="0.25,0.5,0.75,1.0,1.25,1.5")
    ap.add_argument("--check-ties", action="store_true")
    ap.add_argument("--tie-lambda-max", type=float, default=1.0)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--floor-eps", type=float, default=1e-12)
    ap.add_argument("--out", default="results/attack3_product_mode_mean_probe_2026_02_19.json")
    args = ap.parse_args()

    t0 = time.time()
    lambda_grid = parse_lambda_grid(args.lambda_grid)

    print(
        "attack3 product probe: "
        f"max_tree_n={args.max_tree_n}, max_factors={args.max_factors}, "
        f"max_total_vertices={args.max_total_vertices}",
        flush=True,
    )
    print(f"lambda grid: {lambda_grid}", flush=True)
    print(f"check ties: {args.check_ties} (tie_lambda_max={args.tie_lambda_max})", flush=True)
    print("-" * 96, flush=True)

    pool = collect_unique_tree_polys(args.max_tree_n, args.geng)
    print(f"unique factor pool size: {len(pool)}", flush=True)

    lambda_stats: dict[str, dict[str, Any]] = {}
    for lam in lambda_grid:
        key = f"{lam:.12g}"
        lambda_stats[key] = {
            "lambda": lam,
            "tested": 0,
            "fail": 0,
            "min_margin_mode_minus_floor_mu": None,
            "worst_witness": None,
        }

    tie_stats: dict[str, Any] = {
        "enabled": args.check_ties,
        "tie_lambda_max": args.tie_lambda_max,
        "tested": 0,
        "fail": 0,
        "min_margin_mode_minus_floor_mu": None,
        "worst_witness": None,
    }

    combo_count = 0
    for combo in iter_combos(pool, args.max_factors, args.max_total_vertices):
        combo_count += 1
        factors = [pool[i] for i in combo]
        prod = [1]
        total_vertices = 0
        for f in factors:
            prod = _polymul(prod, f["poly"])
            total_vertices += int(f["n"])

        mode_p = mode_index_leftmost(prod)
        factor_ns = [int(f["n"]) for f in factors]
        factor_reps = [str(f["rep_g6"]) for f in factors]

        for lam in lambda_grid:
            key = f"{lam:.12g}"
            ls = lambda_stats[key]
            mu = mean_at_lambda(prod, lam)
            floor_mu = math.floor(mu + args.floor_eps)
            margin = mode_p - floor_mu
            ls["tested"] += 1

            if ls["min_margin_mode_minus_floor_mu"] is None or margin < ls["min_margin_mode_minus_floor_mu"]:
                ls["min_margin_mode_minus_floor_mu"] = margin
                ls["worst_witness"] = {
                    "factor_ns": factor_ns,
                    "factor_reps": factor_reps,
                    "total_vertices": total_vertices,
                    "mode_P": mode_p,
                    "mu_lambda": mu,
                    "floor_mu": floor_mu,
                    "margin": margin,
                    "lambda": lam,
                    "poly": prod,
                }

            if margin < 0:
                ls["fail"] += 1

        if args.check_ties:
            for k in range(1, len(prod)):
                if prod[k - 1] <= 0 or prod[k] <= 0:
                    continue
                lam_tie = prod[k - 1] / prod[k]
                if lam_tie > args.tie_lambda_max + args.tol:
                    continue
                mu = mean_at_lambda(prod, lam_tie)
                floor_mu = math.floor(mu + args.floor_eps)
                margin = mode_p - floor_mu
                tie_stats["tested"] += 1

                if (
                    tie_stats["min_margin_mode_minus_floor_mu"] is None
                    or margin < tie_stats["min_margin_mode_minus_floor_mu"]
                ):
                    tie_stats["min_margin_mode_minus_floor_mu"] = margin
                    tie_stats["worst_witness"] = {
                        "factor_ns": factor_ns,
                        "factor_reps": factor_reps,
                        "total_vertices": total_vertices,
                        "mode_P": mode_p,
                        "k": k,
                        "lambda_tie": lam_tie,
                        "mu_lambda": mu,
                        "floor_mu": floor_mu,
                        "margin": margin,
                        "poly": prod,
                    }

                if margin < 0:
                    tie_stats["fail"] += 1

    summary = {
        "params": {
            "max_tree_n": args.max_tree_n,
            "max_factors": args.max_factors,
            "max_total_vertices": args.max_total_vertices,
            "lambda_grid": lambda_grid,
            "check_ties": args.check_ties,
            "tie_lambda_max": args.tie_lambda_max,
            "tol": args.tol,
            "floor_eps": args.floor_eps,
        },
        "factor_pool_size": len(pool),
        "combo_count": combo_count,
        "lambda_stats": lambda_stats,
        "tie_stats": tie_stats,
        "wall_s": time.time() - t0,
    }

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"combos tested: {combo_count:,}", flush=True)
    for key in sorted(lambda_stats, key=lambda s: float(s)):
        ls = lambda_stats[key]
        print(
            f"lambda={ls['lambda']:.6g}: tested={ls['tested']:,} fail={ls['fail']:,} "
            f"min_margin={ls['min_margin_mode_minus_floor_mu']}",
            flush=True,
        )
    if args.check_ties:
        print(
            f"ties<= {args.tie_lambda_max:.6g}: tested={tie_stats['tested']:,} "
            f"fail={tie_stats['fail']:,} "
            f"min_margin={tie_stats['min_margin_mode_minus_floor_mu']}",
            flush=True,
        )
    print(f"results written to: {args.out}", flush=True)


if __name__ == "__main__":
    main()
