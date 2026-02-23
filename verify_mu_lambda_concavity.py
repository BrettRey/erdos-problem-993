#!/usr/bin/env python3
"""Check discrete concavity/convexity of mu(lambda) for tree IS polynomials.

Given I(x)=sum_k c_k x^k and mu(lambda)=lambda I'(lambda)/I(lambda), this script
samples mu on a lambda-grid and checks second finite differences:
  d2_i = mu_{i+1} - 2*mu_i + mu_{i-1}

Concavity (on the sampled grid): d2_i <= tol for all i
Convexity (on the sampled grid): d2_i >= -tol for all i
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from typing import Any

from graph6 import parse_graph6
from indpoly import independence_poly


def mean_at_lambda(poly: list[int], lam: float) -> float:
    z = 0.0
    mu_num = 0.0
    p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p
        z += w
        mu_num += k * w
        p *= lam
    return mu_num / z if z else 0.0


def profile_d2(poly: list[int], lams: list[float]) -> tuple[float, float]:
    mus = [mean_at_lambda(poly, lam) for lam in lams]
    min_d2 = float("inf")
    max_d2 = float("-inf")
    for i in range(1, len(mus) - 1):
        d2 = mus[i + 1] - 2.0 * mus[i] + mus[i - 1]
        if d2 < min_d2:
            min_d2 = d2
        if d2 > max_d2:
            max_d2 = d2
    return min_d2, max_d2


def fresh_stats() -> dict[str, Any]:
    return {
        "seen": 0,
        "checked": 0,
        "concavity_fail": 0,
        "convexity_fail": 0,
        "max_d2": None,
        "max_d2_witness": None,
        "min_d2": None,
        "min_d2_witness": None,
        "wall_s": 0.0,
    }


def maybe_update_max(
    stats: dict[str, Any],
    key: str,
    witness_key: str,
    value: float,
    witness: dict[str, Any],
) -> None:
    cur = stats[key]
    if cur is None or value > cur:
        stats[key] = value
        stats[witness_key] = witness


def maybe_update_min(
    stats: dict[str, Any],
    key: str,
    witness_key: str,
    value: float,
    witness: dict[str, Any],
) -> None:
    cur = stats[key]
    if cur is None or value < cur:
        stats[key] = value
        stats[witness_key] = witness


def main() -> None:
    ap = argparse.ArgumentParser(description="Concavity/convexity sampler for mu(lambda) on tree IS polynomials.")
    ap.add_argument("--min-n", type=int, default=1)
    ap.add_argument("--max-n", type=int, default=16)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--lam-min", type=float, default=0.1)
    ap.add_argument("--lam-max", type=float, default=1.9)
    ap.add_argument("--lam-steps", type=int, default=19)
    ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--out", default="results/mu_lambda_concavity_scan.json")
    args = ap.parse_args()

    if args.lam_steps < 3:
        raise ValueError("--lam-steps must be >= 3")
    if args.lam_min <= 0 or args.lam_max <= args.lam_min:
        raise ValueError("Require 0 < lam-min < lam-max")

    lams = [
        args.lam_min + (args.lam_max - args.lam_min) * i / (args.lam_steps - 1)
        for i in range(args.lam_steps)
    ]

    per_n: dict[str, dict[str, Any]] = {}
    t_all = time.time()
    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        stats = fresh_stats()

        proc = subprocess.Popen(
            [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"],
            stdout=subprocess.PIPE,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            stats["seen"] += 1
            nn, adj = parse_graph6(line)
            poly = independence_poly(nn, adj)
            if len(poly) <= 1:
                continue
            stats["checked"] += 1

            min_d2, max_d2 = profile_d2(poly, lams)
            g6 = line.decode("ascii").strip()
            witness = {
                "n": nn,
                "g6": g6,
                "min_d2": min_d2,
                "max_d2": max_d2,
            }
            maybe_update_max(stats, "max_d2", "max_d2_witness", max_d2, witness)
            maybe_update_min(stats, "min_d2", "min_d2_witness", min_d2, witness)

            if max_d2 > args.tol:
                stats["concavity_fail"] += 1
            if min_d2 < -args.tol:
                stats["convexity_fail"] += 1

        proc.wait()
        stats["wall_s"] = time.time() - t0
        per_n[str(n)] = stats
        print(
            f"n={n:2d}: seen={stats['seen']:9d} checked={stats['checked']:9d} "
            f"conc_fail={stats['concavity_fail']:5d} conv_fail={stats['convexity_fail']:5d} "
            f"max_d2={stats['max_d2']} min_d2={stats['min_d2']} ({stats['wall_s']:.1f}s)",
            flush=True,
        )

    # aggregate
    tot = fresh_stats()
    for n_key in sorted(per_n, key=lambda s: int(s)):
        s = per_n[n_key]
        for k in ["seen", "checked", "concavity_fail", "convexity_fail"]:
            tot[k] += int(s[k])
        if s["max_d2"] is not None:
            maybe_update_max(tot, "max_d2", "max_d2_witness", float(s["max_d2"]), s["max_d2_witness"])
        if s["min_d2"] is not None:
            maybe_update_min(tot, "min_d2", "min_d2_witness", float(s["min_d2"]), s["min_d2_witness"])
    tot.pop("wall_s", None)

    payload = {
        "params": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "lam_min": args.lam_min,
            "lam_max": args.lam_max,
            "lam_steps": args.lam_steps,
            "tol": args.tol,
        },
        "summary": tot,
        "per_n": {k: per_n[k] for k in sorted(per_n, key=lambda s: int(s))},
    }

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print("-" * 96, flush=True)
    print(
        f"TOTAL seen={tot['seen']:,} checked={tot['checked']:,} "
        f"conc_fail={tot['concavity_fail']:,} conv_fail={tot['convexity_fail']:,} "
        f"max_d2={tot['max_d2']} min_d2={tot['min_d2']} wall={time.time()-t_all:.1f}s",
        flush=True,
    )
    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
