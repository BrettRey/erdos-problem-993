#!/usr/bin/env python3
"""Scan tie-point mean inequalities for tree IS distributions.

For polynomial I_T(x) = sum_k i_k x^k and each k>=1 with i_{k-1}, i_k > 0,
define tie fugacity:

    lambda_k = i_{k-1} / i_k,

so the weighted coefficients satisfy i_{k-1} lambda_k^{k-1} = i_k lambda_k^k.

At this lambda, define the IS-size mean:

    mu(lambda_k) = sum_j j * i_j lambda_k^j / sum_j i_j lambda_k^j.

This script checks (on d_leaf<=1 trees):
  (A) lower tie bound: mu(lambda_k) >= k-1
  (B) upper tie bound: mu(lambda_k) <= k

These are natural candidates for a direct mode-vs-mean proof without LC.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import time
from typing import Any

from conjecture_a_hall_subset_scan import is_dleaf_le_1, parse_graph6
from indpoly import independence_poly


def mode_index(poly: list[int]) -> int:
    return max(range(len(poly)), key=lambda i: poly[i])


def mean_at_lambda(poly: list[int], lam: float) -> float:
    # Small degrees here, but compute stably via recursive powers.
    z = 0.0
    m = 0.0
    p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p
        z += w
        m += k * w
        p *= lam
    return m / z if z > 0.0 else float("nan")


def tie_mode_set(poly: list[int], lam: float, tol: float) -> list[int]:
    # Determine maximizers of i_k * lam^k with tolerance.
    vals: list[float] = []
    p = 1.0
    for ck in poly:
        vals.append(ck * p)
        p *= lam
    vmax = max(vals)
    return [k for k, v in enumerate(vals) if vmax - v <= tol * max(1.0, abs(vmax))]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Tie-point mean scan for d_leaf<=1 trees."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=23)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    total_seen = 0
    total_considered = 0
    tie_points = 0
    jump_like_points = 0

    lower_fail = 0
    upper_fail = 0
    mode_ceil_fail = 0
    darroch_fail = 0  # mode not in {floor(mu), ceil(mu)} at lambda=1
    mode_tie_checked = 0
    mode_tie_lower_fail = 0

    min_lower_margin = math.inf
    min_upper_margin = math.inf
    min_mode_ceil_margin = math.inf
    min_mode_tie_lower_margin = math.inf

    worst_lower: dict[str, Any] | None = None
    worst_upper: dict[str, Any] | None = None
    worst_mode_ceil: dict[str, Any] | None = None
    first_darroch_fail: dict[str, Any] | None = None
    worst_mode_tie_lower: dict[str, Any] | None = None

    per_n: dict[str, dict[str, Any]] = {}
    t_all = time.time()

    print("Tie-point mean scan", flush=True)
    print("Checks: mu(lambda_k) in [k-1, k], mode<=ceil(mu) at lambda=1", flush=True)
    print("-" * 88, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_ties = 0
        n_jump_like = 0
        n_lower_fail = 0
        n_upper_fail = 0
        n_mode_ceil_fail = 0
        n_darroch_fail = 0
        n_mode_tie_checked = 0
        n_mode_tie_lower_fail = 0

        assert proc.stdout is not None
        for line in proc.stdout:
            n0, adj = parse_graph6(line)
            n_seen += 1
            total_seen += 1
            if not is_dleaf_le_1(n0, adj):
                continue

            n_considered += 1
            total_considered += 1
            g6 = line.decode("ascii").strip()

            poly = independence_poly(n0, adj)
            z1 = float(sum(poly))
            mu1 = sum(i * poly[i] for i in range(len(poly))) / z1
            m1 = mode_index(poly)
            ceil_mu1 = math.ceil(mu1 - 1e-15)
            floor_mu1 = math.floor(mu1 + 1e-15)

            m_margin = ceil_mu1 - m1
            if m_margin < min_mode_ceil_margin:
                min_mode_ceil_margin = m_margin
                worst_mode_ceil = {
                    "n": n0,
                    "g6": g6,
                    "mode": m1,
                    "mu": mu1,
                    "ceil_mu": ceil_mu1,
                    "margin": m_margin,
                }
            if m1 > ceil_mu1:
                mode_ceil_fail += 1
                n_mode_ceil_fail += 1

            if m1 not in (floor_mu1, ceil_mu1):
                darroch_fail += 1
                n_darroch_fail += 1
                if first_darroch_fail is None:
                    first_darroch_fail = {
                        "n": n0,
                        "g6": g6,
                        "mode": m1,
                        "mu": mu1,
                        "floor_mu": floor_mu1,
                        "ceil_mu": ceil_mu1,
                        "poly": poly,
                    }

            # Focused tie at k = mode(lambda=1): candidate bridge near lambda=1.
            if m1 >= 1 and poly[m1 - 1] > 0 and poly[m1] > 0:
                n_mode_tie_checked += 1
                mode_tie_checked += 1
                lam_m = poly[m1 - 1] / poly[m1]
                mu_m = mean_at_lambda(poly, lam_m)
                mode_tie_margin = mu_m - (m1 - 1)
                if mode_tie_margin < min_mode_tie_lower_margin:
                    min_mode_tie_lower_margin = mode_tie_margin
                    worst_mode_tie_lower = {
                        "n": n0,
                        "g6": g6,
                        "mode": m1,
                        "lambda_mode": lam_m,
                        "mu_lambda_mode": mu_m,
                        "margin": mode_tie_margin,
                    }
                if mode_tie_margin < -args.tol:
                    n_mode_tie_lower_fail += 1
                    mode_tie_lower_fail += 1

            for k in range(1, len(poly)):
                if poly[k - 1] <= 0 or poly[k] <= 0:
                    continue
                lam = poly[k - 1] / poly[k]
                mu_lam = mean_at_lambda(poly, lam)
                n_ties += 1
                tie_points += 1

                # "Jump-like" tie: both adjacent indices are (near-)maximizers.
                max_set = tie_mode_set(poly, lam, 1e-10)
                if (k - 1) in max_set and k in max_set:
                    n_jump_like += 1
                    jump_like_points += 1

                lower_margin = mu_lam - (k - 1)
                upper_margin = k - mu_lam
                if lower_margin < min_lower_margin:
                    min_lower_margin = lower_margin
                    worst_lower = {
                        "n": n0,
                        "g6": g6,
                        "k": k,
                        "lambda": lam,
                        "mu_lambda": mu_lam,
                        "margin": lower_margin,
                        "jump_like": ((k - 1) in max_set and k in max_set),
                        "mode_set": max_set,
                    }
                if upper_margin < min_upper_margin:
                    min_upper_margin = upper_margin
                    worst_upper = {
                        "n": n0,
                        "g6": g6,
                        "k": k,
                        "lambda": lam,
                        "mu_lambda": mu_lam,
                        "margin": upper_margin,
                        "jump_like": ((k - 1) in max_set and k in max_set),
                        "mode_set": max_set,
                    }

                if lower_margin < -args.tol:
                    lower_fail += 1
                    n_lower_fail += 1
                if upper_margin < -args.tol:
                    upper_fail += 1
                    n_upper_fail += 1

        proc.wait()

        per_n[str(n)] = {
            "seen": n_seen,
            "considered": n_considered,
            "tie_points": n_ties,
            "jump_like_points": n_jump_like,
            "lower_fail": n_lower_fail,
            "upper_fail": n_upper_fail,
            "mode_ceil_fail": n_mode_ceil_fail,
            "darroch_fail": n_darroch_fail,
            "mode_tie_checked": n_mode_tie_checked,
            "mode_tie_lower_fail": n_mode_tie_lower_fail,
            "wall_s": time.time() - t0,
        }
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} ties={n_ties:9d} "
            f"jump_like={n_jump_like:9d} lower_fail={n_lower_fail:4d} upper_fail={n_upper_fail:4d} "
            f"mode<=ceil fail={n_mode_ceil_fail:4d} darroch_fail={n_darroch_fail:4d} "
            f"mode_tie_fail={n_mode_tie_lower_fail:4d} "
            f"({time.time()-t0:.1f}s)",
            flush=True,
        )

    print("-" * 88, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} ties={tie_points:,} "
        f"jump_like={jump_like_points:,} lower_fail={lower_fail:,} upper_fail={upper_fail:,} "
        f"mode_ceil_fail={mode_ceil_fail:,} darroch_fail={darroch_fail:,} "
        f"mode_tie_checked={mode_tie_checked:,} mode_tie_lower_fail={mode_tie_lower_fail:,} "
        f"wall={time.time()-t_all:.1f}s",
        flush=True,
    )
    print(f"Minimum lower tie margin mu-(k-1): {min_lower_margin}", flush=True)
    print(f"Minimum upper tie margin k-mu: {min_upper_margin}", flush=True)
    print(f"Minimum mode<=ceil(mu) margin ceil(mu)-mode: {min_mode_ceil_margin}", flush=True)
    print(
        f"Minimum focused mode-tie lower margin mu(lambda_mode)-(mode-1): {min_mode_tie_lower_margin}",
        flush=True,
    )
    print(f"Worst lower witness: {worst_lower}", flush=True)
    print(f"Worst upper witness: {worst_upper}", flush=True)
    print(f"Worst mode<=ceil witness: {worst_mode_ceil}", flush=True)
    print(f"Worst focused mode-tie witness: {worst_mode_tie_lower}", flush=True)
    if first_darroch_fail:
        print(f"First Darroch failure: {first_darroch_fail}", flush=True)

    payload = {
        "params": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "tol": args.tol,
        },
        "summary": {
            "seen": total_seen,
            "considered": total_considered,
            "tie_points": tie_points,
            "jump_like_points": jump_like_points,
            "lower_fail": lower_fail,
            "upper_fail": upper_fail,
            "mode_ceil_fail": mode_ceil_fail,
            "darroch_fail": darroch_fail,
            "mode_tie_checked": mode_tie_checked,
            "mode_tie_lower_fail": mode_tie_lower_fail,
            "minimum_lower_tie_margin": min_lower_margin,
            "minimum_upper_tie_margin": min_upper_margin,
            "minimum_mode_ceil_margin": min_mode_ceil_margin,
            "minimum_mode_tie_lower_margin": min_mode_tie_lower_margin,
            "worst_lower": worst_lower,
            "worst_upper": worst_upper,
            "worst_mode_ceil": worst_mode_ceil,
            "worst_mode_tie_lower": worst_mode_tie_lower,
            "first_darroch_fail": first_darroch_fail,
            "wall_s": time.time() - t_all,
        },
        "per_n": per_n,
    }

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
