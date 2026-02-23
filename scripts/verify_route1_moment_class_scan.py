#!/usr/bin/env python3
"""Verify moment-class LP-dual bound against Route-1 tail threshold.

This script computes, for each canonical degree-2 bridge decomposition:
  - Threshold = (m-2) - (p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m)/P(lambda)
  - B_max(moment class) via subset-DP + one-step LP dual lower bounds.

The one-step LP dual bound uses factorial moments up to order `s` (default 4):
  q(i) = c0 + c1*i + c2*i(i-1) + ... + cs*i(i-1)...(i-s+1) <= beta_i
  maximize q-moment objective under these constraints.

If B_max < Threshold for any checked tree, the class fails on that tree.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from functools import lru_cache
from typing import Any

import numpy as np
from scipy.optimize import linprog

# Ensure repository root is on sys.path when run from scripts/.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from attack4_common import bridge_decomposition
from conjecture_a_hall_subset_scan import is_dleaf_le_1
from indpoly import _polymul
from trees import trees_geng_raw


def poly_eval(poly: list[int], lam: float) -> float:
    return sum(c * (lam ** i) for i, c in enumerate(poly))


def coeff(poly: list[int], k: int) -> int:
    return poly[k] if 0 <= k < len(poly) else 0


def threshold(decomp: Any, lam: float) -> float:
    m = decomp.m_t
    p_poly = decomp.p_poly
    p_lam = poly_eval(p_poly, lam)
    p_m1 = coeff(p_poly, m - 1)
    p_m = coeff(p_poly, m)
    return (m - 2) - (p_m1 * (lam ** (m - 1)) + 2 * p_m * (lam ** m)) / p_lam


def factorial_basis(i: int, s: int) -> list[float]:
    vals = [1.0]
    cur = 1.0
    for t in range(1, s + 1):
        cur *= (i - (t - 1))
        vals.append(cur)
    return vals


def step_beta(f_poly: list[int], m: int, lam: float, deg_a: int) -> np.ndarray:
    f_lam = poly_eval(f_poly, lam)
    beta = np.zeros(deg_a + 1, dtype=float)
    if m <= 2:
        return beta

    # alpha_r(F) = tail_{>=r+1}(F at lambda) / F(lambda), r=0..m-3
    alphas: list[float] = []
    for r in range(0, m - 2):
        tail = 0.0
        for j in range(r + 1, len(f_poly)):
            tail += f_poly[j] * (lam ** j)
        alphas.append(tail / f_lam)

    for i in range(deg_a + 1):
        s_val = 0.0
        for r, alpha in enumerate(alphas):
            if i <= (m - 3 - r):
                s_val += alpha
        beta[i] = s_val

    return beta


def factorial_moments(a_poly: list[int], lam: float, s: int) -> np.ndarray:
    a_lam = poly_eval(a_poly, lam)
    out = [1.0]
    for t in range(1, s + 1):
        num = 0.0
        for i, ai in enumerate(a_poly):
            if ai == 0:
                continue
            ff = 1.0
            for q in range(t):
                ff *= (i - q)
            num += ai * (lam ** i) * ff
        out.append(num / a_lam)
    return np.array(out, dtype=float)


def lp_step_lower_bound(
    a_poly: list[int], f_poly: list[int], m: int, lam: float, s: int
) -> float:
    deg_a = len(a_poly) - 1
    beta = step_beta(f_poly, m, lam, deg_a)
    mat = np.array([factorial_basis(i, s) for i in range(deg_a + 1)], dtype=float)
    moms = factorial_moments(a_poly, lam, s)

    # maximize c·moms subject to mat*c <= beta
    # linprog minimizes; so minimize -c·moms
    res = linprog(
        c=-moms,
        A_ub=mat,
        b_ub=beta,
        bounds=[(None, None)] * (s + 1),
        method="highs",
    )
    if res.status != 0:
        raise RuntimeError(f"LP failure: status={res.status} msg={res.message}")
    return float(-res.fun)


def bmax_for_decomp(decomp: Any, lam: float, s: int) -> float:
    factors = [list(f) for f in decomp.f_list]
    d = len(factors)

    @lru_cache(maxsize=None)
    def poly_of(mask: int) -> tuple[int, ...]:
        out = [1]
        for j in range(d):
            if mask & (1 << j):
                out = _polymul(out, factors[j])
        return tuple(out)

    @lru_cache(maxsize=None)
    def best(mask: int) -> float:
        if mask == (1 << d) - 1:
            return 0.0
        a_poly = list(poly_of(mask))
        val = -1e100
        for j in range(d):
            if mask & (1 << j):
                continue
            lb = lp_step_lower_bound(a_poly, factors[j], decomp.m_t, lam, s)
            cand = lb + best(mask | (1 << j))
            if cand > val:
                val = cand
        return val

    return best(0)


def scan_range(min_n: int, max_n: int, s: int, progress_every: int) -> dict[str, Any]:
    started = time.time()
    overall_checked = 0
    failures: list[dict[str, Any]] = []
    per_n: list[dict[str, Any]] = []

    for n in range(min_n, max_n + 1):
        t0 = time.time()
        total_n = 0
        checked_n = 0
        min_gap_n: float | None = None
        min_rec_n: dict[str, Any] | None = None

        for nn, adj, raw in trees_geng_raw(n):
            total_n += 1
            if progress_every > 0 and total_n % progress_every == 0:
                print(
                    f"progress n={n} total={total_n} checked={checked_n} "
                    f"elapsed={time.time() - t0:.2f}s",
                    flush=True,
                )

            if not is_dleaf_le_1(nn, adj):
                continue

            g6 = raw.decode("ascii").strip()
            decomp = bridge_decomposition(nn, adj, g6, require_dleaf=True)
            if decomp is None:
                continue

            m = decomp.m_t
            i_m1 = coeff(decomp.poly_t, m - 1)
            i_m = coeff(decomp.poly_t, m)
            lam = i_m1 / i_m
            th = threshold(decomp, lam)
            bmax = bmax_for_decomp(decomp, lam, s=s)
            gap = bmax - th

            checked_n += 1
            overall_checked += 1

            if min_gap_n is None or gap < min_gap_n:
                min_gap_n = gap
                min_rec_n = {
                    "g6": g6,
                    "m": m,
                    "lambda": lam,
                    "threshold": th,
                    "b_max": bmax,
                    "gap": gap,
                }

            if gap < -1e-9:
                fail = {
                    "n": nn,
                    "g6": g6,
                    "m": m,
                    "lambda": lam,
                    "threshold": th,
                    "b_max": bmax,
                    "gap": gap,
                }
                failures.append(fail)
                print(
                    f"FAIL n={nn} g6={g6} gap={gap:.12g} "
                    f"threshold={th:.12g} b_max={bmax:.12g}",
                    flush=True,
                )
                break

        rec = {
            "n": n,
            "total_trees": total_n,
            "checked": checked_n,
            "elapsed_sec": time.time() - t0,
            "min_gap": min_gap_n,
            "min_witness": min_rec_n,
        }
        per_n.append(rec)

        if checked_n == 0:
            print(
                f"n={n:2d} total={total_n:7d} checked={checked_n:6d} "
                f"time={rec['elapsed_sec']:.2f}s",
                flush=True,
            )
        else:
            assert min_rec_n is not None
            print(
                f"n={n:2d} total={total_n:7d} checked={checked_n:6d} "
                f"min_gap={min_gap_n:.9f} min_g6={min_rec_n['g6']} "
                f"time={rec['elapsed_sec']:.2f}s",
                flush=True,
            )

        if failures:
            break

    out = {
        "class": f"C_deg_mu_to_mu{s}",
        "moment_order": s,
        "min_n": min_n,
        "max_n": max_n,
        "checked_total": overall_checked,
        "elapsed_sec": time.time() - started,
        "failures": failures,
        "per_n": per_n,
    }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Scan Route-1 moment class LP-dual bound.")
    ap.add_argument("--min-n", type=int, default=1)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--moment-order", type=int, default=4)
    ap.add_argument("--progress-every", type=int, default=0)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    if args.moment_order < 0:
        raise ValueError("moment-order must be nonnegative")

    payload = scan_range(
        min_n=args.min_n,
        max_n=args.max_n,
        s=args.moment_order,
        progress_every=args.progress_every,
    )

    print(
        f"done checked={payload['checked_total']} failures={len(payload['failures'])} "
        f"elapsed={payload['elapsed_sec']:.2f}s",
        flush=True,
    )

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
