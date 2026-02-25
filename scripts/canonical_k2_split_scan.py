#!/usr/bin/env python3
"""Scan canonical trees for K2-key collisions that split i1/N.

K2 key:
  (d, m, lambda, mu1, mu2)
where
  d     = deg(P)
  m     = leftmost mode index of I(T)
  lambda= i_{m-1}/i_m
  mu1   = lambda * P'(lambda) / P(lambda)
  mu2   = lambda^2 * P''(lambda) / P(lambda)

Canonical gate:
  - is_dleaf_le_1(...)
  - bridge_decomposition(..., require_dleaf=True) exists

Split condition:
  same K2 key but different i1 (equiv. different N = [x]P = i1 - 3).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction
from typing import Any

# Ensure repository root is on sys.path when run from scripts/.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from attack4_common import bridge_decomposition
from conjecture_a_hall_subset_scan import is_dleaf_le_1
from graph6 import parse_graph6
from trees import trees_geng_raw


def coeff(poly: list[int], k: int) -> int:
    return poly[k] if 0 <= k < len(poly) else 0


def poly_eval_frac(poly: list[int], x: Fraction) -> Fraction:
    out = Fraction(0, 1)
    for c in reversed(poly):
        out = out * x + c
    return out


def deriv1(poly: list[int]) -> list[int]:
    return [i * c for i, c in enumerate(poly)][1:]


def deriv2(poly: list[int]) -> list[int]:
    return [i * (i - 1) * c for i, c in enumerate(poly)][2:]


def frac_pair(x: Fraction) -> list[int]:
    return [x.numerator, x.denominator]


def key_tuple(decomp: Any) -> tuple[int, ...] | None:
    p_poly = decomp.p_poly
    i_poly = decomp.poly_t
    m = decomp.m_t

    i_m = coeff(i_poly, m)
    if i_m == 0:
        return None
    i_m1 = coeff(i_poly, m - 1)
    lam = Fraction(i_m1, i_m)

    p_lam = poly_eval_frac(p_poly, lam)
    if p_lam == 0:
        return None
    mu1 = lam * poly_eval_frac(deriv1(p_poly), lam) / p_lam
    mu2 = lam * lam * poly_eval_frac(deriv2(p_poly), lam) / p_lam

    return (
        len(p_poly) - 1,
        m,
        lam.numerator,
        lam.denominator,
        mu1.numerator,
        mu1.denominator,
        mu2.numerator,
        mu2.denominator,
    )


def build_record(decomp: Any) -> dict[str, Any]:
    p_poly = decomp.p_poly
    q_poly = decomp.q_poly
    i_poly = decomp.poly_t
    m = decomp.m_t

    i_m = coeff(i_poly, m)
    i_m1 = coeff(i_poly, m - 1)
    lam = Fraction(i_m1, i_m)
    p_lam = poly_eval_frac(p_poly, lam)
    mu1 = lam * poly_eval_frac(deriv1(p_poly), lam) / p_lam
    mu2 = lam * lam * poly_eval_frac(deriv2(p_poly), lam) / p_lam

    return {
        "n": decomp.n,
        "g6": decomp.g6,
        "leaf": decomp.leaf,
        "support": decomp.support,
        "u": decomp.u,
        "d": len(p_poly) - 1,
        "m": m,
        "lambda": frac_pair(lam),
        "mu1": frac_pair(mu1),
        "mu2": frac_pair(mu2),
        "i1": coeff(i_poly, 1),
        "N": coeff(p_poly, 1),
        "P": p_poly,
        "Q": q_poly,
        "I": i_poly,
    }


def recompute_record_from_g6(n: int, g6: str) -> dict[str, Any]:
    nn, adj = parse_graph6(g6.encode("ascii"))
    if nn != n:
        raise RuntimeError(f"n mismatch for {g6}: expected {n}, got {nn}")
    decomp = bridge_decomposition(nn, adj, g6, require_dleaf=True)
    if decomp is None:
        raise RuntimeError(f"could not rebuild canonical decomposition for {g6}")
    return build_record(decomp)


def scan(
    min_n: int,
    max_n: int,
    progress_every: int,
    within_n_only: bool,
) -> dict[str, Any]:
    started = time.time()

    # Store minimal data; rebuild full record only for split witnesses.
    # Value: (n, g6, i1, N)
    seen: dict[tuple[int, ...], tuple[int, str, int, int]] = {}

    checked_total = 0
    collisions = 0
    skipped_dleaf = 0
    skipped_decomp = 0
    split: dict[str, Any] | None = None
    per_n: list[dict[str, Any]] = []

    for n in range(min_n, max_n + 1):
        if within_n_only:
            seen.clear()

        n_started = time.time()
        total_n = 0
        checked_n = 0

        for nn, adj, raw in trees_geng_raw(n):
            total_n += 1
            if progress_every > 0 and total_n % progress_every == 0:
                print(
                    f"progress n={n} total={total_n} checked={checked_n} "
                    f"unique={len(seen)} collisions={collisions}",
                    flush=True,
                )

            if not is_dleaf_le_1(nn, adj):
                skipped_dleaf += 1
                continue

            g6 = raw.decode("ascii").strip()
            decomp = bridge_decomposition(nn, adj, g6, require_dleaf=True)
            if decomp is None:
                skipped_decomp += 1
                continue

            key = key_tuple(decomp)
            if key is None:
                continue

            i1 = coeff(decomp.poly_t, 1)
            n_coeff = coeff(decomp.p_poly, 1)

            prev = seen.get(key)
            if prev is None:
                seen[key] = (nn, g6, i1, n_coeff)
            else:
                collisions += 1
                prev_n, prev_g6, prev_i1, prev_n_coeff = prev
                if prev_i1 != i1 or prev_n_coeff != n_coeff:
                    split = {
                        "A": recompute_record_from_g6(prev_n, prev_g6),
                        "B": build_record(decomp),
                    }
                    break

            checked_total += 1
            checked_n += 1

        elapsed_n = time.time() - n_started
        per_n.append(
            {
                "n": n,
                "total_trees": total_n,
                "checked": checked_n,
                "unique_keys": len(seen),
                "collisions": collisions,
                "elapsed_sec": elapsed_n,
            }
        )
        print(
            f"n={n:2d} total={total_n:8d} checked={checked_n:7d} "
            f"unique={len(seen):7d} collisions={collisions:7d} "
            f"time={elapsed_n:.2f}s",
            flush=True,
        )

        if split is not None:
            break

    return {
        "scan": "canonical_K2_split_search_exact",
        "key": "(d,m,lambda,mu1,mu2)",
        "within_n_only": within_n_only,
        "min_n": min_n,
        "max_n": per_n[-1]["n"] if per_n else min_n,
        "checked_total": checked_total,
        "unique_keys": len(seen),
        "collisions": collisions,
        "skipped_dleaf": skipped_dleaf,
        "skipped_decomp": skipped_decomp,
        "split_found": split is not None,
        "split": split,
        "per_n": per_n,
        "elapsed_sec": time.time() - started,
    }


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Exact canonical split scan for K2=(d,m,lambda,mu1,mu2)."
    )
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=22)
    ap.add_argument(
        "--within-n-only",
        action="store_true",
        help="Reset keys at each n (intra-layer collisions only).",
    )
    ap.add_argument(
        "--progress-every",
        type=int,
        default=0,
        help="Emit progress every N raw trees within each n (0 disables).",
    )
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    payload = scan(
        min_n=args.min_n,
        max_n=args.max_n,
        progress_every=args.progress_every,
        within_n_only=args.within_n_only,
    )

    print("---", flush=True)
    print(
        f"done checked={payload['checked_total']} unique={payload['unique_keys']} "
        f"collisions={payload['collisions']} split_found={payload['split_found']} "
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
