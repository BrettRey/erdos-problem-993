#!/usr/bin/env python3
"""Find canonical-valid pairs with same P but different derived lambda.

Canonical gate:
  - is_dleaf_le_1(...)
  - bridge_decomposition(..., require_dleaf=True) exists

Pair condition:
  - same P = dp0[u]
  - different canonical derived lambda = i_{m-1}/i_m from I(T)
  - optionally require same m (default True)
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
from trees import trees_geng_raw


def coeff(poly: list[int], k: int) -> int:
    return poly[k] if 0 <= k < len(poly) else 0


def frac_pair(x: Fraction) -> list[int]:
    return [x.numerator, x.denominator]


def record_from_decomp(decomp: Any) -> dict[str, Any] | None:
    i_poly = decomp.poly_t
    m = decomp.m_t
    i_m = coeff(i_poly, m)
    if i_m == 0:
        return None
    i_m1 = coeff(i_poly, m - 1)
    lam = Fraction(i_m1, i_m)

    return {
        "n": decomp.n,
        "g6": decomp.g6,
        "leaf": decomp.leaf,
        "support": decomp.support,
        "u": decomp.u,
        "P": decomp.p_poly,
        "Q": decomp.q_poly,
        "I": i_poly,
        "m": m,
        "lambda": frac_pair(lam),
    }


def scan(
    min_n: int,
    max_n: int,
    require_same_m: bool,
    progress_every: int,
) -> dict[str, Any]:
    started = time.time()

    # key: tuple(P), value: first record seen for this P
    seen_by_p: dict[tuple[int, ...], dict[str, Any]] = {}
    checked_total = 0
    collisions_same_p = 0
    skipped_dleaf = 0
    skipped_decomp = 0
    skipped_zero_im = 0
    witness: dict[str, Any] | None = None
    per_n: list[dict[str, Any]] = []

    for n in range(min_n, max_n + 1):
        n_started = time.time()
        total_n = 0
        checked_n = 0

        for nn, adj, raw in trees_geng_raw(n):
            total_n += 1
            if progress_every > 0 and total_n % progress_every == 0:
                print(
                    f"progress n={n} total={total_n} checked={checked_n} "
                    f"sameP_collisions={collisions_same_p}",
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

            rec = record_from_decomp(decomp)
            if rec is None:
                skipped_zero_im += 1
                continue

            p_key = tuple(rec["P"])
            prev = seen_by_p.get(p_key)
            if prev is None:
                seen_by_p[p_key] = rec
            else:
                collisions_same_p += 1
                same_m = prev["m"] == rec["m"]
                same_lam = prev["lambda"] == rec["lambda"]
                if (not same_lam) and ((not require_same_m) or same_m):
                    witness = {"A": prev, "B": rec}
                    break

            checked_total += 1
            checked_n += 1

        elapsed_n = time.time() - n_started
        per_n.append(
            {
                "n": n,
                "total_trees": total_n,
                "checked": checked_n,
                "unique_P": len(seen_by_p),
                "sameP_collisions": collisions_same_p,
                "elapsed_sec": elapsed_n,
            }
        )
        print(
            f"n={n:2d} total={total_n:8d} checked={checked_n:7d} "
            f"unique_P={len(seen_by_p):7d} sameP_collisions={collisions_same_p:7d} "
            f"time={elapsed_n:.2f}s",
            flush=True,
        )

        if witness is not None:
            break

    return {
        "scan": "canonical_same_P_diff_lambda",
        "min_n": min_n,
        "max_n": per_n[-1]["n"] if per_n else min_n,
        "require_same_m": require_same_m,
        "checked_total": checked_total,
        "unique_P": len(seen_by_p),
        "sameP_collisions": collisions_same_p,
        "skipped_dleaf": skipped_dleaf,
        "skipped_decomp": skipped_decomp,
        "skipped_zero_im": skipped_zero_im,
        "witness_found": witness is not None,
        "witness": witness,
        "per_n": per_n,
        "elapsed_sec": time.time() - started,
    }


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Find canonical pair with same P but different canonical lambda."
    )
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=22)
    ap.add_argument(
        "--allow-different-m",
        action="store_true",
        help="Allow witness with different m (default requires same m).",
    )
    ap.add_argument("--progress-every", type=int, default=0)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    payload = scan(
        min_n=args.min_n,
        max_n=args.max_n,
        require_same_m=not args.allow_different_m,
        progress_every=args.progress_every,
    )

    print("---", flush=True)
    print(
        f"done checked={payload['checked_total']} unique_P={payload['unique_P']} "
        f"sameP_collisions={payload['sameP_collisions']} "
        f"witness_found={payload['witness_found']} elapsed={payload['elapsed_sec']:.2f}s",
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
