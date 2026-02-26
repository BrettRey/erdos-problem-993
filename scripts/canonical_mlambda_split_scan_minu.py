#!/usr/bin/env python3
"""Exhaustive canonical split scan for the projected key (m, lambda).

Uses the same explicit canonical triplet rule as canonical_kstar_split_scan_minu.py:
  choose admissible (leaf,support,u) by minimum (u, support, leaf).

Key tested:
  (m, lambda) where
    m = leftmost mode index of I(T)
    lambda = i_{m-1}/i_m

Split condition:
  same (m,lambda) but different N=[x]P (equiv. different i1-3).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Any

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from trees import trees_geng_raw

from canonical_kstar_split_scan_minu import (
    build_record_from_triplet,
    choose_triplet_min_u,
)


def key_from_record(rec: dict[str, Any]) -> tuple[int, int, int]:
    lam_n, lam_d = rec["lambda"]
    return (int(rec["m"]), int(lam_n), int(lam_d))


def scan(
    min_n: int,
    max_n: int,
    m_min: int,
    progress_every: int,
    within_n_only: bool,
) -> dict[str, Any]:
    started = time.time()

    # key -> (n, g6, N)
    seen: dict[tuple[int, int, int], tuple[int, str, int]] = {}
    checked_total = 0
    skipped_dleaf = 0
    skipped_no_triplet = 0
    skipped_bad = 0
    skipped_m = 0
    collisions = 0
    split: dict[str, Any] | None = None
    per_n: list[dict[str, Any]] = []

    for n in range(min_n, max_n + 1):
        if within_n_only:
            seen.clear()

        total_n = 0
        checked_n = 0
        t0 = time.time()

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

            trip = choose_triplet_min_u(adj)
            if trip is None:
                skipped_no_triplet += 1
                continue

            g6 = raw.decode("ascii").strip()
            rec = build_record_from_triplet(nn, adj, g6, trip)
            if rec is None:
                skipped_bad += 1
                continue
            if rec["m"] < m_min:
                skipped_m += 1
                continue

            key = key_from_record(rec)
            prev = seen.get(key)
            if prev is None:
                seen[key] = (nn, g6, int(rec["N"]))
            else:
                collisions += 1
                pn, pg6, pN = prev
                if pN != int(rec["N"]):
                    split = {
                        "A": {"n": pn, "g6": pg6, "N": pN, "m_lambda": list(key)},
                        "B": {
                            "n": nn,
                            "g6": g6,
                            "N": int(rec["N"]),
                            "m_lambda": list(key),
                        },
                    }
                    break

            checked_total += 1
            checked_n += 1

        per_n.append(
            {
                "n": n,
                "total_trees": total_n,
                "checked": checked_n,
                "unique_keys": len(seen),
                "collisions": collisions,
                "elapsed_sec": time.time() - t0,
            }
        )
        print(
            f"n={n:2d} total={total_n:8d} checked={checked_n:7d} "
            f"unique={len(seen):7d} collisions={collisions:7d} "
            f"time={time.time()-t0:.2f}s",
            flush=True,
        )
        if split is not None:
            break

    return {
        "scan": "canonical_mlambda_split_scan_minu",
        "triplet_rule": "min_u_then_support_then_leaf",
        "key": "(m,lambda)",
        "min_n": min_n,
        "max_n": per_n[-1]["n"] if per_n else min_n,
        "m_min": m_min,
        "within_n_only": within_n_only,
        "checked_total": checked_total,
        "unique_keys": len(seen),
        "collisions": collisions,
        "skipped_dleaf": skipped_dleaf,
        "skipped_no_triplet": skipped_no_triplet,
        "skipped_bad": skipped_bad,
        "skipped_m": skipped_m,
        "split_found": split is not None,
        "split": split,
        "per_n": per_n,
        "elapsed_sec": time.time() - started,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-n", type=int, default=3)
    ap.add_argument("--max-n", type=int, default=24)
    ap.add_argument("--m-min", type=int, default=4)
    ap.add_argument("--progress-every", type=int, default=0)
    ap.add_argument(
        "--within-n-only",
        action="store_true",
        help="Reset key table per n (intra-layer collisions only).",
    )
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    payload = scan(
        min_n=args.min_n,
        max_n=args.max_n,
        m_min=args.m_min,
        progress_every=args.progress_every,
        within_n_only=args.within_n_only,
    )
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

