#!/usr/bin/env python3
"""Profile same-(m,lambda) overlaps under min-u canonical gating.

This ignores rho and groups accepted trees by key (m,lambda). For keys that have
multiple N values, it reports overlap structure and, in particular, checks the
odd/even adjacent condition:
  if N_even = N_odd + 1, does d_even == N_even/2 ?
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
from scripts.canonical_kstar_split_scan_minu import (
    build_record_from_triplet,
    choose_triplet_min_u,
)


def mlambda_key(rec: dict[str, Any]) -> tuple[int, int, int]:
    lam = rec["lambda"]
    if isinstance(lam, dict):
        return (int(rec["m"]), int(lam["num"]), int(lam["den"]))
    return (int(rec["m"]), int(lam[0]), int(lam[1]))


def compact(rec: dict[str, Any], n: int, g6: str) -> dict[str, Any]:
    lam = rec["lambda"]
    if isinstance(lam, dict):
        lam_num = int(lam["num"])
        lam_den = int(lam["den"])
    else:
        lam_num = int(lam[0])
        lam_den = int(lam[1])
    return {
        "n": int(n),
        "g6": g6,
        "N": int(rec["N"]),
        "d": int(rec["d"]),
        "m": int(rec["m"]),
        "lambda": {"num": lam_num, "den": lam_den},
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-n", type=int, default=3)
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--m-min", type=int, default=4)
    ap.add_argument("--max-examples", type=int, default=10)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    t0 = time.time()

    checked = 0
    skipped_dleaf = 0
    skipped_no_triplet = 0
    skipped_bad = 0
    skipped_m = 0

    # (m,lam_num,lam_den) -> N -> representative
    seen: dict[tuple[int, int, int], dict[int, dict[str, Any]]] = {}

    per_n: list[dict[str, Any]] = []

    for n in range(args.min_n, args.max_n + 1):
        total_n = 0
        checked_n = 0
        t_n = time.time()
        for nn, adj, raw in trees_geng_raw(n):
            total_n += 1
            if not is_dleaf_le_1(nn, adj):
                skipped_dleaf += 1
                continue
            trip = choose_triplet_min_u(adj)
            if trip is None:
                skipped_no_triplet += 1
                continue
            rec = build_record_from_triplet(nn, adj, raw.decode("ascii").strip(), trip)
            if rec is None:
                skipped_bad += 1
                continue
            if int(rec["m"]) < args.m_min:
                skipped_m += 1
                continue

            checked += 1
            checked_n += 1

            key = mlambda_key(rec)
            bucket = seen.get(key)
            if bucket is None:
                bucket = {}
                seen[key] = bucket

            n_val = int(rec["N"])
            if n_val not in bucket:
                bucket[n_val] = compact(rec, nn, raw.decode("ascii").strip())

        per_n.append(
            {
                "n": n,
                "total_trees": total_n,
                "checked": checked_n,
                "elapsed_sec": time.time() - t_n,
            }
        )
        print(
            f"n={n:2d} total={total_n:8d} checked={checked_n:7d} keys={len(seen):7d} time={time.time()-t_n:.2f}s",
            flush=True,
        )

    multi_n_keys = 0
    adjacent_pairs = 0
    adjacent_pairs_even_deg_half = 0
    adjacent_examples: list[dict[str, Any]] = []

    for key, bucket in seen.items():
        Ns = sorted(bucket.keys())
        if len(Ns) > 1:
            multi_n_keys += 1
        for i in range(len(Ns) - 1):
            aN = Ns[i]
            bN = Ns[i + 1]
            if bN != aN + 1:
                continue
            adjacent_pairs += 1
            evenN = bN if (bN % 2 == 0) else aN
            oddN = aN if (aN % 2 == 1) else bN
            even_rec = bucket[evenN]
            cond = (int(even_rec["d"]) * 2 == int(evenN))
            if cond:
                adjacent_pairs_even_deg_half += 1
            if len(adjacent_examples) < args.max_examples:
                adjacent_examples.append(
                    {
                        "key": {
                            "m": int(key[0]),
                            "lambda": {"num": int(key[1]), "den": int(key[2])},
                        },
                        "odd": bucket[oddN],
                        "even": bucket[evenN],
                        "condition_even_d_eq_halfN": cond,
                    }
                )

    payload = {
        "scan": "mlambda_overlap_profile_minu",
        "triplet_rule": "min_u_then_support_then_leaf",
        "min_n": args.min_n,
        "max_n": args.max_n,
        "m_min": args.m_min,
        "checked_total": checked,
        "keys_total": len(seen),
        "keys_with_multiple_N": multi_n_keys,
        "adjacent_pairs_total": adjacent_pairs,
        "adjacent_pairs_even_d_eq_halfN": adjacent_pairs_even_deg_half,
        "adjacent_pairs_even_d_eq_halfN_ratio": (
            (adjacent_pairs_even_deg_half / adjacent_pairs) if adjacent_pairs else None
        ),
        "examples": adjacent_examples,
        "skipped_dleaf": skipped_dleaf,
        "skipped_no_triplet": skipped_no_triplet,
        "skipped_bad": skipped_bad,
        "skipped_m": skipped_m,
        "per_n": per_n,
        "elapsed_sec_total": time.time() - t0,
    }

    print(
        "done "
        f"checked={checked} keys={len(seen)} multiN={multi_n_keys} "
        f"adjacent_pairs={adjacent_pairs} "
        f"adjacent_even_deg_halfN={adjacent_pairs_even_deg_half}",
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
