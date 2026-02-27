#!/usr/bin/env python3
"""Exact adjacent split scan for key (m,lambda,rho) under min-u canonical gating.

Searches for any pair of gated canonical trees with:
  same (m, lambda, rho), but |N_A - N_B| = 1,
where N = [x]P from the canonical triplet decomposition.

Pipeline and gate match canonical_projection_battery_minu:
  - is_dleaf_le_1
  - admissible triplet exists
  - canonical triplet by min (u,s,leaf)
  - optional m >= m_min filter

Outputs totals and (if found) full rebuilt witness records including
canonical triplets and P/Q/I coefficient lists.
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

from scripts.canonical_projection_battery_minu import scan_full_keys
from scripts.canonical_kstar_split_scan_minu import rebuild_record


# full key layout from canonical_kstar_split_scan_minu.key_from_record:
# (d,m,lam_n,lam_d,mu1_n,mu1_d,mu2_n,mu2_d,rho_n,rho_d,sig_n,sig_d)


def projection_key_from_full(full_key: tuple[int, ...]) -> tuple[int, ...]:
    return (
        int(full_key[1]),  # m
        int(full_key[2]),
        int(full_key[3]),  # lambda
        int(full_key[8]),
        int(full_key[9]),  # rho
    )


def unpack_projection_key(pkey: tuple[int, ...]) -> dict[str, Any]:
    return {
        "m": int(pkey[0]),
        "lambda": [int(pkey[1]), int(pkey[2])],
        "rho": [int(pkey[3]), int(pkey[4])],
    }


def unpack_full_key(fk: tuple[int, ...]) -> dict[str, Any]:
    return {
        "d": int(fk[0]),
        "m": int(fk[1]),
        "lambda": [int(fk[2]), int(fk[3])],
        "mu1": [int(fk[4]), int(fk[5])],
        "mu2": [int(fk[6]), int(fk[7])],
        "rho": [int(fk[8]), int(fk[9])],
        "sigma": [int(fk[10]), int(fk[11])],
    }


def compact_info(info: dict[str, Any], fk: tuple[int, ...]) -> dict[str, Any]:
    return {
        "n": int(info["n"]),
        "g6": str(info["g6"]),
        "N": int(info["N"]),
        "full_key": list(fk),
        "full_invariants": unpack_full_key(fk),
    }


def rebuild_pair(a: dict[str, Any], b: dict[str, Any]) -> dict[str, Any]:
    rec_a = rebuild_record(int(a["n"]), str(a["g6"]))
    rec_b = rebuild_record(int(b["n"]), str(b["g6"]))
    return {"A": rec_a, "B": rec_b}


def run_projection_adjacent_scan(
    full_seen: dict[tuple[int, ...], dict[str, Any]],
    checked_total: int,
    rebuild_witness: bool,
) -> dict[str, Any]:
    started = time.time()

    # pkey -> map N -> representative full key
    buckets: dict[tuple[int, ...], dict[int, tuple[int, ...]]] = {}

    first_split: dict[str, Any] | None = None
    first_adjacent_split: dict[str, Any] | None = None

    for fk, info in full_seen.items():
        pkey = projection_key_from_full(fk)
        n_val = int(info["N"])

        bucket = buckets.get(pkey)
        if bucket is None:
            bucket = {}
            buckets[pkey] = bucket
        else:
            # any split on (m,lambda,rho)
            if first_split is None and n_val not in bucket:
                n0 = next(iter(bucket.keys()))
                fk0 = bucket[n0]
                info0 = full_seen[fk0]
                a = compact_info(info0, fk0)
                b = compact_info(info, fk)
                first_split = {
                    "projection_key": unpack_projection_key(pkey),
                    "A": a,
                    "B": b,
                }
                if rebuild_witness:
                    first_split["rebuilt"] = rebuild_pair(a, b)

            # adjacent split detection
            if first_adjacent_split is None:
                neigh_fk = None
                neigh_n = None
                if (n_val - 1) in bucket:
                    neigh_n = n_val - 1
                    neigh_fk = bucket[neigh_n]
                elif (n_val + 1) in bucket:
                    neigh_n = n_val + 1
                    neigh_fk = bucket[neigh_n]

                if neigh_fk is not None and neigh_n is not None:
                    info0 = full_seen[neigh_fk]
                    a = compact_info(info0, neigh_fk)
                    b = compact_info(info, fk)
                    first_adjacent_split = {
                        "projection_key": unpack_projection_key(pkey),
                        "delta_N": abs(int(a["N"]) - int(b["N"])),
                        "A": a,
                        "B": b,
                    }
                    if rebuild_witness:
                        first_adjacent_split["rebuilt"] = rebuild_pair(a, b)

        if n_val not in bucket:
            bucket[n_val] = fk

    unique_keys = len(buckets)
    collisions = int(checked_total) - unique_keys

    return {
        "name": "m_lambda_rho",
        "fields": ["m", "lambda", "rho"],
        "checked": int(checked_total),
        "unique_keys": unique_keys,
        "collisions": collisions,
        "split_found": first_split is not None,
        "adjacent_split_found": first_adjacent_split is not None,
        "first_split": first_split,
        "first_adjacent_split": first_adjacent_split,
        "elapsed_sec": time.time() - started,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-n", type=int, default=3)
    ap.add_argument("--max-n", type=int, default=27)
    ap.add_argument("--m-min", type=int, default=4)
    ap.add_argument("--progress-every", type=int, default=0)
    ap.add_argument(
        "--no-rebuild-witness",
        action="store_true",
        help="Do not rebuild full records for split witnesses.",
    )
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    t0 = time.time()
    scan_payload = scan_full_keys(
        min_n=args.min_n,
        max_n=args.max_n,
        m_min=args.m_min,
        progress_every=args.progress_every,
    )
    full_seen = scan_payload.pop("full_seen")

    projection_result = run_projection_adjacent_scan(
        full_seen=full_seen,
        checked_total=int(scan_payload["checked_total"]),
        rebuild_witness=not bool(args.no_rebuild_witness),
    )

    payload = {
        "scan": "adjacent_rho_split_scan_minu",
        "triplet_rule": "min_u_then_support_then_leaf",
        "min_n": args.min_n,
        "max_n": args.max_n,
        "m_min": args.m_min,
        **scan_payload,
        "projection_result": projection_result,
        "elapsed_sec_total": time.time() - t0,
    }

    print(
        f"done checked={payload['checked_total']} full_unique={payload['full_unique_keys']} "
        f"full_collisions={payload['full_collisions']} split={projection_result['split_found']} "
        f"adjacent_split={projection_result['adjacent_split_found']} "
        f"elapsed={payload['elapsed_sec_total']:.2f}s",
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
