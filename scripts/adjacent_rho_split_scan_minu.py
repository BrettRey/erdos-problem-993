#!/usr/bin/env python3
"""Exact adjacent split scan for key (m,lambda,rho) under min-u canonical gating.

Search target:
  same (m, lambda, rho), with |N_A - N_B| = 1, where N=[x]P.

Pipeline:
  - d_leaf gate via is_dleaf_le_1
  - canonical admissible triplet via min (u,s,leaf)
  - canonical record from build_record_from_triplet
  - exact rational invariants (num/den integers)

Unlike the projection-battery path, this script scans in one pass and can:
  - early-abort immediately after the first adjacent split,
  - emit checkpoint JSON after each n-layer.
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
    key_from_record,
    rebuild_record,
)


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


def make_payload(
    *,
    args: argparse.Namespace,
    checked_total: int,
    full_unique_keys: int,
    full_collisions: int,
    skipped_dleaf: int,
    skipped_no_triplet: int,
    skipped_bad: int,
    skipped_m: int,
    per_n: list[dict[str, Any]],
    projection_result: dict[str, Any],
    elapsed_scan: float,
    elapsed_total: float,
    done: bool,
) -> dict[str, Any]:
    return {
        "scan": "adjacent_rho_split_scan_minu",
        "triplet_rule": "min_u_then_support_then_leaf",
        "min_n": args.min_n,
        "max_n": args.max_n,
        "m_min": args.m_min,
        "progress_every": args.progress_every,
        "done": bool(done),
        "checked_total": int(checked_total),
        "full_unique_keys": int(full_unique_keys),
        "full_collisions": int(full_collisions),
        "skipped_dleaf": int(skipped_dleaf),
        "skipped_no_triplet": int(skipped_no_triplet),
        "skipped_bad": int(skipped_bad),
        "skipped_m": int(skipped_m),
        "per_n": per_n,
        "projection_result": projection_result,
        "elapsed_sec_scan": float(elapsed_scan),
        "elapsed_sec_total": float(elapsed_total),
    }


def maybe_write_json(path: str, payload: dict[str, Any]) -> None:
    if not path:
        return
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-n", type=int, default=3)
    ap.add_argument("--max-n", type=int, default=27)
    ap.add_argument("--m-min", type=int, default=4)
    ap.add_argument("--progress-every", type=int, default=2_000_000)
    ap.add_argument(
        "--no-rebuild-witness",
        action="store_true",
        help="Do not rebuild full records for split witnesses.",
    )
    ap.add_argument(
        "--no-early-abort",
        action="store_true",
        help="Continue full scan even after first adjacent split is found.",
    )
    ap.add_argument(
        "--checkpoint-path",
        default="",
        help="Optional JSON checkpoint path written after each n-layer.",
    )
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    t0 = time.time()
    checked_total = 0
    skipped_dleaf = 0
    skipped_no_triplet = 0
    skipped_bad = 0
    skipped_m = 0
    full_collisions = 0
    per_n: list[dict[str, Any]] = []

    # full K* key -> first representative info
    full_seen: dict[tuple[int, ...], dict[str, Any]] = {}
    # projection key (m,lambda,rho) -> map N -> first representative
    proj_seen: dict[tuple[int, ...], dict[int, dict[str, Any]]] = {}

    first_split: dict[str, Any] | None = None
    first_adjacent_split: dict[str, Any] | None = None
    early_abort = False

    for n in range(args.min_n, args.max_n + 1):
        total_n = 0
        checked_n = 0
        t_n = time.time()
        for nn, adj, raw in trees_geng_raw(n):
            total_n += 1
            if args.progress_every > 0 and total_n % args.progress_every == 0:
                print(
                    f"progress n={n} total={total_n} checked={checked_n} "
                    f"unique_full={len(full_seen)} full_collisions={full_collisions}",
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
            if rec["m"] < args.m_min:
                skipped_m += 1
                continue

            checked_total += 1
            checked_n += 1

            fk = key_from_record(rec)
            if fk not in full_seen:
                full_seen[fk] = {
                    "n": int(nn),
                    "g6": g6,
                    "N": int(rec["N"]),
                }
            else:
                full_collisions += 1

            pkey = projection_key_from_full(fk)
            n_val = int(rec["N"])
            bucket = proj_seen.get(pkey)
            if bucket is None:
                bucket = {}
                proj_seen[pkey] = bucket

            if bucket:
                # any split on projection key
                if first_split is None and n_val not in bucket:
                    n0 = next(iter(bucket.keys()))
                    a = bucket[n0]
                    b = compact_info({"n": nn, "g6": g6, "N": n_val}, fk)
                    first_split = {
                        "projection_key": unpack_projection_key(pkey),
                        "A": a,
                        "B": b,
                    }
                    if not args.no_rebuild_witness:
                        first_split["rebuilt"] = rebuild_pair(a, b)

                # adjacent split check
                if first_adjacent_split is None:
                    neigh = None
                    if (n_val - 1) in bucket:
                        neigh = bucket[n_val - 1]
                    elif (n_val + 1) in bucket:
                        neigh = bucket[n_val + 1]
                    if neigh is not None:
                        b = compact_info({"n": nn, "g6": g6, "N": n_val}, fk)
                        first_adjacent_split = {
                            "projection_key": unpack_projection_key(pkey),
                            "delta_N": abs(int(neigh["N"]) - int(b["N"])),
                            "A": neigh,
                            "B": b,
                        }
                        if not args.no_rebuild_witness:
                            first_adjacent_split["rebuilt"] = rebuild_pair(neigh, b)
                        if not args.no_early_abort:
                            early_abort = True
                            if n_val not in bucket:
                                bucket[n_val] = b
                            break

            if n_val not in bucket:
                bucket[n_val] = compact_info({"n": nn, "g6": g6, "N": n_val}, fk)

        per_n.append(
            {
                "n": n,
                "total_trees": total_n,
                "checked": checked_n,
                "unique_full_keys": len(full_seen),
                "full_collisions": full_collisions,
                "elapsed_sec": time.time() - t_n,
            }
        )
        print(
            f"n={n:2d} total={total_n:8d} checked={checked_n:7d} "
            f"unique_full={len(full_seen):7d} full_collisions={full_collisions:7d} "
            f"time={time.time()-t_n:.2f}s",
            flush=True,
        )

        projection_result = {
            "name": "m_lambda_rho",
            "fields": ["m", "lambda", "rho"],
            "checked": int(checked_total),
            "unique_keys": len(proj_seen),
            "collisions": int(checked_total) - len(proj_seen),
            "split_found": first_split is not None,
            "adjacent_split_found": first_adjacent_split is not None,
            "first_split": first_split,
            "first_adjacent_split": first_adjacent_split,
            "elapsed_sec": time.time() - t0,
        }

        if args.checkpoint_path:
            cp_payload = make_payload(
                args=args,
                checked_total=checked_total,
                full_unique_keys=len(full_seen),
                full_collisions=full_collisions,
                skipped_dleaf=skipped_dleaf,
                skipped_no_triplet=skipped_no_triplet,
                skipped_bad=skipped_bad,
                skipped_m=skipped_m,
                per_n=per_n,
                projection_result=projection_result,
                elapsed_scan=time.time() - t0,
                elapsed_total=time.time() - t0,
                done=False,
            )
            maybe_write_json(args.checkpoint_path, cp_payload)

        if early_abort:
            break

    projection_result = {
        "name": "m_lambda_rho",
        "fields": ["m", "lambda", "rho"],
        "checked": int(checked_total),
        "unique_keys": len(proj_seen),
        "collisions": int(checked_total) - len(proj_seen),
        "split_found": first_split is not None,
        "adjacent_split_found": first_adjacent_split is not None,
        "first_split": first_split,
        "first_adjacent_split": first_adjacent_split,
        "elapsed_sec": time.time() - t0,
    }

    payload = make_payload(
        args=args,
        checked_total=checked_total,
        full_unique_keys=len(full_seen),
        full_collisions=full_collisions,
        skipped_dleaf=skipped_dleaf,
        skipped_no_triplet=skipped_no_triplet,
        skipped_bad=skipped_bad,
        skipped_m=skipped_m,
        per_n=per_n,
        projection_result=projection_result,
        elapsed_scan=time.time() - t0,
        elapsed_total=time.time() - t0,
        done=True,
    )

    print(
        f"done checked={payload['checked_total']} full_unique={payload['full_unique_keys']} "
        f"full_collisions={payload['full_collisions']} split={projection_result['split_found']} "
        f"adjacent_split={projection_result['adjacent_split_found']} "
        f"elapsed={payload['elapsed_sec_total']:.2f}s",
        flush=True,
    )

    if args.out:
        maybe_write_json(args.out, payload)
        print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
