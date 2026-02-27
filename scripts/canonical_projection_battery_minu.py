#!/usr/bin/env python3
"""Single-pass canonical projection battery scan under min-u tie-break.

Workflow:
1) Enumerate gated canonical trees once (same pipeline as canonical_kstar_split_scan_minu.py),
   storing first occurrence of each full K* key.
2) Evaluate multiple projected keys against N using only the unique full-key table.

This avoids re-running the expensive tree enumeration per projection.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction
from math import comb
from typing import Any

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from trees import trees_geng_raw

from canonical_kstar_split_scan_minu import (
    build_record_from_triplet,
    choose_triplet_min_u,
    key_from_record,
)

VALID_FIELDS = ("d", "m", "lambda", "mu1", "mu2", "rho", "sigma")
DEFAULT_PROJECTIONS = (
    "m,lambda,d",
    "m,lambda,mu1",
    "m,lambda,rho",
    "m,lambda,sigma",
    "m,lambda,rho,sigma",
    "d,m,lambda,mu1,mu2",
    "d,m,lambda,rho,sigma",
)

# full key layout from canonical_kstar_split_scan_minu.key_from_record:
# (d,m,lam_n,lam_d,mu1_n,mu1_d,mu2_n,mu2_d,rho_n,rho_d,sig_n,sig_d)
FIELD_POS = {
    "d": (0,),
    "m": (1,),
    "lambda": (2, 3),
    "mu1": (4, 5),
    "mu2": (6, 7),
    "rho": (8, 9),
    "sigma": (10, 11),
}


def parse_projection_list(raw: str) -> list[tuple[str, tuple[str, ...]]]:
    out: list[tuple[str, tuple[str, ...]]] = []
    for block in raw.split(";"):
        block = block.strip()
        if not block:
            continue
        fields: list[str] = []
        for f in (x.strip() for x in block.split(",")):
            if not f:
                continue
            if f not in VALID_FIELDS:
                raise ValueError(f"invalid field in projection: {f}")
            if f not in fields:
                fields.append(f)
        if not fields:
            continue
        name = "_".join(fields)
        out.append((name, tuple(fields)))
    if not out:
        raise ValueError("no projections parsed")
    return out


def project_full_key(full_key: tuple[int, ...], fields: tuple[str, ...]) -> tuple[int, ...]:
    out: list[int] = []
    for f in fields:
        for idx in FIELD_POS[f]:
            out.append(int(full_key[idx]))
    return tuple(out)


def unpack_full_key(full_key: tuple[int, ...]) -> dict[str, Any]:
    return {
        "d": full_key[0],
        "m": full_key[1],
        "lambda": [full_key[2], full_key[3]],
        "mu1": [full_key[4], full_key[5]],
        "mu2": [full_key[6], full_key[7]],
        "rho": [full_key[8], full_key[9]],
        "sigma": [full_key[10], full_key[11]],
    }


def scan_full_keys(
    min_n: int,
    max_n: int,
    m_min: int,
    progress_every: int,
) -> dict[str, Any]:
    started = time.time()
    checked_total = 0
    skipped_dleaf = 0
    skipped_no_triplet = 0
    skipped_bad = 0
    skipped_m = 0
    full_collisions = 0
    per_n: list[dict[str, Any]] = []

    # full K* key -> first representative info
    full_seen: dict[tuple[int, ...], dict[str, Any]] = {}

    for n in range(min_n, max_n + 1):
        total_n = 0
        checked_n = 0
        t0 = time.time()
        for nn, adj, raw in trees_geng_raw(n):
            total_n += 1
            if progress_every > 0 and total_n % progress_every == 0:
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
            if rec["m"] < m_min:
                skipped_m += 1
                continue

            full_key = key_from_record(rec)
            if full_key not in full_seen:
                full_seen[full_key] = {
                    "n": int(nn),
                    "g6": g6,
                    "N": int(rec["N"]),
                }
            else:
                full_collisions += 1

            checked_total += 1
            checked_n += 1

        per_n.append(
            {
                "n": n,
                "total_trees": total_n,
                "checked": checked_n,
                "unique_full_keys": len(full_seen),
                "full_collisions": full_collisions,
                "elapsed_sec": time.time() - t0,
            }
        )
        print(
            f"n={n:2d} total={total_n:8d} checked={checked_n:7d} "
            f"unique_full={len(full_seen):7d} full_collisions={full_collisions:7d} "
            f"time={time.time()-t0:.2f}s",
            flush=True,
        )

    return {
        "checked_total": checked_total,
        "full_unique_keys": len(full_seen),
        "full_collisions": full_collisions,
        "skipped_dleaf": skipped_dleaf,
        "skipped_no_triplet": skipped_no_triplet,
        "skipped_bad": skipped_bad,
        "skipped_m": skipped_m,
        "per_n": per_n,
        "elapsed_sec_scan": time.time() - started,
        "full_seen": full_seen,
    }


def run_projection_stage(
    full_seen: dict[tuple[int, ...], dict[str, Any]],
    checked_total: int,
    projections: list[tuple[str, tuple[str, ...]]],
    analyze_c1c2: bool = False,
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for name, fields in projections:
        t0 = time.time()
        seen: dict[tuple[int, ...], tuple[int, tuple[int, ...]]] = {}
        split: dict[str, Any] | None = None
        for fk, info in full_seen.items():
            pkey = project_full_key(fk, fields)
            prev = seen.get(pkey)
            if prev is None:
                seen[pkey] = (int(info["N"]), fk)
            else:
                n0, fk0 = prev
                if int(info["N"]) != n0 and split is None:
                    info0 = full_seen[fk0]
                    split = {
                        "A": {
                            "n": int(info0["n"]),
                            "g6": str(info0["g6"]),
                            "N": int(info0["N"]),
                            "full_key": list(fk0),
                        },
                        "B": {
                            "n": int(info["n"]),
                            "g6": str(info["g6"]),
                            "N": int(info["N"]),
                            "full_key": list(fk),
                        },
                        "A_invariants": unpack_full_key(fk0),
                        "B_invariants": unpack_full_key(fk),
                        "projection_key": list(pkey),
                    }
        unique = len(seen)
        collisions = checked_total - unique
        rec = {
            "name": name,
            "fields": list(fields),
            "checked": checked_total,
            "unique_keys": unique,
            "collisions": collisions,
            "split_found": split is not None,
            "first_split": split,
            "elapsed_sec": time.time() - t0,
        }

        # Optional theorem-coverage analysis for projection (m, lambda, rho):
        # certify singleton interval via (C1,C2) criterion and compare to observed N.
        if analyze_c1c2 and fields == ("m", "lambda", "rho"):
            stats = {
                "keys_total": unique,
                "c1_true": 0,
                "c2_true": 0,
                "c1c2_true": 0,
                "c1c2_and_NeqL": 0,
                "c1c2_and_NneL": 0,
            }

            for pkey, (n_obs, _fk) in seen.items():
                # pkey layout for fields ("m","lambda","rho"):
                # (m, lam_n, lam_d, rho_n, rho_d)
                m = int(pkey[0])
                lam = Fraction(int(pkey[1]), int(pkey[2]))
                rho = Fraction(int(pkey[3]), int(pkey[4]))

                # L1 = ceil(m - 4 + m/lambda)
                l1_val = Fraction(m - 4, 1) + Fraction(m, 1) / lam
                l1 = (l1_val.numerator + l1_val.denominator - 1) // l1_val.denominator

                # r_min = least integer r s.t. (1+lambda)^r >= lambda/rho
                r_min = 0
                lhs = Fraction(1, 1)
                base = Fraction(1, 1) + lam
                target = lam / rho
                while lhs < target:
                    lhs *= base
                    r_min += 1
                l2 = max(0, 2 * r_min - 1)
                l = max(l1, l2)

                # E = smallest even integer strictly greater than L
                e = 2 * ((l + 2) // 2)

                # C1: ((E+5)(E+4))/((E+5-m)(E+4-m)) <= 1+lambda
                den1 = (e + 5 - m) * (e + 4 - m)
                if den1 <= 0:
                    c1 = False
                else:
                    r_e = Fraction((e + 5) * (e + 4), den1)
                    c1 = (r_e <= (Fraction(1, 1) + lam))

                # C2: A*(1+lambda)^(E/2) > binom(E+3,m)
                # A = (1-lambda)*((1+2lambda)+(1+lambda)rho)
                a_const = (Fraction(1, 1) - lam) * (
                    (Fraction(1, 1) + 2 * lam) + (Fraction(1, 1) + lam) * rho
                )
                c2 = (a_const * (base ** (e // 2)) > Fraction(comb(e + 3, m), 1))

                if c1:
                    stats["c1_true"] += 1
                if c2:
                    stats["c2_true"] += 1
                if c1 and c2:
                    stats["c1c2_true"] += 1
                    if int(n_obs) == int(l):
                        stats["c1c2_and_NeqL"] += 1
                    else:
                        stats["c1c2_and_NneL"] += 1

            rec["c1c2_analysis"] = stats

        out.append(rec)
        print(
            f"projection={name:28s} unique={unique:8d} collisions={collisions:8d} "
            f"split={rec['split_found']} time={rec['elapsed_sec']:.2f}s",
            flush=True,
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-n", type=int, default=3)
    ap.add_argument("--max-n", type=int, default=24)
    ap.add_argument("--m-min", type=int, default=4)
    ap.add_argument("--progress-every", type=int, default=0)
    ap.add_argument(
        "--projections",
        default=";".join(DEFAULT_PROJECTIONS),
        help="Semicolon-separated list of comma-separated field lists.",
    )
    ap.add_argument(
        "--analyze-c1c2",
        action="store_true",
        help="For projection (m,lambda,rho), report coverage of the (C1,C2) singleton criterion.",
    )
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    projections = parse_projection_list(args.projections)

    t0 = time.time()
    scan_payload = scan_full_keys(
        min_n=args.min_n,
        max_n=args.max_n,
        m_min=args.m_min,
        progress_every=args.progress_every,
    )
    full_seen = scan_payload.pop("full_seen")
    projection_results = run_projection_stage(
        full_seen=full_seen,
        checked_total=int(scan_payload["checked_total"]),
        projections=projections,
        analyze_c1c2=bool(args.analyze_c1c2),
    )
    total_elapsed = time.time() - t0

    payload = {
        "scan": "canonical_projection_battery_minu",
        "triplet_rule": "min_u_then_support_then_leaf",
        "min_n": args.min_n,
        "max_n": args.max_n,
        "m_min": args.m_min,
        "projection_specs": [list(p[1]) for p in projections],
        **scan_payload,
        "projection_results": projection_results,
        "elapsed_sec_total": total_elapsed,
    }

    print(
        f"done checked={payload['checked_total']} full_unique={payload['full_unique_keys']} "
        f"full_collisions={payload['full_collisions']} elapsed={payload['elapsed_sec_total']:.2f}s",
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
