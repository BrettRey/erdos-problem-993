#!/usr/bin/env python3
"""Scan K*->N split behavior on triplet-invariant trees.

A tree is triplet-invariant if all admissible bridge triplets (leaf,support,u)
produce the same anchored key
  K* = (d,m,lambda,mu1,mu2,rho,sigma).

This script:
1) enumerates nonisomorphic trees in [min_n, max_n]
2) filters by d_leaf<=1
3) computes K* for every admissible triplet
4) keeps trees where K* is triplet-invariant
5) checks for same-K* / different-N collisions among that subclass

Optional:
  --m-min 3 to restrict to the only nontrivial regime after the m=2 identity.
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
    admissible_triplets,
    build_record_from_triplet,
    key_from_record,
)


def scan(
    min_n: int,
    max_n: int,
    m_min: int,
    progress_every: int,
    within_n_only: bool,
) -> dict[str, Any]:
    started = time.time()

    seen: dict[tuple[int, ...], tuple[int, str, int]] = {}
    checked_total = 0
    dleaf_total = 0
    gated_total = 0
    invariant_total = 0
    collisions = 0
    split: dict[str, Any] | None = None

    per_n: list[dict[str, Any]] = []

    for n in range(min_n, max_n + 1):
        if within_n_only:
            seen.clear()

        total_n = 0
        dleaf_n = 0
        gated_n = 0
        invariant_n = 0
        checked_n = 0
        t0 = time.time()

        for nn, adj, raw in trees_geng_raw(n):
            total_n += 1
            if progress_every > 0 and total_n % progress_every == 0:
                print(
                    f"progress n={n} total={total_n} dleaf={dleaf_n} "
                    f"gated={gated_n} invariant={invariant_n} "
                    f"checked={checked_n} collisions={collisions}",
                    flush=True,
                )

            if not is_dleaf_le_1(nn, adj):
                continue
            dleaf_total += 1
            dleaf_n += 1

            g6 = raw.decode("ascii").strip()
            trips = admissible_triplets(adj)
            if not trips:
                continue
            gated_total += 1
            gated_n += 1

            recs: list[dict[str, Any]] = []
            keyset: set[tuple[int, ...]] = set()
            bad = False
            for trip in trips:
                rec = build_record_from_triplet(nn, adj, g6, trip)
                if rec is None:
                    bad = True
                    break
                keyset.add(key_from_record(rec))
                recs.append(rec)
            if bad or not recs:
                continue

            if len(keyset) != 1:
                continue
            invariant_total += 1
            invariant_n += 1

            # Any rec works because the key is invariant across triplets.
            rep = recs[0]
            if rep["m"] < m_min:
                continue

            key = key_from_record(rep)
            prev = seen.get(key)
            if prev is None:
                seen[key] = (nn, g6, rep["N"])
            else:
                collisions += 1
                pn, pg6, pN = prev
                if pN != rep["N"]:
                    # Emit full records for both trees (all triplets) for auditability.
                    def all_recs(n0: int, adj0: list[list[int]], g60: str) -> list[dict[str, Any]]:
                        out: list[dict[str, Any]] = []
                        for t in admissible_triplets(adj0):
                            r = build_record_from_triplet(n0, adj0, g60, t)
                            if r is not None:
                                out.append(r)
                        return out

                    pnn, padj, praw = next(
                        (xn, xa, xr)
                        for xn, xa, xr in trees_geng_raw(pn)
                        if xr.decode("ascii").strip() == pg6
                    )
                    split = {
                        "A": {
                            "n": pn,
                            "g6": pg6,
                            "N": pN,
                            "triplets": all_recs(pnn, padj, pg6),
                        },
                        "B": {
                            "n": nn,
                            "g6": g6,
                            "N": rep["N"],
                            "triplets": recs,
                        },
                    }
                    break

            checked_total += 1
            checked_n += 1

        per_n.append(
            {
                "n": n,
                "total_trees": total_n,
                "dleaf_trees": dleaf_n,
                "gated_trees": gated_n,
                "triplet_invariant_trees": invariant_n,
                "checked": checked_n,
                "unique_keys": len(seen),
                "collisions": collisions,
                "elapsed_sec": time.time() - t0,
            }
        )
        print(
            f"n={n:2d} total={total_n:8d} dleaf={dleaf_n:7d} gated={gated_n:7d} "
            f"invariant={invariant_n:6d} checked={checked_n:6d} "
            f"unique={len(seen):7d} collisions={collisions:6d} "
            f"time={time.time()-t0:.2f}s",
            flush=True,
        )
        if split is not None:
            break

    return {
        "scan": "canonical_kstar_triplet_invariance_scan",
        "key": "(d,m,lambda,mu1,mu2,rho,sigma)",
        "min_n": min_n,
        "max_n": per_n[-1]["n"] if per_n else min_n,
        "m_min": m_min,
        "within_n_only": within_n_only,
        "dleaf_total": dleaf_total,
        "gated_total": gated_total,
        "triplet_invariant_total": invariant_total,
        "checked_total": checked_total,
        "unique_keys": len(seen),
        "collisions": collisions,
        "split_found": split is not None,
        "split": split,
        "per_n": per_n,
        "elapsed_sec": time.time() - started,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-n", type=int, default=3)
    ap.add_argument("--max-n", type=int, default=22)
    ap.add_argument("--m-min", type=int, default=0)
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
        f"done dleaf={payload['dleaf_total']} gated={payload['gated_total']} "
        f"invariant={payload['triplet_invariant_total']} checked={payload['checked_total']} "
        f"unique={payload['unique_keys']} collisions={payload['collisions']} "
        f"split_found={payload['split_found']} elapsed={payload['elapsed_sec']:.2f}s",
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

