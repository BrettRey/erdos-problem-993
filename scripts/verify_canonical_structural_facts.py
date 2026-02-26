#!/usr/bin/env python3
"""Verify structural facts about canonical bridge decomposition on d_leaf<=1 trees.

Checks (exact, exhaustive over nonisomorphic trees in a size range):
1) Canonical lambda = i_{m-1}/i_m is always strictly < 1.
2) Number of admissible degree-2 bridge triplets is never exactly 1
   under d_leaf<=1.
3) For m=2 cases, the universal tree identity holds:
     i2 = (i1-1)(i1-2)/2
   (equivalently lambda = i1/i2 determines i1).
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

from attack4_common import bridge_decomposition
from conjecture_a_hall_subset_scan import is_dleaf_le_1
from trees import trees_geng_raw


def count_admissible_triplets(adj: list[list[int]]) -> int:
    """Count all admissible (leaf, support, u) triplets in one tree."""
    n = len(adj)
    out = 0
    for leaf in range(n):
        if len(adj[leaf]) != 1:
            continue
        support = adj[leaf][0]
        if len(adj[support]) != 2:
            continue
        out += 1
    return out


def run(min_n: int, max_n: int, progress_every: int) -> dict[str, Any]:
    started = time.time()

    checked_total = 0
    dleaf_total = 0
    gated_total = 0

    one_triplet_count = 0
    one_triplet_examples: list[dict[str, Any]] = []

    min_lambda: float | None = None
    max_lambda = -1.0
    max_lambda_witness: dict[str, Any] | None = None

    m2_cases = 0
    m2_identity_failures = 0
    m2_failure_example: dict[str, Any] | None = None

    per_n: list[dict[str, Any]] = []

    for n in range(min_n, max_n + 1):
        checked_n = 0
        dleaf_n = 0
        gated_n = 0
        one_triplet_n = 0
        t0 = time.time()

        for nn, adj, raw in trees_geng_raw(n):
            checked_total += 1
            checked_n += 1
            if progress_every > 0 and checked_n % progress_every == 0:
                print(
                    f"n={n} checked_n={checked_n} dleaf_n={dleaf_n} gated_n={gated_n}",
                    flush=True,
                )

            if not is_dleaf_le_1(nn, adj):
                continue
            dleaf_total += 1
            dleaf_n += 1

            tcnt = count_admissible_triplets(adj)
            if tcnt == 1:
                one_triplet_count += 1
                one_triplet_n += 1
                if len(one_triplet_examples) < 10:
                    one_triplet_examples.append(
                        {"n": n, "g6": raw.decode("ascii").strip()}
                    )

            g6 = raw.decode("ascii").strip()
            decomp = bridge_decomposition(nn, adj, g6, require_dleaf=True)
            if decomp is None:
                continue

            gated_total += 1
            gated_n += 1

            m = decomp.m_t
            i = decomp.poly_t
            i_m = i[m]
            i_m1 = i[m - 1] if m - 1 >= 0 else 0
            lam = i_m1 / i_m if i_m else 0.0

            if min_lambda is None or lam < min_lambda:
                min_lambda = lam
            if lam > max_lambda:
                max_lambda = lam
                max_lambda_witness = {
                    "n": n,
                    "g6": g6,
                    "m": m,
                    "i_m_minus_1": i_m1,
                    "i_m": i_m,
                }

            if m == 2:
                m2_cases += 1
                i1 = i[1] if len(i) > 1 else 0
                i2 = i[2] if len(i) > 2 else 0
                rhs = (i1 - 1) * (i1 - 2) // 2
                if i2 != rhs:
                    m2_identity_failures += 1
                    if m2_failure_example is None:
                        m2_failure_example = {
                            "n": n,
                            "g6": g6,
                            "i1": i1,
                            "i2": i2,
                            "rhs": rhs,
                        }

        per_n.append(
            {
                "n": n,
                "checked": checked_n,
                "dleaf": dleaf_n,
                "gated": gated_n,
                "one_triplet": one_triplet_n,
                "elapsed_sec": time.time() - t0,
            }
        )
        print(
            f"done n={n} checked={checked_n} dleaf={dleaf_n} gated={gated_n} "
            f"one_triplet={one_triplet_n}",
            flush=True,
        )

    payload = {
        "scan": "verify_canonical_structural_facts",
        "min_n": min_n,
        "max_n": max_n,
        "checked_total": checked_total,
        "dleaf_total": dleaf_total,
        "gated_total": gated_total,
        "one_triplet_count": one_triplet_count,
        "one_triplet_examples": one_triplet_examples,
        "lambda_lt_1_holds": bool(max_lambda < 1.0),
        "min_lambda": min_lambda,
        "max_lambda": max_lambda,
        "max_lambda_witness": max_lambda_witness,
        "m2_cases": m2_cases,
        "m2_identity_failures": m2_identity_failures,
        "m2_failure_example": m2_failure_example,
        "per_n": per_n,
        "elapsed_sec": time.time() - started,
    }
    return payload


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-n", type=int, default=3)
    ap.add_argument("--max-n", type=int, default=22)
    ap.add_argument("--progress-every", type=int, default=0)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    payload = run(args.min_n, args.max_n, args.progress_every)
    print(
        f"done checked={payload['checked_total']} dleaf={payload['dleaf_total']} "
        f"gated={payload['gated_total']} one_triplet={payload['one_triplet_count']} "
        f"lambda_lt_1={payload['lambda_lt_1_holds']} m2_failures={payload['m2_identity_failures']}",
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

