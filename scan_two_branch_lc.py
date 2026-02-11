#!/usr/bin/env python3
"""Exhaustive scan of C2 trees (<= 2 branch vertices) for LC/SLC/unimodality.

C2 := {T tree : #{v in V(T) : deg(v) >= 3} <= 2}.

Output schema matches `results/two_branch_lc_n24.json`.
"""

from __future__ import annotations

import argparse
import json
import time
from datetime import datetime, timezone

from indpoly import (
    independence_poly,
    is_log_concave,
    is_unimodal,
    log_concavity_ratio,
    near_miss_ratio,
)
from trees import trees


def count_branch_vertices(adj: list[list[int]]) -> int:
    return sum(1 for nbrs in adj if len(nbrs) >= 3)


def is_strict_log_concave(seq: list[int]) -> bool:
    for k in range(1, len(seq) - 1):
        if seq[k] * seq[k] <= seq[k - 1] * seq[k + 1]:
            return False
    return True


def failure_payload(
    n: int,
    branch_vertices: int,
    poly: list[int],
    lc_ratio: float,
    lc_pos: int,
    nm_ratio: float,
    nm_pos: int,
) -> dict:
    return {
        "n": n,
        "branch_vertices": branch_vertices,
        "lc_ratio": lc_ratio,
        "lc_pos": lc_pos,
        "nm_ratio": nm_ratio,
        "nm_pos": nm_pos,
        "poly": poly,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Scan C2 trees (<=2 branch vertices) for LC/SLC/unimodality"
    )
    parser.add_argument("--min-n", type=int, default=1, help="Minimum n (default: 1)")
    parser.add_argument("--max-n", type=int, default=24, help="Maximum n (default: 24)")
    parser.add_argument(
        "--backend",
        type=str,
        default="auto",
        choices=["auto", "geng", "networkx"],
        help="Tree enumeration backend (default: auto)",
    )
    parser.add_argument(
        "--out",
        type=str,
        default="results/two_branch_lc_n24.json",
        help="Output JSON path",
    )
    args = parser.parse_args()

    if args.min_n < 1:
        raise ValueError("min-n must be >= 1")
    if args.max_n < args.min_n:
        raise ValueError("max-n must be >= min-n")

    t0_global = time.time()

    by_n: list[dict] = []
    lc_failures: list[dict] = []
    strict_lc_failures: list[dict] = []
    unimodality_failures: list[dict] = []

    totals = {
        "total_trees": 0,
        "c2_trees": 0,
        "lc_failures": 0,
        "strict_lc_failures": 0,
        "unimodality_failures": 0,
    }

    worst_lc_ratio = {
        "ratio": 0.0,
        "n": None,
        "branch_vertices": None,
        "lc_pos": None,
        "nm_ratio": None,
        "nm_pos": None,
        "poly": None,
    }
    worst_near_miss_ratio = {
        "ratio": 0.0,
        "n": None,
        "branch_vertices": None,
        "nm_pos": None,
        "lc_ratio": None,
        "lc_pos": None,
        "poly": None,
    }

    for n in range(args.min_n, args.max_n + 1):
        t0_n = time.time()
        total_n = 0
        c2_n = 0
        lc_fail_n = 0
        strict_lc_fail_n = 0
        unimodal_fail_n = 0

        for _, adj in trees(n, backend=args.backend):
            total_n += 1
            branch_vertices = count_branch_vertices(adj)
            if branch_vertices > 2:
                continue

            c2_n += 1
            poly = independence_poly(n, adj)

            lc = is_log_concave(poly)
            strict_lc = is_strict_log_concave(poly)
            unimodal = is_unimodal(poly)
            lc_ratio, lc_pos = log_concavity_ratio(poly)
            nm_ratio, nm_pos = near_miss_ratio(poly)

            if lc_ratio > worst_lc_ratio["ratio"]:
                worst_lc_ratio = {
                    "ratio": lc_ratio,
                    "n": n,
                    "branch_vertices": branch_vertices,
                    "lc_pos": lc_pos,
                    "nm_ratio": nm_ratio,
                    "nm_pos": nm_pos,
                    "poly": poly,
                }

            if nm_ratio > worst_near_miss_ratio["ratio"]:
                worst_near_miss_ratio = {
                    "ratio": nm_ratio,
                    "n": n,
                    "branch_vertices": branch_vertices,
                    "nm_pos": nm_pos,
                    "lc_ratio": lc_ratio,
                    "lc_pos": lc_pos,
                    "poly": poly,
                }

            if not lc:
                lc_fail_n += 1
                totals["lc_failures"] += 1
                lc_failures.append(
                    failure_payload(
                        n,
                        branch_vertices,
                        poly,
                        lc_ratio,
                        lc_pos,
                        nm_ratio,
                        nm_pos,
                    )
                )

            if not strict_lc:
                strict_lc_fail_n += 1
                totals["strict_lc_failures"] += 1
                strict_lc_failures.append(
                    failure_payload(
                        n,
                        branch_vertices,
                        poly,
                        lc_ratio,
                        lc_pos,
                        nm_ratio,
                        nm_pos,
                    )
                )

            if not unimodal:
                unimodal_fail_n += 1
                totals["unimodality_failures"] += 1
                unimodality_failures.append(
                    failure_payload(
                        n,
                        branch_vertices,
                        poly,
                        lc_ratio,
                        lc_pos,
                        nm_ratio,
                        nm_pos,
                    )
                )

        totals["total_trees"] += total_n
        totals["c2_trees"] += c2_n

        row = {
            "n": n,
            "total_trees": total_n,
            "c2_trees": c2_n,
            "lc_failures": lc_fail_n,
            "strict_lc_failures": strict_lc_fail_n,
            "unimodality_failures": unimodal_fail_n,
            "elapsed_s": round(time.time() - t0_n, 3),
        }
        by_n.append(row)

        print(
            f"n={n:2d} total={total_n:9d} c2={c2_n:9d} "
            f"lc_fail={lc_fail_n:6d} strict_lc_fail={strict_lc_fail_n:6d} "
            f"uni_fail={unimodal_fail_n:6d} t={row['elapsed_s']:.2f}s",
            flush=True,
        )

    out = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "class_definition": "C2 = trees with at most two branch vertices (degree >= 3)",
        "range": {"min_n": args.min_n, "max_n": args.max_n},
        "backend": args.backend,
        "by_n": by_n,
        "totals": totals,
        "worst_lc_ratio": worst_lc_ratio,
        "worst_near_miss_ratio": worst_near_miss_ratio,
        "lc_failures": lc_failures,
        "strict_lc_failures": strict_lc_failures,
        "unimodality_failures": unimodality_failures,
        "elapsed_total_s": round(time.time() - t0_global, 3),
    }

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)

    print("\nTotals:")
    print(json.dumps(totals, indent=2))
    print(f"Worst LC ratio: {worst_lc_ratio['ratio']:.12f}")
    print(f"Worst near-miss ratio: {worst_near_miss_ratio['ratio']:.12f}")
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
