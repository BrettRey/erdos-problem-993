#!/usr/bin/env python3
"""Stress-test conditional-index localization for signed PB laws."""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path
from typing import Any

from analyze_signed_conditionals import analyze_row, compact_analysis
from probe_signed_pb_reserve import (
    add_finite_skellam_rows,
    add_random_rows,
    add_two_binomial_rows,
)


def best_by_cutoff(
    analyses: list[dict[str, Any]],
    cutoffs: list[float],
    key: str,
    *,
    reverse: bool,
) -> list[dict[str, Any]]:
    out = []
    for cutoff in cutoffs:
        candidates = [row for row in analyses if row["variance"] >= cutoff]
        best = sorted(candidates, key=lambda row: row[key], reverse=reverse)
        out.append(
            {
                "variance_cutoff": cutoff,
                "best": compact_analysis(best[0]) if best else None,
            }
        )
    return out


def finite_rows(rows: list[dict[str, Any]], key: str) -> list[dict[str, Any]]:
    return [row for row in rows if math.isfinite(row[key])]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--random-samples", type=int, default=2000)
    parser.add_argument("--max-groups", type=int, default=6)
    parser.add_argument("--binomial-grid-size", type=int, default=31)
    parser.add_argument("--candidate-bound", type=float, default=3.0)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/signed_pb_conditional_index_probe_2026-07-03.json"),
    )
    args = parser.parse_args()

    rng = random.Random(args.seed)
    rows: list[dict[str, Any]] = []
    add_random_rows(
        rows,
        rng=rng,
        samples=args.random_samples,
        max_groups=args.max_groups,
    )
    add_two_binomial_rows(rows, grid_size=args.binomial_grid_size)
    add_finite_skellam_rows(rows)

    analyses = []
    for index, row in enumerate(rows):
        analysis = analyze_row(row, f"generated_rank={index + 1}")
        if analysis is not None:
            analyses.append(analysis)

    x_key = "conditional_x_index_plus_one_over_variance"
    y_key = "conditional_y_index_plus_one_over_variance"
    x_rows = finite_rows(analyses, x_key)
    y_rows = finite_rows(analyses, y_key)
    x_failures = [
        row
        for row in x_rows
        if row["variance"] >= 1.0 and row[x_key] > args.candidate_bound
    ]
    y_failures = [
        row
        for row in y_rows
        if row["variance"] >= 1.0 and row[y_key] > args.candidate_bound
    ]
    identity_errors = [
        row["max_identity_error"]
        for row in analyses
        if math.isfinite(row["max_identity_error"])
    ]
    summary = {
        "source": {
            "kind": "signed_pb_conditional_index_probe",
            "seed": args.seed,
            "random_samples": args.random_samples,
            "max_groups": args.max_groups,
            "binomial_grid_size": args.binomial_grid_size,
            "candidate_bound": args.candidate_bound,
        },
        "processed": len(rows),
        "analyzed": len(analyses),
        "max_identity_error": max(identity_errors, default=math.nan),
        "candidate_bound_failures": {
            "x_index_over_variance": len(x_failures),
            "y_index_over_variance": len(y_failures),
        },
        "max_x_index_by_variance_cutoff": best_by_cutoff(
            x_rows,
            [1, 2, 5, 10, 20, 50],
            x_key,
            reverse=True,
        ),
        "max_y_index_by_variance_cutoff": best_by_cutoff(
            y_rows,
            [1, 2, 5, 10, 20, 50],
            y_key,
            reverse=True,
        ),
        "largest_x_index_over_variance": [
            compact_analysis(row)
            for row in sorted(x_rows, key=lambda row: row[x_key], reverse=True)[
                : args.top
            ]
        ],
        "largest_y_index_over_variance": [
            compact_analysis(row)
            for row in sorted(y_rows, key=lambda row: row[y_key], reverse=True)[
                : args.top
            ]
        ],
        "x_index_bound_failures": [
            compact_analysis(row)
            for row in sorted(x_failures, key=lambda row: row[x_key], reverse=True)[
                : args.top
            ]
        ],
        "y_index_bound_failures": [
            compact_analysis(row)
            for row in sorted(y_failures, key=lambda row: row[y_key], reverse=True)[
                : args.top
            ]
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Wrote {args.out}; processed={len(rows)}, analyzed={len(analyses)}, "
        f"x_failures={len(x_failures)}, y_failures={len(y_failures)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
