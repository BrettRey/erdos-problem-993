#!/usr/bin/env python3
"""Probe effective signed ratio drop at first descent.

For signed coefficients c_z, the useful sufficient diagnostic is

    1 - (c_{D+1}/c_D) / (c_D/c_{D-1}),

where D is the first strict descent.  Since c_D/c_{D-1} < 1, a lower bound on
this effective ratio drop implies the same lower bound for the reserve
1 - c_{D+1}/c_D.  This script treats that as a falsification target only.
"""

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


def best_by_cutoff(rows: list[dict[str, Any]], cutoffs: list[float]) -> list[dict[str, Any]]:
    out = []
    for cutoff in cutoffs:
        candidates = [row for row in rows if row["variance"] >= cutoff]
        best = min(
            candidates,
            key=lambda row: row["variance_times_effective_ratio_drop"],
            default=None,
        )
        out.append({"variance_cutoff": cutoff, "best": compact_analysis(best) if best else None})
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--random-samples", type=int, default=2000)
    parser.add_argument("--max-groups", type=int, default=6)
    parser.add_argument("--binomial-grid-size", type=int, default=31)
    parser.add_argument("--candidate-constant", type=float, default=0.25)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/signed_pb_ratio_drop_probe_2026-07-03.json"),
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

    feasible = [
        row
        for row in analyses
        if row["variance"] >= 1.0
        and math.isfinite(row["variance_times_effective_ratio_drop"])
    ]
    failures = [
        row
        for row in feasible
        if row["variance_times_effective_ratio_drop"] < args.candidate_constant
    ]
    identity_errors = [
        row["max_identity_error"]
        for row in analyses
        if math.isfinite(row["max_identity_error"])
    ]
    summary = {
        "source": {
            "kind": "signed_pb_effective_ratio_drop_probe",
            "seed": args.seed,
            "random_samples": args.random_samples,
            "max_groups": args.max_groups,
            "binomial_grid_size": args.binomial_grid_size,
            "candidate_constant": args.candidate_constant,
        },
        "processed": len(rows),
        "analyzed": len(analyses),
        "variance_ge_1": len(feasible),
        "candidate_failures": len(failures),
        "max_identity_error": max(identity_errors, default=math.nan),
        "best_by_variance_cutoff": best_by_cutoff(
            feasible,
            [1, 2, 5, 10, 20, 50],
        ),
        "smallest_effective_ratio_drops": [
            compact_analysis(row)
            for row in sorted(
                feasible,
                key=lambda row: row["variance_times_effective_ratio_drop"],
            )[: args.top]
        ],
        "candidate_failures_rows": [
            compact_analysis(row)
            for row in sorted(
                failures,
                key=lambda row: row["variance_times_effective_ratio_drop"],
            )[: args.top]
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Wrote {args.out}; processed={len(rows)}, analyzed={len(analyses)}, "
        f"candidate_failures={len(failures)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
