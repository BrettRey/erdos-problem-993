#!/usr/bin/env python3
"""Adversarial breaker for the signed PB effective-drop route.

This is a falsification harness for issue #5.  It attacks signed laws

    Z = X - Y,

where X and Y are independent low-probability Poisson-binomial sums.  The
targets are the sufficient effective-drop diagnostic

    Delta_eff = 1 - (c_{D+1}/c_D) / (c_D/c_{D-1})

and the raw reserve

    1 - c_{D+1}/c_D

at the first strict descent D in signed support coordinates.  Search output is
only computational evidence unless a strict counterexample certificate appears.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import differential_evolution, minimize

from analyze_signed_conditionals import analyze_row, compact_analysis
from probe_signed_pb_reserve import Block, grouped_pmf, signed_metric


TOL = 1e-10


def clean_blocks(blocks: list[Block]) -> list[Block]:
    clean = [(int(count), min(max(float(p), 0.0), 0.5)) for count, p in blocks if count > 0]
    clean.sort(key=lambda item: item[1])
    merged: list[Block] = []
    for count, p in clean:
        if merged and abs(merged[-1][1] - p) <= 1e-12:
            old_count, old_p = merged[-1]
            merged[-1] = (old_count + count, old_p)
        else:
            merged.append((count, p))
    return merged


def block_key(blocks: list[Block]) -> tuple[tuple[int, int], ...]:
    return tuple((count, round(p * 10**12)) for count, p in clean_blocks(blocks))


def row_key(x_blocks: list[Block], y_blocks: list[Block]) -> tuple[Any, ...]:
    return (block_key(x_blocks), block_key(y_blocks))


def finite_poisson_blocks(
    mean: float,
    *,
    scale: int,
    max_side_n: int,
    min_support: int = 1,
) -> list[Block]:
    if mean <= 0.0:
        return [(min_support, 0.0)]
    n = max(min_support, int(math.ceil(scale * max(mean, 0.05))), int(math.ceil(2.0 * mean)))
    n = min(max_side_n, n)
    n = max(n, int(math.ceil(2.0 * mean)))
    return [(n, min(0.5, mean / n))]


def metric_with_analysis(
    x_blocks: list[Block],
    y_blocks: list[Block],
    *,
    source: str,
    details: dict[str, Any] | None = None,
) -> dict[str, Any] | None:
    x_blocks = clean_blocks(x_blocks)
    y_blocks = clean_blocks(y_blocks)
    metric = signed_metric(x_blocks=x_blocks, y_blocks=y_blocks, kind=source, details=details)
    if metric is None:
        return None
    analysis = analyze_row(metric, source)
    if analysis is None:
        return None
    analysis["x_blocks"] = metric["x_blocks"]
    analysis["y_blocks"] = metric["y_blocks"]
    analysis["x_mean"] = metric["x_mean"]
    analysis["y_mean"] = metric["y_mean"]
    analysis["search_family"] = source
    if details:
        analysis.update(details)
    return analysis


def compact_breaker(row: dict[str, Any]) -> dict[str, Any]:
    out = compact_analysis(row)
    for key in [
        "search_family",
        "x_blocks",
        "y_blocks",
        "x_mean",
        "y_mean",
        "lambda",
        "eta",
        "total_variance_parameter",
        "split",
        "scale",
        "pattern",
        "target",
        "optimizer_success",
        "optimizer_fun",
        "side_balance_floor",
    ]:
        if key in row:
            out[key] = row[key]
    return out


def failure_certificate(row: dict[str, Any]) -> dict[str, Any]:
    x_blocks = [(int(count), float(p)) for count, p in row["x_blocks"]]
    y_blocks = [(int(count), float(p)) for count, p in row["y_blocks"]]
    x_pmf = grouped_pmf(x_blocks)
    y_pmf = grouped_pmf(y_blocks)
    signed = np.convolve(x_pmf, y_pmf[::-1])
    support_start = -(len(y_pmf) - 1)
    return {
        **compact_breaker(row),
        "support_start": support_start,
        "x_pmf": [float(v) for v in x_pmf],
        "y_pmf": [float(v) for v in y_pmf],
        "signed_pmf": [
            {"z": support_start + index, "probability": float(prob)}
            for index, prob in enumerate(signed)
        ],
    }


def add_candidate(
    rows: list[dict[str, Any]],
    seen: set[tuple[Any, ...]],
    counters: dict[str, int],
    x_blocks: list[Block],
    y_blocks: list[Block],
    *,
    source: str,
    details: dict[str, Any] | None = None,
) -> None:
    counters["attempted"] += 1
    key = row_key(x_blocks, y_blocks)
    if key in seen:
        counters["duplicates"] += 1
        return
    seen.add(key)
    row = metric_with_analysis(x_blocks, y_blocks, source=source, details=details)
    if row is None:
        counters["not_analyzed"] += 1
        return
    rows.append(row)
    counters["analyzed"] += 1


def add_finite_skellam_boundary(
    rows: list[dict[str, Any]],
    seen: set[tuple[Any, ...]],
    counters: dict[str, int],
    *,
    max_side_n: int,
) -> None:
    totals = [
        1.0,
        1.001,
        1.005,
        1.01,
        1.025,
        1.05,
        1.1,
        1.25,
        1.5,
        2.0,
        3.0,
        5.0,
        10.0,
        20.0,
        50.0,
    ]
    splits = [
        0.0,
        1e-8,
        1e-7,
        1e-6,
        1e-5,
        1e-4,
        1e-3,
        0.003,
        0.01,
        0.03,
        0.05,
        0.1,
        0.25,
        0.5,
        0.75,
        0.9,
        0.95,
        0.99,
        0.999,
        0.999999,
    ]
    scales = [20, 50, 100, 250, 500, 1000]
    for total in totals:
        for split in splits:
            lam = total * (1.0 - split)
            eta = total * split
            for scale in scales:
                x_blocks = finite_poisson_blocks(lam, scale=scale, max_side_n=max_side_n)
                y_blocks = finite_poisson_blocks(eta, scale=scale, max_side_n=max_side_n)
                details = {
                    "total_variance_parameter": total,
                    "split": split,
                    "lambda": lam,
                    "eta": eta,
                    "scale": scale,
                }
                add_candidate(
                    rows,
                    seen,
                    counters,
                    x_blocks,
                    y_blocks,
                    source="finite_skellam_boundary",
                    details=details,
                )


def add_near_one_sided_perturbations(
    rows: list[dict[str, Any]],
    seen: set[tuple[Any, ...]],
    counters: dict[str, int],
    *,
    max_side_n: int,
) -> None:
    dominant_means = [0.95, 0.99, 1.0, 1.005, 1.01, 1.025, 1.05, 1.1, 1.25, 1.5, 2.0]
    perturb_means = [0.0, 1e-7, 1e-5, 1e-4, 1e-3, 0.003, 0.01, 0.03, 0.1, 0.25]
    half_counts = [0, 1, 2, 3, 4, 5, 8]
    scales = [100, 250, 500, 1000]
    for lam in dominant_means:
        for eta in perturb_means:
            for half_count in half_counts:
                for scale in scales:
                    x_sparse = finite_poisson_blocks(lam, scale=scale, max_side_n=max_side_n)
                    y_sparse = finite_poisson_blocks(eta, scale=scale, max_side_n=max_side_n)
                    x_blocks = x_sparse
                    y_blocks = ([(half_count, 0.5)] if half_count else []) + y_sparse
                    details = {
                        "lambda": lam,
                        "eta": eta,
                        "scale": scale,
                        "half_count_perturbation": half_count,
                    }
                    add_candidate(
                        rows,
                        seen,
                        counters,
                        x_blocks,
                        y_blocks,
                        source="near_one_sided_x_dominant",
                        details=details,
                    )
                    add_candidate(
                        rows,
                        seen,
                        counters,
                        y_blocks,
                        x_blocks,
                        source="near_one_sided_y_dominant",
                        details=details,
                    )


def add_boundary_heavy_grid(
    rows: list[dict[str, Any]],
    seen: set[tuple[Any, ...]],
    counters: dict[str, int],
    *,
    max_side_n: int,
) -> None:
    sparse_means = [0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.01, 1.25, 1.5, 2.0, 3.0]
    half_counts = range(0, 9)
    scales = [50, 200, 800]
    for hx in half_counts:
        for hy in half_counts:
            for lam in sparse_means:
                for eta in sparse_means:
                    for scale in scales:
                        x_blocks = ([(hx, 0.5)] if hx else []) + finite_poisson_blocks(
                            lam,
                            scale=scale,
                            max_side_n=max_side_n,
                        )
                        y_blocks = ([(hy, 0.5)] if hy else []) + finite_poisson_blocks(
                            eta,
                            scale=scale,
                            max_side_n=max_side_n,
                        )
                        add_candidate(
                            rows,
                            seen,
                            counters,
                            x_blocks,
                            y_blocks,
                            source="boundary_heavy_grid",
                            details={
                                "x_half_count": hx,
                                "y_half_count": hy,
                                "lambda": lam,
                                "eta": eta,
                                "scale": scale,
                            },
                        )


def random_blocks(
    rng: random.Random,
    *,
    n: int,
    groups: int,
    family: str,
) -> list[Block]:
    groups = min(groups, n)
    cuts = sorted(rng.sample(range(1, n), groups - 1)) if groups > 1 else []
    counts = [cuts[0], *[cuts[i] - cuts[i - 1] for i in range(1, len(cuts))], n - cuts[-1]] if cuts else [n]
    probs = []
    for _ in counts:
        if family == "sparse":
            probs.append(0.5 * rng.betavariate(0.25, 6.0))
        elif family == "near_half":
            probs.append(0.5 * (1.0 - rng.random() ** 5))
        elif family == "edge":
            probs.append(0.5 * rng.betavariate(0.25, 0.25))
        elif family == "mixed":
            probs.append(rng.choice([0.5 * rng.random(), 0.5 * rng.betavariate(0.25, 6.0), 0.5]))
        else:
            probs.append(0.5 * rng.random())
    return clean_blocks(list(zip(counts, probs)))


def add_random_large_grouped(
    rows: list[dict[str, Any]],
    seen: set[tuple[Any, ...]],
    counters: dict[str, int],
    *,
    rng: random.Random,
    samples: int,
    max_groups: int,
) -> None:
    n_values = [20, 50, 100, 200, 400, 800, 1200]
    families = ["sparse", "near_half", "edge", "mixed", "uniform"]
    for sample in range(samples):
        x_n = rng.choice(n_values)
        y_n = rng.choice(n_values)
        x_groups = rng.randint(1, min(max_groups, x_n))
        y_groups = rng.randint(1, min(max_groups, y_n))
        add_candidate(
            rows,
            seen,
            counters,
            random_blocks(rng, n=x_n, groups=x_groups, family=rng.choice(families)),
            random_blocks(rng, n=y_n, groups=y_groups, family=rng.choice(families)),
            source="random_large_grouped",
            details={"sample": sample},
        )


def evaluate_pattern_values(
    x_counts: list[int],
    y_counts: list[int],
    values: list[float],
    *,
    source: str,
    details: dict[str, Any],
) -> dict[str, Any] | None:
    x_values = values[: len(x_counts)]
    y_values = values[len(x_counts) :]
    x_blocks = list(zip(x_counts, x_values))
    y_blocks = list(zip(y_counts, y_values))
    return metric_with_analysis(x_blocks, y_blocks, source=source, details=details)


def target_value(row: dict[str, Any], target: str) -> float:
    if target == "reserve":
        return float(row["variance_times_reserve"])
    return float(row["variance_times_effective_ratio_drop"])


def optimizer_penalty(row: dict[str, Any] | None, *, target: str, variance_cutoff: float, side_floor: float) -> float:
    if row is None:
        return 1_000.0
    penalty = 0.0
    if row["variance"] < variance_cutoff:
        penalty += 200.0 + 50.0 * (variance_cutoff - row["variance"])
    side_required = side_floor * row["variance"]
    side_shortfall = max(0.0, side_required - row["x_variance"]) + max(
        0.0,
        side_required - row["y_variance"],
    )
    if side_shortfall > 0.0:
        penalty += 200.0 + 50.0 * side_shortfall
    return target_value(row, target) + penalty


def add_optimizer_rows(
    rows: list[dict[str, Any]],
    seen: set[tuple[Any, ...]],
    counters: dict[str, int],
    *,
    rng: random.Random,
    de_maxiter: int,
    de_restarts: int,
    de_popsize: int,
    variance_cutoff: float,
) -> None:
    patterns = [
        ([200], [50]),
        ([500], [50]),
        ([1000], [20]),
        ([1200], [20]),
        ([200], [200]),
        ([500], [500]),
        ([100, 100], [50]),
        ([200, 50], [50, 10]),
        ([400, 80, 8], [80, 8]),
        ([800, 100, 10], [50, 10, 2]),
        ([200, 20, 5, 1], [50, 10, 2, 1]),
        ([100, 100, 50, 10], [100, 50, 10, 2]),
    ]
    side_floors = [0.0, 0.02, 0.05, 0.1]
    targets = ["effective", "reserve"]
    for pattern_index, (x_counts, y_counts) in enumerate(patterns):
        dimension = len(x_counts) + len(y_counts)
        bounds = [(0.0, 0.5)] * dimension
        for target in targets:
            for side_floor in side_floors:
                for restart in range(de_restarts):
                    seed = rng.randrange(1, 2**31 - 1)

                    def objective(values: Any) -> float:
                        row = evaluate_pattern_values(
                            x_counts,
                            y_counts,
                            [float(v) for v in values],
                            source="differential_evolution_breaker",
                            details={
                                "pattern": [x_counts, y_counts],
                                "target": target,
                                "side_balance_floor": side_floor,
                            },
                        )
                        return optimizer_penalty(
                            row,
                            target=target,
                            variance_cutoff=variance_cutoff,
                            side_floor=side_floor,
                        )

                    result = differential_evolution(
                        objective,
                        bounds,
                        seed=seed,
                        maxiter=de_maxiter,
                        popsize=de_popsize,
                        polish=False,
                        updating="immediate",
                        workers=1,
                        tol=1e-7,
                    )
                    starts = [list(float(v) for v in result.x)]
                    starts.extend(seed_vectors_for_pattern(x_counts, y_counts, rng))
                    for start_index, start in enumerate(starts[:4]):
                        local = minimize(
                            objective,
                            start,
                            method="Nelder-Mead",
                            options={"maxiter": 140, "xatol": 1e-8, "fatol": 1e-8},
                        )
                        values = [float(min(max(v, 0.0), 0.5)) for v in local.x]
                        details = {
                            "pattern": [x_counts, y_counts],
                            "target": target,
                            "side_balance_floor": side_floor,
                            "optimizer_success": bool(local.success),
                            "optimizer_fun": float(local.fun),
                            "de_fun": float(result.fun),
                            "de_seed": seed,
                            "pattern_index": pattern_index,
                            "local_start_index": start_index,
                        }
                        x_blocks = list(zip(x_counts, values[: len(x_counts)]))
                        y_blocks = list(zip(y_counts, values[len(x_counts) :]))
                        add_candidate(
                            rows,
                            seen,
                            counters,
                            x_blocks,
                            y_blocks,
                            source="optimizer_breaker",
                            details=details,
                        )


def seed_vectors_for_pattern(
    x_counts: list[int],
    y_counts: list[int],
    rng: random.Random,
) -> list[list[float]]:
    dimension = len(x_counts) + len(y_counts)
    seeds: list[list[float]] = []
    seeds.append([0.005] * len(x_counts) + [0.0001] * len(y_counts))
    seeds.append([1.0 / max(2.0, sum(x_counts))] * len(x_counts) + [1e-5] * len(y_counts))
    seeds.append([0.5] * dimension)
    seeds.append([0.005] * dimension)
    for _ in range(3):
        seeds.append([0.5 * rng.betavariate(0.25, 6.0) for _ in range(dimension)])
    return seeds


def best_by_cutoff(
    rows: list[dict[str, Any]],
    cutoffs: list[float],
    *,
    key: str,
) -> list[dict[str, Any]]:
    out = []
    for cutoff in cutoffs:
        candidates = [row for row in rows if row["variance"] >= cutoff and math.isfinite(row[key])]
        best = min(candidates, key=lambda row: row[key], default=None)
        out.append({"variance_cutoff": cutoff, "best": compact_breaker(best) if best else None})
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--candidate-constant", type=float, default=0.25)
    parser.add_argument("--variance-cutoff", type=float, default=1.0)
    parser.add_argument("--max-side-n", type=int, default=1200)
    parser.add_argument("--random-samples", type=int, default=1200)
    parser.add_argument("--max-groups", type=int, default=8)
    parser.add_argument("--de-maxiter", type=int, default=12)
    parser.add_argument("--de-restarts", type=int, default=1)
    parser.add_argument("--de-popsize", type=int, default=4)
    parser.add_argument("--top", type=int, default=25)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/signed_ratio_drop_breaker_2026-07-04.json"),
    )
    args = parser.parse_args()

    rng = random.Random(args.seed)
    rows: list[dict[str, Any]] = []
    seen: set[tuple[Any, ...]] = set()
    counters = {"attempted": 0, "duplicates": 0, "not_analyzed": 0, "analyzed": 0}

    add_finite_skellam_boundary(rows, seen, counters, max_side_n=args.max_side_n)
    add_near_one_sided_perturbations(rows, seen, counters, max_side_n=args.max_side_n)
    add_boundary_heavy_grid(rows, seen, counters, max_side_n=args.max_side_n)
    add_random_large_grouped(
        rows,
        seen,
        counters,
        rng=rng,
        samples=args.random_samples,
        max_groups=args.max_groups,
    )
    add_optimizer_rows(
        rows,
        seen,
        counters,
        rng=rng,
        de_maxiter=args.de_maxiter,
        de_restarts=args.de_restarts,
        de_popsize=args.de_popsize,
        variance_cutoff=args.variance_cutoff,
    )

    feasible = [
        row
        for row in rows
        if row["variance"] >= args.variance_cutoff
        and math.isfinite(row["variance_times_effective_ratio_drop"])
        and math.isfinite(row["variance_times_reserve"])
    ]
    effective_failures = [
        row
        for row in feasible
        if row["variance_times_effective_ratio_drop"] < args.candidate_constant - TOL
    ]
    reserve_failures = [
        row
        for row in feasible
        if row["variance_times_reserve"] < args.candidate_constant - TOL
    ]
    side_floors = [0.0, 0.02, 0.05, 0.1, 0.25]
    summary = {
        "source": {
            "kind": "signed_ratio_drop_breaker",
            "seed": args.seed,
            "candidate_constant": args.candidate_constant,
            "variance_cutoff": args.variance_cutoff,
            "max_side_n": args.max_side_n,
            "random_samples": args.random_samples,
            "max_groups": args.max_groups,
            "de_maxiter": args.de_maxiter,
            "de_restarts": args.de_restarts,
            "de_popsize": args.de_popsize,
        },
        "counts": {
            **counters,
            "feasible": len(feasible),
            "effective_candidate_failures": len(effective_failures),
            "reserve_candidate_failures": len(reserve_failures),
        },
        "best_by_effective_ratio_drop": [
            compact_breaker(row)
            for row in sorted(feasible, key=lambda row: row["variance_times_effective_ratio_drop"])[
                : args.top
            ]
        ],
        "best_by_raw_reserve": [
            compact_breaker(row)
            for row in sorted(feasible, key=lambda row: row["variance_times_reserve"])[: args.top]
        ],
        "best_effective_by_variance_cutoff": best_by_cutoff(
            feasible,
            [1, 1.25, 1.5, 2, 5, 10, 20, 50],
            key="variance_times_effective_ratio_drop",
        ),
        "best_reserve_by_variance_cutoff": best_by_cutoff(
            feasible,
            [1, 1.25, 1.5, 2, 5, 10, 20, 50],
            key="variance_times_reserve",
        ),
        "best_effective_by_side_balance": [
            {
                "side_variance_fraction_floor": floor,
                "best": compact_breaker(
                    min(
                        (
                            row
                            for row in feasible
                            if row["side_variance_fraction"] >= floor
                        ),
                        key=lambda row: row["variance_times_effective_ratio_drop"],
                        default=None,
                    )
                )
                if any(row["side_variance_fraction"] >= floor for row in feasible)
                else None,
            }
            for floor in side_floors
        ],
        "effective_failure_certificates": [
            failure_certificate(row)
            for row in sorted(effective_failures, key=lambda row: row["variance_times_effective_ratio_drop"])[
                : args.top
            ]
        ],
        "reserve_failure_certificates": [
            failure_certificate(row)
            for row in sorted(reserve_failures, key=lambda row: row["variance_times_reserve"])[
                : args.top
            ]
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Wrote {args.out}; attempted={counters['attempted']}, "
        f"analyzed={counters['analyzed']}, feasible={len(feasible)}, "
        f"effective_failures={len(effective_failures)}, "
        f"reserve_failures={len(reserve_failures)}",
        flush=True,
    )
    return 1 if effective_failures or reserve_failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
