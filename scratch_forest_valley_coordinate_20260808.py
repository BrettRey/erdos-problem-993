#!/usr/bin/env python3
"""Coordinate search around the strongest forest-product near valleys.

The semigroup beam search supplies multisets of non-log-concave tree
polynomials whose product has a shallow descent followed by a later rebound.
This program performs exact-factor coordinate moves (replace, add, or delete
one component) while ranking with long-double convolution.  Any apparent
non-unimodality is immediately replayed with arbitrary-precision integers.

A witness is a forest counterexample to the Alavi--Malde--Schwenk--Erdos
unimodality conjecture.  A negative run is only computational evidence.
"""

from __future__ import annotations

import argparse
import json
import random
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

from scratch_forest_valley_semigroup_20260807 import (
    State,
    bush_library,
    exact_record,
    normalize,
    rebound_score,
)


def product_profile(
    factors: tuple[int, ...], profiles: list[np.ndarray],
) -> np.ndarray:
    out = np.asarray([1.0], dtype=np.longdouble)
    for factor in factors:
        out = normalize(np.convolve(out, profiles[factor]))
    return out


def score_profile(
    profile: np.ndarray, threshold: float, objective: str,
    min_separation: int,
) -> tuple[float, int, int, bool]:
    return rebound_score(
        profile, threshold, objective, min_separation,
    )


def seed_factor_sets(paths: list[Path], limit: int) -> list[tuple[int, ...]]:
    rows: list[tuple[float, tuple[int, ...]]] = []
    for path in paths:
        data = json.loads(path.read_text())
        for record in data.get("top_exact", []):
            factors = tuple(sorted(int(value) for value in record["factors"]))
            score = float(record.get("thresholded_rebound_ratio", 0.0))
            rows.append((score, factors))
    rows.sort(reverse=True)
    distinct: list[tuple[int, ...]] = []
    seen: set[tuple[int, ...]] = set()
    for _, factors in rows:
        if factors in seen:
            continue
        seen.add(factors)
        distinct.append(factors)
        if len(distinct) >= limit:
            break
    if not distinct:
        raise ValueError("no top_exact factor sets found in the inputs")
    return distinct


def other_profiles(
    factors: tuple[int, ...], profiles: list[np.ndarray],
) -> list[np.ndarray]:
    prefix = [np.asarray([1.0], dtype=np.longdouble)]
    for factor in factors:
        prefix.append(normalize(np.convolve(prefix[-1], profiles[factor])))
    suffix = [np.asarray([1.0], dtype=np.longdouble)] * (len(factors) + 1)
    for index in range(len(factors) - 1, -1, -1):
        suffix[index] = normalize(np.convolve(
            profiles[factors[index]], suffix[index + 1],
        ))
    return [
        normalize(np.convolve(prefix[index], suffix[index + 1]))
        for index in range(len(factors))
    ]


def encoded_state(
    factors: tuple[int, ...], profile: np.ndarray, score_data,
    seeds,
) -> dict[str, object]:
    score, dip, rise, witness = score_data
    return {
        "factors": list(factors),
        "labels": [seeds[index].label for index in factors],
        "component_orders": [seeds[index].order for index in factors],
        "order": sum(seeds[index].order for index in factors),
        "factor_count": len(factors),
        "float_score": score,
        "float_dip": dip,
        "float_rise": rise,
        "float_witness": witness,
        "profile_length": len(profile),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--inputs", type=Path, nargs="+",
        default=[Path("results/forest_valley_semigroup_20260807.json")],
    )
    parser.add_argument("--starts", type=int, default=8)
    parser.add_argument("--sweeps", type=int, default=20)
    parser.add_argument("--threshold", type=float, default=0.001)
    parser.add_argument(
        "--objective",
        choices=("suffix_value", "later_slope", "slope_rebound"),
        default="suffix_value",
    )
    parser.add_argument("--min-separation", type=int, default=1)
    parser.add_argument("--min-factors", type=int, default=2)
    parser.add_argument("--max-factors", type=int, default=40)
    parser.add_argument(
        "--library-sample", type=int, default=0,
        help="0 uses all seeds; otherwise sample this many per sweep",
    )
    parser.add_argument("--seed", type=int, default=202608086)
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/forest_valley_coordinate_20260808.json"),
    )
    args = parser.parse_args()
    started = time.time()
    rng = random.Random(args.seed)
    exact_threshold = Fraction(str(args.threshold))
    seeds = bush_library()
    profiles = [normalize(seed.poly) for seed in seeds]
    starts = seed_factor_sets(args.inputs, args.starts)
    print(json.dumps({
        "event": "library",
        "seeds": len(seeds),
        "starts": len(starts),
        "elapsed_seconds": time.time() - started,
    }), flush=True)

    leaders: list[dict[str, object]] = []
    witness = None
    evaluated = 0
    for restart, initial in enumerate(starts):
        factors = initial
        profile = product_profile(factors, profiles)
        current = score_profile(
            profile, args.threshold, args.objective, args.min_separation,
        )
        for sweep in range(args.sweeps):
            if args.library_sample and args.library_sample < len(seeds):
                library = rng.sample(range(len(seeds)), args.library_sample)
                library.extend(factors)
                library = sorted(set(library))
            else:
                library = range(len(seeds))

            best_score = current[0]
            best_factors = factors
            best_profile = profile
            best_data = current
            seen: set[tuple[int, ...]] = {factors}

            # Replacement coordinates use prefix/suffix products, avoiding a
            # full product rebuild for each of the 4,445 seed choices.
            for position, other in enumerate(other_profiles(factors, profiles)):
                for candidate in library:
                    if candidate == factors[position]:
                        continue
                    moved = tuple(sorted(
                        factors[:position] + factors[position + 1:] + (candidate,)
                    ))
                    if moved in seen:
                        continue
                    seen.add(moved)
                    candidate_profile = normalize(np.convolve(other, profiles[candidate]))
                    data = score_profile(
                        candidate_profile, args.threshold, args.objective,
                        args.min_separation,
                    )
                    evaluated += 1
                    if data[3]:
                        exact = exact_record(
                            moved, seeds, exact_threshold, args.objective,
                            args.min_separation,
                        )
                        if not exact["unimodal"]:
                            witness = exact
                            best_factors, best_profile, best_data = moved, candidate_profile, data
                            break
                    if data[0] > best_score + 1e-18:
                        best_score = data[0]
                        best_factors, best_profile, best_data = moved, candidate_profile, data
                if witness is not None:
                    break
            if witness is not None:
                factors, profile, current = best_factors, best_profile, best_data
                break

            # Addition and deletion let the optimizer find its own product
            # length instead of being trapped at the beam seed's depth.
            if len(factors) < args.max_factors:
                for candidate in library:
                    moved = tuple(sorted(factors + (candidate,)))
                    candidate_profile = normalize(np.convolve(profile, profiles[candidate]))
                    data = score_profile(
                        candidate_profile, args.threshold, args.objective,
                        args.min_separation,
                    )
                    evaluated += 1
                    if data[3]:
                        exact = exact_record(
                            moved, seeds, exact_threshold, args.objective,
                            args.min_separation,
                        )
                        if not exact["unimodal"]:
                            witness = exact
                            best_factors, best_profile, best_data = moved, candidate_profile, data
                            break
                    if data[0] > best_score + 1e-18:
                        best_score = data[0]
                        best_factors, best_profile, best_data = moved, candidate_profile, data
            if witness is None and len(factors) > args.min_factors:
                for position, other in enumerate(other_profiles(factors, profiles)):
                    moved = factors[:position] + factors[position + 1:]
                    if moved in seen:
                        continue
                    data = score_profile(
                        other, args.threshold, args.objective,
                        args.min_separation,
                    )
                    evaluated += 1
                    if data[0] > best_score + 1e-18:
                        best_score = data[0]
                        best_factors, best_profile, best_data = moved, other, data

            improved = best_factors != factors
            factors, profile, current = best_factors, best_profile, best_data
            print(json.dumps({
                "event": "sweep",
                "restart": restart,
                "sweep": sweep,
                "evaluated": evaluated,
                "factor_count": len(factors),
                "score": current[0],
                "dip": current[1],
                "rise": current[2],
                "improved": improved,
                "witness": witness is not None,
                "elapsed_seconds": time.time() - started,
            }), flush=True)
            if witness is not None or not improved:
                break

        exact = exact_record(
            factors, seeds, exact_threshold, args.objective,
            args.min_separation,
        )
        leaders.append({
            **encoded_state(factors, profile, current, seeds),
            "exact": exact,
        })
        if witness is not None:
            break

    leaders.sort(
        key=lambda row: row["exact"]["objective_ratio"], reverse=True,
    )
    report = {
        "claim_scope": "exact-verified computational search; no witness is not a proof",
        "inputs": [str(path) for path in args.inputs],
        "seed": args.seed,
        "threshold": str(exact_threshold),
        "objective": args.objective,
        "min_separation": args.min_separation,
        "library_size": len(seeds),
        "starts": len(starts),
        "sweeps": args.sweeps,
        "library_sample": args.library_sample,
        "evaluated": evaluated,
        "counterexample": witness,
        "leaders": leaders,
        "elapsed_seconds": time.time() - started,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "counterexample": witness is not None,
        "evaluated": evaluated,
        "best_exact_score": leaders[0]["exact"]["objective_ratio"],
        "out": str(args.out),
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    raise SystemExit(1 if witness is not None else 0)


if __name__ == "__main__":
    main()
