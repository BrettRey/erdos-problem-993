#!/usr/bin/env python3
"""Direct-valley semigroup search over all 4,445 Blair bush seeds.

The public Blair sweep generated every non-log-concave bush seed but mixed only
the 80 most severe seeds.  Severity of a log-concavity defect is not a direct
proxy for non-unimodality.  This script samples the full pair surface, retains a
diverse beam under the actual post-first-descent slope, and grows products to
larger forests.  Floating point is used only to rank candidates.  Every apparent
ascent and every reported finalist is replayed with Python integers.
"""

from __future__ import annotations

import argparse
import heapq
import itertools
import json
import math
import random
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

from scripts.valley_search import _polymul, bouquet_poly


@dataclass(frozen=True)
class Seed:
    label: str
    order: int
    poly: tuple[int, ...]


@dataclass
class State:
    score: float
    factors: tuple[int, ...]
    profile: np.ndarray
    dip_pos: int
    rise_pos: int


def is_lc(poly: tuple[int, ...]) -> bool:
    return all(poly[k] * poly[k] >= poly[k - 1] * poly[k + 1]
               for k in range(1, len(poly) - 1))


def is_unimodal(poly: list[int] | tuple[int, ...]) -> bool:
    down = False
    for left, right in zip(poly, poly[1:]):
        if right < left:
            down = True
        elif down and right > left:
            return False
    return True


def bush_library(order_cap: int = 60) -> list[Seed]:
    btypes = [
        (kind, children, pendant)
        for children in range(2, 7)
        for pendant in range(1, 4)
        for kind in ("u", "s")
    ]

    distinct: dict[tuple[int, ...], Seed] = {}
    for root_degree in range(2, 6):
        for combo in itertools.combinations_with_replacement(range(len(btypes)), root_degree):
            gadgets: list[tuple[int, ...]] = []
            tokens: list[str] = []
            order = 1
            for index in combo:
                kind, children, pendant = btypes[index]
                legs = [pendant + 1] * children
                if kind == "s":
                    legs[0] += 1
                gadgets.append(tuple(sorted(legs)))
                tokens.append(f"{kind}{children}p{pendant}")
                order += 1 + sum(legs)
            if order > order_cap:
                continue
            poly = tuple(bouquet_poly(gadgets))
            if is_lc(poly):
                continue
            label = f"bush({','.join(tokens)})|n={order}"
            distinct.setdefault(poly, Seed(label, order, poly))
    return list(distinct.values())


def normalize(poly: tuple[int, ...] | np.ndarray) -> np.ndarray:
    arr = np.asarray(poly, dtype=np.longdouble)
    maximum = np.max(arr)
    return arr / maximum if maximum else arr


def rebound_score(
    profile: np.ndarray, threshold: float, objective: str, min_separation: int
) -> tuple[float, int, int, bool]:
    """Thresholded suffix rebound plus threshold-free nonunimodality flag."""
    down = False
    witness = False
    for left, right in zip(profile, profile[1:]):
        if right < left:
            down = True
        elif down and right > left * (1.0 + 1e-14):
            witness = True

    prefix = np.maximum.accumulate(profile)
    if objective == "slope_rebound":
        first_descent = next(
            (index for index in range(len(profile) - 1)
             if profile[index + 1] < profile[index]),
            -1,
        )
        if first_descent < 0:
            return 0.0, -1, -1, witness
        alpha = len(profile) - 1
        tail_start = math.ceil((2 * alpha - 1) / 3)
        min_num = profile[first_descent + 1]
        min_den = profile[first_descent]
        min_pos = first_descent
        best = 0.0
        best_b = -1
        best_c = -1
        for c in range(first_descent + 1, min(tail_start, len(profile) - 1)):
            numerator = profile[c + 1]
            denominator = profile[c]
            if (
                c - min_pos >= min_separation
                and numerator * min_den > min_num * denominator * (1.0 + 1e-14)
            ):
                ratio = float(numerator / denominator)
                if ratio > best:
                    best, best_b, best_c = ratio, min_pos, c
            if numerator * min_den < min_num * denominator:
                min_num, min_den, min_pos = numerator, denominator, c
        return best, best_b, best_c, witness

    if objective == "later_slope":
        eligible = [False] * len(profile)
        for b in range(1, len(profile) - 1):
            # Repeated normalized convolutions can underflow the extreme
            # low-degree tail to zero.  The exact polynomial is positive
            # there, so a numerical 0 >= (1+t) * 0 is not a descent.
            eligible[b] = (
                profile[b] > 0
                and prefix[b - 1] >= (1.0 + threshold) * profile[b]
            )
        alpha = len(profile) - 1
        tail_start = math.ceil((2 * alpha - 1) / 3)
        best = 0.0
        best_b = -1
        best_c = -1
        latest_eligible = -1
        for c in range(1, min(tail_start, len(profile) - 1)):
            add_b = c - min_separation
            if add_b >= 1 and eligible[add_b]:
                latest_eligible = add_b
            if latest_eligible < 0 or profile[c] == 0:
                continue
            ratio = float(profile[c + 1] / profile[c])
            if ratio > best:
                best = ratio
                best_b = latest_eligible
                best_c = c
        return best, best_b, best_c, witness

    suffix = np.maximum.accumulate(profile[::-1])[::-1]
    best = 0.0
    best_b = -1
    best_c = -1
    for b in range(1, len(profile) - 1):
        if profile[b] == 0 or prefix[b - 1] < (1.0 + threshold) * profile[b]:
            continue
        ratio = float(suffix[b + 1] / profile[b])
        if ratio > best:
            best = ratio
            best_b = b
            target = suffix[b + 1]
            best_c = next(c for c in range(b + 1, len(profile)) if profile[c] == target)
    return best, best_b, best_c, witness


def exact_product(factors: tuple[int, ...], seeds: list[Seed]) -> list[int]:
    out = [1]
    for index in factors:
        out = _polymul(out, list(seeds[index].poly))
    return out


def exact_record(
    factors: tuple[int, ...], seeds: list[Seed], threshold: Fraction,
    objective: str, min_separation: int,
) -> dict:
    poly = exact_product(factors, seeds)
    first = next((k for k in range(len(poly) - 1) if poly[k + 1] < poly[k]), -1)
    best_num = 0
    best_den = 1
    best_k = -1
    witness_k = -1
    if first >= 0:
        for k in range(first, len(poly) - 1):
            if poly[k + 1] * best_den > best_num * poly[k]:
                best_num, best_den, best_k = poly[k + 1], poly[k], k
            if poly[k + 1] > poly[k] and witness_k < 0:
                witness_k = k
    prefix = poly[0]
    threshold_num = threshold.numerator
    threshold_den = threshold.denominator
    threshold_best_num = 0
    threshold_best_den = 1
    threshold_b = -1
    threshold_c = -1
    suffix = [0] * len(poly)
    running = 0
    for k in range(len(poly) - 1, -1, -1):
        suffix[k] = running
        running = max(running, poly[k])
    for b in range(1, len(poly) - 1):
        if prefix * threshold_den >= (threshold_den + threshold_num) * poly[b]:
            if suffix[b] * threshold_best_den > threshold_best_num * poly[b]:
                threshold_best_num = suffix[b]
                threshold_best_den = poly[b]
                threshold_b = b
                threshold_c = next(c for c in range(b + 1, len(poly)) if poly[c] == suffix[b])
        prefix = max(prefix, poly[b])
    eligible = [False] * len(poly)
    prefix = poly[0]
    for b in range(1, len(poly) - 1):
        eligible[b] = prefix * threshold_den >= (threshold_den + threshold_num) * poly[b]
        prefix = max(prefix, poly[b])
    alpha = len(poly) - 1
    tail_start = math.ceil((2 * alpha - 1) / 3)
    slope_num = 0
    slope_den = 1
    slope_b = -1
    slope_c = -1
    latest_eligible = -1
    for c in range(1, min(tail_start, len(poly) - 1)):
        add_b = c - min_separation
        if add_b >= 1 and eligible[add_b]:
            latest_eligible = add_b
        if latest_eligible >= 0 and poly[c + 1] * slope_den > slope_num * poly[c]:
            slope_num, slope_den = poly[c + 1], poly[c]
            slope_b, slope_c = latest_eligible, c
    rebound_num = 0
    rebound_den = 1
    rebound_b = -1
    rebound_c = -1
    if first >= 0:
        min_num = poly[first + 1]
        min_den = poly[first]
        min_pos = first
        for c in range(first + 1, min(tail_start, len(poly) - 1)):
            numerator = poly[c + 1]
            denominator = poly[c]
            if (
                c - min_pos >= min_separation
                and numerator * min_den > min_num * denominator
                and numerator * rebound_den > rebound_num * denominator
            ):
                rebound_num, rebound_den = numerator, denominator
                rebound_b, rebound_c = min_pos, c
            if numerator * min_den < min_num * denominator:
                min_num, min_den, min_pos = numerator, denominator, c
    if objective == "slope_rebound":
        objective_num, objective_den = rebound_num, rebound_den
        objective_b, objective_c = rebound_b, rebound_c
    elif objective == "later_slope":
        objective_num, objective_den = slope_num, slope_den
        objective_b, objective_c = slope_b, slope_c
    else:
        objective_num, objective_den = threshold_best_num, threshold_best_den
        objective_b, objective_c = threshold_b, threshold_c
    return {
        "factors": list(factors),
        "labels": [seeds[i].label for i in factors],
        "component_orders": [seeds[i].order for i in factors],
        "order": sum(seeds[i].order for i in factors),
        "alpha": len(poly) - 1,
        "unimodal": is_unimodal(poly),
        "first_descent_transition": first,
        "best_post_descent_transition": best_k,
        "best_ratio_numerator": best_num,
        "best_ratio_denominator": best_den,
        "best_ratio": best_num / best_den if best_den else 0.0,
        "witness_transition": witness_k,
        "threshold": str(threshold),
        "thresholded_rebound_dip": threshold_b,
        "thresholded_rebound_rise": threshold_c,
        "thresholded_rebound_ratio_numerator": threshold_best_num,
        "thresholded_rebound_ratio_denominator": threshold_best_den,
        "thresholded_rebound_ratio": threshold_best_num / threshold_best_den,
        "later_slope_dip": slope_b,
        "later_slope_transition": slope_c,
        "later_slope_ratio_numerator": slope_num,
        "later_slope_ratio_denominator": slope_den,
        "later_slope_ratio": slope_num / slope_den,
        "slope_rebound_dip": rebound_b,
        "slope_rebound_transition": rebound_c,
        "slope_rebound_ratio_numerator": rebound_num,
        "slope_rebound_ratio_denominator": rebound_den,
        "slope_rebound_ratio": rebound_num / rebound_den,
        "objective": objective,
        "objective_dip": objective_b,
        "objective_position": objective_c,
        "objective_ratio_numerator": objective_num,
        "objective_ratio_denominator": objective_den,
        "objective_ratio": objective_num / objective_den,
        "polynomial": poly if witness_k >= 0 else None,
    }


def retain(states: list[State], beam: int) -> list[State]:
    """Keep global leaders plus leaders in coarse location/length bins."""
    by_key: dict[tuple[int, int, int], State] = {}
    for state in states:
        key = (len(state.factors), len(state.profile) // 8,
               state.rise_pos - state.dip_pos)
        old = by_key.get(key)
        if old is None or state.score > old.score:
            by_key[key] = state
    diverse = sorted(by_key.values(), key=lambda state: state.score, reverse=True)
    global_best = heapq.nlargest(beam, states, key=lambda state: state.score)
    merged: dict[tuple[int, ...], State] = {}
    for state in global_best + diverse[:beam]:
        old = merged.get(state.factors)
        if old is None or state.score > old.score:
            merged[state.factors] = state
    return sorted(merged.values(), key=lambda state: state.score, reverse=True)[:beam]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-samples", type=int, default=5_000_000)
    parser.add_argument(
        "--start-singles", action="store_true",
        help="seed the beam with individual non-LC bushes instead of random pairs",
    )
    parser.add_argument("--beam", type=int, default=1200)
    parser.add_argument("--children-per-state", type=int, default=240)
    parser.add_argument("--max-factors", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260807)
    parser.add_argument("--threshold", type=float, default=0.001)
    parser.add_argument("--objective", choices=(
        "suffix_value", "later_slope", "slope_rebound",
    ),
                        default="later_slope")
    parser.add_argument("--min-separation", type=int, default=2)
    parser.add_argument("--out", type=Path,
                        default=Path("results/forest_valley_semigroup_20260807.json"))
    args = parser.parse_args()

    rng = random.Random(args.seed)
    exact_threshold = Fraction(str(args.threshold))
    start = time.time()
    seeds = bush_library()
    profiles = [normalize(seed.poly) for seed in seeds]
    print(json.dumps({"event": "library", "count": len(seeds),
                      "seconds": time.time() - start}), flush=True)

    candidate_heap: list[tuple[float, int, State]] = []
    pair_bins: dict[tuple[int, int], State] = {}
    serial = 0
    exact_witness: dict | None = None
    for sample in range(1, args.pair_samples + 1):
        left = rng.randrange(len(seeds))
        right = rng.randrange(len(seeds))
        factors = tuple(sorted((left, right)))
        profile = normalize(np.convolve(profiles[left], profiles[right]))
        score, dip, rise, float_witness = rebound_score(
            profile, args.threshold, args.objective, args.min_separation)
        state = State(score, factors, profile, dip, rise)
        serial += 1
        entry = (score, serial, state)
        if len(candidate_heap) < args.beam * 4:
            heapq.heappush(candidate_heap, entry)
        elif score > candidate_heap[0][0]:
            heapq.heapreplace(candidate_heap, entry)
        bin_key = (len(profile) // 8, rise - dip)
        old = pair_bins.get(bin_key)
        if old is None or score > old.score:
            pair_bins[bin_key] = state
        if float_witness:
            record = exact_record(factors, seeds, exact_threshold,
                                  args.objective, args.min_separation)
            if not record["unimodal"]:
                exact_witness = record
                break
        if sample % 1_000_000 == 0:
            candidates = [entry[2] for entry in candidate_heap] + list(pair_bins.values())
            leaders = retain(candidates, args.beam * 2)
            print(json.dumps({"event": "pairs", "sample": sample,
                              "best": leaders[0].score,
                              "elapsed": time.time() - start}), flush=True)

    if args.start_singles:
        candidates = []
        for index, profile in enumerate(profiles):
            score, dip, rise, _ = rebound_score(
                profile, args.threshold, args.objective, args.min_separation,
            )
            candidates.append(State(score, (index,), profile, dip, rise))
        depth_start = 2
    else:
        candidates = [entry[2] for entry in candidate_heap] + list(pair_bins.values())
        depth_start = 3
    beam_states = retain(candidates, args.beam)
    depth_summaries = []
    if exact_witness is None:
        for depth in range(depth_start, args.max_factors + 1):
            grown: list[State] = []
            seen: set[tuple[int, ...]] = set()
            for parent in beam_states:
                for _ in range(args.children_per_state):
                    child = rng.randrange(len(seeds))
                    factors = tuple(sorted(parent.factors + (child,)))
                    if factors in seen:
                        continue
                    seen.add(factors)
                    profile = normalize(np.convolve(parent.profile, profiles[child]))
                    score, dip, rise, float_witness = rebound_score(
                        profile, args.threshold, args.objective, args.min_separation)
                    state = State(score, factors, profile, dip, rise)
                    grown.append(state)
                    if float_witness:
                        record = exact_record(factors, seeds, exact_threshold,
                                              args.objective, args.min_separation)
                        if not record["unimodal"]:
                            exact_witness = record
                            break
                if exact_witness is not None:
                    break
            if exact_witness is not None:
                break
            beam_states = retain(grown, args.beam)
            summary = {"depth": depth, "tested": len(grown),
                       "best_float_score": beam_states[0].score,
                       "elapsed": time.time() - start}
            depth_summaries.append(summary)
            print(json.dumps({"event": "depth", **summary}), flush=True)

    finalists = []
    if exact_witness is None:
        finalists = [exact_record(state.factors, seeds, exact_threshold,
                                   args.objective, args.min_separation)
                      for state in beam_states[:100]]
        finalists.sort(key=lambda row: row["objective_ratio"], reverse=True)
        exact_witness = next((row for row in finalists if not row["unimodal"]), None)

    payload = {
        "claim_scope": "computational search; not a proof if no witness is found",
        "seed": args.seed,
        "descent_threshold": args.threshold,
        "objective": args.objective,
        "min_separation": args.min_separation,
        "library_size": len(seeds),
        "pair_samples_requested": args.pair_samples,
        "start_singles": args.start_singles,
        "beam": args.beam,
        "children_per_state": args.children_per_state,
        "max_factors": args.max_factors,
        "elapsed_seconds": time.time() - start,
        "witness": exact_witness,
        "depth_summaries": depth_summaries,
        "top_exact": finalists[:25],
    }
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"event": "done", "witness": exact_witness is not None,
                      "best": finalists[0]["objective_ratio"] if finalists else None,
                      "out": str(args.out), "elapsed": time.time() - start}), flush=True)


if __name__ == "__main__":
    main()
