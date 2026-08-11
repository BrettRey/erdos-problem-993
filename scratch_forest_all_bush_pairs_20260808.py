#!/usr/bin/env python3
"""Exhaustively screen every unordered pair of Blair bush polynomials.

Every one of the 4,445 distinct inputs is the independence polynomial of a
tree and is non-log-concave.  Their products are independence polynomials of
two-component forests.  A floating FFT screens all unordered pairs; every
apparent coefficient valley and every numerical leader is replayed by exact
integer convolution.

The FFT screen deliberately inspects only coefficients at least ``floor``
times the row maximum.  This avoids treating roundoff in the extreme tails as
signal.  Thus a negative run is an exhaustive screen at the stated floor, not
a proof about all exact coefficients.  Reversing the inputs gives the same
absolute scale, so it does not cure FFT cancellation in those tails.
"""

from __future__ import annotations

import argparse
import heapq
import json
import time
from pathlib import Path

import numpy as np

from scratch_forest_valley_semigroup_20260807 import (
    bush_library,
    exact_record,
    is_unimodal,
)
from scripts.valley_search import _polymul


def exact_pair(left: int, right: int, seeds) -> dict[str, object]:
    poly = _polymul(list(seeds[left].poly), list(seeds[right].poly))
    first = next(
        (k for k in range(len(poly) - 1) if poly[k + 1] < poly[k]), -1,
    )
    best_num, best_den, best_k = 0, 1, -1
    witness = -1
    if first >= 0:
        for k in range(first + 1, len(poly) - 1):
            if poly[k + 1] * best_den > best_num * poly[k]:
                best_num, best_den, best_k = poly[k + 1], poly[k], k
            if poly[k + 1] > poly[k] and witness < 0:
                witness = k
    return {
        "factors": [left, right],
        "labels": [seeds[left].label, seeds[right].label],
        "component_orders": [seeds[left].order, seeds[right].order],
        "order": seeds[left].order + seeds[right].order,
        "unimodal": is_unimodal(poly),
        "first_descent_transition": first,
        "best_post_descent_transition": best_k,
        "best_ratio_numerator": best_num,
        "best_ratio_denominator": best_den,
        "best_ratio": best_num / best_den,
        "witness_transition": witness,
        "polynomial": poly if witness >= 0 else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--floor", type=float, default=1e-11)
    parser.add_argument("--tolerance", type=float, default=2e-10)
    parser.add_argument("--chunk", type=int, default=512)
    parser.add_argument("--leaders", type=int, default=200)
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/forest_all_bush_pairs_20260808.json"),
    )
    args = parser.parse_args()
    started = time.time()
    seeds = bush_library()
    count = len(seeds)
    max_coefficients = max(len(seed.poly) for seed in seeds)
    convolution_length = 2 * max_coefficients - 1
    fft_length = 1 << (convolution_length - 1).bit_length()

    padded = np.zeros((count, fft_length), dtype=np.float64)
    lengths = np.empty(count, dtype=np.int16)
    for index, seed in enumerate(seeds):
        profile = np.asarray(seed.poly, dtype=np.float64)
        profile /= profile.max()
        padded[index, : len(profile)] = profile
        lengths[index] = len(profile)
    transforms = np.fft.rfft(padded, axis=1)
    del padded

    # Min-heap entries are (score, serial, left, right, float positions).
    leaders: list[tuple[float, int, int, int, int, int]] = []
    serial = 0
    tested = 0
    numerical_witnesses = 0
    exact_witness = None

    for left in range(count):
        for start in range(left, count, args.chunk):
            stop = min(count, start + args.chunk)
            rows = np.fft.irfft(
                transforms[start:stop] * transforms[left],
                n=fft_length,
                axis=1,
            )[:, :convolution_length]
            row_max = rows.max(axis=1, keepdims=True)
            valid = rows >= args.floor * row_max
            denominator = rows[:, :-1]
            numerator = rows[:, 1:]
            valid_edge = valid[:, :-1] & valid[:, 1:] & (denominator > 0)
            ratios = np.full_like(denominator, -np.inf)
            np.divide(numerator, denominator, out=ratios, where=valid_edge)

            for offset, right in enumerate(range(start, stop)):
                degree = int(lengths[left] + lengths[right] - 2)
                row_ratios = ratios[offset, :degree]
                descent_positions = np.flatnonzero(
                    np.isfinite(row_ratios)
                    & (row_ratios < 1.0 - args.tolerance),
                )
                tested += 1
                if not len(descent_positions):
                    continue
                first = int(descent_positions[0])
                later = row_ratios[first + 1:]
                if not len(later):
                    continue
                relative = int(np.argmax(later))
                position = first + 1 + relative
                score = float(later[relative])
                if score > 1.0 + args.tolerance:
                    numerical_witnesses += 1
                    exact = exact_pair(left, right, seeds)
                    if not exact["unimodal"]:
                        exact_witness = exact
                        break
                serial += 1
                entry = (score, serial, left, right, first, position)
                if len(leaders) < args.leaders:
                    heapq.heappush(leaders, entry)
                elif score > leaders[0][0]:
                    heapq.heapreplace(leaders, entry)
            if exact_witness is not None:
                break
        if exact_witness is not None:
            break
        if (left + 1) % 250 == 0 or left + 1 == count:
            print(json.dumps({
                "event": "progress",
                "left": left + 1,
                "library_size": count,
                "tested": tested,
                "best_float_score": max(leaders)[0] if leaders else None,
                "numerical_witnesses": numerical_witnesses,
                "elapsed_seconds": time.time() - started,
            }), flush=True)

    exact_leaders = []
    for _, _, left, right, _, _ in sorted(leaders, reverse=True):
        exact_leaders.append(exact_pair(left, right, seeds))
    exact_leaders.sort(key=lambda row: row["best_ratio"], reverse=True)
    if exact_witness is None:
        exact_witness = next(
            (row for row in exact_leaders if not row["unimodal"]), None,
        )

    report = {
        "claim_scope": (
            "all unordered pairs FFT-screened at the stated relative floor; "
            "all leaders and apparent valleys replayed exactly"
        ),
        "library_size": count,
        "unordered_pairs": count * (count + 1) // 2,
        "tested": tested,
        "relative_floor": args.floor,
        "tolerance": args.tolerance,
        "fft_length": fft_length,
        "numerical_witnesses": numerical_witnesses,
        "counterexample": exact_witness,
        "top_exact": exact_leaders,
        "elapsed_seconds": time.time() - started,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "tested": tested,
        "counterexample": exact_witness is not None,
        "best_exact_ratio": exact_leaders[0]["best_ratio"]
        if exact_leaders else None,
        "elapsed_seconds": report["elapsed_seconds"],
        "out": str(args.out),
    }), flush=True)


if __name__ == "__main__":
    main()
