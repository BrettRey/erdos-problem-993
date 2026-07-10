#!/usr/bin/env python3
"""High-precision verification of the Skellam effective-drop theorem.

This is a falsification/verification harness, not the proof.  It scans a
variance/imbalance grid and refines every observed mode-transition curve.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import mpmath as mp


DEFAULT_VARIANCES = "1,1.000001,1.001,1.01,1.1,1.5,2,3,5,10,30,100"


def parse_values(raw: str) -> list[mp.mpf]:
    return [mp.mpf(part.strip()) for part in raw.split(",") if part.strip()]


def log_skellam_pmf(z: int, variance: mp.mpf, mean: mp.mpf) -> mp.mpf:
    lam = (variance + mean) / 2
    eta = (variance - mean) / 2
    if lam == 0:
        if z > 0:
            return mp.ninf
        count = -z
        return -eta + count * mp.log(eta) - mp.loggamma(count + 1)
    if eta == 0:
        if z < 0:
            return mp.ninf
        return -lam + z * mp.log(lam) - mp.loggamma(z + 1)
    x = 2 * mp.sqrt(lam * eta)
    return (
        -variance
        + mp.mpf(z) * mp.log(lam / eta) / 2
        + mp.log(mp.besseli(abs(z), x))
    )


def largest_mode(variance: mp.mpf, mean: mp.mpf) -> int:
    radius = int(mp.ceil(8 * mp.sqrt(variance) + 12))
    lo = int(mp.floor(mean)) - radius
    hi = int(mp.ceil(mean)) + radius
    values = [(log_skellam_pmf(z, variance, mean), z) for z in range(lo, hi + 1)]
    best = max(value for value, _ in values)
    tolerance = mp.power(10, -(mp.mp.dps - 15))
    tied = [z for value, z in values if abs(value - best) <= tolerance]
    if not tied:
        raise AssertionError("failed to locate a Skellam mode")
    return max(tied)


def effective_drop(
    variance: mp.mpf,
    mean: mp.mpf,
    descent: int,
) -> mp.mpf:
    log_ratio_product = (
        log_skellam_pmf(descent + 1, variance, mean)
        + log_skellam_pmf(descent - 1, variance, mean)
        - 2 * log_skellam_pmf(descent, variance, mean)
    )
    return -mp.expm1(log_ratio_product)


def baricz_lower_bound(variance: mp.mpf, mean: mp.mpf, descent: int) -> mp.mpf:
    # Baricz's theorem assumes x>0, hence both Poisson rates are positive.
    # The two zero-rate boundaries are checked directly by the quarter bound.
    if abs(mean) == variance:
        return mp.mpf(0)
    x = mp.sqrt(variance * variance - mean * mean)
    nu = mp.mpf(descent)
    a = nu + mp.mpf("0.5")
    return (a / (nu + 1)) / mp.sqrt(x * x + a * a)


def analyze_point(variance: mp.mpf, imbalance: mp.mpf) -> dict[str, Any]:
    mean = variance * imbalance
    mode = largest_mode(variance, mean)
    descent = mode + 1
    terminal = not mp.isfinite(log_skellam_pmf(descent, variance, mean))
    drop = mp.mpf(1) if terminal else effective_drop(variance, mean, descent)
    scaled = variance * drop
    lower = (
        mp.mpf(0)
        if terminal
        else baricz_lower_bound(variance, mean, abs(descent))
    )
    tolerance = mp.power(10, -(mp.mp.dps - 20))

    if mode < mean - 1 - tolerance or mode > mean + 1 + tolerance:
        raise AssertionError(
            f"mode outside closed Darroch window: M={mode}, mean={mean}"
        )
    if scaled <= mp.mpf("0.25") - tolerance:
        raise AssertionError(
            f"quarter bound failed: V*Delta={scaled}, V={variance}, mean={mean}"
        )
    if lower and drop <= lower - tolerance:
        raise AssertionError(
            f"Baricz lower bound failed: Delta={drop}, lower={lower}"
        )

    return {
        "variance": variance,
        "imbalance": imbalance,
        "mean": mean,
        "mode": mode,
        "descent": descent,
        "terminal_descent": terminal,
        "effective_drop": drop,
        "variance_times_effective_drop": scaled,
        "baricz_lower_bound": lower,
        "variance_times_baricz_lower_bound": variance * lower,
    }


def bisect_transition(
    variance: mp.mpf,
    left: mp.mpf,
    right: mp.mpf,
    left_mode: int,
    right_mode: int,
) -> mp.mpf:
    if right_mode != left_mode + 1:
        raise AssertionError(
            f"unexpected mode jump {left_mode}->{right_mode} at V={variance}"
        )

    def difference(imbalance: mp.mpf) -> mp.mpf:
        mean = variance * imbalance
        return (
            log_skellam_pmf(right_mode, variance, mean)
            - log_skellam_pmf(left_mode, variance, mean)
        )

    f_left = difference(left)
    f_right = difference(right)
    if f_left > 0 or f_right < 0:
        raise AssertionError("transition bracket does not straddle a mode tie")
    # Two bisection steps per requested decimal digit comfortably resolve
    # the dps/2-sized one-sided transition probes used below.
    for _ in range(mp.mp.dps * 2):
        mid = (left + right) / 2
        f_mid = difference(mid)
        if f_mid < 0:
            left = mid
        else:
            right = mid
    return (left + right) / 2


def compact(row: dict[str, Any]) -> dict[str, Any]:
    return {
        key: (mp.nstr(value, 50) if isinstance(value, mp.mpf) else value)
        for key, value in row.items()
    }


def update_best(
    best: dict[str, Any] | None,
    row: dict[str, Any],
) -> dict[str, Any]:
    if best is None:
        return row
    if (
        row["variance_times_effective_drop"]
        < best["variance_times_effective_drop"]
    ):
        return row
    return best


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--v-values", default=DEFAULT_VARIANCES)
    parser.add_argument("--imbalance-grid", type=int, default=101)
    parser.add_argument("--boundary-exponents", type=int, default=12)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/skellam_effective_drop_verify_2026-07-10.json"),
    )
    args = parser.parse_args()
    if args.dps < 50:
        parser.error("--dps must be at least 50")
    if args.imbalance_grid < 3:
        parser.error("--imbalance-grid must be at least 3")
    mp.mp.dps = args.dps

    variances = parse_values(args.v_values)
    base_grid = [
        -1 + 2 * mp.mpf(index) / (args.imbalance_grid - 1)
        for index in range(args.imbalance_grid)
    ]
    boundary = []
    for exponent in range(1, args.boundary_exponents + 1):
        epsilon = mp.power(10, -exponent)
        boundary.extend(
            [-1 + epsilon, -epsilon, epsilon, 1 - epsilon]
        )
    imbalances = sorted(set([*base_grid, *boundary]))

    processed = 0
    terminal_descents = 0
    best: dict[str, Any] | None = None
    transition_best: dict[str, Any] | None = None
    transitions: list[dict[str, Any]] = []

    for variance in variances:
        previous_u: mp.mpf | None = None
        previous_mode: int | None = None
        for imbalance in imbalances:
            row = analyze_point(variance, imbalance)
            processed += 1
            terminal_descents += int(row["terminal_descent"])
            best = update_best(best, row)
            mode = int(row["mode"])
            if previous_mode is not None and mode != previous_mode:
                if (
                    previous_u == mp.mpf(-1)
                    and previous_mode == 0
                    and mode == -1
                ):
                    previous_u = imbalance
                    previous_mode = mode
                    continue
                if mode < previous_mode:
                    raise AssertionError(
                        f"unexpected decreasing mode jump "
                        f"{previous_mode}->{mode} at V={variance}"
                    )
                for left_mode in range(previous_mode, mode):
                    right_mode = left_mode + 1
                    root = bisect_transition(
                        variance,
                        previous_u,
                        imbalance,
                        left_mode,
                        right_mode,
                    )
                    epsilon = mp.power(10, -(args.dps // 2))
                    samples = [
                        max(mp.mpf(-1), root - epsilon),
                        root,
                        min(mp.mpf(1), root + epsilon),
                    ]
                    transition_rows = []
                    for sample in samples:
                        transition_row = analyze_point(variance, sample)
                        processed += 1
                        terminal_descents += int(
                            transition_row["terminal_descent"]
                        )
                        best = update_best(best, transition_row)
                        transition_best = update_best(
                            transition_best,
                            transition_row,
                        )
                        transition_rows.append(compact(transition_row))
                    transitions.append(
                        {
                            "variance": mp.nstr(variance, 50),
                            "imbalance_root": mp.nstr(root, 50),
                            "left_mode": left_mode,
                            "right_mode": right_mode,
                            "samples": transition_rows,
                        }
                    )
            previous_u = imbalance
            previous_mode = mode

    if best is None:
        raise AssertionError("no rows processed")
    result = {
        "kind": "skellam_effective_drop_high_precision_verification",
        "precision_dps": args.dps,
        "variance_values": [mp.nstr(value, 50) for value in variances],
        "imbalance_grid": args.imbalance_grid,
        "boundary_exponents": args.boundary_exponents,
        "processed": processed,
        "mode_transitions_refined": len(transitions),
        "quarter_bound_failures": 0,
        "baricz_bound_failures": 0,
        "darroch_window_failures": 0,
        "terminal_descents": terminal_descents,
        "best": compact(best),
        "best_transition_sample": (
            compact(transition_best) if transition_best is not None else None
        ),
        "transitions": transitions,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "transitions"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
