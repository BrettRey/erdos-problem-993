"""Probe threshold growth for the fixed-r hub-on certificates.

This script is diagnostic, not a proof certificate.  It asks:

    For a given r, how large does the asymptotic threshold need to be before
    the existing hub-on mode and Route-2 perturbation bounds pass?

The hub-off symbolic reserve is intentionally not checked here; that remains
the bottleneck for r beyond the current sampled certificate range.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import ceil
from pathlib import Path

import sympy as sp

from fixed_r_huboff_certificate import stabilized_shifts
from fixed_r_hubon_mode_certificate import best_witness_term, shifted_positive
from fixed_r_hubon_route2_perturbation import (
    eval_poly_fraction,
    tail_step_decreasing_expr,
    tail_weighted_step_decreasing_expr,
)
from route2_spider_lane_scan import path_polys


t = sp.symbols("t", integer=True, positive=True)


def thresholds_to_try(start: int, stop: int, step: int) -> list[int]:
    out = []
    x = start
    while x <= stop:
        out.append(x)
        x += step
    return out


def hubon_passes(r: int, threshold: int, margin_denom: int, reserve_denom: int) -> dict:
    paths = path_polys(max(2, r))
    pr = paths[r]
    prm1 = paths[r - 1]
    hub_on_sum = sum(prm1)
    shifts = stabilized_shifts(r, threshold, paths)
    mixture_constant = (
        2
        * eval_poly_fraction(prm1, Fraction(2))
        / eval_poly_fraction(pr, Fraction(1, 2))
    )

    residues = []
    ok = True
    for residue, shift in shifts.items():
        a0 = threshold
        while a0 % 3 != residue:
            a0 += 1
        threshold_t = ceil((threshold - residue) / 3)
        m0 = (2 * a0 + shift) // 3
        s, tail_bound = best_witness_term(pr, hub_on_sum, a0, m0)
        n0 = m0 - s

        mode_bound = 2 * tail_bound
        mode_target = Fraction(1, margin_denom * a0)
        mode_threshold_ok = mode_bound < mode_target

        a = 3 * t + residue
        m = sp.simplify((2 * a + shift) / 3)
        n = sp.simplify(m - s)
        mode_decreasing_ok = shifted_positive(
            a * (a + 1) * (a + 2)
            - 2 * (n + 1) * (n + 2) * (a - n + 1),
            threshold_t,
        )

        degree_bound = a0 + r
        lambda_shift_bound = 2 * degree_bound * degree_bound * tail_bound
        mixture_bound = degree_bound * mixture_constant * Fraction(3, 4) ** (a0 - 1)
        perturb_bound = lambda_shift_bound + mixture_bound
        reserve = Fraction(1, reserve_denom * a0)
        perturb_threshold_ok = perturb_bound < reserve
        # The reserve target is proportional to 1/a, so certify decrease of
        # a times each perturbation term.
        perturb_decreasing_ok = shifted_positive(
            tail_weighted_step_decreasing_expr(a, n, a * (a + r) ** 2), threshold_t
        ) and ((a0 + 3) * (a0 + r + 3) * 27 < a0 * (a0 + r) * 64)

        residue_ok = (
            mode_threshold_ok
            and mode_decreasing_ok
            and perturb_threshold_ok
            and perturb_decreasing_ok
        )
        ok = ok and residue_ok
        residues.append(
            {
                "residue": residue,
                "shift": shift,
                "a0": a0,
                "m0": m0,
                "witness_s": s,
                "n0": n0,
                "mode_bound_float": float(mode_bound),
                "mode_target_float": float(mode_target),
                "perturb_bound_float": float(perturb_bound),
                "reserve_float": float(reserve),
                "ok": residue_ok,
            }
        )

    return {
        "r": r,
        "threshold": threshold,
        "margin_denom": margin_denom,
        "reserve_denom": reserve_denom,
        "shifts": shifts,
        "ok": ok,
        "max_mode_bound_float": max(row["mode_bound_float"] for row in residues),
        "max_perturb_bound_float": max(row["perturb_bound_float"] for row in residues),
        "min_reserve_float": min(row["reserve_float"] for row in residues),
        "residues": residues,
    }


def probe_r(r: int, start: int, stop: int, step: int, margin_denom: int, reserve_denom: int) -> dict:
    attempts = []
    for threshold in thresholds_to_try(start, stop, step):
        rec = hubon_passes(r, threshold, margin_denom, reserve_denom)
        attempts.append(
            {
                "threshold": threshold,
                "ok": rec["ok"],
                "max_mode_bound_float": rec["max_mode_bound_float"],
                "max_perturb_bound_float": rec["max_perturb_bound_float"],
                "min_reserve_float": rec["min_reserve_float"],
                "shifts": rec["shifts"],
            }
        )
        if rec["ok"]:
            return {"r": r, "first_ok": rec, "attempts": attempts}
    return {"r": r, "first_ok": None, "attempts": attempts}


def parse_r_values(raw: str) -> list[int]:
    return [int(part) for part in raw.split(",") if part.strip()]


def main() -> None:
    ap = argparse.ArgumentParser(description="Probe fixed-r hub-on threshold growth.")
    ap.add_argument("--r-values", default="80,100,120,160,200")
    ap.add_argument("--start", type=int, default=200)
    ap.add_argument("--stop", type=int, default=600)
    ap.add_argument("--step", type=int, default=50)
    ap.add_argument("--margin-denom", type=int, default=1000)
    ap.add_argument("--reserve-denom", type=int, default=1000)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    records = [
        probe_r(r, args.start, args.stop, args.step, args.margin_denom, args.reserve_denom)
        for r in parse_r_values(args.r_values)
    ]
    for rec in records:
        if rec["first_ok"] is None:
            print(f"r={rec['r']}: no passing threshold in grid")
            continue
        first = rec["first_ok"]
        print(
            f"r={rec['r']}: first_ok threshold={first['threshold']} "
            f"max_perturb={first['max_perturb_bound_float']:.3e} "
            f"min_reserve={first['min_reserve_float']:.3e}"
        )

    payload = {
        "params": {
            "r_values": parse_r_values(args.r_values),
            "start": args.start,
            "stop": args.stop,
            "step": args.step,
            "margin_denom": args.margin_denom,
            "reserve_denom": args.reserve_denom,
        },
        "records": records,
    }
    if args.out:
        path = Path(args.out)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2))
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
