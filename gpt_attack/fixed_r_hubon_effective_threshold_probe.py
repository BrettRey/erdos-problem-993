"""Effective hub-on threshold probe for fixed-r spider lanes.

The older hub-on perturbation certificate compares the perturbation to
``1/(R a)``.  After the hub-off asymptotic work, we have a stronger eventual
reserve target:

    reserve >= eta_reserve/a,  eta_reserve = C_{r,q}/2.

Likewise, the F-mode adjacent margins have explicit first-order constants, so
mode domination can compare against ``eta_mode/a`` rather than ``1/(M a)``.

This diagnostic searches for thresholds where the existing elementary hub-on
tail bounds are small enough for these sharper targets.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import ceil
from pathlib import Path

import sympy as sp

from fixed_r_effective_threshold_probe import canonical_shift
from fixed_r_huboff_asymptotic import raw_moments_at_one, reserve_constant_moment_formula
from fixed_r_hubon_mode_certificate import best_witness_term, shifted_positive
from fixed_r_hubon_route2_perturbation import (
    eval_poly_fraction,
    tail_weighted_step_decreasing_expr,
)
from route2_spider_lane_scan import path_polys


t = sp.symbols("t", integer=True, positive=True)


def parse_r_values(raw: str) -> list[int]:
    return [int(part) for part in raw.split(",") if part.strip()]


def first_order_lambda_constant(path_poly: list[int], shift: int) -> Fraction:
    """Return L(D) where F_{m-1}/F_m = 1 + L(D)/a + O(a^-2)."""
    m1, _, _ = raw_moments_at_one(path_poly)
    return Fraction(3 * shift, 2) - Fraction(9, 2) * m1 - 3


def effective_mode_eta(path_poly: list[int], shift: int) -> Fraction:
    left_const = -first_order_lambda_constant(path_poly, shift)
    right_const = first_order_lambda_constant(path_poly, shift + 3)
    eta = min(left_const, right_const) / 2
    if eta <= 0:
        raise ValueError((shift, left_const, right_const))
    return eta


def effective_reserve_eta(path_poly: list[int], shift: int) -> Fraction:
    reserve_const = reserve_constant_moment_formula(path_poly, shift)
    if reserve_const <= 0:
        raise ValueError((shift, reserve_const))
    return reserve_const / 2


def first_a_in_residue(threshold: int, residue: int) -> int:
    a0 = threshold
    while a0 % 3 != residue:
        a0 += 1
    return a0


def hubon_effective_passes(r: int, threshold: int) -> dict:
    paths = path_polys(max(2, r))
    pr = paths[r]
    prm1 = paths[r - 1]
    hub_on_sum = sum(prm1)
    mixture_constant = (
        2
        * eval_poly_fraction(prm1, Fraction(2))
        / eval_poly_fraction(pr, Fraction(1, 2))
    )

    residues = []
    ok = True
    for residue in [0, 1, 2]:
        shift = canonical_shift(pr, r, residue)
        a0 = first_a_in_residue(threshold, residue)
        threshold_t = ceil((threshold - residue) / 3)
        m0 = (2 * a0 + shift) // 3
        s, tail_bound = best_witness_term(pr, hub_on_sum, a0, m0)
        n0 = m0 - s

        mode_eta = effective_mode_eta(pr, shift)
        reserve_eta = effective_reserve_eta(pr, shift)

        mode_bound = 2 * tail_bound
        mode_target = mode_eta / a0
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
        reserve_target = reserve_eta / a0
        perturb_threshold_ok = perturb_bound < reserve_target
        tail_small_ok = tail_bound < Fraction(1, 2)
        lambda_decreasing_ok = shifted_positive(
            tail_weighted_step_decreasing_expr(a, n, a * (a + r) ** 2),
            threshold_t,
        )
        mixture_decreasing_ok = (a0 + 3) * (a0 + r + 3) * 27 < a0 * (a0 + r) * 64

        residue_ok = (
            mode_threshold_ok
            and mode_decreasing_ok
            and perturb_threshold_ok
            and tail_small_ok
            and lambda_decreasing_ok
            and mixture_decreasing_ok
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
                "mode_eta": str(mode_eta),
                "reserve_eta": str(reserve_eta),
                "mode_bound_float": float(mode_bound),
                "mode_target_float": float(mode_target),
                "perturb_bound_float": float(perturb_bound),
                "reserve_target_float": float(reserve_target),
                "tail_small_ok": tail_small_ok,
                "mode_threshold_ok": mode_threshold_ok,
                "mode_decreasing_ok": mode_decreasing_ok,
                "perturb_threshold_ok": perturb_threshold_ok,
                "lambda_decreasing_ok": lambda_decreasing_ok,
                "mixture_decreasing_ok": mixture_decreasing_ok,
                "ok": residue_ok,
            }
        )

    return {
        "r": r,
        "threshold": threshold,
        "ok": ok,
        "max_mode_bound_float": max(row["mode_bound_float"] for row in residues),
        "min_mode_target_float": min(row["mode_target_float"] for row in residues),
        "max_perturb_bound_float": max(row["perturb_bound_float"] for row in residues),
        "min_reserve_target_float": min(row["reserve_target_float"] for row in residues),
        "residues": residues,
    }


def thresholds_to_try(start: int, stop: int, step: int) -> list[int]:
    out = []
    value = start
    while value <= stop:
        out.append(value)
        value += step
    return out


def probe_r(r: int, start: int, stop: int, step: int) -> dict:
    attempts = []
    for threshold in thresholds_to_try(start, stop, step):
        rec = hubon_effective_passes(r, threshold)
        attempts.append(
            {
                "threshold": threshold,
                "ok": rec["ok"],
                "max_mode_bound_float": rec["max_mode_bound_float"],
                "min_mode_target_float": rec["min_mode_target_float"],
                "max_perturb_bound_float": rec["max_perturb_bound_float"],
                "min_reserve_target_float": rec["min_reserve_target_float"],
            }
        )
        if rec["ok"]:
            return {"r": r, "first_ok": rec, "attempts": attempts}
    return {"r": r, "first_ok": None, "attempts": attempts}


def main() -> None:
    ap = argparse.ArgumentParser(description="Probe effective fixed-r hub-on thresholds.")
    ap.add_argument("--r-values", default="80,120,160,200")
    ap.add_argument("--start", type=int, default=20)
    ap.add_argument("--stop", type=int, default=400)
    ap.add_argument("--step", type=int, default=10)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    records = [
        probe_r(r, args.start, args.stop, args.step)
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
            f"min_reserve_target={first['min_reserve_target_float']:.3e}"
        )

    payload = {
        "params": {
            "r_values": parse_r_values(args.r_values),
            "start": args.start,
            "stop": args.stop,
            "step": args.step,
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
