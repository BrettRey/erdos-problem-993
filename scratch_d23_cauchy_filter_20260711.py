#!/usr/bin/env python3
"""High-precision Cauchy filter for the D23 depth-seven FFT artefact.

The double-precision search first sees a descent where coefficients are only
1e-14 of the mode.  This script evaluates the recursive generating function
directly on a saddle-point circle and extracts a few tilted coefficients by a
high-precision discrete Cauchy integral.  With M much larger than the tilted
standard deviation, aliases are exponentially remote.  Repeating at two M
and precision values makes the numerical-falsification role explicit; this is
not used as an exact counterexample certificate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

from scratch_d23_two_type_phase_search_20260711 import (
    decoration_data,
    specs,
)
from scratch_d17_checkpointed_search_20260711 import append_record, utc_now


def polynomial_value_theta(coefficients: list[int], z: mp.mpc) -> tuple[mp.mpc, mp.mpc]:
    value = mp.mpc(0)
    derivative = mp.mpc(0)
    power = mp.mpc(1)
    for k, coefficient in enumerate(coefficients):
        value += coefficient * power
        derivative += k * coefficient * power
        power *= z
    return value, derivative


def states_value_theta(layers, data, z: mp.mpc):
    excluded = (mp.mpc(1), mp.mpc(0))
    selected = (z, z)
    for layer in layers:
        old_e, old_et = excluded
        total = excluded[0] + selected[0]
        total_t = excluded[1] + selected[1]
        decor = data[layer.decor]
        dp, dpt = polynomial_value_theta(decor.total_exact, z)
        de, det = polynomial_value_theta(decor.excluded_exact, z)
        e = total**layer.branch * dp
        et = (
            layer.branch * total ** (layer.branch - 1) * total_t * dp
            + total**layer.branch * dpt
        )
        core = old_e**layer.branch * de
        core_t = (
            layer.branch * old_e ** (layer.branch - 1) * old_et * de
            + old_e**layer.branch * det
        )
        excluded = (e, et)
        selected = (z * core, z * core + z * core_t)
    return excluded, selected


def candidate_value_theta(z, left_layers, right_layers, data, selected_only: bool):
    le, ls = states_value_theta(left_layers, data, z)
    re, rs = states_value_theta(right_layers, data, z)
    if selected_only:
        return rs
    value = le[0] * re[0] + ls[0] * re[0] + le[0] * rs[0]
    theta = (
        le[1] * re[0] + le[0] * re[1]
        + ls[1] * re[0] + ls[0] * re[1]
        + le[1] * rs[0] + le[0] * rs[1]
    )
    return value, theta


def saddle_radius(k, left_layers, right_layers, data, selected_only: bool):
    low = mp.mpf("0.1")
    high = mp.mpf("1.0")
    for _ in range(100):
        mid = mp.sqrt(low * high)
        value, theta = candidate_value_theta(mid, left_layers, right_layers, data, selected_only)
        mean = mp.re(theta / value)
        if mean < k:
            low = mid
        else:
            high = mid
    return mp.sqrt(low * high)


def cauchy_ratios(
    ranks: list[int],
    left_layers,
    right_layers,
    data,
    selected_only: bool,
    points: int,
):
    centre = ranks[len(ranks) // 2]
    radius = saddle_radius(centre, left_layers, right_layers, data, selected_only)
    base, _ = candidate_value_theta(radius, left_layers, right_layers, data, selected_only)
    omega = mp.e ** (2j * mp.pi / points)
    steps = {k: omega ** (-k) for k in ranks}
    weights = {k: mp.mpc(1) for k in ranks}
    sums = {k: mp.mpc(0) for k in ranks}
    z = mp.mpc(radius)
    for _ in range(points):
        value, _ = candidate_value_theta(z, left_layers, right_layers, data, selected_only)
        normalized = value / base
        for k in ranks:
            sums[k] += normalized * weights[k]
            weights[k] *= steps[k]
        z *= omega
    tilted = {k: mp.re(sums[k] / points) for k in ranks}
    ratios = {
        k: tilted[k + 1] / tilted[k] / radius
        for k in ranks
        if k + 1 in tilted
    }
    return radius, tilted, ratios


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--points", type=int, default=4096)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--selected-only", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default="results/d23_cauchy_filter_20260711.jsonl",
    )
    args = parser.parse_args()
    mp.mp.dps = args.dps

    data = decoration_data()
    profiles = dict(specs(7, data, 500000))
    left = profiles["transient27_o0_d1"]
    right = profiles["critical5_o0_d7"]
    ranks = [114843, 114844, 114845, 115952, 115953, 115954, 116055, 116056, 116057]
    if args.selected_only:
        ranks = [114848, 114849, 114850, 115952, 115953, 115954, 116055, 116056, 116057]
    radius, tilted, ratios = cauchy_ratios(
        ranks, left, right, data, args.selected_only, args.points
    )
    result = {
        "at": utc_now(),
        "kind": "d23_high_precision_cauchy_filter",
        "points": args.points,
        "dps": args.dps,
        "selected_only": args.selected_only,
        "radius": mp.nstr(radius, 30),
        "tilted": {str(k): mp.nstr(value, 30) for k, value in tilted.items()},
        "ratios": {str(k): mp.nstr(value, 30) for k, value in ratios.items()},
    }
    append_record(args.output, result)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
