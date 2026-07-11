#!/usr/bin/env python3
"""Tilted FFT search on recursively nested D25 hard-hub states.

Start with B_0=TG(4,6) rooted at its distinguished root.  Given a rooted
child state B and an integer m, form

    H_m(B) = a new root c joined to G(40,20) and to m copies of B.

Thus a sequence ``m1,m2,...`` means
``H_mk(...H_m2(H_m1(TG))...)``.  The script evaluates ordered log-concavity
and GSB at the last required prefix rank.  Exact base polynomials are tilted
to their saddle point and all large products use floating FFT convolution.
Depth one is audited against the exact D25 certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass

import numpy as np
from scipy.signal import fftconvolve

from scratch_d25_hard_branch_composition_search_20260711 import make_state
from scripts.analyze_prufer_corpus import (
    make_bautista_ramos_tree,
    make_galvin_tree,
)


@dataclass(frozen=True)
class Moment:
    log_value: float
    mean: float
    variance: float


@dataclass(frozen=True)
class MomentState:
    total: Moment
    excluded: Moment
    reserve: Moment


@dataclass(frozen=True)
class DegreeState:
    total: int
    excluded: int
    reserve: int
    order: int
    hubs: int


@dataclass(frozen=True)
class Scaled:
    values: np.ndarray
    log_scale: float


@dataclass(frozen=True)
class ScaledState:
    total: Scaled
    excluded: Scaled
    reserve: Scaled


def moment_product(left: Moment, right: Moment) -> Moment:
    return Moment(
        left.log_value + right.log_value,
        left.mean + right.mean,
        left.variance + right.variance,
    )


def moment_power(item: Moment, exponent: int) -> Moment:
    return Moment(
        exponent * item.log_value,
        exponent * item.mean,
        exponent * item.variance,
    )


def moment_shift(item: Moment, log_z: float) -> Moment:
    return Moment(item.log_value + log_z, item.mean + 1.0, item.variance)


def moment_add(left: Moment, right: Moment) -> Moment:
    top = max(left.log_value, right.log_value)
    left_weight = math.exp(left.log_value - top)
    right_weight = math.exp(right.log_value - top)
    total_weight = left_weight + right_weight
    probability = right_weight / total_weight
    mean = (1.0 - probability) * left.mean + probability * right.mean
    variance = (
        (1.0 - probability)
        * (left.variance + (left.mean - mean) ** 2)
        + probability * (right.variance + (right.mean - mean) ** 2)
    )
    return Moment(top + math.log(total_weight), mean, variance)


def base_moment(coefficients: tuple[int, ...], log_z: float) -> Moment:
    logs = np.array(
        [math.log(value) + degree * log_z if value else -math.inf
         for degree, value in enumerate(coefficients)],
        dtype=float,
    )
    top = float(np.max(logs))
    weights = np.exp(logs - top)
    total = float(np.sum(weights))
    weights /= total
    degrees = np.arange(len(coefficients), dtype=float)
    mean = float(np.dot(degrees, weights))
    variance = float(np.dot((degrees - mean) ** 2, weights))
    return Moment(top + math.log(total), mean, variance)


def base_moment_state(state, log_z: float) -> MomentState:
    return MomentState(*(
        base_moment(poly, log_z)
        for poly in (state.total, state.excluded, state.reserve)
    ))


def compose_moment(
    parent: MomentState, child: MomentState, multiplicity: int, log_z: float
) -> MomentState:
    excluded = moment_product(
        parent.total, moment_power(child.total, multiplicity)
    )
    reserve = moment_product(
        parent.excluded, moment_power(child.excluded, multiplicity)
    )
    total = moment_add(excluded, moment_shift(reserve, log_z))
    return MomentState(total, excluded, reserve)


def trim_degree(coefficients: tuple[int, ...]) -> int:
    return len(coefficients) - 1


def base_degree_state(state, order: int) -> DegreeState:
    return DegreeState(
        trim_degree(state.total),
        trim_degree(state.excluded),
        trim_degree(state.reserve),
        order,
        0,
    )


def compose_degree(
    parent: DegreeState, child: DegreeState, multiplicity: int
) -> DegreeState:
    excluded = parent.total + multiplicity * child.total
    reserve = parent.excluded + multiplicity * child.excluded
    total = max(excluded, 1 + reserve)
    return DegreeState(
        total,
        excluded,
        reserve,
        1 + parent.order + multiplicity * child.order,
        1 + multiplicity * child.hubs,
    )


def normalize(item: Scaled) -> Scaled:
    maximum = float(np.max(item.values))
    if not maximum > 0.0:
        raise ArithmeticError("lost all polynomial mass")
    return Scaled(item.values / maximum, item.log_scale + math.log(maximum))


def base_scaled(coefficients: tuple[int, ...], log_z: float) -> Scaled:
    logs = np.array(
        [math.log(value) + degree * log_z if value else -math.inf
         for degree, value in enumerate(coefficients)],
        dtype=float,
    )
    top = float(np.max(logs))
    values = np.exp(logs - top)
    values[~np.isfinite(logs)] = 0.0
    return normalize(Scaled(values, top))


def base_scaled_state(state, log_z: float) -> ScaledState:
    return ScaledState(*(
        base_scaled(poly, log_z)
        for poly in (state.total, state.excluded, state.reserve)
    ))


def multiply(left: Scaled, right: Scaled) -> Scaled:
    output_length = len(left.values) + len(right.values) - 1
    if min(len(left.values), len(right.values)) < 96 or output_length < 2048:
        values = np.convolve(left.values, right.values)
    else:
        values = fftconvolve(left.values, right.values)
    # FFT roundoff can create tiny negative coefficients.
    values[values < 0.0] = 0.0
    return normalize(Scaled(values, left.log_scale + right.log_scale))


def power(item: Scaled, exponent: int) -> Scaled:
    result = Scaled(np.array([1.0]), 0.0)
    factor = item
    while exponent:
        if exponent & 1:
            result = multiply(result, factor)
        exponent //= 2
        if exponent:
            factor = multiply(factor, factor)
    return result


def shift(item: Scaled, log_z: float) -> Scaled:
    return Scaled(np.pad(item.values, (1, 0)), item.log_scale + log_z)


def weighted(item: Scaled, weight: int) -> Scaled:
    return Scaled(item.values, item.log_scale + math.log(weight))


def add(left: Scaled, right: Scaled) -> Scaled:
    top = max(left.log_scale, right.log_scale)
    length = max(len(left.values), len(right.values))
    values = np.zeros(length, dtype=float)
    values[: len(left.values)] += math.exp(left.log_scale - top) * left.values
    values[: len(right.values)] += math.exp(right.log_scale - top) * right.values
    return normalize(Scaled(values, top))


def compose_scaled(
    parent: ScaledState, child: ScaledState, multiplicity: int, log_z: float
) -> ScaledState:
    excluded = multiply(parent.total, power(child.total, multiplicity))
    reserve = multiply(parent.excluded, power(child.excluded, multiplicity))
    total = add(excluded, shift(reserve, log_z))
    return ScaledState(total, excluded, reserve)


def bottom_inverse_ratio(
    branches: list[tuple[ScaledState, int, int]], rank: int, log_z: float
) -> float:
    """Inverse-degree center ratio for a star of rooted branch states.

    Each item is ``(state, multiplicity, full degree of the branch root)``.
    This is used to check that bottom hard hubs remain locally obstructive
    after they acquire an exterior branch from the next nesting level.
    """
    power_cache: dict[tuple[int, str, int], Scaled] = {}

    def get_power(index: int, field: str, exponent: int) -> Scaled:
        key = (index, field, exponent)
        if key not in power_cache:
            power_cache[key] = power(getattr(branches[index][0], field), exponent)
        return power_cache[key]

    def product_items(items: list[Scaled]) -> Scaled:
        result = Scaled(np.array([1.0]), 0.0)
        for item in items:
            result = multiply(result, item)
        return result

    total_product = product_items([
        get_power(index, "total", count)
        for index, (_, count, _) in enumerate(branches)
    ])
    excluded_product = product_items([
        get_power(index, "excluded", count)
        for index, (_, count, _) in enumerate(branches)
    ])
    whole = add(total_product, shift(excluded_product, log_z))

    a_polynomials: list[Scaled] = []
    joint_terms: list[Scaled] = []
    for index, (state, count, _) in enumerate(branches):
        factors = [state.reserve, get_power(index, "total", count - 1)]
        factors += [
            get_power(other, "total", other_count)
            for other, (_, other_count, _) in enumerate(branches)
            if other != index
        ]
        a_polynomials.append(product_items(factors))
        if count >= 2:
            factors = [
                multiply(state.reserve, state.reserve),
                get_power(index, "total", count - 2),
            ]
            factors += [
                get_power(other, "total", other_count)
                for other, (_, other_count, _) in enumerate(branches)
                if other != index
            ]
            joint_terms.append(weighted(product_items(factors), count * (count - 1)))
        for other in range(index + 1, len(branches)):
            other_state, other_count, _ = branches[other]
            factors = [
                state.reserve,
                other_state.reserve,
                get_power(index, "total", count - 1),
                get_power(other, "total", other_count - 1),
            ]
            factors += [
                get_power(third, "total", third_count)
                for third, (_, third_count, _) in enumerate(branches)
                if third not in (index, other)
            ]
            joint_terms.append(
                weighted(product_items(factors), 2 * count * other_count)
            )
    joint = joint_terms[0]
    for term in joint_terms[1:]:
        joint = add(joint, term)

    def coefficient_log(item: Scaled) -> float:
        value = float(item.values[rank])
        if not value > 0.0:
            raise ArithmeticError("target coefficient vanished")
        return math.log(value) + item.log_scale

    log_n = coefficient_log(whole)
    log_j = coefficient_log(joint)
    log_ac = coefficient_log(excluded_product)
    log_a = [coefficient_log(item) for item in a_polynomials]
    top_a = max(log_a)
    sum_a_scaled = sum(
        count * math.exp(value - top_a)
        for value, (_, count, _) in zip(log_a, branches)
    )
    log_sum_a = top_a + math.log(sum_a_scaled)
    top = max(log_n + log_j, 2 * log_sum_a, log_ac + log_sum_a)
    n_joint = math.exp(log_n + log_j - top)
    sum_a = math.exp(log_sum_a - top / 2)
    sum_diagonal = sum(
        count * math.exp(2 * value - top)
        for value, (_, count, _) in zip(log_a, branches)
    )
    product_ordered = sum_a * sum_a - sum_diagonal
    q_ordered = n_joint - product_ordered
    inverse_budget = sum(
        count * math.exp(2 * value - top) / degree
        for value, (_, count, degree) in zip(log_a, branches)
    ) + math.exp(log_ac + log_sum_a - top)
    return q_ordered / inverse_budget


def parse_sequence(text: str) -> tuple[int, ...]:
    sequence = tuple(int(part) for part in text.split(",") if part)
    if not sequence or any(value < 1 for value in sequence):
        raise argparse.ArgumentTypeError("use comma-separated positive integers")
    return sequence


def evaluate(
    sequence: tuple[int, ...], parent_state, child_state, audit_bottom: bool = False
) -> dict:
    parent_degree = base_degree_state(parent_state, 1641)
    degree = base_degree_state(child_state, 162)
    for multiplicity in sequence:
        degree = compose_degree(parent_degree, degree, multiplicity)
    alpha = degree.total
    limit = math.ceil((2 * alpha - 1) / 3)
    rank = limit - 2

    def nested_moment(log_z: float) -> MomentState:
        parent = base_moment_state(parent_state, log_z)
        child = base_moment_state(child_state, log_z)
        for multiplicity in sequence:
            child = compose_moment(parent, child, multiplicity, log_z)
        return child

    lo, hi = -8.0, 8.0
    for _ in range(70):
        mid = (lo + hi) / 2.0
        if nested_moment(mid).total.mean < rank + 1:
            lo = mid
        else:
            hi = mid
    log_z = (lo + hi) / 2.0
    z = math.exp(log_z)

    parent_scaled = base_scaled_state(parent_state, log_z)
    base_child_scaled = base_scaled_state(child_state, log_z)
    child_scaled = base_child_scaled
    scaled_levels = [child_scaled]
    for multiplicity in sequence:
        child_scaled = compose_scaled(
            parent_scaled, child_scaled, multiplicity, log_z
        )
        scaled_levels.append(child_scaled)
    values = child_scaled.total.values
    if rank + 2 >= len(values):
        raise AssertionError((rank, len(values), degree))
    w0, w1, w2 = map(float, values[rank : rank + 3])
    ordered = (rank + 2) * w0 * w2 / ((rank + 1) * w1 * w1)
    gsb = (rank + 2) * w0 * w2 / (
        (rank + 1) * w1 * w1 + z * w0 * w1
    )
    saddle = nested_moment(log_z).total
    result = {
        "sequence_inner_to_outer": sequence,
        "depth": len(sequence),
        "n": degree.order,
        "alpha": alpha,
        "hubs": degree.hubs,
        "rank": rank,
        "log_z": log_z,
        "z": z,
        "saddle_mean": saddle.mean,
        "saddle_variance": saddle.variance,
        "ordered_lc_ratio": ordered,
        "ordered_lc_margin": 1.0 - ordered,
        "gsb_ratio": gsb,
        "gsb_margin": 1.0 - gsb,
        "coefficient_window": [w0, w1, w2],
    }
    if audit_bottom and len(sequence) >= 2:
        inner_m, outer_m = sequence[0], sequence[1]
        inner_state = scaled_levels[1]
        # The component outside a chosen inner hub is a top hub with one
        # fewer copy of the inner state.
        exterior = compose_scaled(
            parent_scaled, inner_state, outer_m - 1, log_z
        )
        result["bottom_inverse_ratio"] = bottom_inverse_ratio(
            [
                (parent_scaled, 1, 41),
                (base_child_scaled, inner_m, 6),
                (exterior, 1, outer_m + 1),
            ],
            rank,
            log_z,
        )
        result["bottom_inverse_fails"] = result["bottom_inverse_ratio"] > 1.0
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "sequences",
        nargs="*",
        type=parse_sequence,
        default=[(48,), (32, 32), (40, 40), (48, 48)],
    )
    parser.add_argument("--audit-bottom-local", action="store_true")
    args = parser.parse_args()
    parent_adj = make_galvin_tree(40, 20)
    child_adj = make_bautista_ramos_tree(4, 6)
    parent = make_state("G(40,20)", parent_adj, 0)
    child = make_state("TG(4,6)", child_adj, 0)
    for sequence in args.sequences:
        row = evaluate(sequence, parent, child, args.audit_bottom_local)
        if sequence == (48,):
            # Exact certificate: 0.999176513939652145136139...
            assert abs(row["ordered_lc_ratio"] - 0.9991765139396521) < 2e-10
        print(json.dumps(row, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
