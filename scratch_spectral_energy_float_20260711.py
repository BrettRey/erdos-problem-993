#!/usr/bin/env python3
"""Fast symmetric-tree pressure test for the D14 energy inequality.

This evaluates the marked residual-forest recurrence at scalar values of the
vertex marker y.  Gauss--Legendre quadrature then uses

    sum_C h(C)^2/q(C) = integral_0^1 sum_C h(C)^2 y^(q(C)-1) dy.

It is a floating-point falsification tool for large double regular trees.  An
exact failure must be replayed by ``scratch_spectral_energy_dp_20260711.py``.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import product
from math import exp, log

import numpy as np
from numpy.polynomial.legendre import leggauss


@dataclass(frozen=True)
class MP:
    # Rows are count, sum h, sum h^2, and sum S2 by selected-set rank.
    values: np.ndarray
    log_scale: float


def normalized(values: np.ndarray, log_scale: float) -> MP | None:
    size = float(np.max(np.abs(values[0]))) if values.size else 0.0
    if size == 0.0:
        return None
    return MP(values / size, log_scale + log(size))


ONE = MP(np.array([[1.0], [0.0], [0.0], [0.0]]), 0.0)


def mp_add(items: list[MP | None]) -> MP | None:
    nonzero = [item for item in items if item is not None]
    if not nonzero:
        return None
    scale = max(item.log_scale for item in nonzero)
    length = max(item.values.shape[1] for item in nonzero)
    values = np.zeros((4, length))
    for item in nonzero:
        weight = exp(item.log_scale - scale)
        values[:, : item.values.shape[1]] += weight * item.values
    return normalized(values, scale)


def mp_mul(a: MP | None, b: MP | None) -> MP | None:
    if a is None or b is None:
        return None
    an, ah, ah2, as2 = a.values
    bn, bh, bh2, bs2 = b.values
    count = np.convolve(an, bn)
    edge = np.convolve(ah, bn) + np.convolve(an, bh)
    edge_square = (
        np.convolve(ah2, bn)
        + 2 * np.convolve(ah, bh)
        + np.convolve(an, bh2)
    )
    degree_square = np.convolve(as2, bn) + np.convolve(an, bs2)
    return normalized(
        np.vstack((count, edge, edge_square, degree_square)),
        a.log_scale + b.log_scale,
    )


def mp_product(items: tuple[MP | None, ...]) -> MP | None:
    result: MP | None = ONE
    for item in items:
        result = mp_mul(result, item)
    return result


def add_local(
    item: MP | None,
    *,
    y: float,
    selected: int,
    addable: int,
    parent_edge: int,
    degree_square: int,
) -> MP | None:
    if item is None:
        return None
    count, edge, edge2, s2 = item.values
    values = np.vstack(
        (
            count,
            edge + parent_edge * count,
            edge2 + 2 * parent_edge * edge + parent_edge * count,
            s2 + degree_square * count,
        )
    )
    if addable:
        values *= y
    if selected:
        values = np.pad(values, ((0, 0), (1, 0)))
    return normalized(values, item.log_scale)


State = tuple[
    dict[tuple[int, int], MP | None],
    dict[tuple[int, int], MP | None],
    dict[tuple[int, int], MP | None],
]
ENVIRONMENTS = ((0, 0), (0, 1), (1, 0))


def vertex_state(children: list[State], y: float) -> State:
    selected: dict[tuple[int, int], MP | None] = {}
    no_selected: dict[tuple[int, int], MP | None] = {}
    with_selected: dict[tuple[int, int], MP | None] = {}

    for parent_selected, parent_addable in ENVIRONMENTS:
        if parent_selected:
            selected[(parent_selected, parent_addable)] = None
        else:
            terms = []
            for types in product((1, 2), repeat=len(children)):
                terms.append(
                    mp_product(
                        tuple(
                            children[index][child_type][(1, 0)]
                            for index, child_type in enumerate(types)
                        )
                    )
                )
            selected[(parent_selected, parent_addable)] = add_local(
                mp_add(terms),
                y=y,
                selected=1,
                addable=0,
                parent_edge=0,
                degree_square=0,
            )

        vertex_addable = 1 - parent_selected
        terms = []
        for types in product((1, 2), repeat=len(children)):
            child_product = mp_product(
                tuple(
                    children[index][child_type][
                        (0, vertex_addable)
                    ]
                    for index, child_type in enumerate(types)
                )
            )
            addable_children = sum(child_type == 1 for child_type in types)
            residual_degree = (
                parent_addable + addable_children
                if vertex_addable else 0
            )
            terms.append(
                add_local(
                    child_product,
                    y=y,
                    selected=0,
                    addable=vertex_addable,
                    parent_edge=vertex_addable * parent_addable,
                    degree_square=residual_degree * residual_degree,
                )
            )
        no_selected[(parent_selected, parent_addable)] = mp_add(terms)

        terms = []
        for types in product((0, 1, 2), repeat=len(children)):
            if 0 not in types:
                continue
            terms.append(
                mp_product(
                    tuple(
                        children[index][child_type][(0, 0)]
                        for index, child_type in enumerate(types)
                    )
                )
            )
        with_selected[(parent_selected, parent_addable)] = mp_add(terms)

    return selected, no_selected, with_selected


def double_regular_distribution(
    branching: int, depth: int, y: float
) -> MP:
    levels = [vertex_state([], y)]
    for _ in range(depth):
        levels.append(vertex_state([levels[-1]] * branching, y))
    full_root = vertex_state(
        [levels[depth - 1]] * branching + [levels[depth]], y
    )
    result = mp_add(
        [state[(0, 0)] for state in full_root]
    )
    assert result is not None
    return result


def double_spherical_distribution(
    branching: tuple[int, ...], y: float
) -> MP:
    """Distribution for two identical rooted spherical trees joined at roots.

    ``branching[j]`` is the number of children at distance ``j`` from each
    central root.
    """
    subtree = vertex_state([], y)
    below_root: State | None = None
    for level, factor in enumerate(reversed(branching)):
        subtree = vertex_state([subtree] * factor, y)
        if level == len(branching) - 2:
            below_root = subtree
    if len(branching) == 1:
        below_root = vertex_state([], y)
    assert below_root is not None
    full_root = vertex_state(
        [below_root] * branching[0] + [subtree], y
    )
    result = mp_add([state[(0, 0)] for state in full_root])
    assert result is not None
    return result


def energy_ratio(
    branching: int, depth: int, rank: int, quadrature: int
) -> float:
    at_one = double_regular_distribution(branching, depth, 1.0)
    count, edge, _, degree_square = at_one.values
    denominator = (rank + 1) * count[rank + 1] + edge[rank]
    degree_term = degree_square[rank - 1]

    nodes, weights = leggauss(quadrature)
    integral_scaled = 0.0
    for node, weight in zip((nodes + 1) / 2, weights / 2):
        marked = double_regular_distribution(branching, depth, float(node))
        coefficient = marked.values[2, rank - 1]
        integral_scaled += (
            weight
            * exp(marked.log_scale - at_one.log_scale)
            * coefficient
            / node
        )
    return float((degree_term - 4 * integral_scaled) / denominator)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--branching", type=int, default=2)
    parser.add_argument("--depth", type=int, default=6)
    parser.add_argument("--rank", type=int)
    parser.add_argument("--quadrature", type=int, default=96)
    args = parser.parse_args()
    # For the balanced binary family alpha is available cheaply from the
    # count support.  The last GSB-prefix rank is the pressure point.
    at_one = double_regular_distribution(
        args.branching, args.depth, 1.0
    )
    alpha = len(at_one.values[0]) - 1
    rank = args.rank
    if rank is None:
        rank = (2 * alpha + 1) // 3 - 2
    ratio = energy_ratio(
        args.branching, args.depth, rank, args.quadrature
    )
    print(
        {
            "branching": args.branching,
            "depth": args.depth,
            "alpha": alpha,
            "rank": rank,
            "quadrature": args.quadrature,
            "energy_ratio": ratio,
        }
    )


if __name__ == "__main__":
    main()
