#!/usr/bin/env python3
"""Exact bipartition-imbalance Rayleigh quotients for the D14 down-up chain.

For each independent set, the DP tracks its selected black-minus-white
imbalance, together with the size and black-minus-white imbalance of its
addable-vertex forest.  This is enough to evaluate exactly the variance of
the selected imbalance and its down-up Dirichlet form.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction

from indpoly import independence_poly
from scratch_spectral_energy_dp_20260711 import make_double_regular_tree


@dataclass(frozen=True)
class Moment:
    count: int
    add_sum: int
    add_square_sum: int
    selected_sum: int
    selected_square_sum: int


Poly = dict[tuple[int, int], Moment]
ONE: Poly = {(0, 0): Moment(1, 0, 0, 0, 0)}


def poly_add(a: Poly, b: Poly, sign: int = 1) -> Poly:
    out: dict[tuple[int, int], list[int]] = {
        key: [
            value.count,
            value.add_sum,
            value.add_square_sum,
            value.selected_sum,
            value.selected_square_sum,
        ]
        for key, value in a.items()
    }
    for key, value in b.items():
        target = out.setdefault(key, [0, 0, 0, 0, 0])
        target[0] += sign * value.count
        target[1] += sign * value.add_sum
        target[2] += sign * value.add_square_sum
        target[3] += sign * value.selected_sum
        target[4] += sign * value.selected_square_sum
    return {
        key: Moment(*value)
        for key, value in out.items()
        if any(value)
    }


def poly_mul(a: Poly, b: Poly) -> Poly:
    out: dict[tuple[int, int], list[int]] = defaultdict(
        lambda: [0, 0, 0, 0, 0]
    )
    for (rank_a, q_a), x in a.items():
        for (rank_b, q_b), y in b.items():
            target = out[(rank_a + rank_b, q_a + q_b)]
            target[0] += x.count * y.count
            target[1] += x.add_sum * y.count + x.count * y.add_sum
            target[2] += (
                x.add_square_sum * y.count
                + 2 * x.add_sum * y.add_sum
                + x.count * y.add_square_sum
            )
            target[3] += (
                x.selected_sum * y.count + x.count * y.selected_sum
            )
            target[4] += (
                x.selected_square_sum * y.count
                + 2 * x.selected_sum * y.selected_sum
                + x.count * y.selected_square_sum
            )
    return {key: Moment(*value) for key, value in out.items()}


def shift(
    polynomial: Poly,
    *,
    rank: int = 0,
    q: int = 0,
    add: int = 0,
    selected: int = 0,
) -> Poly:
    out: Poly = {}
    for (old_rank, old_q), value in polynomial.items():
        out[(old_rank + rank, old_q + q)] = Moment(
            value.count,
            value.add_sum + add * value.count,
            value.add_square_sum
            + 2 * add * value.add_sum
            + add * add * value.count,
            value.selected_sum + selected * value.count,
            value.selected_square_sum
            + 2 * selected * value.selected_sum
            + selected * selected * value.count,
        )
    return out


def color_distribution(adj: list[list[int]], root: int = 0) -> Poly:
    n = len(adj)
    parent = [-2] * n
    color = [0] * n
    parent[root] = -1
    color[root] = 1
    order = [root]
    for vertex in order:
        for neighbor in adj[vertex]:
            if parent[neighbor] == -2:
                parent[neighbor] = vertex
                color[neighbor] = -color[vertex]
                order.append(neighbor)

    selected: list[dict[tuple[int, int], Poly]] = [{} for _ in adj]
    no_selected_child: list[dict[tuple[int, int], Poly]] = [{} for _ in adj]
    with_selected_child: list[dict[tuple[int, int], Poly]] = [{} for _ in adj]
    environments = ((0, 0), (0, 1), (1, 0))

    for vertex in reversed(order):
        children = [u for u in adj[vertex] if parent[u] == vertex]
        for parent_selected, parent_addable in environments:
            product = ONE
            for child in children:
                product = poly_mul(
                    product,
                    poly_add(
                        no_selected_child[child][(1, 0)],
                        with_selected_child[child][(1, 0)],
                    ),
                )
            selected[vertex][(parent_selected, parent_addable)] = (
                {}
                if parent_selected
                else shift(product, rank=1, selected=color[vertex])
            )

            vertex_addable = 1 - parent_selected
            no_product = ONE
            for child in children:
                no_product = poly_mul(
                    no_product,
                    poly_add(
                        shift(
                            no_selected_child[child][
                                (0, vertex_addable)
                            ],
                            add=0,
                        ),
                        with_selected_child[child][
                            (0, vertex_addable)
                        ],
                    ),
                )
            no_selected_child[vertex][
                (parent_selected, parent_addable)
            ] = shift(
                no_product,
                q=vertex_addable,
                add=color[vertex] if vertex_addable else 0,
            )

            arbitrary = ONE
            without_selected = ONE
            for child in children:
                arbitrary = poly_mul(
                    arbitrary,
                    poly_add(
                        selected[child][(0, 0)],
                        poly_add(
                            no_selected_child[child][(0, 0)],
                            with_selected_child[child][(0, 0)],
                        ),
                    ),
                )
                without_selected = poly_mul(
                    without_selected,
                    poly_add(
                        no_selected_child[child][(0, 0)],
                        with_selected_child[child][(0, 0)],
                    ),
                )
            with_selected_child[vertex][
                (parent_selected, parent_addable)
            ] = poly_add(arbitrary, without_selected, -1)

    return poly_add(
        selected[root][(0, 0)],
        poly_add(
            no_selected_child[root][(0, 0)],
            with_selected_child[root][(0, 0)],
        ),
    )


def rayleigh(distribution: Poly, polynomial: list[int], rank: int) -> Fraction:
    count = polynomial[rank]
    selected_sum = 0
    selected_square_sum = 0
    for (set_rank, _), value in distribution.items():
        if set_rank == rank:
            selected_sum += value.selected_sum
            selected_square_sum += value.selected_square_sum
    variance = Fraction(selected_square_sum, count) - Fraction(
        selected_sum * selected_sum, count * count
    )

    numerator = Fraction(0)
    for (set_rank, q), value in distribution.items():
        if set_rank == rank - 1 and q:
            numerator += value.count * q - Fraction(
                value.add_square_sum, q
            )
    dirichlet = numerator / (rank * count)
    return dirichlet / variance if variance else Fraction(0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--branching", type=int, default=2)
    parser.add_argument("--depth", type=int, default=4)
    args = parser.parse_args()
    adj = make_double_regular_tree(args.branching, args.depth)
    distribution = color_distribution(adj)
    polynomial = independence_poly(len(adj), adj)
    counts = [0] * len(polynomial)
    for (rank, _), value in distribution.items():
        counts[rank] += value.count
    assert counts == polynomial
    alpha = len(polynomial) - 1
    limit = (2 * alpha + 1) // 3
    best: tuple[Fraction, int] | None = None
    for rank in range(1, limit - 1):
        quotient = rayleigh(distribution, polynomial, rank)
        scaled = 2 * rank * quotient
        if best is None or scaled < best[0]:
            best = (scaled, rank)
    print(
        {
            "n": len(adj),
            "alpha": alpha,
            "best_2r_rayleigh": best,
        }
    )


if __name__ == "__main__":
    main()
