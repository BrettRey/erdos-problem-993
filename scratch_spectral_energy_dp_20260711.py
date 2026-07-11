#!/usr/bin/env python3
"""Exact marked DP for the D14 aggregate spectral-energy inequality.

For an independent set C, let H_C be the forest induced by vertices which
can be added to C, and put

    q=|V(H_C)|, h=|E(H_C)|, S2=sum_v d_H(v)^2.

At every rank the DP retains the distribution in q together with the sums
of h, h^2, and S2.  It therefore evaluates exactly the proposed aggregate
energy deficit

    sum_C [q(q+h-1)-2h-(r+1)S2+4r h^2/q]

for C of size r-1.  This is a falsification tool, not a proof of the
inequality.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction

from indpoly import independence_poly
from targeted import make_T_m_t_1


@dataclass(frozen=True)
class Moment:
    count: int
    edge_sum: int
    edge_square_sum: int
    degree_square_sum: int


Poly = dict[tuple[int, int], Moment]
ZERO = Moment(0, 0, 0, 0)
ONE: Poly = {(0, 0): Moment(1, 0, 0, 0)}
RANK_LIMIT: int | None = None


def poly_add(a: Poly, b: Poly, sign: int = 1) -> Poly:
    out: dict[tuple[int, int], list[int]] = {
        key: [value.count, value.edge_sum, value.edge_square_sum,
              value.degree_square_sum]
        for key, value in a.items()
    }
    for key, value in b.items():
        target = out.setdefault(key, [0, 0, 0, 0])
        target[0] += sign * value.count
        target[1] += sign * value.edge_sum
        target[2] += sign * value.edge_square_sum
        target[3] += sign * value.degree_square_sum
    return {
        key: Moment(*value)
        for key, value in out.items()
        if any(value)
    }


def poly_mul(a: Poly, b: Poly) -> Poly:
    out: dict[tuple[int, int], list[int]] = defaultdict(
        lambda: [0, 0, 0, 0]
    )
    for (rank_a, q_a), x in a.items():
        for (rank_b, q_b), y in b.items():
            if (
                RANK_LIMIT is not None
                and rank_a + rank_b > RANK_LIMIT
            ):
                continue
            target = out[(rank_a + rank_b, q_a + q_b)]
            target[0] += x.count * y.count
            target[1] += x.edge_sum * y.count + x.count * y.edge_sum
            target[2] += (
                x.edge_square_sum * y.count
                + 2 * x.edge_sum * y.edge_sum
                + x.count * y.edge_square_sum
            )
            target[3] += (
                x.degree_square_sum * y.count
                + x.count * y.degree_square_sum
            )
    return {key: Moment(*value) for key, value in out.items()}


def add_local(
    polynomial: Poly,
    *,
    rank_shift: int,
    q_shift: int,
    edge_shift: int,
    degree_square_shift: int,
) -> Poly:
    out: Poly = {}
    for (rank, q), value in polynomial.items():
        if RANK_LIMIT is not None and rank + rank_shift > RANK_LIMIT:
            continue
        out[(rank + rank_shift, q + q_shift)] = Moment(
            value.count,
            value.edge_sum + edge_shift * value.count,
            value.edge_square_sum
            + 2 * edge_shift * value.edge_sum
            + edge_shift * edge_shift * value.count,
            value.degree_square_sum
            + degree_square_shift * value.count,
        )
    return out


def energy_distribution(
    adj: list[list[int]], root: int = 0, max_rank: int | None = None
) -> Poly:
    """Return exact (rank,q)->(count,sum h,sum h^2,sum S2) data."""
    global RANK_LIMIT
    RANK_LIMIT = max_rank
    n = len(adj)
    parent = [-2] * n
    parent[root] = -1
    order = [root]
    for vertex in order:
        for neighbor in adj[vertex]:
            if parent[neighbor] == -2:
                parent[neighbor] = vertex
                order.append(neighbor)

    # Each state is indexed by the two parent-environment bits
    # (parent selected, parent addable).  Only (0,0), (0,1), and (1,0)
    # are realizable, but retaining a dictionary keeps the recurrence clear.
    selected: list[dict[tuple[int, int], Poly]] = [{} for _ in adj]
    no_selected_child: list[dict[tuple[int, int], Poly]] = [{} for _ in adj]
    with_selected_child: list[dict[tuple[int, int], Poly]] = [{} for _ in adj]
    environments = ((0, 0), (0, 1), (1, 0))
    signatures: list[tuple[object, ...] | None] = [None] * n
    cache: dict[
        tuple[object, ...],
        tuple[
            dict[tuple[int, int], Poly],
            dict[tuple[int, int], Poly],
            dict[tuple[int, int], Poly],
        ],
    ] = {}

    for vertex in reversed(order):
        children = [u for u in adj[vertex] if parent[u] == vertex]
        signature: tuple[object, ...] = tuple(
            sorted((signatures[u] for u in children), key=repr)
        )
        signatures[vertex] = signature
        if signature in cache:
            (
                selected[vertex],
                no_selected_child[vertex],
                with_selected_child[vertex],
            ) = cache[signature]
            continue

        for parent_selected, parent_addable in environments:
            # Vertex selected: every child is unselected, and the selected
            # vertex is neither addable nor adjacent to an addable parent.
            product = ONE
            for child in children:
                options = poly_add(
                    no_selected_child[child][(1, 0)],
                    with_selected_child[child][(1, 0)],
                )
                product = poly_mul(product, options)
            selected[vertex][(parent_selected, parent_addable)] = (
                {} if parent_selected else add_local(
                    product,
                    rank_shift=1,
                    q_shift=0,
                    edge_shift=0,
                    degree_square_shift=0,
                )
            )

            # Vertex unselected with no selected child.  Track the number of
            # addable children so that the local residual degree is exact.
            vertex_addable = 1 - parent_selected
            temp: dict[tuple[int, int, int], Moment] = {
                (0, 0, 0): Moment(1, 0, 0, 0)
            }
            for child in children:
                new: dict[tuple[int, int, int], list[int]] = defaultdict(
                    lambda: [0, 0, 0, 0]
                )
                for child_type, child_addable in (
                    (no_selected_child, 1),
                    (with_selected_child, 0),
                ):
                    child_poly = child_type[child][
                        (0, vertex_addable)
                    ]
                    for (rank, q, addable_children), x in temp.items():
                        for (child_rank, child_q), y in child_poly.items():
                            target = new[
                                (
                                    rank + child_rank,
                                    q + child_q,
                                    addable_children + child_addable,
                                )
                            ]
                            target[0] += x.count * y.count
                            target[1] += (
                                x.edge_sum * y.count
                                + x.count * y.edge_sum
                            )
                            target[2] += (
                                x.edge_square_sum * y.count
                                + 2 * x.edge_sum * y.edge_sum
                                + x.count * y.edge_square_sum
                            )
                            target[3] += (
                                x.degree_square_sum * y.count
                                + x.count * y.degree_square_sum
                            )
                temp = {
                    key: Moment(*value) for key, value in new.items()
                }

            no_poly: Poly = {}
            for (rank, q, addable_children), value in temp.items():
                residual_degree = (
                    parent_addable + addable_children
                    if vertex_addable else 0
                )
                local = add_local(
                    {(rank, q): value},
                    rank_shift=0,
                    q_shift=vertex_addable,
                    edge_shift=vertex_addable * parent_addable,
                    degree_square_shift=residual_degree * residual_degree,
                )
                no_poly = poly_add(no_poly, local)
            no_selected_child[vertex][
                (parent_selected, parent_addable)
            ] = no_poly

            # Vertex unselected with at least one selected child.  Here the
            # vertex is not addable, but a no-selected-child child is addable.
            arbitrary = ONE
            without_selected = ONE
            for child in children:
                all_options = poly_add(
                    selected[child][(0, 0)],
                    poly_add(
                        no_selected_child[child][(0, 0)],
                        with_selected_child[child][(0, 0)],
                    ),
                )
                nonselected_options = poly_add(
                    no_selected_child[child][(0, 0)],
                    with_selected_child[child][(0, 0)],
                )
                arbitrary = poly_mul(arbitrary, all_options)
                without_selected = poly_mul(
                    without_selected, nonselected_options
                )
            with_selected_child[vertex][
                (parent_selected, parent_addable)
            ] = poly_add(arbitrary, without_selected, -1)

        cache[signature] = (
            selected[vertex],
            no_selected_child[vertex],
            with_selected_child[vertex],
        )

    return poly_add(
        selected[root][(0, 0)],
        poly_add(
            no_selected_child[root][(0, 0)],
            with_selected_child[root][(0, 0)],
        ),
    )


def energy_deficit(
    distribution: Poly, rank: int
) -> Fraction:
    """Return the D14 aggregate deficit for C-rank ``rank-1``."""
    total = Fraction(0)
    for (set_rank, q), value in distribution.items():
        if set_rank != rank - 1 or q == 0:
            continue
        total += (
            q * (q - 1) * value.count
            + (q - 2) * value.edge_sum
            - (rank + 1) * value.degree_square_sum
            + Fraction(4 * rank * value.edge_square_sum, q)
        )
    return total


def energy_ratio(distribution: Poly, rank: int) -> Fraction:
    """Return ``r * Dirichlet numerator / extension-energy budget``."""
    dirichlet = Fraction(0)
    budget = Fraction(0)
    for (set_rank, q), value in distribution.items():
        if set_rank != rank - 1 or q == 0:
            continue
        dirichlet += (
            value.degree_square_sum
            - Fraction(4 * value.edge_square_sum, q)
        )
        budget += (
            q * (q - 1) * value.count
            + (q - 2) * value.edge_sum
            - value.degree_square_sum
        )
    if budget == 0:
        return Fraction(0)
    return rank * dirichlet / budget


def mean_only_ratio(
    distribution: Poly, polynomial: list[int], rank: int
) -> Fraction:
    """Return the ratio in the stronger proposal ``Dsum <= (r+1)i[r+1]``."""
    dirichlet = Fraction(0)
    for (set_rank, q), value in distribution.items():
        if set_rank != rank - 1 or q == 0:
            continue
        dirichlet += (
            value.degree_square_sum
            - Fraction(4 * value.edge_square_sum, q)
        )
    denominator = (rank + 1) * polynomial[rank + 1]
    return dirichlet / denominator if denominator else Fraction(0)


def make_double_regular_tree(branching: int, depth: int) -> list[list[int]]:
    adj = [[], []]
    adj[0].append(1)
    adj[1].append(0)
    for root in (0, 1):
        layer = [root]
        for _ in range(depth):
            next_layer: list[int] = []
            for vertex in layer:
                for _ in range(branching):
                    child = len(adj)
                    adj.append([vertex])
                    adj[vertex].append(child)
                    next_layer.append(child)
            layer = next_layer
    return adj


def make_double_spherical_tree(
    branching: tuple[int, ...]
) -> list[list[int]]:
    adj = [[], []]
    adj[0].append(1)
    adj[1].append(0)
    for root in (0, 1):
        layer = [root]
        for factor in branching:
            next_layer: list[int] = []
            for vertex in layer:
                for _ in range(factor):
                    child = len(adj)
                    adj.append([vertex])
                    adj[vertex].append(child)
                    next_layer.append(child)
            layer = next_layer
    return adj


def check_tree(adj: list[list[int]], label: str) -> None:
    polynomial = independence_poly(len(adj), adj)
    alpha = len(polynomial) - 1
    limit = (2 * alpha + 1) // 3
    distribution = energy_distribution(adj, max_rank=limit - 2)
    counts = [0] * len(polynomial)
    for (rank, _), value in distribution.items():
        if rank < len(counts):
            counts[rank] += value.count
    assert counts[: limit - 1] == polynomial[: limit - 1]
    minimum: tuple[Fraction, int] | None = None
    maximum_ratio: tuple[Fraction, int] | None = None
    maximum_mean_only_ratio: tuple[Fraction, int] | None = None
    for rank in range(1, limit - 1):
        deficit = energy_deficit(distribution, rank)
        ratio = energy_ratio(distribution, rank)
        mean_ratio = mean_only_ratio(distribution, polynomial, rank)
        if minimum is None or deficit < minimum[0]:
            minimum = (deficit, rank)
        if maximum_ratio is None or ratio > maximum_ratio[0]:
            maximum_ratio = (ratio, rank)
        if (
            maximum_mean_only_ratio is None
            or mean_ratio > maximum_mean_only_ratio[0]
        ):
            maximum_mean_only_ratio = (mean_ratio, rank)
        if deficit < 0:
            raise AssertionError((label, len(adj), alpha, rank, deficit))
    print(
        {
            "label": label,
            "n": len(adj),
            "alpha": alpha,
            "minimum": minimum,
            "maximum_ratio": maximum_ratio,
            "maximum_mean_only_ratio": maximum_mean_only_ratio,
        },
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--double-branching", type=int, default=2)
    parser.add_argument("--double-depth", type=int, default=3)
    parser.add_argument("--galvin-m", type=int, default=4)
    parser.add_argument("--galvin-t", type=int, default=3)
    parser.add_argument(
        "--spherical",
        help="comma-separated root-outward branching sequence",
    )
    args = parser.parse_args()
    if args.spherical:
        sequence = tuple(int(value) for value in args.spherical.split(","))
        check_tree(
            make_double_spherical_tree(sequence),
            f"double_spherical{sequence}",
        )
        return
    check_tree(
        make_double_regular_tree(
            args.double_branching, args.double_depth
        ),
        f"double_regular({args.double_branching},{args.double_depth})",
    )
    _, galvin = make_T_m_t_1(args.galvin_m, args.galvin_t)
    check_tree(galvin, f"T({args.galvin_m},{args.galvin_t},1)")


if __name__ == "__main__":
    main()
