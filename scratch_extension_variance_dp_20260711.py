#!/usr/bin/env python3
"""Exact rank-wise moments of the one-vertex extension count on a tree.

For an independent set S, let e(S) be the number of vertices that can be
adjoined to S.  This computes, for every rank r,

    z[r]  = # independent r-sets,
    d1[r] = sum_{|S|=r} e(S),
    d2[r] = sum_{|S|=r} e(S)(e(S)-1).

The DP is the y=1 two-jet of sum_S x^|S| y^e(S).  It is intended as a
falsification tool for the prefix conjecture Var(e) <= E[e], equivalently
z[r] * d2[r] <= d1[r]^2.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from math import ceil

from indpoly import _polyadd, _polymul, independence_poly
from trees import trees


Poly = list[int]


@dataclass(frozen=True)
class Jet:
    value: Poly
    first: Poly
    second: Poly


ZERO = Jet([0], [0], [0])
ONE = Jet([1], [0], [0])


def poly_sub(a: Poly, b: Poly) -> Poly:
    out = [0] * max(len(a), len(b))
    for i, value in enumerate(a):
        out[i] += value
    for i, value in enumerate(b):
        out[i] -= value
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def jet_add(a: Jet, b: Jet) -> Jet:
    return Jet(
        _polyadd(a.value, b.value),
        _polyadd(a.first, b.first),
        _polyadd(a.second, b.second),
    )


def jet_sub(a: Jet, b: Jet) -> Jet:
    return Jet(
        poly_sub(a.value, b.value),
        poly_sub(a.first, b.first),
        poly_sub(a.second, b.second),
    )


def jet_mul(a: Jet, b: Jet) -> Jet:
    value = _polymul(a.value, b.value)
    first = _polyadd(
        _polymul(a.first, b.value),
        _polymul(a.value, b.first),
    )
    second = _polyadd(
        _polyadd(
            _polymul(a.second, b.value),
            [2 * coefficient for coefficient in _polymul(a.first, b.first)],
        ),
        _polymul(a.value, b.second),
    )
    return Jet(value, first, second)


def jet_shift_x(a: Jet) -> Jet:
    return Jet([0] + a.value, [0] + a.first, [0] + a.second)


def jet_mul_y(a: Jet) -> Jet:
    # At y=1: (yf)'=f+f' and (yf)''=2f'+f''.
    return Jet(
        a.value,
        _polyadd(a.value, a.first),
        _polyadd([2 * coefficient for coefficient in a.first], a.second),
    )


def extension_moment_jets(adj: list[list[int]], root: int = 0) -> Jet:
    """Return the rank generating polynomials for 1, e, and e(e-1)."""
    n = len(adj)
    if n == 0:
        return ONE

    parent = [-2] * n
    parent[root] = -1
    order = [root]
    for vertex in order:
        for neighbor in adj[vertex]:
            if parent[neighbor] == -2:
                parent[neighbor] = vertex
                order.append(neighbor)

    selected = [ZERO] * n
    excluded_no_selected_child = [ZERO] * n
    excluded_with_selected_child = [ZERO] * n

    for vertex in reversed(order):
        no_selected_child = ONE
        arbitrary_children = ONE
        selected_children_blocked = ONE

        for child in adj[vertex]:
            if parent[child] != vertex:
                continue
            child_excluded_parent_free = jet_add(
                jet_mul_y(excluded_no_selected_child[child]),
                excluded_with_selected_child[child],
            )
            child_excluded_parent_selected = jet_add(
                excluded_no_selected_child[child],
                excluded_with_selected_child[child],
            )
            no_selected_child = jet_mul(
                no_selected_child, child_excluded_parent_free
            )
            arbitrary_children = jet_mul(
                arbitrary_children,
                jet_add(selected[child], child_excluded_parent_free),
            )
            selected_children_blocked = jet_mul(
                selected_children_blocked, child_excluded_parent_selected
            )

        excluded_no_selected_child[vertex] = no_selected_child
        excluded_with_selected_child[vertex] = jet_sub(
            arbitrary_children, no_selected_child
        )
        selected[vertex] = jet_shift_x(selected_children_blocked)

    return jet_add(
        selected[root],
        jet_add(
            jet_mul_y(excluded_no_selected_child[root]),
            excluded_with_selected_child[root],
        ),
    )


def first_prefix_failure(
    adj: list[list[int]],
) -> tuple[int, int, int, int, int] | None:
    return first_prefix_failure_from_jet(extension_moment_jets(adj))


def first_prefix_failure_from_jet(
    jet: Jet,
) -> tuple[int, int, int, int, int] | None:
    alpha = len(jet.value) - 1
    limit = ceil((2 * alpha - 1) / 3)
    for rank in range(limit - 1):
        lhs = jet.value[rank] * jet.second[rank]
        rhs = jet.first[rank] * jet.first[rank]
        if lhs > rhs:
            return rank, alpha, limit, lhs, rhs
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=1)
    parser.add_argument("--max-n", type=int, default=18)
    parser.add_argument("--backend", choices=("auto", "geng", "networkx"), default="auto")
    args = parser.parse_args()

    for n in range(args.min_n, args.max_n + 1):
        count = 0
        for _, adj in trees(n, backend=args.backend):
            jet = extension_moment_jets(adj)
            expected = independence_poly(n, adj)
            if jet.value != expected:
                raise AssertionError((n, jet.value, expected, adj))
            failure = first_prefix_failure_from_jet(jet)
            if failure is not None:
                print({"n": n, "failure": failure, "adj": adj})
                raise SystemExit(1)
            count += 1
        print({"n": n, "trees": count, "prefix_failures": 0}, flush=True)


if __name__ == "__main__":
    main()
