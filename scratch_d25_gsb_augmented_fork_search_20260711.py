#!/usr/bin/env python3
"""Exact search for a GSB-augmented inverse-degree fork allocation.

For each center c, let

    Q_c = 2 [N J_c - sum_{u<v in N(c)} a_u a_v],
    B_c = sum_{u in N(c)} a_u^2/deg(u) + a_c sum_{u in N(c)} a_u,
    S_c = sum_{u in N(c)} a_u/deg(u).

The candidate one-parameter local inequality is

    Q_c <= B_c + N [lambda a_c + (1-lambda) S_c].

Both allocations sum to N*d1, so summing over c gives exactly the GSB
budget.  This driver intersects the exact rational lambda constraints.
"""

from __future__ import annotations

import argparse
import math
import random
from fractions import Fraction

from indpoly import _polyadd, independence_poly
from scratch_d18_covariance_search_20260711 import (
    cavity_messages,
    coeff,
    product,
    prufer_tree,
)
from trees import trees


def constraints(adj: list[list[int]], prefix_only: bool = True):
    poly = independence_poly(len(adj), adj)
    alpha = len(poly) - 1
    limit = math.ceil((2 * alpha - 1) / 3)
    total, excluded = cavity_messages(adj)
    one = [product([excluded[(neighbor, vertex)] for neighbor in adj[vertex]])
           for vertex in range(len(adj))]
    joint_by_center = [[0] for _ in adj]
    for center, neighbors in enumerate(adj):
        for left_index, left in enumerate(neighbors):
            for right in neighbors[left_index + 1 :]:
                factors = [
                    excluded[(other, left)]
                    for other in adj[left]
                    if other != center
                ]
                factors += [
                    excluded[(other, right)]
                    for other in adj[right]
                    if other != center
                ]
                factors += [
                    total[(other, center)]
                    for other in neighbors
                    if other not in (left, right)
                ]
                joint_by_center[center] = _polyadd(
                    joint_by_center[center], product(factors)
                )

    for rank, n_r in enumerate(poly):
        if prefix_only and not (rank >= 4 and rank <= limit - 2):
            continue
        a = [coeff(item, rank) for item in one]
        for center, neighbors in enumerate(adj):
            product_pairs = sum(
                a[left] * a[right]
                for index, left in enumerate(neighbors)
                for right in neighbors[index + 1 :]
            )
            q_ordered = 2 * (
                n_r * coeff(joint_by_center[center], rank) - product_pairs
            )
            neighbor_sum = sum(a[neighbor] for neighbor in neighbors)
            base = sum(
                Fraction(a[neighbor] ** 2, len(adj[neighbor]))
                for neighbor in neighbors
            ) + a[center] * neighbor_sum
            s_center = sum(
                Fraction(a[neighbor], len(adj[neighbor]))
                for neighbor in neighbors
            )
            constant_budget = base + n_r * s_center
            slope = n_r * (a[center] - s_center)
            residual = q_ordered - constant_budget
            if slope > 0:
                kind, bound = "lower", residual / slope
            elif slope < 0:
                kind, bound = "upper", residual / slope
            elif residual <= 0:
                kind, bound = "none", Fraction(0)
            else:
                kind, bound = "impossible", Fraction(0)
            yield {
                "alpha": alpha,
                "limit": limit,
                "rank": rank,
                "center": center,
                "N": n_r,
                "a_c": a[center],
                "S_c": s_center,
                "q_ordered": q_ordered,
                "base": base,
                "kind": kind,
                "bound": bound,
                "gap_lambda_half": (
                    base + n_r * (a[center] + s_center) / 2 - q_ordered
                ),
            }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exhaustive-max-n", type=int, default=15)
    parser.add_argument("--random", type=int, default=0)
    parser.add_argument("--random-max-n", type=int, default=100)
    parser.add_argument("--seed", type=int, default=25993)
    parser.add_argument("--all-ranks", action="store_true")
    args = parser.parse_args()
    lower, upper = Fraction(0), Fraction(1)
    lower_witness = upper_witness = None
    minimum_half = None
    tested_rows = 0

    def inspect(label: str, adj: list[list[int]]) -> None:
        nonlocal lower, upper, lower_witness, upper_witness, minimum_half, tested_rows
        for row in constraints(adj, prefix_only=not args.all_ranks):
            tested_rows += 1
            payload = (label, len(adj), row)
            if minimum_half is None or row["gap_lambda_half"] < minimum_half[0]:
                minimum_half = (row["gap_lambda_half"], payload)
            if row["kind"] == "lower" and row["bound"] > lower:
                lower, lower_witness = row["bound"], payload
            elif row["kind"] == "upper" and row["bound"] < upper:
                upper, upper_witness = row["bound"], payload
            elif row["kind"] == "impossible":
                lower, upper = Fraction(1), Fraction(0)
                lower_witness = upper_witness = payload
            if lower > upper:
                raise RuntimeError({
                    "empty_interval": (lower, upper),
                    "lower_witness": lower_witness,
                    "upper_witness": upper_witness,
                })

    for n in range(1, args.exhaustive_max_n + 1):
        count = 0
        for index, (_, adj) in enumerate(trees(n)):
            inspect(f"tree:{n}:{index}", adj)
            count += 1
        print({
            "n": n,
            "trees": count,
            "rows": tested_rows,
            "lambda_interval": (lower, upper),
            "float_interval": (float(lower), float(upper)),
            "min_half": minimum_half,
        }, flush=True)

    rng = random.Random(args.seed)
    for trial in range(args.random):
        n = rng.randint(2, args.random_max_n)
        inspect(f"prufer:{trial}", prufer_tree(n, rng))
    print({
        "complete": True,
        "rows": tested_rows,
        "lambda_interval": (lower, upper),
        "float_interval": (float(lower), float(upper)),
        "lower_witness": lower_witness,
        "upper_witness": upper_witness,
        "min_half": minimum_half,
    }, flush=True)


if __name__ == "__main__":
    main()
