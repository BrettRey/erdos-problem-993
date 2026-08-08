#!/usr/bin/env python3
"""Exact certificate for bounded pendant cores in the two-hub class.

For fixed pendant-arm multisets on two hubs, the connector polynomials obey

    F_t = F_{t-1} + x F_{t-2}.

The accompanying note proves that log-concavity for every connector length
follows once the adjacent-hub polynomial B and contracted-hub polynomial C
are log-concave and satisfy C ~_p B and xC ~_p B.  This script exhaustively
checks those hypotheses for every unordered pair of integer partitions whose
total weight is at most the requested bound.  All arithmetic is integral.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Iterator

from indpoly import independence_poly, is_log_concave, is_unimodal
from verify_connector_partial_sync_route_20260808 import path_poly, shift
from verify_double_star_log_concavity_20260808 import poly_add, poly_mul
from verify_two_hub_ratio_order_obstructions_20260808 import (
    partial_synchronization_failure,
)

Poly = list[int]
ArmVector = tuple[int, ...]
State = tuple[ArmVector, Poly, Poly]


def integer_partitions(total: int, minimum: int = 1) -> Iterator[ArmVector]:
    """Yield nondecreasing positive-part partitions of ``total``."""
    if total == 0:
        yield ()
        return
    for first in range(minimum, total + 1):
        for rest in integer_partitions(total - first, first):
            yield (first,) + rest


def poly_product(polynomials: Iterator[Poly]) -> Poly:
    out = [1]
    for polynomial in polynomials:
        out = poly_mul(out, polynomial)
    return out


def rooted_arm_state(arms: ArmVector) -> tuple[Poly, Poly]:
    """Return (Q,R): hub excluded and included arm factors.

    An arm of length a contributes P_a when the hub is excluded and P_{a-1}
    when the hub is included.
    """
    excluded = poly_product(path_poly(length) for length in arms)
    included = poly_product(path_poly(length - 1) for length in arms)
    return excluded, included


def adjacent_and_contracted(
    left_state: State,
    right_state: State,
) -> tuple[Poly, Poly]:
    """Return B=F_0 and C for the contracted pair of hubs."""
    _, q_left, r_left = left_state
    _, q_right, r_right = right_state
    q_product = poly_mul(q_left, q_right)
    mixed = poly_add(
        poly_mul(r_left, q_right),
        poly_mul(q_left, r_right),
    )
    included_product = poly_mul(r_left, r_right)
    adjacent = poly_add(q_product, shift(mixed))
    contracted = poly_add(q_product, shift(included_product))
    return adjacent, contracted


def connector_formula(
    left_state: State,
    right_state: State,
    internal_vertices: int,
) -> Poly:
    """Return the exact independence polynomial for the two-hub tree."""
    _, q_left, r_left = left_state
    _, q_right, r_right = right_state
    q_product = poly_mul(q_left, q_right)
    mixed = poly_add(
        poly_mul(r_left, q_right),
        poly_mul(q_left, r_right),
    )
    included_product = poly_mul(r_left, r_right)
    return poly_add(
        poly_add(
            poly_mul(q_product, path_poly(internal_vertices)),
            shift(
                poly_mul(mixed, path_poly(internal_vertices - 1)),
            ),
        ),
        shift(
            poly_mul(
                included_product,
                path_poly(internal_vertices - 2),
            ),
            2,
        ),
    )


def connector_adjacency(
    left_arms: ArmVector,
    right_arms: ArmVector,
    internal_vertices: int,
) -> list[list[int]]:
    """Build the corresponding tree for independent DP replay."""
    order = 2 + internal_vertices + sum(left_arms) + sum(right_arms)
    adjacency = [[] for _ in range(order)]

    def add_edge(left: int, right: int) -> None:
        adjacency[left].append(right)
        adjacency[right].append(left)

    right_hub = internal_vertices + 1
    for vertex in range(right_hub):
        add_edge(vertex, vertex + 1)

    next_vertex = right_hub + 1
    for hub, arms in ((0, left_arms), (right_hub, right_arms)):
        for length in arms:
            previous = hub
            for _ in range(length):
                add_edge(previous, next_vertex)
                previous = next_vertex
                next_vertex += 1
    assert next_vertex == order
    return adjacency


def states_by_weight(maximum_weight: int) -> dict[int, list[State]]:
    by_weight: dict[int, list[State]] = defaultdict(list)
    for weight in range(maximum_weight + 1):
        for arms in integer_partitions(weight):
            excluded, included = rooted_arm_state(arms)
            by_weight[weight].append((arms, excluded, included))
    return dict(by_weight)


def verify(
    maximum_pendant_vertices: int,
    dp_maximum_pendant_vertices: int,
    dp_maximum_internal_vertices: int,
) -> dict[str, object]:
    by_weight = states_by_weight(maximum_pendant_vertices)
    digest = hashlib.sha256()
    pairs_by_total_weight = [0] * (maximum_pendant_vertices + 1)
    pair_count = 0

    for left_weight in range(maximum_pendant_vertices + 1):
        for right_weight in range(
            left_weight,
            maximum_pendant_vertices + 1 - left_weight,
        ):
            left_states = by_weight[left_weight]
            right_states = by_weight[right_weight]
            for left_index, left_state in enumerate(left_states):
                right_start = left_index if left_weight == right_weight else 0
                for right_state in right_states[right_start:]:
                    adjacent, contracted = adjacent_and_contracted(
                        left_state,
                        right_state,
                    )
                    assert is_log_concave(adjacent)
                    assert is_unimodal(adjacent)
                    assert is_log_concave(contracted)
                    assert is_unimodal(contracted)
                    assert partial_synchronization_failure(
                        contracted,
                        adjacent,
                    ) is None
                    assert partial_synchronization_failure(
                        shift(contracted),
                        adjacent,
                    ) is None

                    canonical_record = [
                        list(left_state[0]),
                        list(right_state[0]),
                        adjacent,
                        contracted,
                    ]
                    digest.update(
                        json.dumps(
                            canonical_record,
                            separators=(",", ":"),
                        ).encode("ascii")
                    )
                    digest.update(b"\n")
                    total_weight = left_weight + right_weight
                    pairs_by_total_weight[total_weight] += 1
                    pair_count += 1

    # Independent tree-DP and recurrence audit on a smaller complete box.
    dp_bound = min(dp_maximum_pendant_vertices, maximum_pendant_vertices)
    dp_states = states_by_weight(dp_bound)
    dp_cases = 0
    recurrence_cases = 0
    for left_weight in range(dp_bound + 1):
        for right_weight in range(left_weight, dp_bound + 1 - left_weight):
            for left_index, left_state in enumerate(dp_states[left_weight]):
                right_start = left_index if left_weight == right_weight else 0
                for right_state in dp_states[right_weight][right_start:]:
                    previous_previous: Poly | None = None
                    previous: Poly | None = None
                    for internal in range(dp_maximum_internal_vertices + 1):
                        formula = connector_formula(
                            left_state,
                            right_state,
                            internal,
                        )
                        adjacency = connector_adjacency(
                            left_state[0],
                            right_state[0],
                            internal,
                        )
                        assert independence_poly(len(adjacency), adjacency) == formula
                        assert is_log_concave(formula)
                        assert is_unimodal(formula)
                        dp_cases += 1
                        if internal >= 2:
                            assert previous_previous is not None
                            assert previous is not None
                            assert formula == poly_add(
                                previous,
                                shift(previous_previous),
                            )
                            recurrence_cases += 1
                        previous_previous, previous = previous, formula

    return {
        "claim_scope": (
            "Computer-assisted theorem: every two-hub tree whose two "
            "pendant-arm multisets have total weight at most "
            f"{maximum_pendant_vertices} has a log-concave independence "
            "polynomial for every number of internal connector vertices."
        ),
        "theorem_vs_computation": {
            "algebraic_part": (
                "The connector formula, Fibonacci recurrence, and partial-"
                "synchronization induction cover every connector length."
            ),
            "exact_computational_part": (
                "The finite integer-partition enumeration verifies the two "
                "base polynomials and two base relations for every bounded "
                "pendant core."
            ),
            "not_claimed": (
                "No conclusion is made for pendant cores above the stated "
                "bound, for trees with three or more branch vertices, or "
                "for Erdos Problem 993 in full."
            ),
        },
        "maximum_total_pendant_vertices": maximum_pendant_vertices,
        "partition_vectors": sum(len(states) for states in by_weight.values()),
        "unordered_arm_vector_pairs": pair_count,
        "pairs_by_total_pendant_weight": pairs_by_total_weight,
        "canonical_enumeration_sha256": digest.hexdigest(),
        "base_checks_per_pair": [
            "B is log-concave and unimodal",
            "C is log-concave and unimodal",
            "C is partially synchronized with B",
            "xC is partially synchronized with B",
        ],
        "dp_audit": {
            "maximum_total_pendant_vertices": dp_bound,
            "maximum_internal_connector_vertices": (
                dp_maximum_internal_vertices
            ),
            "tree_dp_formula_cases": dp_cases,
            "recurrence_cases": recurrence_cases,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-pendant", type=int, default=24)
    parser.add_argument("--dp-max-pendant", type=int, default=8)
    parser.add_argument("--dp-max-internal", type=int, default=4)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/c2_bounded_pendant_core_20260808.json"),
    )
    args = parser.parse_args()
    report = verify(
        args.max_pendant,
        args.dp_max_pendant,
        args.dp_max_internal,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"event": "verified", "out": str(args.out), **report}))


if __name__ == "__main__":
    main()
