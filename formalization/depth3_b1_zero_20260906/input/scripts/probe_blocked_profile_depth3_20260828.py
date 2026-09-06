#!/usr/bin/env python3
"""Bounded exact probe of the blocked/extendable depth-3 decomposition."""

from __future__ import annotations

import argparse
import json
import sys
from itertools import product
from pathlib import Path

import networkx as nx

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from indpoly import independence_poly


DEFAULT_OUTPUT = Path("results/blocked_profile_depth3_probe_20260828.json")


def matching_bags(tree: nx.Graph) -> list[tuple[int, ...]]:
    matching = nx.max_weight_matching(tree, maxcardinality=True)
    mate: dict[int, int] = {}
    for left, right in matching:
        mate[left] = right
        mate[right] = left

    bags: list[tuple[int, ...]] = []
    seen: set[int] = set()
    for vertex in tree.nodes():
        if vertex in seen:
            continue
        bag = tuple(sorted((vertex, mate[vertex]))) if vertex in mate else (vertex,)
        bags.append(bag)
        seen.update(bag)
    return bags


def maximum_assignment_code(tree: nx.Graph, bags: list[tuple[int, ...]]) -> tuple[int, ...]:
    matched_indices = [index for index, bag in enumerate(bags) if len(bag) == 2]
    code: list[int] = []
    for choices in product((0, 1), repeat=len(matched_indices)):
        selected = {bag[0] for bag in bags if len(bag) == 1}
        word = 0
        for index, choice in zip(matched_indices, choices):
            selected.add(bags[index][choice])
            word |= choice << index
        if all(not (left in selected and right in selected) for left, right in tree.edges()):
            code.append(word)
    return tuple(code)


def projection_profile(code: tuple[int, ...], bag_count: int) -> list[int]:
    profile = [0] * (bag_count + 1)
    for mask in range(1 << bag_count):
        profile[mask.bit_count()] += len({word & mask for word in code})
    return profile


def brute_force_extendable_profile(tree: nx.Graph, max_defect: int = 4) -> list[int]:
    """Count independent sets contained in a maximum one, without bag data."""
    order = tree.number_of_nodes()
    all_bits = (1 << order) - 1
    edge_masks = [(1 << left) | (1 << right) for left, right in tree.edges()]
    independent = [
        mask
        for mask in range(1 << order)
        if all(mask & edge_mask != edge_mask for edge_mask in edge_masks)
    ]
    alpha = max(mask.bit_count() for mask in independent)
    maximum = [mask for mask in independent if mask.bit_count() == alpha]

    profile = []
    for defect in range(max_defect + 1):
        target_size = alpha - defect
        profile.append(
            sum(
                mask.bit_count() == target_size
                and any(mask & (all_bits ^ maximum_mask) == 0 for maximum_mask in maximum)
                for mask in independent
            )
        )
    return profile


def witness_row(
    tree: nx.Graph,
    bag_count: int,
    extendable: list[int],
    blocked: list[int],
    full: list[int],
    extendable_turan: int,
    blocked_turan: int,
    cross: int,
    relative_normalized_lc_margin: int,
) -> dict[str, object]:
    return {
        "n": tree.number_of_nodes(),
        "alpha_bag_count": bag_count,
        "graph6": nx.to_graph6_bytes(tree, header=False).decode().strip(),
        "e_0_to_4": extendable,
        "b_0_to_4": blocked,
        "s_0_to_4": full,
        "extendable_turan": extendable_turan,
        "blocked_turan": blocked_turan,
        "cross_term": cross,
        "relative_normalized_lc_margin": relative_normalized_lc_margin,
        "cross_reserve_margin": extendable_turan + cross,
        "full_margin": extendable_turan + blocked_turan + cross,
    }


def run_probe(max_n: int) -> dict[str, object]:
    counters = {
        "trees_checked": 0,
        "bruteforce_extendable_profile_checks": 0,
        "decomposition_identity_failures": 0,
        "blocked_profile_depth3_failures": 0,
        "cross_reserve_depth3_failures": 0,
        "cross_term_negative": 0,
        "negative_cross_with_b2_zero": 0,
        "negative_cross_with_b2_positive_b4_zero": 0,
        "negative_cross_with_b2_b4_positive": 0,
        "negative_cross_actual_reserve_without_positive_b3_failures": 0,
        "negative_cross_pascal_reserve_failures": 0,
        "negative_cross_pascal_reserve_without_positive_b3_failures": 0,
        "relative_normalized_lc_failures": 0,
        "negative_cross_with_relative_normalized_lc_failure": 0,
        "relative_normalized_lc_failures_with_nonnegative_cross": 0,
        "blocked_b3_positive": 0,
        "discard_positive_b3_bound_failures": 0,
        "full_depth3_failures": 0,
    }
    first: dict[str, dict[str, object]] = {}
    smallest_cross_reserve: dict[str, object] | None = None
    smallest_pascal_cross_reserve: dict[str, object] | None = None
    worst_pascal_load: dict[str, object] | None = None

    for order in range(2, max_n + 1):
        for tree in nx.generators.nonisomorphic_trees(order):
            bags = matching_bags(tree)
            bag_count = len(bags)
            if bag_count < 4:
                continue
            code = maximum_assignment_code(tree, bags)
            projection = projection_profile(code, bag_count)
            extendable = [projection[bag_count - defect] for defect in range(5)]

            adjacency = [list(tree.neighbors(vertex)) for vertex in range(order)]
            polynomial = independence_poly(order, adjacency)
            if len(polynomial) - 1 != bag_count:
                raise AssertionError("bag count does not equal independence number")
            full = [polynomial[bag_count - defect] for defect in range(5)]
            blocked = [full[index] - extendable[index] for index in range(5)]
            if any(value < 0 for value in blocked):
                raise AssertionError("negative blocked count")
            if order <= 10:
                brute_force_extendable = brute_force_extendable_profile(tree)
                if extendable != brute_force_extendable:
                    graph6 = nx.to_graph6_bytes(tree, header=False).decode().strip()
                    raise AssertionError(
                        "matching-bag projection disagrees with brute-force "
                        f"extendability for {graph6}: {extendable} != "
                        f"{brute_force_extendable}"
                    )
                counters["bruteforce_extendable_profile_checks"] += 1

            extendable_turan = extendable[3] ** 2 - extendable[2] * extendable[4]
            blocked_turan = blocked[3] ** 2 - blocked[2] * blocked[4]
            cross = (
                2 * extendable[3] * blocked[3]
                - extendable[2] * blocked[4]
                - extendable[4] * blocked[2]
            )
            relative_normalized_lc_margin = (
                blocked[3] ** 2 * extendable[2] * extendable[4]
                - blocked[2] * blocked[4] * extendable[3] ** 2
            )
            full_margin = full[3] ** 2 - full[2] * full[4]
            row = witness_row(
                tree,
                bag_count,
                extendable,
                blocked,
                full,
                extendable_turan,
                blocked_turan,
                cross,
                relative_normalized_lc_margin,
            )

            if full_margin != extendable_turan + blocked_turan + cross:
                counters["decomposition_identity_failures"] += 1
                first.setdefault("decomposition_identity_failure", row)
            if blocked_turan < 0:
                counters["blocked_profile_depth3_failures"] += 1
                first.setdefault("blocked_profile_depth3_failure", row)
            if extendable_turan + cross < 0:
                counters["cross_reserve_depth3_failures"] += 1
                first.setdefault("cross_reserve_depth3_failure", row)
            if cross < 0:
                counters["cross_term_negative"] += 1
                first.setdefault("cross_term_negative", row)
                if blocked[2] == 0:
                    counters["negative_cross_with_b2_zero"] += 1
                elif blocked[4] == 0:
                    counters["negative_cross_with_b2_positive_b4_zero"] += 1
                else:
                    counters["negative_cross_with_b2_b4_positive"] += 1
                    first.setdefault("negative_cross_with_b2_b4_positive", row)
                negative_endpoints = (
                    extendable[2] * blocked[4]
                    + extendable[4] * blocked[2]
                )
                actual_reserve_without_positive_b3 = (
                    extendable_turan - negative_endpoints
                )
                pascal_cross_reserve_numerator = (
                    (5 * bag_count + 17) * extendable[2] * extendable[4]
                    + 27 * (bag_count - 3) * cross
                )
                pascal_without_positive_b3_numerator = (
                    (5 * bag_count + 17) * extendable[2] * extendable[4]
                    - 27 * (bag_count - 3) * negative_endpoints
                )
                if actual_reserve_without_positive_b3 < 0:
                    counters[
                        "negative_cross_actual_reserve_without_positive_b3_failures"
                    ] += 1
                if pascal_cross_reserve_numerator < 0:
                    counters["negative_cross_pascal_reserve_failures"] += 1
                    first.setdefault("negative_cross_pascal_reserve_failure", row)
                if pascal_without_positive_b3_numerator < 0:
                    counters[
                        "negative_cross_pascal_reserve_without_positive_b3_failures"
                    ] += 1

                reserve_row = dict(row)
                reserve_row.update(
                    {
                        "negative_endpoint_terms": negative_endpoints,
                        "actual_reserve_without_positive_b3": (
                            actual_reserve_without_positive_b3
                        ),
                        "pascal_cross_reserve_numerator": (
                            pascal_cross_reserve_numerator
                        ),
                        "pascal_without_positive_b3_numerator": (
                            pascal_without_positive_b3_numerator
                        ),
                        "pascal_load_numerator": (
                            -27 * (bag_count - 3) * cross
                        ),
                        "pascal_load_denominator": (
                            (5 * bag_count + 17)
                            * extendable[2]
                            * extendable[4]
                        ),
                    }
                )
                if (
                    smallest_pascal_cross_reserve is None
                    or pascal_cross_reserve_numerator
                    < smallest_pascal_cross_reserve[
                        "pascal_cross_reserve_numerator"
                    ]
                ):
                    smallest_pascal_cross_reserve = reserve_row
                if worst_pascal_load is None or (
                    reserve_row["pascal_load_numerator"]
                    * worst_pascal_load["pascal_load_denominator"]
                    > worst_pascal_load["pascal_load_numerator"]
                    * reserve_row["pascal_load_denominator"]
                ):
                    worst_pascal_load = reserve_row
            if relative_normalized_lc_margin < 0:
                counters["relative_normalized_lc_failures"] += 1
                first.setdefault("relative_normalized_lc_failure", row)
                if cross >= 0:
                    counters[
                        "relative_normalized_lc_failures_with_nonnegative_cross"
                    ] += 1
            if cross < 0 and relative_normalized_lc_margin < 0:
                counters[
                    "negative_cross_with_relative_normalized_lc_failure"
                ] += 1
                first.setdefault(
                    "negative_cross_with_relative_normalized_lc_failure", row
                )
            if blocked[3] > 0:
                counters["blocked_b3_positive"] += 1
            negative_without_b3 = (
                extendable[2] * blocked[4]
                + extendable[4] * blocked[2]
                + blocked[2] * blocked[4]
            )
            if extendable_turan < negative_without_b3:
                counters["discard_positive_b3_bound_failures"] += 1
                first.setdefault("discard_positive_b3_bound_failure", row)
            if full_margin < 0:
                counters["full_depth3_failures"] += 1
                first.setdefault("full_depth3_failure", row)
            if (
                smallest_cross_reserve is None
                or row["cross_reserve_margin"] < smallest_cross_reserve["cross_reserve_margin"]
            ):
                smallest_cross_reserve = row
            counters["trees_checked"] += 1

    return {
        "certificate_date": "2026-08-28",
        "kind": "bounded_blocked_profile_depth3_probe",
        "claim_status": "bounded computational evidence, not a theorem",
        "scope": {
            "tree_orders": [2, max_n],
            "minimum_alpha_bag_count": 4,
            "enumeration": "networkx nonisomorphic trees",
            "arithmetic": "exact integers",
            "independent_cross_check": "brute-force containment in a maximum independent set through n=10",
        },
        "identities": {
            "split": "s_d = e_d + b_d",
            "margin": "s_3^2-s_2*s_4 = (e_3^2-e_2*e_4) + (b_3^2-b_2*b_4) + (2*e_3*b_3-e_2*b_4-e_4*b_2)",
            "relative_normalized_lc_margin": "b_3^2*e_2*e_4-b_2*b_4*e_3^2",
        },
        "abstract_kill_test": {
            "scope": (
                "Positive numerical profiles showing that the proved extendable "
                "reserve and blocked log-concavity alone do not imply R>=0 "
                "when C<0; these are not asserted to be tree profiles."
            ),
            "alpha": 5,
            "e_2_to_4": [10, 20, 20],
            "b_2_to_4": [1, 10, 100],
            "extendable_reserve_lhs": 21600,
            "extendable_reserve_rhs": 19200,
            "blocked_lc_margin": 0,
            "cross_term": -620,
            "relative_normalized_lc_margin": -20000,
        },
        "last_extended": "2026-09-04",
        "counters": counters,
        "first_witnesses": first,
        "smallest_cross_reserve_margin": smallest_cross_reserve,
        "smallest_pascal_cross_reserve_numerator": (
            smallest_pascal_cross_reserve
        ),
        "worst_pascal_reserve_load": worst_pascal_load,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=15)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.max_n < 5:
        raise SystemExit("--max-n must be at least 5")
    result = run_probe(args.max_n)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    counters = result["counters"]
    print(
        "PASS: "
        f"{counters['trees_checked']} trees; "
        f"blocked_LC_failures={counters['blocked_profile_depth3_failures']}; "
        f"cross_reserve_failures={counters['cross_reserve_depth3_failures']}; "
        f"full_depth3_failures={counters['full_depth3_failures']}"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
