#!/usr/bin/env python3
"""Probe closure and depth-3 ratios of canonical obstruction classes."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from itertools import combinations, product
from pathlib import Path

import networkx as nx

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from indpoly import independence_poly
from scripts.probe_blocked_profile_depth3_20260828 import (
    matching_bags,
    maximum_assignment_code,
    projection_profile,
)
from scripts.probe_first_violation_partition_20260904 import (
    code_value,
    forbidden_certificates,
    independent_mask,
    quotient_skeleton,
    violated_certificates,
)


DEFAULT_OUTPUT = Path("results/obstruction_class_structure_probe_20260904.json")


Certificate = tuple[object, ...]
Assignment = tuple[int, ...]


def certificate_key(certificate: Certificate, skeleton: nx.Graph) -> tuple[object, ...]:
    """Put unary obstructions first, then pair obstructions by path length."""
    if certificate[0] == 0:
        _, index, value = certificate
        return (0, index, value)
    _, left, left_value, right, right_value = certificate
    return (
        1,
        nx.shortest_path_length(skeleton, left, right),
        left,
        left_value,
        right,
        right_value,
    )


def selected_mask(
    assignment: Assignment, bags: list[tuple[int, ...]]
) -> int:
    mask = 0
    for index, value in enumerate(assignment):
        if value >= 0:
            mask |= 1 << bags[index][value]
    return mask


def induced_independence_poly(
    tree: nx.Graph, removed: set[int]
) -> list[int]:
    remaining = [vertex for vertex in tree.nodes() if vertex not in removed]
    relabel = {vertex: index for index, vertex in enumerate(remaining)}
    adjacency = [
        [relabel[neighbor] for neighbor in tree.neighbors(vertex) if neighbor in relabel]
        for vertex in remaining
    ]
    return independence_poly(len(remaining), adjacency)


def coefficient(polynomial: list[int], degree: int) -> int:
    if degree < 0 or degree >= len(polynomial):
        return 0
    return polynomial[degree]


def profile_margins(profile: list[int], extendable: list[int]) -> tuple[int, int]:
    class_lc = profile[3] ** 2 - profile[2] * profile[4]
    relative = (
        profile[3] ** 2 * extendable[2] * extendable[4]
        - profile[2] * profile[4] * extendable[3] ** 2
    )
    return class_lc, relative


def witness_row(
    *,
    tree: nx.Graph,
    bags: list[tuple[int, ...]],
    certificate: Certificate | None,
    profile: list[int],
    extendable: list[int],
    blocked: list[int],
    cross: int,
    class_lc: int,
    relative: int,
    prefix_index: int | None = None,
) -> dict[str, object]:
    row: dict[str, object] = {
        "n": tree.number_of_nodes(),
        "graph6": nx.to_graph6_bytes(tree, header=False).decode().strip(),
        "bags": [list(bag) for bag in bags],
        "profile_2_to_4": profile[2:5],
        "extendable_2_to_4": extendable[2:5],
        "blocked_2_to_4": blocked[2:5],
        "full_cross_term": cross,
        "class_lc_margin": class_lc,
        "relative_ratio_margin": relative,
    }
    if certificate is not None:
        row["certificate"] = list(certificate)
    if prefix_index is not None:
        row["prefix_index"] = prefix_index
    return row


def run_probe(max_n: int) -> dict[str, object]:
    counters: Counter[str] = Counter()
    for key in (
        "blocked_fillings_missing",
        "prefix_fill_failures",
        "same_class_fill_failures",
    ):
        counters[key] = 0
    by_defect = {
        str(defect): Counter(
            feasible=0,
            extendable=0,
            blocked=0,
            feasible_fillings=0,
            same_class_fill_failures=0,
            prefix_fill_failures=0,
        )
        for defect in range(2, 5)
    }
    assigned_certificate_kind: Counter[str] = Counter()
    assigned_certificate_kind_on_negative_cross_trees: Counter[str] = Counter()
    assigned_pair_path_length: Counter[str] = Counter()
    class_tests: Counter[str] = Counter()
    prefix_tests: Counter[str] = Counter()
    pairwise_class_tests: Counter[str] = Counter()
    atomic_group_tests: Counter[str] = Counter()
    coarse_group_tests: Counter[str] = Counter()
    atomic_group_pair_tests: Counter[str] = Counter()
    coarse_group_pair_tests: Counter[str] = Counter()
    negative_cross_group_signatures: Counter[str] = Counter()
    formula_checks: Counter[str] = Counter()
    for key in (
        "pair_group_formula_failures",
        "unary_class_formula_failures",
        "unary_group_formula_failures",
    ):
        formula_checks[key] = 0
    for test_counter in (
        class_tests,
        prefix_tests,
        pairwise_class_tests,
        atomic_group_tests,
        coarse_group_tests,
        atomic_group_pair_tests,
        coarse_group_pair_tests,
    ):
        for key in (
            "lc_failures",
            "relative_ratio_failures",
            "relative_ratio_failures_on_negative_cross_trees",
            "two_sided_profiles",
            "all_levels_positive_profiles",
        ):
            test_counter[key] = 0
    first: dict[str, object] = {}

    for order in range(2, max_n + 1):
        for tree in nx.generators.nonisomorphic_trees(order):
            bags = matching_bags(tree)
            alpha = len(bags)
            if alpha < 4:
                continue

            code = maximum_assignment_code(tree, bags)
            projection = projection_profile(code, alpha)
            extendable = [projection[alpha - defect] for defect in range(5)]
            polynomial = independence_poly(
                order, [list(tree.neighbors(vertex)) for vertex in range(order)]
            )
            full = [polynomial[alpha - defect] for defect in range(5)]
            blocked = [full[index] - extendable[index] for index in range(5)]
            cross = (
                2 * extendable[3] * blocked[3]
                - extendable[2] * blocked[4]
                - extendable[4] * blocked[2]
            )
            full_relative = (
                blocked[3] ** 2 * extendable[2] * extendable[4]
                - blocked[2] * blocked[4] * extendable[3] ** 2
            )
            adverse = cross < 0
            counters["negative_cross_trees"] += adverse
            counters["full_relative_ratio_failures"] += full_relative < 0
            counters[
                "negative_cross_with_full_relative_ratio_failure"
            ] += adverse and full_relative < 0

            skeleton = quotient_skeleton(tree, bags)
            certificates = tuple(
                sorted(
                    forbidden_certificates(code, bags),
                    key=lambda certificate: certificate_key(certificate, skeleton),
                )
            )
            certificate_rank = {
                certificate: rank for rank, certificate in enumerate(certificates)
            }
            unary_certificates = [
                certificate for certificate in certificates if certificate[0] == 0
            ]
            forbidden_vertices = [
                bags[index][value]
                for _, index, value in unary_certificates
            ]
            avoiding_polynomial = induced_independence_poly(
                tree, set(forbidden_vertices)
            )
            avoiding_profile = [
                coefficient(avoiding_polynomial, alpha - defect)
                for defect in range(5)
            ]
            adjacency_masks = [
                sum(1 << neighbor for neighbor in tree.neighbors(vertex))
                for vertex in range(order)
            ]

            records: dict[int, dict[Assignment, int]] = {
                defect: {} for defect in range(2, 5)
            }
            class_profiles: dict[Certificate, list[int]] = defaultdict(
                lambda: [0] * 5
            )
            observed = {
                defect: Counter(feasible=0, extendable=0, blocked=0)
                for defect in range(2, 5)
            }

            for defect in range(2, 5):
                if defect > alpha:
                    continue
                for empty_tuple in combinations(range(alpha), defect):
                    empty = set(empty_tuple)
                    occupied = [index for index in range(alpha) if index not in empty]
                    for chosen_values in product(
                        *(range(len(bags[index])) for index in occupied)
                    ):
                        chosen_vertices = tuple(
                            bags[index][value]
                            for index, value in zip(occupied, chosen_values)
                        )
                        if not independent_mask(chosen_vertices, adjacency_masks):
                            continue

                        observed[defect]["feasible"] += 1
                        by_defect[str(defect)]["feasible"] += 1
                        counters["feasible_assignments"] += 1
                        partial = dict(zip(occupied, chosen_values))
                        if any(
                            all(
                                code_value(word, bags, index) == value
                                for index, value in partial.items()
                            )
                            for word in code
                        ):
                            observed[defect]["extendable"] += 1
                            by_defect[str(defect)]["extendable"] += 1
                            counters["extendable_assignments"] += 1
                            continue

                        violations = violated_certificates(partial, certificates)
                        if not violations:
                            graph6 = nx.to_graph6_bytes(
                                tree, header=False
                            ).decode().strip()
                            raise AssertionError(
                                f"blocked assignment has no certificate: {graph6}"
                            )
                        first_certificate = min(
                            violations, key=certificate_rank.__getitem__
                        )
                        rank = certificate_rank[first_certificate]
                        assignment = tuple(partial.get(index, -1) for index in range(alpha))
                        records[defect][assignment] = rank
                        class_profiles[first_certificate][defect] += 1
                        observed[defect]["blocked"] += 1
                        by_defect[str(defect)]["blocked"] += 1
                        counters["blocked_assignments"] += 1
                        if first_certificate[0] == 0:
                            assigned_certificate_kind["unary"] += 1
                            if adverse:
                                assigned_certificate_kind_on_negative_cross_trees[
                                    "unary"
                                ] += 1
                        else:
                            assigned_certificate_kind["pair"] += 1
                            if adverse:
                                assigned_certificate_kind_on_negative_cross_trees[
                                    "pair"
                                ] += 1
                            _, left, _, right, _ = first_certificate
                            path_length = nx.shortest_path_length(
                                skeleton, left, right
                            )
                            assigned_pair_path_length[str(path_length)] += 1

            for defect in range(2, 5):
                if observed[defect]["feasible"] != full[defect]:
                    raise AssertionError("full-profile reconstruction failed")
                if observed[defect]["extendable"] != extendable[defect]:
                    raise AssertionError("extendable-profile reconstruction failed")
                if observed[defect]["blocked"] != blocked[defect]:
                    raise AssertionError("blocked-profile reconstruction failed")

            for defect in (3, 4):
                for assignment, rank in records[defect].items():
                    mask = selected_mask(assignment, bags)
                    has_feasible_filling = False
                    has_same_class_failure = False
                    for index, old_value in enumerate(assignment):
                        if old_value >= 0:
                            continue
                        for value, vertex in enumerate(bags[index]):
                            if adjacency_masks[vertex] & mask:
                                continue
                            has_feasible_filling = True
                            counters["feasible_one_bag_fillings"] += 1
                            by_defect[str(defect)]["feasible_fillings"] += 1
                            filled = list(assignment)
                            filled[index] = value
                            filled_tuple = tuple(filled)
                            target_rank = records[defect - 1].get(filled_tuple)
                            if target_rank is None:
                                first.setdefault(
                                    "blocked_filling_missing",
                                    {
                                        "n": order,
                                        "graph6": nx.to_graph6_bytes(
                                            tree, header=False
                                        ).decode().strip(),
                                        "defect": defect,
                                        "assignment": list(assignment),
                                        "filled_assignment": filled,
                                        "source_rank": rank,
                                    },
                                )
                                counters["blocked_fillings_missing"] += 1
                                continue
                            if target_rank != rank:
                                counters["same_class_fill_failures"] += 1
                                by_defect[str(defect)][
                                    "same_class_fill_failures"
                                ] += 1
                                has_same_class_failure = True
                                first.setdefault(
                                    "same_class_fill_failure",
                                    {
                                        "n": order,
                                        "graph6": nx.to_graph6_bytes(
                                            tree, header=False
                                        ).decode().strip(),
                                        "defect": defect,
                                        "assignment": list(assignment),
                                        "filled_assignment": filled,
                                        "source_rank": rank,
                                        "target_rank": target_rank,
                                        "source_certificate": list(certificates[rank]),
                                        "target_certificate": list(
                                            certificates[target_rank]
                                        ),
                                    },
                                )
                            if target_rank > rank:
                                counters["prefix_fill_failures"] += 1
                                by_defect[str(defect)][
                                    "prefix_fill_failures"
                                ] += 1
                                first.setdefault(
                                    "prefix_fill_failure",
                                    {
                                        "n": order,
                                        "graph6": nx.to_graph6_bytes(
                                            tree, header=False
                                        ).decode().strip(),
                                        "defect": defect,
                                        "source_rank": rank,
                                        "target_rank": target_rank,
                                    },
                                )
                    if not has_feasible_filling:
                        counters["assignments_without_feasible_filling"] += 1
                    if has_same_class_failure:
                        counters[
                            "assignments_with_same_class_fill_failure"
                        ] += 1

            ordered_profiles = [
                (certificate, class_profiles[certificate])
                for certificate in certificates
                if certificate in class_profiles
            ]
            counters["tree_local_classes"] += len(ordered_profiles)
            cumulative = [0] * 5
            atomic_groups: dict[str, list[int]] = defaultdict(lambda: [0] * 5)
            coarse_groups: dict[str, list[int]] = defaultdict(lambda: [0] * 5)
            for prefix_index, (certificate, profile) in enumerate(ordered_profiles):
                class_lc, relative = profile_margins(profile, extendable)
                class_tests["profiles"] += 1
                class_tests["two_sided_profiles"] += (
                    profile[2] > 0 and profile[4] > 0
                )
                class_tests["all_levels_positive_profiles"] += all(
                    profile[defect] > 0 for defect in range(2, 5)
                )
                class_tests["lc_failures"] += class_lc < 0
                class_tests["relative_ratio_failures"] += relative < 0
                class_tests["profiles_on_negative_cross_trees"] += adverse
                class_tests[
                    "lc_failures_on_negative_cross_trees"
                ] += adverse and class_lc < 0
                class_tests[
                    "relative_ratio_failures_on_negative_cross_trees"
                ] += adverse and relative < 0
                if class_lc < 0:
                    first.setdefault(
                        "class_lc_failure",
                        witness_row(
                            tree=tree,
                            bags=bags,
                            certificate=certificate,
                            profile=profile,
                            extendable=extendable,
                            blocked=blocked,
                            cross=cross,
                            class_lc=class_lc,
                            relative=relative,
                        ),
                    )
                if relative < 0:
                    first.setdefault(
                        "class_relative_ratio_failure",
                        witness_row(
                            tree=tree,
                            bags=bags,
                            certificate=certificate,
                            profile=profile,
                            extendable=extendable,
                            blocked=blocked,
                            cross=cross,
                            class_lc=class_lc,
                            relative=relative,
                        ),
                    )
                    if adverse:
                        first.setdefault(
                            "class_relative_ratio_failure_on_negative_cross_tree",
                            witness_row(
                                tree=tree,
                                bags=bags,
                                certificate=certificate,
                                profile=profile,
                                extendable=extendable,
                                blocked=blocked,
                                cross=cross,
                                class_lc=class_lc,
                                relative=relative,
                            ),
                        )

                for defect in range(2, 5):
                    cumulative[defect] += profile[defect]
                cumulative_lc, cumulative_relative = profile_margins(
                    cumulative, extendable
                )
                prefix_tests["profiles"] += 1
                prefix_tests["two_sided_profiles"] += (
                    cumulative[2] > 0 and cumulative[4] > 0
                )
                prefix_tests["all_levels_positive_profiles"] += all(
                    cumulative[defect] > 0 for defect in range(2, 5)
                )
                prefix_tests["lc_failures"] += cumulative_lc < 0
                prefix_tests[
                    "relative_ratio_failures"
                ] += cumulative_relative < 0
                prefix_tests[
                    "profiles_on_negative_cross_trees"
                ] += adverse
                prefix_tests[
                    "lc_failures_on_negative_cross_trees"
                ] += adverse and cumulative_lc < 0
                prefix_tests[
                    "relative_ratio_failures_on_negative_cross_trees"
                ] += adverse and cumulative_relative < 0
                if cumulative_relative < 0 and adverse:
                    first.setdefault(
                        "prefix_relative_ratio_failure_on_negative_cross_tree",
                        witness_row(
                            tree=tree,
                            bags=bags,
                            certificate=certificate,
                            profile=cumulative.copy(),
                            extendable=extendable,
                            blocked=blocked,
                            cross=cross,
                            class_lc=cumulative_lc,
                            relative=cumulative_relative,
                            prefix_index=prefix_index,
                        ),
                    )
                if cumulative_relative < 0:
                    first.setdefault(
                        "prefix_relative_ratio_failure",
                        witness_row(
                            tree=tree,
                            bags=bags,
                            certificate=certificate,
                            profile=cumulative.copy(),
                            extendable=extendable,
                            blocked=blocked,
                            cross=cross,
                            class_lc=cumulative_lc,
                            relative=cumulative_relative,
                            prefix_index=prefix_index,
                        ),
                    )

                if certificate[0] == 0:
                    atomic_key = "unary"
                    coarse_key = "unary"
                else:
                    _, left, _, right, _ = certificate
                    path_length = nx.shortest_path_length(skeleton, left, right)
                    atomic_key = f"pair_path_{path_length}"
                    coarse_key = "pair"
                for defect in range(2, 5):
                    atomic_groups[atomic_key][defect] += profile[defect]
                    coarse_groups[coarse_key][defect] += profile[defect]

            if cumulative[2:5] != blocked[2:5]:
                graph6 = nx.to_graph6_bytes(tree, header=False).decode().strip()
                raise AssertionError(
                    "class profiles do not reconstruct blocked profile: "
                    f"{graph6}, {cumulative[2:5]} != {blocked[2:5]}"
                )

            earlier_forbidden: set[int] = set()
            for certificate in unary_certificates:
                _, index, value = certificate
                vertex = bags[index][value]
                removed = set(tree.neighbors(vertex)) | {vertex} | earlier_forbidden
                residual_polynomial = induced_independence_poly(tree, removed)
                predicted_profile = [
                    coefficient(residual_polynomial, alpha - defect - 1)
                    for defect in range(5)
                ]
                observed_profile = class_profiles.get(certificate, [0] * 5)
                formula_checks["unary_classes_checked"] += 1
                if predicted_profile[2:5] != observed_profile[2:5]:
                    formula_checks["unary_class_formula_failures"] += 1
                    first.setdefault(
                        "unary_class_formula_failure",
                        {
                            "n": order,
                            "graph6": nx.to_graph6_bytes(
                                tree, header=False
                            ).decode().strip(),
                            "certificate": list(certificate),
                            "predicted_2_to_4": predicted_profile[2:5],
                            "observed_2_to_4": observed_profile[2:5],
                        },
                    )
                earlier_forbidden.add(vertex)

            predicted_unary_group = [
                full[defect] - avoiding_profile[defect] for defect in range(5)
            ]
            predicted_pair_group = [
                avoiding_profile[defect] - extendable[defect]
                for defect in range(5)
            ]
            observed_unary_group = coarse_groups.get("unary", [0] * 5)
            observed_pair_group = coarse_groups.get("pair", [0] * 5)
            formula_checks["coarse_groups_checked"] += 2
            if predicted_unary_group[2:5] != observed_unary_group[2:5]:
                formula_checks["unary_group_formula_failures"] += 1
                raise AssertionError("unary group formula failed")
            if predicted_pair_group[2:5] != observed_pair_group[2:5]:
                formula_checks["pair_group_formula_failures"] += 1
                raise AssertionError("pair group formula failed")

            decomposed_relative = sum(
                profile_margins(profile, extendable)[1]
                for _, profile in ordered_profiles
            )
            for left_index, (left_certificate, left_profile) in enumerate(
                ordered_profiles
            ):
                for right_certificate, right_profile in ordered_profiles[
                    left_index + 1 :
                ]:
                    pair_margin = (
                        2
                        * left_profile[3]
                        * right_profile[3]
                        * extendable[2]
                        * extendable[4]
                        - (
                            left_profile[2] * right_profile[4]
                            + right_profile[2] * left_profile[4]
                        )
                        * extendable[3] ** 2
                    )
                    decomposed_relative += pair_margin
                    pairwise_class_tests["profiles"] += 1
                    pairwise_class_tests["relative_ratio_failures"] += (
                        pair_margin < 0
                    )
                    pairwise_class_tests[
                        "profiles_on_negative_cross_trees"
                    ] += adverse
                    pairwise_class_tests[
                        "relative_ratio_failures_on_negative_cross_trees"
                    ] += adverse and pair_margin < 0
                    if pair_margin < 0:
                        pair_witness = {
                            "n": order,
                            "graph6": nx.to_graph6_bytes(
                                tree, header=False
                            ).decode().strip(),
                            "bags": [list(bag) for bag in bags],
                            "left_certificate": list(left_certificate),
                            "left_profile_2_to_4": left_profile[2:5],
                            "right_certificate": list(right_certificate),
                            "right_profile_2_to_4": right_profile[2:5],
                            "extendable_2_to_4": extendable[2:5],
                            "blocked_2_to_4": blocked[2:5],
                            "full_cross_term": cross,
                            "pair_relative_ratio_margin": pair_margin,
                        }
                        first.setdefault(
                            "pairwise_class_relative_ratio_failure", pair_witness
                        )
                        if adverse:
                            first.setdefault(
                                "pairwise_class_relative_ratio_failure_on_negative_cross_tree",
                                pair_witness,
                            )
            if decomposed_relative != full_relative:
                raise AssertionError(
                    "class-pair decomposition of relative margin failed"
                )

            for test_counter, groups in (
                (atomic_group_tests, atomic_groups),
                (coarse_group_tests, coarse_groups),
            ):
                for group_key, profile in groups.items():
                    group_lc, group_relative = profile_margins(profile, extendable)
                    test_counter["profiles"] += 1
                    test_counter["two_sided_profiles"] += (
                        profile[2] > 0 and profile[4] > 0
                    )
                    test_counter["all_levels_positive_profiles"] += all(
                        profile[defect] > 0 for defect in range(2, 5)
                    )
                    test_counter["lc_failures"] += group_lc < 0
                    test_counter["relative_ratio_failures"] += group_relative < 0
                    test_counter["profiles_on_negative_cross_trees"] += adverse
                    test_counter[
                        "lc_failures_on_negative_cross_trees"
                    ] += adverse and group_lc < 0
                    test_counter[
                        "relative_ratio_failures_on_negative_cross_trees"
                    ] += adverse and group_relative < 0
                    if group_relative < 0 and adverse:
                        first.setdefault(
                            f"{group_key}_relative_ratio_failure_on_negative_cross_tree",
                            witness_row(
                                tree=tree,
                                bags=bags,
                                certificate=None,
                                profile=profile,
                                extendable=extendable,
                                blocked=blocked,
                                cross=cross,
                                class_lc=group_lc,
                                relative=group_relative,
                            ),
                        )

            for test_counter, groups in (
                (atomic_group_pair_tests, atomic_groups),
                (coarse_group_pair_tests, coarse_groups),
            ):
                group_items = sorted(groups.items())
                for left_index, (left_key, left_profile) in enumerate(group_items):
                    for right_key, right_profile in group_items[left_index + 1 :]:
                        pair_margin = (
                            2
                            * left_profile[3]
                            * right_profile[3]
                            * extendable[2]
                            * extendable[4]
                            - (
                                left_profile[2] * right_profile[4]
                                + right_profile[2] * left_profile[4]
                            )
                            * extendable[3] ** 2
                        )
                        test_counter["profiles"] += 1
                        test_counter["relative_ratio_failures"] += pair_margin < 0
                        test_counter[
                            "profiles_on_negative_cross_trees"
                        ] += adverse
                        test_counter[
                            "relative_ratio_failures_on_negative_cross_trees"
                        ] += adverse and pair_margin < 0
                        if pair_margin < 0 and adverse:
                            witness_key = (
                                f"{left_key}_with_{right_key}_group_pair_"
                                "failure_on_negative_cross_tree"
                            )
                            first.setdefault(
                                witness_key,
                                {
                                    "n": order,
                                    "graph6": nx.to_graph6_bytes(
                                        tree, header=False
                                    ).decode().strip(),
                                    "left_group": left_key,
                                    "left_profile_2_to_4": left_profile[2:5],
                                    "right_group": right_key,
                                    "right_profile_2_to_4": right_profile[2:5],
                                    "extendable_2_to_4": extendable[2:5],
                                    "full_cross_term": cross,
                                    "pair_relative_ratio_margin": pair_margin,
                                },
                            )

            if adverse:
                coarse_signature = "+".join(sorted(coarse_groups))
                atomic_signature = "+".join(sorted(atomic_groups))
                negative_cross_group_signatures[
                    f"coarse:{coarse_signature}"
                ] += 1
                negative_cross_group_signatures[
                    f"atomic:{atomic_signature}"
                ] += 1

            counters["trees_checked"] += 1

    return {
        "certificate_date": "2026-09-04",
        "kind": "bounded_obstruction_class_structure_probe",
        "claim_status": "exact bounded computation, not a theorem",
        "scope": {
            "tree_orders": [2, max_n],
            "defects": [2, 3, 4],
            "minimum_alpha_bag_count": 4,
            "enumeration": "networkx nonisomorphic trees",
            "arithmetic": "exact integers",
            "matching_choice": "networkx maximum-cardinality matching",
            "certificate_order": (
                "unary lexicographic first, then pair comparisons by increasing "
                "quotient-skeleton path length and lexicographic order"
            ),
        },
        "closure_tests": {
            "individual_class": (
                "Every feasible filling of one erased bag retains the exact same "
                "first certificate."
            ),
            "cumulative_prefix": (
                "Every feasible filling of one erased bag has first-certificate "
                "rank no later than the source assignment."
            ),
        },
        "ratio_test": (
            "For profile c, test c_3^2*e_2*e_4 >= "
            "c_2*c_4*e_3^2."
        ),
        "exact_formulas": {
            "forbidden_vertex_set": (
                "U is the set of vertices contained in no maximum independent set."
            ),
            "unary_group": (
                "u_d = i_{alpha-d}(T) - i_{alpha-d}(T-U)"
            ),
            "pair_group": (
                "p_d = i_{alpha-d}(T-U) - e_d"
            ),
            "ordered_unary_class": (
                "For the i-th forbidden vertex v_i, c_{i,d} is the coefficient "
                "of x^(alpha-d-1) in I(T-N[v_i]-{v_j:j<i};x)."
            ),
        },
        "counters": dict(sorted(counters.items())),
        "by_defect": {
            defect: dict(sorted(values.items()))
            for defect, values in by_defect.items()
        },
        "assigned_certificate_kind": dict(sorted(assigned_certificate_kind.items())),
        "assigned_certificate_kind_on_negative_cross_trees": dict(
            sorted(assigned_certificate_kind_on_negative_cross_trees.items())
        ),
        "assigned_pair_path_length": dict(
            sorted(assigned_pair_path_length.items(), key=lambda item: int(item[0]))
        ),
        "individual_class_tests": dict(sorted(class_tests.items())),
        "cumulative_prefix_tests": dict(sorted(prefix_tests.items())),
        "pairwise_class_tests": dict(sorted(pairwise_class_tests.items())),
        "atomic_group_tests": dict(sorted(atomic_group_tests.items())),
        "coarse_group_tests": dict(sorted(coarse_group_tests.items())),
        "atomic_group_pair_tests": dict(sorted(atomic_group_pair_tests.items())),
        "coarse_group_pair_tests": dict(sorted(coarse_group_pair_tests.items())),
        "negative_cross_group_signatures": dict(
            sorted(negative_cross_group_signatures.items())
        ),
        "formula_checks": dict(sorted(formula_checks.items())),
        "first_witnesses": first,
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
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    counters = result["counters"]
    prefix = result["cumulative_prefix_tests"]
    print(
        "PASS: "
        f"{counters['trees_checked']} trees; "
        f"prefix_fill_failures={counters.get('prefix_fill_failures', 0)}; "
        "negative-cross prefix-ratio failures="
        f"{prefix.get('relative_ratio_failures_on_negative_cross_trees', 0)}"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
