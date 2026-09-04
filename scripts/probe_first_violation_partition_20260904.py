#!/usr/bin/env python3
"""Exact first-violation partition for blocked matching-bag assignments."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
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


DEFAULT_OUTPUT = Path("results/first_violation_partition_probe_20260904.json")


def quotient_skeleton(
    tree: nx.Graph, bags: list[tuple[int, ...]]
) -> nx.Graph:
    bag_of = {
        vertex: bag_index
        for bag_index, bag in enumerate(bags)
        for vertex in bag
    }
    skeleton = nx.Graph()
    skeleton.add_nodes_from(range(len(bags)))
    for left, right in tree.edges():
        left_bag = bag_of[left]
        right_bag = bag_of[right]
        if left_bag != right_bag:
            skeleton.add_edge(left_bag, right_bag)
    if not nx.is_tree(skeleton):
        raise AssertionError("matching-bag quotient is not a tree")
    return skeleton


def code_value(word: int, bags: list[tuple[int, ...]], index: int) -> int:
    if len(bags[index]) == 1:
        return 0
    return (word >> index) & 1


def forbidden_certificates(
    code: tuple[int, ...], bags: list[tuple[int, ...]]
) -> tuple[tuple[object, ...], ...]:
    """Return all minimal unary/pair projection obstructions in fixed order."""
    allowed: list[set[int]] = []
    for index, bag in enumerate(bags):
        values = {code_value(word, bags, index) for word in code}
        allowed.append(values)

    certificates: list[tuple[object, ...]] = []
    for index, bag in enumerate(bags):
        for value in range(len(bag)):
            if value not in allowed[index]:
                certificates.append((0, index, value))

    for left, right in combinations(range(len(bags)), 2):
        for left_value in sorted(allowed[left]):
            for right_value in sorted(allowed[right]):
                if not any(
                    code_value(word, bags, left) == left_value
                    and code_value(word, bags, right) == right_value
                    for word in code
                ):
                    certificates.append(
                        (1, left, left_value, right, right_value)
                    )
    return tuple(sorted(certificates))


def violated_certificates(
    assignment: dict[int, int], certificates: tuple[tuple[object, ...], ...]
) -> list[tuple[object, ...]]:
    violations: list[tuple[object, ...]] = []
    for certificate in certificates:
        if certificate[0] == 0:
            _, index, value = certificate
            if assignment.get(index) == value:
                violations.append(certificate)
        else:
            _, left, left_value, right, right_value = certificate
            if (
                assignment.get(left) == left_value
                and assignment.get(right) == right_value
            ):
                violations.append(certificate)
    return violations


def extends_code(
    assignment: dict[int, int],
    code: tuple[int, ...],
    bags: list[tuple[int, ...]],
) -> bool:
    return any(
        all(code_value(word, bags, index) == value for index, value in assignment.items())
        for word in code
    )


def independent_mask(chosen: tuple[int, ...], adjacency_masks: list[int]) -> bool:
    mask = 0
    for vertex in chosen:
        if adjacency_masks[vertex] & mask:
            return False
        mask |= 1 << vertex
    return True


def run_probe(max_n: int) -> dict[str, object]:
    counters: Counter[str] = Counter()
    by_defect = {
        str(defect): Counter(
            candidates=0, feasible=0, extendable=0, blocked=0, partitioned=0
        )
        for defect in range(2, 5)
    }
    assigned_certificate_kind: Counter[str] = Counter()
    assigned_path_length: Counter[str] = Counter()
    first: dict[str, object] = {}

    for order in range(2, max_n + 1):
        for tree in nx.generators.nonisomorphic_trees(order):
            bags = matching_bags(tree)
            alpha = len(bags)
            if alpha < 4:
                continue

            code = maximum_assignment_code(tree, bags)
            projection = projection_profile(code, alpha)
            extendable_profile = [
                projection[alpha - defect] for defect in range(5)
            ]
            polynomial = independence_poly(
                order, [list(tree.neighbors(vertex)) for vertex in range(order)]
            )
            full_profile = [polynomial[alpha - defect] for defect in range(5)]
            blocked_profile = [
                full_profile[defect] - extendable_profile[defect]
                for defect in range(5)
            ]

            skeleton = quotient_skeleton(tree, bags)
            certificates = forbidden_certificates(code, bags)
            adjacency_masks = [
                sum(1 << neighbor for neighbor in tree.neighbors(vertex))
                for vertex in range(order)
            ]
            graph6 = nx.to_graph6_bytes(tree, header=False).decode().strip()

            observed = {
                defect: Counter(feasible=0, extendable=0, blocked=0, partitioned=0)
                for defect in range(2, 5)
            }
            class_counts: Counter[tuple[object, ...]] = Counter()

            for defect in range(2, 5):
                if defect > alpha:
                    continue
                for empty_tuple in combinations(range(alpha), defect):
                    empty = set(empty_tuple)
                    occupied = [index for index in range(alpha) if index not in empty]
                    for chosen_values in product(*(range(len(bags[index])) for index in occupied)):
                        counters["candidate_assignments"] += 1
                        by_defect[str(defect)]["candidates"] += 1
                        chosen_vertices = tuple(
                            bags[index][value]
                            for index, value in zip(occupied, chosen_values)
                        )
                        if not independent_mask(chosen_vertices, adjacency_masks):
                            continue

                        assignment = dict(zip(occupied, chosen_values))
                        observed[defect]["feasible"] += 1
                        by_defect[str(defect)]["feasible"] += 1
                        counters["feasible_assignments"] += 1

                        extendable = extends_code(assignment, code, bags)
                        violations = violated_certificates(assignment, certificates)
                        if extendable:
                            observed[defect]["extendable"] += 1
                            by_defect[str(defect)]["extendable"] += 1
                            counters["extendable_assignments"] += 1
                            if violations:
                                raise AssertionError(
                                    f"extendable assignment violates certificate: {graph6}"
                                )
                            continue

                        observed[defect]["blocked"] += 1
                        by_defect[str(defect)]["blocked"] += 1
                        counters["blocked_assignments"] += 1
                        if not violations:
                            raise AssertionError(
                                f"blocked assignment has no unary/pair violation: {graph6}"
                            )

                        first_violation = violations[0]
                        class_counts[first_violation] += 1
                        observed[defect]["partitioned"] += 1
                        by_defect[str(defect)]["partitioned"] += 1
                        counters["partitioned_assignments"] += 1
                        counters["total_violations_on_blocked_assignments"] += len(
                            violations
                        )
                        if len(violations) > 1:
                            counters["multiply_certified_blocked_assignments"] += 1
                            first.setdefault(
                                "multiply_certified_blocked_assignment",
                                {
                                    "n": order,
                                    "graph6": graph6,
                                    "defect": defect,
                                    "empty_bags": list(empty_tuple),
                                    "assignment": {
                                        str(index): value
                                        for index, value in sorted(assignment.items())
                                    },
                                    "violations": [list(item) for item in violations],
                                },
                            )

                        if first_violation[0] == 0:
                            assigned_certificate_kind["unary"] += 1
                        else:
                            assigned_certificate_kind["pair"] += 1
                            _, left, _, right, _ = first_violation
                            path_length = nx.shortest_path_length(skeleton, left, right)
                            assigned_path_length[str(path_length)] += 1
                            if path_length < 2:
                                raise AssertionError(
                                    "feasible assignment violates an adjacent-bag constraint"
                                )

            for defect in range(2, 5):
                if observed[defect]["feasible"] != full_profile[defect]:
                    raise AssertionError(
                        f"full profile mismatch for {graph6}, defect {defect}"
                    )
                if observed[defect]["extendable"] != extendable_profile[defect]:
                    raise AssertionError(
                        f"extendable profile mismatch for {graph6}, defect {defect}"
                    )
                if observed[defect]["blocked"] != blocked_profile[defect]:
                    raise AssertionError(
                        f"blocked profile mismatch for {graph6}, defect {defect}"
                    )
                if observed[defect]["partitioned"] != blocked_profile[defect]:
                    raise AssertionError(
                        f"partition mismatch for {graph6}, defect {defect}"
                    )
            if sum(class_counts.values()) != sum(blocked_profile[2:5]):
                raise AssertionError(f"class sum mismatch for {graph6}")
            counters["trees_checked"] += 1
            counters["distinct_tree_local_classes"] += len(class_counts)

    return {
        "certificate_date": "2026-09-04",
        "kind": "bounded_first_violation_partition_probe",
        "claim_status": "exact bounded computation, not a theorem",
        "scope": {
            "tree_orders": [2, max_n],
            "defects": [2, 3, 4],
            "minimum_alpha_bag_count": 4,
            "enumeration": "networkx nonisomorphic trees",
            "arithmetic": "exact integers",
            "matching_choice": "networkx maximum-cardinality matching",
            "certificate_order": "unary before pair, then bag/value lexicographic order",
        },
        "partition_rule": (
            "Assign each blocked feasible partial bag assignment to the first "
            "unary or pair projection obstruction of the maximum-assignment code."
        ),
        "counters": dict(sorted(counters.items())),
        "by_defect": {
            defect: dict(sorted(values.items()))
            for defect, values in by_defect.items()
        },
        "assigned_certificate_kind": dict(sorted(assigned_certificate_kind.items())),
        "assigned_pair_path_length": dict(sorted(assigned_path_length.items())),
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
    print(
        "PASS: "
        f"{counters['trees_checked']} trees; "
        f"{counters['feasible_assignments']} feasible assignments; "
        f"{counters['blocked_assignments']} blocked assignments partitioned"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
