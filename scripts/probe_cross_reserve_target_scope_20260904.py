#!/usr/bin/env python3
"""Random port-labeled bag-tree probe in the actual depth-3 target scope."""

from __future__ import annotations

import argparse
import json
import random
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

import networkx as nx

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from indpoly import independence_poly, is_unimodal
from scripts.probe_blocked_profile_depth3_20260828 import (
    matching_bags,
    maximum_assignment_code,
)


DEFAULT_OUTPUT = Path("results/cross_reserve_target_scope_probe_20260904.json")
DEFAULT_SEED = 993_20260904


def random_bag_tree(
    alpha: int, deficiency: int, rng: random.Random
) -> nx.Graph:
    """Expand a random bag skeleton with a proposed maximum matching."""
    singleton_indices = set(rng.sample(range(alpha), deficiency))
    skeleton = nx.from_prufer_sequence(
        [rng.randrange(alpha) for _ in range(alpha - 2)]
    )

    bags: list[tuple[int, ...]] = []
    next_vertex = 0
    for index in range(alpha):
        size = 1 if index in singleton_indices else 2
        bags.append(tuple(range(next_vertex, next_vertex + size)))
        next_vertex += size

    tree = nx.Graph()
    tree.add_nodes_from(range(next_vertex))
    for bag in bags:
        if len(bag) == 2:
            tree.add_edge(*bag)
    for left, right in skeleton.edges():
        tree.add_edge(rng.choice(bags[left]), rng.choice(bags[right]))
    if not nx.is_tree(tree):
        raise AssertionError("expanded bag skeleton is not a tree")
    return tree


def projection_defect_profile(
    code: tuple[int, ...], bag_count: int, max_defect: int = 4
) -> list[int]:
    """Compute only the top projection layers needed by the depth-3 test."""
    full_mask = (1 << bag_count) - 1
    profile: list[int] = []
    for defect in range(max_defect + 1):
        total = 0
        for empty in combinations(range(bag_count), defect):
            empty_mask = 0
            for index in empty:
                empty_mask |= 1 << index
            retained_mask = full_mask ^ empty_mask
            total += len({word & retained_mask for word in code})
        profile.append(total)
    return profile


def result_row(
    tree: nx.Graph,
    alpha: int,
    deficiency: int,
    extendable: list[int],
    blocked: list[int],
    full: list[int],
) -> dict[str, object]:
    extendable_turan = extendable[3] ** 2 - extendable[2] * extendable[4]
    blocked_turan = blocked[3] ** 2 - blocked[2] * blocked[4]
    cross = (
        2 * extendable[3] * blocked[3]
        - extendable[2] * blocked[4]
        - extendable[4] * blocked[2]
    )
    negative_endpoints = (
        extendable[2] * blocked[4] + extendable[4] * blocked[2]
    )
    reserve_denominator = (
        (5 * alpha + 17) * extendable[2] * extendable[4]
    )
    endpoint_load_numerator = 27 * (alpha - 3) * negative_endpoints
    cross_load_numerator = -27 * (alpha - 3) * cross
    return {
        "n": tree.number_of_nodes(),
        "alpha": alpha,
        "deficiency": deficiency,
        "graph6": nx.to_graph6_bytes(tree, header=False).decode().strip(),
        "e_0_to_4": extendable,
        "b_0_to_4": blocked,
        "s_0_to_4": full,
        "extendable_turan": extendable_turan,
        "blocked_turan": blocked_turan,
        "cross_term": cross,
        "cross_reserve_margin": extendable_turan + cross,
        "full_margin": extendable_turan + blocked_turan + cross,
        "negative_endpoint_terms": negative_endpoints,
        "pascal_reserve_denominator": reserve_denominator,
        "endpoint_load_numerator": endpoint_load_numerator,
        "cross_load_numerator": cross_load_numerator,
        "pascal_endpoint_margin_numerator": (
            reserve_denominator - endpoint_load_numerator
        ),
        "pascal_cross_margin_numerator": (
            reserve_denominator - cross_load_numerator
        ),
    }


def ratio_is_larger(
    candidate: dict[str, object],
    current: dict[str, object] | None,
    numerator_key: str,
) -> bool:
    if current is None:
        return True
    return (
        candidate[numerator_key] * current["pascal_reserve_denominator"]
        > current[numerator_key] * candidate["pascal_reserve_denominator"]
    )


def run_probe(
    accepted_per_cell: int, max_trials_per_cell: int, seed: int
) -> dict[str, object]:
    rng = random.Random(seed)
    counters: Counter[str] = Counter()
    by_cell: dict[str, Counter[str]] = {}
    first: dict[str, object] = {}
    worst_endpoint_load: dict[str, object] | None = None
    worst_cross_load: dict[str, object] | None = None

    cells = [
        (alpha, deficiency)
        for alpha in range(17, 20)
        for deficiency in range(6)
        if 33 <= 2 * alpha - deficiency <= 38
    ]
    for alpha, proposed_deficiency in cells:
        cell_key = f"alpha_{alpha}_deficiency_{proposed_deficiency}"
        cell = Counter()
        by_cell[cell_key] = cell
        while (
            cell["accepted"] < accepted_per_cell
            and cell["trials"] < max_trials_per_cell
        ):
            cell["trials"] += 1
            counters["trials"] += 1
            tree = random_bag_tree(alpha, proposed_deficiency, rng)
            bags = matching_bags(tree)
            actual_alpha = len(bags)
            actual_deficiency = 2 * actual_alpha - tree.number_of_nodes()
            if actual_alpha != alpha or actual_deficiency != proposed_deficiency:
                cell["rejected_nonmaximum_proposed_matching"] += 1
                counters["rejected_nonmaximum_proposed_matching"] += 1
                continue

            cell["accepted"] += 1
            counters["accepted"] += 1
            code = maximum_assignment_code(tree, bags)
            extendable = projection_defect_profile(code, alpha)
            polynomial = independence_poly(
                tree.number_of_nodes(),
                [list(tree.neighbors(vertex)) for vertex in tree.nodes()],
            )
            if len(polynomial) - 1 != alpha:
                raise AssertionError("matching bags disagree with independence number")
            full = [polynomial[alpha - defect] for defect in range(5)]
            blocked = [full[index] - extendable[index] for index in range(5)]
            if any(value < 0 for value in blocked):
                raise AssertionError("negative blocked coefficient")
            row = result_row(
                tree, alpha, actual_deficiency, extendable, blocked, full
            )

            if not is_unimodal(polynomial):
                counters["unimodality_failures"] += 1
                first.setdefault("unimodality_failure", row)
            if row["blocked_turan"] < 0:
                counters["blocked_lc_failures"] += 1
                first.setdefault("blocked_lc_failure", row)
            if row["cross_reserve_margin"] < 0:
                counters["cross_reserve_failures"] += 1
                first.setdefault("cross_reserve_failure", row)
            if row["full_margin"] < 0:
                counters["full_depth3_failures"] += 1
                first.setdefault("full_depth3_failure", row)

            if row["cross_term"] < 0:
                counters["negative_cross"] += 1
                cell["negative_cross"] += 1
                if row["pascal_endpoint_margin_numerator"] < 0:
                    counters["pascal_endpoint_bound_failures"] += 1
                    cell["pascal_endpoint_bound_failures"] += 1
                    first.setdefault("pascal_endpoint_bound_failure", row)
                if row["pascal_cross_margin_numerator"] < 0:
                    counters["pascal_cross_bound_failures"] += 1
                    cell["pascal_cross_bound_failures"] += 1
                    first.setdefault("pascal_cross_bound_failure", row)
                if ratio_is_larger(
                    row,
                    worst_endpoint_load,
                    "endpoint_load_numerator",
                ):
                    worst_endpoint_load = row
                if ratio_is_larger(
                    row,
                    worst_cross_load,
                    "cross_load_numerator",
                ):
                    worst_cross_load = row

    for key in (
        "blocked_lc_failures",
        "cross_reserve_failures",
        "full_depth3_failures",
        "negative_cross",
        "pascal_cross_bound_failures",
        "pascal_endpoint_bound_failures",
        "unimodality_failures",
    ):
        counters[key] += 0
    return {
        "certificate_date": "2026-09-04",
        "kind": "random_target_scope_cross_reserve_probe",
        "claim_status": "bounded random evidence, not a theorem",
        "scope": {
            "alpha_values": [17, 18, 19],
            "maximum_deficiency": 5,
            "tree_orders": [33, 38],
            "cells": [
                {"alpha": alpha, "deficiency": deficiency}
                for alpha, deficiency in cells
            ],
            "sampling": (
                "Random labeled quotient trees, random singleton locations, "
                "and random endpoint ports, conditioned on the proposed "
                "matching being maximum."
            ),
            "accepted_per_cell": accepted_per_cell,
            "max_trials_per_cell": max_trials_per_cell,
            "seed": seed,
            "arithmetic": "exact integers",
        },
        "candidate": (
            "When C<0, (5*alpha+17)*e_2*e_4 >= "
            "27*(alpha-3)*(e_2*b_4+e_4*b_2)."
        ),
        "counters": dict(sorted(counters.items())),
        "by_cell": {
            key: dict(sorted(value.items())) for key, value in by_cell.items()
        },
        "worst_endpoint_load": worst_endpoint_load,
        "worst_cross_load": worst_cross_load,
        "first_witnesses": first,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accepted-per-cell", type=int, default=10)
    parser.add_argument("--max-trials-per-cell", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.accepted_per_cell < 1 or args.max_trials_per_cell < 1:
        raise SystemExit("sample and trial counts must be positive")
    result = run_probe(
        args.accepted_per_cell, args.max_trials_per_cell, args.seed
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    counters = result["counters"]
    print(
        "PASS: "
        f"{counters['accepted']} accepted target-scope trees; "
        f"negative_cross={counters['negative_cross']}; "
        "endpoint_bound_failures="
        f"{counters['pascal_endpoint_bound_failures']}; "
        f"unimodality_failures={counters['unimodality_failures']}"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
