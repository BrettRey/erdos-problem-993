#!/usr/bin/env python3
"""Mine small-tree adverse-correction seeds for each target deficiency."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction
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


DEFAULT_OUTPUT = Path("results/negative_correction_seeds_n18_20260904.json")
TARGET_CELLS = tuple(
    (alpha, deficiency)
    for alpha in range(17, 20)
    for deficiency in range(6)
    if 33 <= 2 * alpha - deficiency <= 38
)


def seed_row(
    tree: nx.Graph,
    alpha: int,
    deficiency: int,
    extendable: list[int],
    blocked: list[int],
) -> dict[str, object]:
    extendable_turan = extendable[3] ** 2 - extendable[2] * extendable[4]
    blocked_turan = blocked[3] ** 2 - blocked[2] * blocked[4]
    cross = (
        2 * extendable[3] * blocked[3]
        - extendable[2] * blocked[4]
        - extendable[4] * blocked[2]
    )
    correction = blocked_turan + cross
    combined_endpoints = (
        extendable[2] * blocked[4]
        + extendable[4] * blocked[2]
        + blocked[2] * blocked[4]
    )
    denominator = (5 * alpha + 17) * extendable[2] * extendable[4]
    numerator = 27 * (alpha - 3) * combined_endpoints
    return {
        "n": tree.number_of_nodes(),
        "alpha": alpha,
        "deficiency": deficiency,
        "graph6": nx.to_graph6_bytes(tree, header=False).decode().strip(),
        "e_0_to_4": extendable,
        "b_0_to_4": blocked,
        "blocked_turan": blocked_turan,
        "cross_term": cross,
        "combined_correction": correction,
        "combined_endpoint_load_numerator": numerator,
        "pascal_reserve_denominator": denominator,
        "combined_endpoint_load_ratio": float(Fraction(numerator, denominator)),
        "full_margin": extendable_turan + correction,
    }


def run_mining(max_n: int, seeds_per_cell: int) -> dict[str, object]:
    counters: Counter[str] = Counter()
    by_deficiency: Counter[int] = Counter()
    candidates: dict[str, list[dict[str, object]]] = defaultdict(list)

    for order in range(2, max_n + 1):
        for tree in nx.nonisomorphic_trees(order):
            bags = matching_bags(tree)
            alpha = len(bags)
            deficiency = 2 * alpha - order
            if alpha < 4 or deficiency > 5:
                continue
            counters["trees_evaluated"] += 1
            code = maximum_assignment_code(tree, bags)
            projection = projection_profile(code, alpha)
            extendable = [projection[alpha - defect] for defect in range(5)]
            polynomial = independence_poly(
                order,
                [list(tree.neighbors(vertex)) for vertex in range(order)],
            )
            full = [polynomial[alpha - defect] for defect in range(5)]
            blocked = [full[index] - extendable[index] for index in range(5)]
            row = seed_row(tree, alpha, deficiency, extendable, blocked)
            if row["combined_correction"] >= 0:
                continue

            counters["negative_combined_correction"] += 1
            by_deficiency[deficiency] += 1
            for target_alpha, target_deficiency in TARGET_CELLS:
                if deficiency != target_deficiency or alpha > target_alpha:
                    continue
                key = f"alpha_{target_alpha}_deficiency_{deficiency}"
                candidates[key].append(row)

    top_by_cell: dict[str, list[dict[str, object]]] = {}
    for target_alpha, deficiency in TARGET_CELLS:
        key = f"alpha_{target_alpha}_deficiency_{deficiency}"
        ranked = sorted(
            candidates[key],
            key=lambda row: Fraction(
                row["combined_endpoint_load_numerator"],
                row["pascal_reserve_denominator"],
            ),
            reverse=True,
        )
        top_by_cell[key] = ranked[:seeds_per_cell]

    return {
        "certificate_date": "2026-09-04",
        "kind": "negative_combined_correction_seed_mining",
        "claim_status": "seed selection only, not theorem evidence",
        "scope": {
            "tree_orders": [2, max_n],
            "maximum_deficiency": 5,
            "target_cells": [
                {"alpha": alpha, "deficiency": deficiency}
                for alpha, deficiency in TARGET_CELLS
            ],
            "seeds_per_cell": seeds_per_cell,
            "enumeration": "networkx nonisomorphic trees",
            "arithmetic": "exact integers",
        },
        "selection": (
            "For each target cell, retain same-deficiency trees with adverse "
            "combined correction, ranked by combined endpoint load divided "
            "by the guaranteed Pascal reserve."
        ),
        "counters": dict(sorted(counters.items())),
        "negative_combined_correction_by_deficiency": {
            str(key): value for key, value in sorted(by_deficiency.items())
        },
        "top_seeds_by_target_cell": top_by_cell,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=18)
    parser.add_argument("--seeds-per-cell", type=int, default=5)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.max_n < 5 or args.seeds_per_cell < 1:
        raise SystemExit("invalid positive scope")
    result = run_mining(args.max_n, args.seeds_per_cell)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    counts = result["counters"]
    covered = sum(bool(rows) for rows in result["top_seeds_by_target_cell"].values())
    print(
        "PASS: "
        f"{counts['trees_evaluated']} trees; "
        f"negative_combined={counts.get('negative_combined_correction', 0)}; "
        f"covered_target_cells={covered}/{len(TARGET_CELLS)}"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
