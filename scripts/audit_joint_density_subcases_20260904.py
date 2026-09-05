#!/usr/bin/env python3
"""Audit proof-sized subcases of the conditional joint density bound."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
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


DEFAULT_OUTPUT = Path("results/joint_density_subcases_n18_20260904.json")


def row(
    tree: nx.Graph,
    alpha: int,
    deficiency: int,
    extendable: list[int],
    blocked: list[int],
) -> dict[str, object]:
    correction = (
        blocked[3] ** 2
        - blocked[2] * blocked[4]
        + 2 * extendable[3] * blocked[3]
        - extendable[2] * blocked[4]
        - extendable[4] * blocked[2]
    )
    load = (
        extendable[2] * blocked[4]
        + extendable[4] * blocked[2]
        + blocked[2] * blocked[4]
    )
    reserve = (5 * alpha + 17) * extendable[2] * extendable[4]
    scaled_load = 27 * (alpha - 3) * load
    return {
        "n": tree.number_of_nodes(),
        "alpha": alpha,
        "deficiency": deficiency,
        "graph6": nx.to_graph6_bytes(tree, header=False).decode().strip(),
        "e_0_to_4": extendable,
        "b_0_to_4": blocked,
        "combined_correction": correction,
        "combined_endpoint_load_numerator": scaled_load,
        "pascal_reserve_denominator": reserve,
        "combined_endpoint_load_ratio": float(Fraction(scaled_load, reserve)),
        "b2_density": float(Fraction(blocked[2], extendable[2])),
        "b4_density": float(Fraction(blocked[4], extendable[4])),
    }


def larger_fraction(
    candidate: dict[str, object],
    incumbent: dict[str, object] | None,
    numerator: int,
    denominator: int,
    numerator_key: str,
    denominator_key: str,
) -> bool:
    if incumbent is None:
        return True
    return (
        numerator * int(incumbent[denominator_key])
        > int(incumbent[numerator_key]) * denominator
    )


def run_audit(max_n: int) -> dict[str, object]:
    counters: Counter[str] = Counter()
    adverse_by_deficiency: Counter[int] = Counter()
    first: dict[str, dict[str, object]] = {}
    extrema: dict[str, dict[str, object]] = {}

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
            result = row(tree, alpha, deficiency, extendable, blocked)
            correction = int(result["combined_correction"])

            if alpha >= 5 and blocked[2] == 0:
                counters["alpha_ge_5_b2_zero"] += 1
                if blocked[3] > 0:
                    counters["b2_zero_b3_positive"] += 1
                if blocked[4] > 0:
                    counters["b2_zero_b4_positive"] += 1
                if 5 * blocked[3] > extendable[3]:
                    counters["b2_zero_five_b3_le_e3_failures"] += 1
                    first.setdefault("b2_zero_five_b3_le_e3_failure", result)
                if 5 * blocked[4] > extendable[4]:
                    counters["b2_zero_five_b4_le_e4_failures"] += 1
                    first.setdefault("b2_zero_five_b4_le_e4_failure", result)
                if 5 * blocked[4] == extendable[4]:
                    counters["b2_zero_five_b4_le_e4_equalities"] += 1
                    first.setdefault("b2_zero_five_b4_le_e4_equality", result)
                if result["combined_endpoint_load_numerator"] > result[
                    "pascal_reserve_denominator"
                ]:
                    counters["b2_zero_conditional_endpoint_failures"] += 1
                    first.setdefault("b2_zero_conditional_endpoint_failure", result)
                candidate = {
                    **result,
                    "ratio_numerator": blocked[4],
                    "ratio_denominator": extendable[4],
                }
                current = extrema.get("largest_b4_over_e4_when_b2_zero")
                if larger_fraction(
                    candidate,
                    current,
                    blocked[4],
                    extendable[4],
                    "ratio_numerator",
                    "ratio_denominator",
                ):
                    extrema["largest_b4_over_e4_when_b2_zero"] = candidate

            if correction >= 0:
                continue
            counters["negative_combined_correction"] += 1
            adverse_by_deficiency[deficiency] += 1
            if blocked[2] == 0:
                counters["negative_combined_with_b2_zero"] += 1
            else:
                counters["negative_combined_with_b2_positive"] += 1
            if deficiency <= 1 and alpha >= 5:
                counters["low_deficiency_alpha_ge_5_negative_combined"] += 1
                first.setdefault("low_deficiency_alpha_ge_5_negative_combined", result)
            if result["combined_endpoint_load_numerator"] > result[
                "pascal_reserve_denominator"
            ]:
                counters["conditional_endpoint_failures"] += 1
                first.setdefault("conditional_endpoint_failure", result)

            for name, index in (("b2", 2), ("b4", 4)):
                candidate = {
                    **result,
                    "ratio_numerator": blocked[index],
                    "ratio_denominator": extendable[index],
                }
                key = f"largest_{name}_density_on_negative_combined"
                current = extrema.get(key)
                if larger_fraction(
                    candidate,
                    current,
                    blocked[index],
                    extendable[index],
                    "ratio_numerator",
                    "ratio_denominator",
                ):
                    extrema[key] = candidate

    for key in (
        "alpha_ge_5_b2_zero",
        "b2_zero_b3_positive",
        "b2_zero_b4_positive",
        "b2_zero_conditional_endpoint_failures",
        "b2_zero_five_b3_le_e3_failures",
        "b2_zero_five_b4_le_e4_failures",
        "b2_zero_five_b4_le_e4_equalities",
        "conditional_endpoint_failures",
        "low_deficiency_alpha_ge_5_negative_combined",
        "negative_combined_correction",
        "negative_combined_with_b2_positive",
        "negative_combined_with_b2_zero",
    ):
        counters[key] += 0
    return {
        "certificate_date": "2026-09-04",
        "kind": "joint_density_analytic_subcase_audit",
        "claim_status": "bounded exact evidence and kill tests, not a theorem",
        "scope": {
            "tree_orders": [2, max_n],
            "minimum_alpha": 4,
            "maximum_deficiency": 5,
            "enumeration": "networkx nonisomorphic trees",
            "arithmetic": "exact integers",
        },
        "tested_candidates": {
            "b2_zero_five_b3_bound": "alpha>=5 and b_2=0 implies 5*b_3<=e_3",
            "b2_zero_five_bound": "alpha>=5 and b_2=0 implies 5*b_4<=e_4",
            "low_deficiency_nonadversity": (
                "alpha>=5 and deficiency<=1 implies D>=0"
            ),
            "conditional_joint_bound": (
                "D<0 implies 27(alpha-3)(e_2b_4+e_4b_2+b_2b_4) "
                "<= (5alpha+17)e_2e_4"
            ),
        },
        "counters": dict(sorted(counters.items())),
        "negative_combined_by_deficiency": {
            str(key): value for key, value in sorted(adverse_by_deficiency.items())
        },
        "extrema": extrema,
        "first_witnesses": first,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=18)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.max_n < 5:
        raise SystemExit("--max-n must be at least 5")
    result = run_audit(args.max_n)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    counters = result["counters"]
    print(
        "PASS: "
        f"{counters['trees_evaluated']} trees; "
        f"b2_zero_five_bound_failures="
        f"{counters['b2_zero_five_b4_le_e4_failures']}; "
        f"low_deficiency_adverse="
        f"{counters['low_deficiency_alpha_ge_5_negative_combined']}; "
        f"conditional_failures={counters['conditional_endpoint_failures']}"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
