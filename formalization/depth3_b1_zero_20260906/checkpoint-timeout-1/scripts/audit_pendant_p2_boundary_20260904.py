#!/usr/bin/env python3
"""Audit root-conditioned inequalities for the pendant-P2 recurrence."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

import networkx as nx


DEFAULT_OUTPUT = Path("results/pendant_p2_boundary_n14_20260904.json")


def independent_masks(tree: nx.Graph) -> list[int]:
    adjacency = [0] * tree.number_of_nodes()
    for left, right in tree.edges():
        adjacency[left] |= 1 << right
        adjacency[right] |= 1 << left

    masks: list[int] = []
    for mask in range(1 << tree.number_of_nodes()):
        remaining = mask
        while remaining:
            bit = remaining & -remaining
            vertex = bit.bit_length() - 1
            remaining -= bit
            if adjacency[vertex] & remaining:
                break
        else:
            masks.append(mask)
    return masks


def contained_in_any(mask: int, containers: list[int]) -> bool:
    return any(mask & ~container == 0 for container in containers)


def profile(
    independent: list[int], maximum: list[int], alpha: int
) -> tuple[list[int], list[int]]:
    extendable: list[int] = []
    blocked: list[int] = []
    for defect in range(5):
        layer = [mask for mask in independent if mask.bit_count() == alpha - defect]
        extendable_count = sum(contained_in_any(mask, maximum) for mask in layer)
        extendable.append(extendable_count)
        blocked.append(len(layer) - extendable_count)
    return extendable, blocked


def rooted_profile(
    independent: list[int], maximum: list[int], alpha: int, root: int
) -> tuple[list[int], list[int]]:
    maximum_avoiding_root = [mask for mask in maximum if not mask & (1 << root)]
    generated: list[int] = []
    residual: list[int] = []
    for defect in range(5):
        layer = [
            mask
            for mask in independent
            if mask.bit_count() == alpha - defect and not mask & (1 << root)
        ]
        generated_count = sum(
            contained_in_any(mask, maximum_avoiding_root) for mask in layer
        )
        generated.append(generated_count)
        residual.append(len(layer) - generated_count)
    return generated, residual


def row(
    tree: nx.Graph,
    root: int,
    alpha: int,
    extendable: list[int],
    blocked: list[int],
    generated: list[int],
    residual: list[int],
) -> dict[str, object]:
    lifted_extendable = [
        (extendable[defect - 1] if defect else 0)
        + extendable[defect]
        + generated[defect]
        for defect in range(5)
    ]
    lifted_blocked = [
        (blocked[defect - 1] if defect else 0)
        + blocked[defect]
        + residual[defect]
        for defect in range(5)
    ]
    return {
        "n": tree.number_of_nodes(),
        "graph6": nx.to_graph6_bytes(tree, header=False).decode().strip(),
        "root": root,
        "alpha": alpha,
        "e_0_to_4": extendable,
        "b_0_to_4": blocked,
        "q_0_to_4": generated,
        "h_0_to_4": residual,
        "lifted_e_0_to_4": lifted_extendable,
        "lifted_b_0_to_4": lifted_blocked,
        "h3_boundary_lhs": 5 * residual[3],
        "h3_boundary_rhs": extendable[2] + generated[3],
        "h4_boundary_lhs": 5 * residual[4],
        "h4_boundary_rhs": generated[4],
    }


def update_extremum(
    extrema: dict[str, dict[str, object]],
    name: str,
    candidate: dict[str, object],
    numerator_key: str,
    denominator_key: str,
) -> None:
    denominator = int(candidate[denominator_key])
    if denominator == 0:
        return
    numerator = int(candidate[numerator_key])
    incumbent = extrema.get(name)
    if incumbent is None or (
        numerator * int(incumbent[denominator_key])
        > int(incumbent[numerator_key]) * denominator
    ):
        extrema[name] = {
            **candidate,
            "ratio": float(Fraction(numerator, denominator)),
        }


def run_audit(max_n: int) -> dict[str, object]:
    counters: Counter[str] = Counter()
    first: dict[str, dict[str, object]] = {}
    extrema: dict[str, dict[str, object]] = {}

    for order in range(2, max_n + 1):
        for tree in nx.nonisomorphic_trees(order):
            counters["trees_evaluated"] += 1
            independent = independent_masks(tree)
            alpha = max(mask.bit_count() for mask in independent)
            if alpha < 4:
                continue
            maximum = [mask for mask in independent if mask.bit_count() == alpha]
            extendable, blocked = profile(independent, maximum, alpha)
            if blocked[1] or blocked[2]:
                continue
            counters["base_trees_with_b1_b2_zero"] += 1

            for root in tree.nodes():
                counters["roots_tested_on_candidate_bases"] += 1
                generated, residual = rooted_profile(
                    independent, maximum, alpha, root
                )
                if residual[2]:
                    continue
                counters["eligible_rooted_extensions"] += 1
                candidate = row(
                    tree,
                    root,
                    alpha,
                    extendable,
                    blocked,
                    generated,
                    residual,
                )
                if candidate["lifted_b_0_to_4"][2] != 0:
                    raise AssertionError("eligible extension does not have b_2^+=0")

                if candidate["h3_boundary_lhs"] > candidate["h3_boundary_rhs"]:
                    counters["h3_boundary_failures"] += 1
                    first.setdefault("h3_boundary_failure", candidate)
                if candidate["h4_boundary_lhs"] > candidate["h4_boundary_rhs"]:
                    counters["h4_boundary_failures"] += 1
                    first.setdefault("h4_boundary_failure", candidate)
                if (
                    5 * candidate["lifted_b_0_to_4"][3]
                    > candidate["lifted_e_0_to_4"][3]
                ):
                    counters["lifted_five_b3_le_e3_failures"] += 1
                    first.setdefault("lifted_five_b3_le_e3_failure", candidate)
                if (
                    5 * candidate["lifted_b_0_to_4"][4]
                    > candidate["lifted_e_0_to_4"][4]
                ):
                    counters["lifted_five_b4_le_e4_failures"] += 1
                    first.setdefault("lifted_five_b4_le_e4_failure", candidate)

                update_extremum(
                    extrema,
                    "largest_h3_boundary_ratio",
                    candidate,
                    "h3_boundary_lhs",
                    "h3_boundary_rhs",
                )
                update_extremum(
                    extrema,
                    "largest_h4_boundary_ratio",
                    candidate,
                    "h4_boundary_lhs",
                    "h4_boundary_rhs",
                )

    for key in (
        "base_trees_with_b1_b2_zero",
        "eligible_rooted_extensions",
        "h3_boundary_failures",
        "h4_boundary_failures",
        "lifted_five_b3_le_e3_failures",
        "lifted_five_b4_le_e4_failures",
        "roots_tested_on_candidate_bases",
        "trees_evaluated",
    ):
        counters[key] += 0
    return {
        "certificate_date": "2026-09-04",
        "kind": "pendant_p2_root_conditioned_boundary_audit",
        "claim_status": "bounded exact evidence and kill tests, not a theorem",
        "scope": {
            "base_tree_orders": [2, max_n],
            "minimum_base_alpha": 4,
            "enumeration": "networkx nonisomorphic trees and every root",
            "arithmetic": "exact integers",
        },
        "recurrence": {
            "extendable_polynomial": "E_plus=(1+x)E+xQ",
            "blocked_polynomial": "B_plus=(1+x)B+xH",
            "rooted_remainder": "H=I(T-root)-Q",
        },
        "eligibility": "b_1=b_2=h_2=0, equivalently b_2^+=0",
        "tested_candidates": {
            "h3_boundary": "5*h_3<=e_2+q_3",
            "h4_boundary": "5*h_4<=q_4",
            "lifted_b3": "5*b_3^+<=e_3^+",
            "lifted_b4": "5*b_4^+<=e_4^+",
        },
        "counters": dict(sorted(counters.items())),
        "extrema": extrema,
        "first_witnesses": first,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=14)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
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
        f"{counters['eligible_rooted_extensions']} eligible rooted extensions; "
        f"h3_failures={counters['h3_boundary_failures']}; "
        f"h4_failures={counters['h4_boundary_failures']}"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
