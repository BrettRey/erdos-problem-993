#!/usr/bin/env python3
"""Targeted matched-pair lifts of the first exact negative-cross witness."""

from __future__ import annotations

import argparse
import json
import random
import sys
from collections import Counter
from fractions import Fraction
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
from scripts.probe_cross_reserve_target_scope_20260904 import (
    projection_defect_profile,
    result_row,
)


BASE_GRAPH6 = "PpG`A?@O??g??@O???I??O??"
DEFAULT_OUTPUT = Path("results/cross_reserve_witness_lifts_20260904.json")
DEFAULT_SEED = 993_20260904
LIFT_PAIRS = 8


def canonical_tree_code(tree: nx.Graph) -> str:
    """Return the exact AHU isomorphism code of an unrooted tree."""

    def rooted_code(vertex: int, parent: int | None) -> str:
        children = sorted(
            rooted_code(child, vertex)
            for child in tree.neighbors(vertex)
            if child != parent
        )
        return "(" + "".join(children) + ")"

    centers = nx.center(tree)
    if len(centers) == 1:
        return rooted_code(centers[0], None)
    if len(centers) != 2:
        raise AssertionError("a tree must have one or two centers")
    left, right = centers
    halves = sorted((rooted_code(left, right), rooted_code(right, left)))
    return "(" + "".join(halves) + ")"


def lift_from_plan(
    base: nx.Graph,
    base_bags: list[tuple[int, ...]],
    plan: tuple[int, ...],
) -> tuple[nx.Graph, list[tuple[int, ...]]]:
    """Append one pendant matched pair for every attachment in ``plan``."""
    tree = base.copy()
    bags = list(base_bags)
    for attachment in plan:
        if attachment not in tree:
            raise ValueError(f"attachment {attachment} does not yet exist")
        inner = tree.number_of_nodes()
        leaf = inner + 1
        tree.add_edge(inner, leaf)
        tree.add_edge(attachment, inner)
        bags.append((inner, leaf))
    if not nx.is_tree(tree):
        raise AssertionError("matched-pair lift is not a tree")
    return tree, bags


def evaluate_lift(
    tree: nx.Graph,
    bags: list[tuple[int, ...]],
) -> dict[str, object]:
    """Evaluate the exact depth-3 decomposition for a prescribed matching."""
    alpha = len(bags)
    prescribed_matching_size = sum(len(bag) == 2 for bag in bags)
    computed_matching_size = len(
        nx.max_weight_matching(tree, maxcardinality=True)
    )
    if computed_matching_size != prescribed_matching_size:
        raise AssertionError("prescribed matching is not maximum")

    code = maximum_assignment_code(tree, bags)
    extendable = projection_defect_profile(code, alpha)
    polynomial = independence_poly(
        tree.number_of_nodes(),
        [list(tree.neighbors(vertex)) for vertex in tree.nodes()],
    )
    if len(polynomial) - 1 != alpha:
        raise AssertionError("prescribed matching gives the wrong alpha")
    full = [polynomial[alpha - defect] for defect in range(5)]
    blocked = [full[index] - extendable[index] for index in range(5)]
    if any(value < 0 for value in blocked):
        raise AssertionError("negative blocked coefficient")

    row = result_row(
        tree,
        alpha,
        2 * alpha - tree.number_of_nodes(),
        extendable,
        blocked,
        full,
    )
    combined_endpoint_load = (
        row["negative_endpoint_terms"] + blocked[2] * blocked[4]
    )
    row["combined_endpoint_terms"] = combined_endpoint_load
    row["combined_endpoint_load_numerator"] = (
        27 * (alpha - 3) * combined_endpoint_load
    )
    row["pascal_combined_endpoint_margin_numerator"] = (
        row["pascal_reserve_denominator"]
        - row["combined_endpoint_load_numerator"]
    )
    row["combined_correction"] = row["cross_term"] + row["blocked_turan"]
    row["unimodal"] = is_unimodal(polynomial)
    row["matching_size"] = computed_matching_size
    return row


def deterministic_plans(base_order: int) -> list[tuple[str, tuple[int, ...]]]:
    """Cover concentrated fans and the two port choices along a chain."""
    plans: list[tuple[str, tuple[int, ...]]] = []
    for anchor in range(base_order):
        plans.append(("fan", (anchor,) * LIFT_PAIRS))

        leaf_chain: list[int] = []
        root_chain: list[int] = []
        leaf_attachment = anchor
        root_attachment = anchor
        for step in range(LIFT_PAIRS):
            inner = base_order + 2 * step
            leaf = inner + 1
            leaf_chain.append(leaf_attachment)
            root_chain.append(root_attachment)
            leaf_attachment = leaf
            root_attachment = inner
        plans.append(("leaf_chain", tuple(leaf_chain)))
        plans.append(("root_chain", tuple(root_chain)))
    return plans


def random_plan(base_order: int, rng: random.Random) -> tuple[str, tuple[int, ...]]:
    """Sample a lift plan, biased toward concentrated and near-chain shapes."""
    strategy = rng.choice(
        ("fan_noise", "fan_noise", "two_anchor", "base_uniform", "chain_noise", "all_uniform")
    )
    anchor = rng.randrange(base_order)
    second_anchor = rng.randrange(base_order)
    attachments: list[int] = []
    previous_inner: int | None = None
    previous_leaf: int | None = None
    for step in range(LIFT_PAIRS):
        existing_order = base_order + 2 * step
        if strategy == "fan_noise":
            attachment = (
                anchor if rng.random() < 0.72 else rng.randrange(existing_order)
            )
        elif strategy == "two_anchor":
            attachment = rng.choice((anchor, second_anchor))
        elif strategy == "base_uniform":
            attachment = rng.randrange(base_order)
        elif strategy == "chain_noise" and previous_inner is not None:
            if rng.random() < 0.72:
                attachment = rng.choice((previous_inner, previous_leaf))
            else:
                attachment = rng.randrange(existing_order)
        else:
            attachment = rng.randrange(existing_order)
        attachments.append(attachment)
        previous_inner = existing_order
        previous_leaf = existing_order + 1
    return strategy, tuple(attachments)


def mutate_plan(
    plan: tuple[int, ...], base_order: int, rng: random.Random
) -> tuple[int, ...]:
    """Change one or two valid attachment choices in a lift plan."""
    child = list(plan)
    for step in rng.sample(range(LIFT_PAIRS), rng.choice((1, 1, 1, 2))):
        child[step] = rng.randrange(base_order + 2 * step)
    return tuple(child)


def ratio(row: dict[str, object], numerator_key: str) -> Fraction:
    return Fraction(
        int(row[numerator_key]), int(row["pascal_reserve_denominator"])
    )


def compact_candidate(candidate: dict[str, object]) -> dict[str, object]:
    row = candidate["row"]
    return {
        "strategy": candidate["strategy"],
        "attachment_plan": list(candidate["plan"]),
        "endpoint_load_ratio": float(ratio(row, "endpoint_load_numerator")),
        "combined_endpoint_load_ratio": float(
            ratio(row, "combined_endpoint_load_numerator")
        ),
        "cross_load_ratio": float(ratio(row, "cross_load_numerator")),
        **row,
    }


def run_probe(target: int, exploration: int, seed: int) -> dict[str, object]:
    if exploration < 1 or exploration > target:
        raise ValueError("exploration must lie between 1 and target")

    rng = random.Random(seed)
    base = nx.from_graph6_bytes(BASE_GRAPH6.encode("ascii"))
    base_bags = matching_bags(base)
    base_row = evaluate_lift(base, base_bags)
    if (
        base_row["n"] != 17
        or base_row["alpha"] != 11
        or base_row["deficiency"] != 5
        or base_row["cross_term"] != -5409
    ):
        raise AssertionError("base witness did not replay exactly")

    counters: Counter[str] = Counter()
    strategy_counts: Counter[str] = Counter()
    negative_strategy_counts: Counter[str] = Counter()
    first: dict[str, object] = {}
    evaluated: list[dict[str, object]] = []
    seen_codes: set[str] = set()
    queued = deterministic_plans(base.number_of_nodes())

    attempts = 0
    max_attempts = 100 * target
    while len(evaluated) < target and attempts < max_attempts:
        attempts += 1
        if queued:
            strategy, plan = queued.pop(0)
        elif len(evaluated) < exploration:
            strategy, plan = random_plan(base.number_of_nodes(), rng)
        else:
            negative = [
                candidate
                for candidate in evaluated
                if candidate["row"]["cross_term"] < 0
            ]
            endpoint_ranked = sorted(
                evaluated,
                key=lambda candidate: ratio(
                    candidate["row"], "endpoint_load_numerator"
                ),
                reverse=True,
            )[:12]
            cross_ranked = sorted(
                evaluated,
                key=lambda candidate: ratio(
                    candidate["row"], "cross_load_numerator"
                ),
                reverse=True,
            )[:12]
            negative_ranked = sorted(
                negative,
                key=lambda candidate: ratio(
                    candidate["row"], "endpoint_load_numerator"
                ),
                reverse=True,
            )[:12]
            parent = rng.choice(negative_ranked + endpoint_ranked + cross_ranked)
            plan = mutate_plan(parent["plan"], base.number_of_nodes(), rng)
            strategy = "mutated_" + parent["strategy"]

        tree, bags = lift_from_plan(base, base_bags, plan)
        code = canonical_tree_code(tree)
        if code in seen_codes:
            counters["isomorphic_duplicates"] += 1
            continue
        seen_codes.add(code)

        row = evaluate_lift(tree, bags)
        candidate = {"strategy": strategy, "plan": plan, "row": row}
        evaluated.append(candidate)
        strategy_counts[strategy] += 1
        counters["evaluated_nonisomorphic_lifts"] += 1

        if row["cross_term"] < 0:
            counters["negative_cross"] += 1
            negative_strategy_counts[strategy] += 1
            if row["pascal_endpoint_margin_numerator"] < 0:
                counters["conditional_endpoint_bound_failures"] += 1
                first.setdefault(
                    "conditional_endpoint_bound_failure",
                    compact_candidate(candidate),
                )
        if row["combined_correction"] < 0:
            counters["negative_combined_correction"] += 1
            if row["pascal_combined_endpoint_margin_numerator"] < 0:
                counters[
                    "strong_combined_bound_failures_on_negative_correction"
                ] += 1
                first.setdefault(
                    "strong_combined_bound_failure_on_negative_correction",
                    compact_candidate(candidate),
                )
        if row["pascal_combined_endpoint_margin_numerator"] < 0:
            counters["strong_combined_bound_failures"] += 1
            first.setdefault(
                "strong_combined_bound_failure", compact_candidate(candidate)
            )
        if row["blocked_turan"] < 0:
            counters["blocked_lc_failures"] += 1
            first.setdefault("blocked_lc_failure", compact_candidate(candidate))
        if row["cross_reserve_margin"] < 0:
            counters["cross_reserve_failures"] += 1
            first.setdefault("cross_reserve_failure", compact_candidate(candidate))
        if row["full_margin"] < 0:
            counters["full_depth3_failures"] += 1
            first.setdefault("full_depth3_failure", compact_candidate(candidate))
        if not row["unimodal"]:
            counters["unimodality_failures"] += 1
            first.setdefault("unimodality_failure", compact_candidate(candidate))

    if len(evaluated) != target:
        raise RuntimeError(
            f"generated only {len(evaluated)} nonisomorphic lifts in {attempts} attempts"
        )

    for key in (
        "blocked_lc_failures",
        "conditional_endpoint_bound_failures",
        "cross_reserve_failures",
        "full_depth3_failures",
        "negative_cross",
        "negative_combined_correction",
        "strong_combined_bound_failures",
        "strong_combined_bound_failures_on_negative_correction",
        "unimodality_failures",
    ):
        counters[key] += 0

    negative = [
        candidate
        for candidate in evaluated
        if candidate["row"]["cross_term"] < 0
    ]
    top_negative_endpoint = sorted(
        negative,
        key=lambda candidate: ratio(
            candidate["row"], "endpoint_load_numerator"
        ),
        reverse=True,
    )[:10]
    top_negative_cross = sorted(
        negative,
        key=lambda candidate: ratio(candidate["row"], "cross_load_numerator"),
        reverse=True,
    )[:10]
    most_negative_cross = sorted(
        negative, key=lambda candidate: candidate["row"]["cross_term"]
    )[:10]
    negative_combined = [
        candidate
        for candidate in evaluated
        if candidate["row"]["combined_correction"] < 0
    ]
    top_negative_combined_load = sorted(
        negative_combined,
        key=lambda candidate: ratio(
            candidate["row"], "combined_endpoint_load_numerator"
        ),
        reverse=True,
    )[:10]

    return {
        "certificate_date": "2026-09-04",
        "kind": "targeted_negative_cross_witness_lift_probe",
        "claim_status": "bounded targeted evidence, not a theorem",
        "base_witness": base_row,
        "construction": {
            "description": (
                "Append eight pendant P2 pairs to the n=17 negative-cross "
                "witness while preserving its chosen maximum matching."
            ),
            "matching_fact": (
                "Appending a pendant P2 increases both order by two and "
                "maximum matching size by one, so eight lifts take "
                "(n,alpha,deficiency)=(17,11,5) to (33,19,5)."
            ),
            "target_cell": {"n": 33, "alpha": 19, "deficiency": 5},
            "lift_pairs": LIFT_PAIRS,
            "prescribed_matching_rechecked_on_every_lift": True,
        },
        "search": {
            "target_nonisomorphic_lifts": target,
            "exploration_lifts_before_adaptive_mutation": exploration,
            "seed": seed,
            "attempts": attempts,
            "arithmetic": "exact integers",
            "deduplication": "exact AHU unrooted-tree isomorphism code",
        },
        "candidates": {
            "conditional_cross_endpoint": (
                "When C<0, (5*alpha+17)*e_2*e_4 >= "
                "27*(alpha-3)*(e_2*b_4+e_4*b_2)."
            ),
            "strong_combined_endpoint": (
                "When C+(b_3^2-b_2*b_4)<0, "
                "(5*alpha+17)*e_2*e_4 >= 27*(alpha-3)*"
                "(e_2*b_4+e_4*b_2+b_2*b_4)."
            ),
        },
        "counters": dict(sorted(counters.items())),
        "strategy_counts": dict(sorted(strategy_counts.items())),
        "negative_cross_strategy_counts": dict(
            sorted(negative_strategy_counts.items())
        ),
        "top_negative_cross_by_endpoint_load": [
            compact_candidate(candidate) for candidate in top_negative_endpoint
        ],
        "top_negative_cross_by_cross_load": [
            compact_candidate(candidate) for candidate in top_negative_cross
        ],
        "most_negative_cross": [
            compact_candidate(candidate) for candidate in most_negative_cross
        ],
        "top_negative_correction_by_combined_endpoint_load": [
            compact_candidate(candidate)
            for candidate in top_negative_combined_load
        ],
        "first_witnesses": first,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=int, default=500)
    parser.add_argument("--exploration", type=int, default=350)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.target < 1:
        raise SystemExit("target must be positive")
    result = run_probe(args.target, args.exploration, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    counters = result["counters"]
    print(
        "PASS: "
        f"{counters['evaluated_nonisomorphic_lifts']} nonisomorphic lifts; "
        f"negative_cross={counters['negative_cross']}; "
        "endpoint_bound_failures="
        f"{counters['conditional_endpoint_bound_failures']}; "
        "combined_bound_failures="
        f"{counters['strong_combined_bound_failures_on_negative_correction']}; "
        f"unimodality_failures={counters['unimodality_failures']}"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
