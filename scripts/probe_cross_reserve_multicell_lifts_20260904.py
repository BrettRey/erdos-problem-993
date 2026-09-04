#!/usr/bin/env python3
"""Target adverse-correction lifts across every actual depth-3 target cell."""

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

from scripts.probe_blocked_profile_depth3_20260828 import matching_bags
from scripts.probe_cross_reserve_witness_lifts_20260904 import (
    canonical_tree_code,
    evaluate_lift,
    lift_from_plan,
)


DEFAULT_INPUT = Path("results/negative_correction_seeds_n18_20260904.json")
DEFAULT_OUTPUT = Path("results/cross_reserve_multicell_lifts_20260904.json")
DEFAULT_RANDOM_SEED = 993_20260904
TARGET_CELLS = tuple(
    (alpha, deficiency)
    for alpha in range(17, 20)
    for deficiency in range(6)
    if 33 <= 2 * alpha - deficiency <= 38
)


def fraction(row: dict[str, object], numerator_key: str) -> Fraction:
    return Fraction(
        int(row[numerator_key]), int(row["pascal_reserve_denominator"])
    )


def seed_from_row(row: dict[str, object]) -> dict[str, object]:
    tree = nx.from_graph6_bytes(str(row["graph6"]).encode("ascii"))
    bags = matching_bags(tree)
    replay = evaluate_lift(tree, bags)
    if replay["combined_correction"] >= 0:
        raise AssertionError("mined adverse seed did not replay")
    return {
        "tree": tree,
        "bags": bags,
        "row": replay,
        "origin": "enumerated_adverse_seed",
        "origin_graph6": row["graph6"],
    }


def deficiency_one_seeds(source: dict[str, object]) -> list[dict[str, object]]:
    """Make valid deficiency-one seeds from the retained perfect-match seeds."""
    rows: list[dict[str, object]] = []
    seen_source_codes: set[str] = set()
    seen_lift_codes: set[str] = set()
    for key, candidates in source["top_seeds_by_target_cell"].items():
        if not key.endswith("deficiency_0"):
            continue
        for source_row in candidates:
            source_graph6 = str(source_row["graph6"])
            if source_graph6 in seen_source_codes:
                continue
            seen_source_codes.add(source_graph6)
            base = nx.from_graph6_bytes(source_graph6.encode("ascii"))
            base_bags = matching_bags(base)
            for anchor in range(base.number_of_nodes()):
                tree = base.copy()
                singleton = tree.number_of_nodes()
                tree.add_edge(anchor, singleton)
                bags = list(base_bags) + [(singleton,)]
                try:
                    row = evaluate_lift(tree, bags)
                except AssertionError:
                    continue
                if row["deficiency"] != 1:
                    raise AssertionError("singleton lift has wrong deficiency")
                tree_code = canonical_tree_code(tree)
                if tree_code in seen_lift_codes:
                    continue
                seen_lift_codes.add(tree_code)
                rows.append(
                    {
                        "tree": tree,
                        "bags": bags,
                        "row": row,
                        "origin": "singleton_leaf_on_deficiency_zero_seed",
                        "origin_graph6": source_graph6,
                        "singleton_anchor": anchor,
                    }
                )

    # No adverse deficiency-one seed occurs through n=18.  Interleave the
    # closest-to-adverse seeds with the seeds having the heaviest endpoint load.
    closest = sorted(
        rows,
        key=lambda seed: (
            int(seed["row"]["combined_correction"]),
            -fraction(seed["row"], "combined_endpoint_load_numerator"),
        ),
    )
    heaviest = sorted(
        rows,
        key=lambda seed: fraction(
            seed["row"], "combined_endpoint_load_numerator"
        ),
        reverse=True,
    )
    selected: list[dict[str, object]] = []
    selected_codes: set[str] = set()
    for seed in closest[:6] + heaviest[:6]:
        code = canonical_tree_code(seed["tree"])
        if code not in selected_codes:
            selected.append(seed)
            selected_codes.add(code)
    return selected


def seeds_for_cell(
    source: dict[str, object],
    alpha: int,
    deficiency: int,
    delta_one: list[dict[str, object]],
) -> list[dict[str, object]]:
    if deficiency == 1:
        seeds = delta_one
    else:
        key = f"alpha_{alpha}_deficiency_{deficiency}"
        seeds = [seed_from_row(row) for row in source["top_seeds_by_target_cell"][key]]
    usable = [seed for seed in seeds if int(seed["row"]["alpha"]) <= alpha]
    if not usable:
        raise AssertionError(f"no usable seeds for {(alpha, deficiency)}")
    return usable


def chain_plan(
    base_order: int, pair_count: int, anchor: int, use_leaf: bool
) -> tuple[int, ...]:
    attachments: list[int] = []
    current = anchor
    for step in range(pair_count):
        attachments.append(current)
        inner = base_order + 2 * step
        current = inner + 1 if use_leaf else inner
    return tuple(attachments)


def deterministic_tasks(
    seeds: list[dict[str, object]], target_alpha: int
) -> list[tuple[int, str, tuple[int, ...]]]:
    tasks: list[tuple[int, str, tuple[int, ...]]] = []
    max_order = max(seed["tree"].number_of_nodes() for seed in seeds)
    for anchor in range(max_order):
        for seed_index, seed in enumerate(seeds):
            base_order = seed["tree"].number_of_nodes()
            if anchor >= base_order:
                continue
            pair_count = target_alpha - int(seed["row"]["alpha"])
            tasks.extend(
                (
                    (seed_index, "fan", (anchor,) * pair_count),
                    (
                        seed_index,
                        "leaf_chain",
                        chain_plan(base_order, pair_count, anchor, True),
                    ),
                    (
                        seed_index,
                        "inner_chain",
                        chain_plan(base_order, pair_count, anchor, False),
                    ),
                )
            )
    return tasks


def random_plan(
    base_order: int, pair_count: int, rng: random.Random
) -> tuple[str, tuple[int, ...]]:
    strategy = rng.choice(
        (
            "fan_noise",
            "fan_noise",
            "two_anchor",
            "base_uniform",
            "chain_noise",
            "all_uniform",
        )
    )
    anchor = rng.randrange(base_order)
    second_anchor = rng.randrange(base_order)
    attachments: list[int] = []
    previous_inner: int | None = None
    previous_leaf: int | None = None
    for step in range(pair_count):
        existing_order = base_order + 2 * step
        if strategy == "fan_noise":
            attachment = anchor if rng.random() < 0.72 else rng.randrange(existing_order)
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
    child = list(plan)
    edit_count = min(len(child), rng.choice((1, 1, 1, 2)))
    for step in rng.sample(range(len(child)), edit_count):
        child[step] = rng.randrange(base_order + 2 * step)
    return tuple(child)


def compact_candidate(candidate: dict[str, object]) -> dict[str, object]:
    seed = candidate["seed"]
    row = candidate["row"]
    return {
        "seed_origin": seed["origin"],
        "seed_graph6": seed["row"]["graph6"],
        "seed_combined_correction": seed["row"]["combined_correction"],
        "seed_adverse": seed["row"]["combined_correction"] < 0,
        "strategy": candidate["strategy"],
        "attachment_plan": list(candidate["plan"]),
        "combined_endpoint_load_ratio": float(
            fraction(row, "combined_endpoint_load_numerator")
        ),
        **row,
    }


def evaluate_cell(
    source: dict[str, object],
    alpha: int,
    deficiency: int,
    target: int,
    exploration: int,
    rng: random.Random,
    delta_one: list[dict[str, object]],
) -> dict[str, object]:
    seeds = seeds_for_cell(source, alpha, deficiency, delta_one)
    queued = deterministic_tasks(seeds, alpha)
    counters: Counter[str] = Counter()
    evaluated: list[dict[str, object]] = []
    seen_codes: set[str] = set()
    attempts = 0
    max_attempts = 200 * target

    while len(evaluated) < target and attempts < max_attempts:
        attempts += 1
        if queued and len(evaluated) < exploration:
            seed_index, strategy, plan = queued.pop(0)
        elif len(evaluated) < exploration:
            seed_index = rng.randrange(len(seeds))
            seed = seeds[seed_index]
            strategy, plan = random_plan(
                seed["tree"].number_of_nodes(),
                alpha - int(seed["row"]["alpha"]),
                rng,
            )
        else:
            adverse = [
                candidate
                for candidate in evaluated
                if candidate["row"]["combined_correction"] < 0
            ]
            high_load = sorted(
                evaluated,
                key=lambda candidate: fraction(
                    candidate["row"], "combined_endpoint_load_numerator"
                ),
                reverse=True,
            )[:12]
            low_correction = sorted(
                evaluated,
                key=lambda candidate: Fraction(
                    int(candidate["row"]["combined_correction"]),
                    max(
                        1,
                        int(candidate["row"]["e_0_to_4"][2])
                        * int(candidate["row"]["e_0_to_4"][4]),
                    ),
                ),
            )[:12]
            adverse_load = sorted(
                adverse,
                key=lambda candidate: fraction(
                    candidate["row"], "combined_endpoint_load_numerator"
                ),
                reverse=True,
            )[:12]
            parent = rng.choice(adverse_load + low_correction + high_load)
            seed_index = int(parent["seed_index"])
            seed = seeds[seed_index]
            plan = mutate_plan(
                parent["plan"], seed["tree"].number_of_nodes(), rng
            )
            strategy = "mutated_" + str(parent["strategy"])

        seed = seeds[seed_index]
        tree, bags = lift_from_plan(seed["tree"], seed["bags"], plan)
        tree_code = canonical_tree_code(tree)
        if tree_code in seen_codes:
            counters["isomorphic_duplicates"] += 1
            continue
        seen_codes.add(tree_code)

        row = evaluate_lift(tree, bags)
        if row["alpha"] != alpha or row["deficiency"] != deficiency:
            raise AssertionError("lift missed its target cell")
        candidate = {
            "seed_index": seed_index,
            "seed": seed,
            "strategy": strategy,
            "plan": plan,
            "row": row,
        }
        evaluated.append(candidate)
        counters["evaluated_nonisomorphic_lifts"] += 1

        if row["combined_correction"] < 0:
            counters["negative_combined_correction"] += 1
            if row["pascal_combined_endpoint_margin_numerator"] < 0:
                counters["joint_conditional_bound_failures"] += 1
        if row["pascal_combined_endpoint_margin_numerator"] < 0:
            counters["unconditional_strong_bound_failures"] += 1
        if row["blocked_turan"] < 0:
            counters["blocked_lc_failures"] += 1
        if row["full_margin"] < 0:
            counters["full_depth3_failures"] += 1
        if not row["unimodal"]:
            counters["unimodality_failures"] += 1

    if len(evaluated) != target:
        raise RuntimeError(
            f"generated only {len(evaluated)} lifts for {(alpha, deficiency)}"
        )
    for key in (
        "blocked_lc_failures",
        "full_depth3_failures",
        "joint_conditional_bound_failures",
        "negative_combined_correction",
        "unconditional_strong_bound_failures",
        "unimodality_failures",
    ):
        counters[key] += 0

    adverse = [
        candidate
        for candidate in evaluated
        if candidate["row"]["combined_correction"] < 0
    ]
    top_adverse = sorted(
        adverse,
        key=lambda candidate: fraction(
            candidate["row"], "combined_endpoint_load_numerator"
        ),
        reverse=True,
    )[:5]
    closest = sorted(
        evaluated,
        key=lambda candidate: Fraction(
            int(candidate["row"]["combined_correction"]),
            max(
                1,
                int(candidate["row"]["e_0_to_4"][2])
                * int(candidate["row"]["e_0_to_4"][4]),
            ),
        ),
    )[:3]
    return {
        "target": {"n": 2 * alpha - deficiency, "alpha": alpha, "deficiency": deficiency},
        "seed_count": len(seeds),
        "adverse_seed_count": sum(
            int(seed["row"]["combined_correction"] < 0) for seed in seeds
        ),
        "attempts": attempts,
        "counters": dict(sorted(counters.items())),
        "top_adverse_by_combined_endpoint_load": [
            compact_candidate(candidate) for candidate in top_adverse
        ],
        "closest_to_adverse": [compact_candidate(candidate) for candidate in closest],
    }


def run_probe(
    source: dict[str, object], target: int, exploration: int, random_seed: int
) -> dict[str, object]:
    if not 1 <= exploration <= target:
        raise ValueError("exploration must lie between 1 and target")
    rng = random.Random(random_seed)
    delta_one = deficiency_one_seeds(source)
    if any(seed["row"]["combined_correction"] < 0 for seed in delta_one):
        raise AssertionError("unexpected adverse deficiency-one seed")

    by_cell: dict[str, object] = {}
    aggregate: Counter[str] = Counter()
    for alpha, deficiency in TARGET_CELLS:
        key = f"alpha_{alpha}_deficiency_{deficiency}"
        cell = evaluate_cell(
            source,
            alpha,
            deficiency,
            target,
            exploration,
            rng,
            delta_one,
        )
        by_cell[key] = cell
        aggregate.update(cell["counters"])
        print(
            f"{key}: {cell['counters']['evaluated_nonisomorphic_lifts']} lifts, "
            f"negative_combined={cell['counters']['negative_combined_correction']}, "
            f"joint_failures={cell['counters']['joint_conditional_bound_failures']}",
            flush=True,
        )

    return {
        "certificate_date": "2026-09-04",
        "kind": "targeted_multicell_negative_correction_lift_probe",
        "claim_status": "bounded targeted evidence, not a theorem",
        "scope": {
            "target_cells": [
                {"alpha": alpha, "deficiency": deficiency}
                for alpha, deficiency in TARGET_CELLS
            ],
            "nonisomorphic_lifts_per_cell": target,
            "exploration_lifts_per_cell": exploration,
            "random_seed": random_seed,
            "arithmetic": "exact integers",
            "deduplication": "exact AHU unrooted-tree isomorphism code within each cell",
        },
        "construction": {
            "description": (
                "Append pendant matched P2 pairs to same-deficiency small seeds; "
                "each pair raises n by two and alpha by one while preserving deficiency."
            ),
            "deficiency_one_caveat": (
                "No negative combined-correction seed was found through n=18. "
                "The deficiency-one cells instead start from the closest and "
                "heaviest valid singleton-leaf transforms of deficiency-zero seeds."
            ),
            "deficiency_one_seed_count": len(delta_one),
            "prescribed_matching_rechecked_on_every_lift": True,
        },
        "candidate": (
            "When D=C+(b_3^2-b_2*b_4)<0, "
            "(5*alpha+17)*e_2*e_4 >= 27*(alpha-3)*"
            "(e_2*b_4+e_4*b_2+b_2*b_4)."
        ),
        "aggregate_counters": dict(sorted(aggregate.items())),
        "by_cell": by_cell,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--target", type=int, default=50)
    parser.add_argument("--exploration", type=int, default=35)
    parser.add_argument("--seed", type=int, default=DEFAULT_RANDOM_SEED)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.target < 1:
        raise SystemExit("target must be positive")
    source = json.loads(args.input.read_text(encoding="utf-8"))
    result = run_probe(source, args.target, args.exploration, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    counters = result["aggregate_counters"]
    adverse_cells = sum(
        cell["counters"]["negative_combined_correction"] > 0
        for cell in result["by_cell"].values()
    )
    print(
        "PASS: "
        f"{counters['evaluated_nonisomorphic_lifts']} lifts; "
        f"adverse_cells={adverse_cells}/{len(TARGET_CELLS)}; "
        f"negative_combined={counters['negative_combined_correction']}; "
        f"joint_failures={counters['joint_conditional_bound_failures']}; "
        f"unimodality_failures={counters['unimodality_failures']}"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
