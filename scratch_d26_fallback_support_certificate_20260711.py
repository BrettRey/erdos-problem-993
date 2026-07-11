#!/usr/bin/env python3
"""Audit the fallback-corrected D26 support inequality.

The bare support exchange is false.  Its one-fallback correction is

    sum_{j in R,q not in R} w(R-j+q) + w(R) >= r w(R),       (FSE)

or equivalently neighbour_weight >= (r-1)w(R).  A stronger local sufficient
condition is the per-deletion cavity inequality

    r * sum_{q not in R} w(R-j+q)
        >= 2(r-1) w(R-j)                                    (FCE)

for every j in R.  Indeed w(R)<=2w(R-j), and summing FCE over j proves FSE.

This script has two separately runnable exact audits:

* ``--mode exhaustive`` checks every prefix support through a requested tree
  order (14 by default);
* ``--mode stress`` uses conditioned block-tree DP on random trees of order
  15--30 and all-double path/comb/star/random block constructions, and also
  replays the exact 20-vertex obstruction to the bare inequality.

Passing either audit is evidence, not a proof of FCE or FSE.
"""

from __future__ import annotations

import argparse
import json
import random
from fractions import Fraction
from pathlib import Path

import networkx as nx

from scratch_block_support_exchange_search_20260711 import (
    matching_blocks,
    support_weights,
)
from scratch_d17_checkpointed_search_20260711 import append_record, utc_now
from scratch_d26_marked_pair_obstruction_certificate_20260711 import prefix_last
from scratch_d26_support_exchange_obstruction_certificate_20260711 import witness
from scratch_d26_support_exchange_stress_20260711 import (
    comb_tree,
    doubled_block_tree,
    exchange_row,
    make_csp,
    maximum_matching_blocks,
    random_tree,
    sampled_support,
)


def serialise(row: dict[str, object]) -> dict[str, object]:
    return {
        key: str(value) if isinstance(value, Fraction) else value
        for key, value in row.items()
    }


def exhaustive(max_order: int) -> dict[str, object]:
    trees = 0
    supports = 0
    deletion_checks = 0
    sharp_fse: tuple[Fraction, dict[str, object]] | None = None
    sharp_fce: tuple[Fraction, dict[str, object]] | None = None

    for n in range(2, max_order + 1):
        for tree_index, tree in enumerate(nx.nonisomorphic_trees(n)):
            trees += 1
            blocks, block_of = matching_blocks(tree)
            weights = support_weights(tree, block_of)
            alpha = len(blocks)
            for support, weight in weights.items():
                rank = support.bit_count()
                if rank < 1 or rank > prefix_last(alpha):
                    continue
                supports += 1
                empty = [q for q in range(alpha) if not (support >> q) & 1]
                neighbour = 0
                for deleted in range(alpha):
                    if not (support >> deleted) & 1:
                        continue
                    deletion_checks += 1
                    cavity_support = support ^ (1 << deleted)
                    cavity_weight = weights[cavity_support]
                    extension_weight = sum(
                        weights[cavity_support | (1 << q)] for q in empty
                    )
                    neighbour += extension_weight
                    cavity_gap = (
                        rank * extension_weight
                        - 2 * (rank - 1) * cavity_weight
                    )
                    assert cavity_gap >= 0, (
                        n,
                        tree_index,
                        sorted(tree.edges()),
                        blocks,
                        support,
                        deleted,
                        extension_weight,
                        cavity_weight,
                    )
                    cavity_slack = Fraction(cavity_gap, rank * cavity_weight)
                    cavity_row = {
                        "n": n,
                        "tree_index": tree_index,
                        "edges": sorted([sorted(edge) for edge in tree.edges()]),
                        "blocks": [list(block) for block in blocks],
                        "alpha": alpha,
                        "rank": rank,
                        "support": support,
                        "deleted": deleted,
                        "extension_weight": extension_weight,
                        "cavity_weight": cavity_weight,
                        "normalised_slack": cavity_slack,
                    }
                    if sharp_fce is None or cavity_slack < sharp_fce[0]:
                        sharp_fce = (cavity_slack, cavity_row)

                fallback_gap = neighbour - (rank - 1) * weight
                assert fallback_gap >= 0
                fallback_slack = Fraction(fallback_gap, weight)
                fallback_row = {
                    "n": n,
                    "tree_index": tree_index,
                    "edges": sorted([sorted(edge) for edge in tree.edges()]),
                    "blocks": [list(block) for block in blocks],
                    "alpha": alpha,
                    "rank": rank,
                    "support": support,
                    "weight": weight,
                    "neighbour_weight": neighbour,
                    "normalised_slack": fallback_slack,
                }
                if sharp_fse is None or fallback_slack < sharp_fse[0]:
                    sharp_fse = (fallback_slack, fallback_row)

    assert sharp_fse is not None and sharp_fce is not None
    return {
        "at": utc_now(),
        "kind": "d26_fallback_support_exhaustive",
        "certificate": "passed",
        "max_order": max_order,
        "trees": trees,
        "prefix_supports": supports,
        "deletion_checks": deletion_checks,
        "sharp_fse": serialise(sharp_fse[1]),
        "sharp_fce": serialise(sharp_fce[1]),
        "statement": "FSE and the stronger per-deletion FCE passed every audited prefix support; this is an exact finite audit, not a proof.",
    }


def stress(trials: int, seed: int) -> dict[str, object]:
    rng = random.Random(seed)
    checks = 0
    deletion_checks = 0
    sharp_fse: tuple[Fraction, dict[str, object]] | None = None
    sharp_fce: tuple[Fraction, dict[str, object]] | None = None

    def audit(
        tree: nx.Graph,
        blocks: tuple[tuple[int, ...], ...],
        support: int,
        label: str,
    ) -> None:
        nonlocal checks, deletion_checks, sharp_fse, sharp_fce
        row = exchange_row(make_csp(tree, blocks), support)
        rank = int(row["rank"])
        weight = int(row["weight"])
        neighbour = int(row["neighbour_weight"])
        checks += 1

        fallback_gap = neighbour - (rank - 1) * weight
        assert fallback_gap >= 0, (label, row, blocks, sorted(tree.edges()))
        fallback_slack = Fraction(fallback_gap, weight)
        fallback_row = {
            "label": label,
            "alpha": row["alpha"],
            "rank": rank,
            "support": support,
            "weight": weight,
            "neighbour_weight": neighbour,
            "normalised_slack": fallback_slack,
        }
        if sharp_fse is None or fallback_slack < sharp_fse[0]:
            sharp_fse = (fallback_slack, fallback_row)

        for deleted, extension_weight, cavity_weight in row["cavity_rows"]:
            deletion_checks += 1
            cavity_gap = (
                rank * extension_weight - 2 * (rank - 1) * cavity_weight
            )
            assert cavity_gap >= 0, (
                label,
                row,
                deleted,
                blocks,
                sorted(tree.edges()),
            )
            cavity_slack = Fraction(cavity_gap, rank * cavity_weight)
            cavity_row = {
                "label": label,
                "alpha": row["alpha"],
                "rank": rank,
                "support": support,
                "deleted": deleted,
                "extension_weight": extension_weight,
                "cavity_weight": cavity_weight,
                "normalised_slack": cavity_slack,
            }
            if sharp_fce is None or cavity_slack < sharp_fce[0]:
                sharp_fce = (cavity_slack, cavity_row)

    obstruction_tree, obstruction_blocks, obstruction_support = witness()
    audit(
        obstruction_tree,
        obstruction_blocks,
        obstruction_support,
        "bare-SE-obstruction",
    )

    for trial in range(trials):
        n = rng.randint(15, 30)
        tree = random_tree(n, rng)
        blocks = maximum_matching_blocks(tree)
        if prefix_last(len(blocks)) >= 1:
            support = sampled_support(len(blocks), rng, trial % 2 == 0)
            audit(tree, blocks, support, f"random-tree-{trial}")

    families = ("path", "comb", "star", "random")
    for trial in range(trials):
        alpha = rng.randint(7, 30)
        family = families[trial % len(families)]
        if family == "path":
            contracted = nx.path_graph(alpha)
        elif family == "comb":
            contracted = comb_tree(alpha)
        elif family == "star":
            contracted = nx.star_graph(alpha - 1)
        else:
            contracted = random_tree(alpha, rng)
        tree, blocks = doubled_block_tree(contracted, rng)
        support = sampled_support(alpha, rng, trial % 3 != 0)
        audit(tree, blocks, support, f"block-{family}-{trial}")

    assert sharp_fse is not None and sharp_fce is not None
    return {
        "at": utc_now(),
        "kind": "d26_fallback_support_stress",
        "certificate": "passed",
        "seed": seed,
        "trials_per_family": trials,
        "support_checks": checks,
        "deletion_checks": deletion_checks,
        "sharp_fse": serialise(sharp_fse[1]),
        "sharp_fce": serialise(sharp_fce[1]),
        "statement": "FSE and FCE passed the conditioned-DP stress corpus, including the exact obstruction to bare SE; this is falsification evidence, not a proof.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("exhaustive", "stress"), required=True)
    parser.add_argument("--max-order", type=int, default=14)
    parser.add_argument("--trials", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=261994)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/d26_fallback_support_certificate_20260711.jsonl"),
    )
    args = parser.parse_args()
    result = (
        exhaustive(args.max_order)
        if args.mode == "exhaustive"
        else stress(args.trials, args.seed)
    )
    append_record(args.output, result)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
