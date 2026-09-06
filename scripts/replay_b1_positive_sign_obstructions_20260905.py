#!/usr/bin/env python3
"""Replay actual-rung refutations of b1-positive nonadversity shortcuts."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import networkx as nx

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from indpoly import is_unimodal
from scripts.audit_b1_zero_forbidden_core_20260905 import (
    audit_graph,
    coefficient,
    polynomial,
    roles,
)

SEED_GRAPH6 = "`pG`A?@O??g??@O???I??O???G???G?O???@?A?????_@?????@?A??????G?O??????C?G???????G????O????@"
MUTATED_GRAPH6 = "`oG`A?@O?_g??@O???I??O???G???G?O???@?A?????_@?????@?A??????G?O??????C?G???????G????O????@"


def profile(graph: nx.Graph) -> dict:
    row = audit_graph(graph)
    alpha = row["alpha"]
    full = polynomial(graph)
    allowed, forced = roles(graph, alpha)
    allowed_poly = polynomial(graph.subgraph(allowed))
    unary = [coefficient(full, alpha - d) - coefficient(allowed_poly, alpha - d)
             for d in range(5)]
    pair = [coefficient(allowed_poly, alpha - d) - row["e_0_to_4"][d]
            for d in range(5)]
    assert all(unary[d] + pair[d] == row["b_0_to_4"][d] for d in range(5))
    full_margin = full[alpha - 3] ** 2 - full[alpha - 2] * full[alpha - 4]
    assert (row["n"], alpha, row["deficiency"]) == (33, 19, 5)
    assert row["combined_correction"] < 0 and full_margin > 0
    assert is_unimodal(full)
    row.update(unary_b_0_to_4=unary, pair_b_0_to_4=pair,
               forbidden_vertices=sorted(set(graph) - allowed),
               forced_vertices=sorted(forced), full_depth3_margin=full_margin,
               independence_polynomial=full, unimodal=True)
    return row


def run() -> dict:
    seed = nx.from_graph6_bytes(SEED_GRAPH6.encode())
    changed = seed.copy()
    assert changed.degree(3) == 1 and changed.has_edge(3, 2)
    changed.remove_edge(3, 2)
    changed.add_edge(3, 10)
    assert nx.is_tree(changed)
    assert nx.to_graph6_bytes(changed, header=False).decode().strip() == MUTATED_GRAPH6
    before, after = profile(seed), profile(changed)
    assert before["unary_b_0_to_4"][1] == 0
    assert before["pair_b_0_to_4"][1] == 8
    assert after["unary_b_0_to_4"][1] == 8
    assert after["pair_b_0_to_4"][1] == 0
    assert after["combined_correction"] == -178212783
    return {"kind": "b1_positive_sign_obstruction_replay", "certificate_date": "2026-09-05",
            "scope": "Actual n=33, alpha=19, deficiency=5 trees; auxiliary refutations only",
            "discovery": "Ninth admissible single-leaf reattachment tested from archived b1-positive adverse seeds",
            "mutation": {"leaf": 3, "old_parent": 2, "new_parent": 10},
            "seed": before, "unary_b1_positive_counterexample": after,
            "conclusion": "Neither pair-only nor unary-only defect-one obstruction makes D nonnegative. Both full polynomials remain unimodal, with positive depth-three margins."}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=REPO / "results/b1_positive_sign_obstructions_20260905.json")
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    row = result["unary_b1_positive_counterexample"]
    print(json.dumps({key: row[key] for key in
                     ("unary_b_0_to_4", "pair_b_0_to_4", "combined_correction",
                      "full_depth3_margin", "unimodal")}, indent=2))


if __name__ == "__main__":
    main()
