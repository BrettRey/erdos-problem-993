#!/usr/bin/env python3
"""Exact replay of three fixed witnesses; no search and no output-file writes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import networkx as nx

import indpoly
from scripts.audit_b1_zero_forbidden_core_20260905 import (
    coefficient,
    extendable_defects,
    polynomial,
    roles,
)

# Make this replay independent of NumPy availability and fixed-width arithmetic.
indpoly._HAS_NUMPY = False
ROOT = Path(__file__).resolve().parent


def require(condition: bool, message: str) -> None:
    """Do not let python -O disable certificate comparisons."""
    if not condition:
        raise AssertionError(message)


def analyze(graph: nx.Graph) -> dict:
    """Recompute intrinsic counts and distinguish the three failure levels."""
    if graph.is_directed() or graph.is_multigraph():
        raise ValueError("expected a finite simple undirected tree")
    if not graph or nx.number_of_selfloops(graph) or not nx.is_tree(graph):
        raise ValueError("expected a nonempty tree")
    graph = nx.convert_node_labels_to_integers(graph)
    poly = polynomial(graph)
    alpha = len(poly) - 1
    if alpha < 4:
        raise ValueError("depth-three replay requires independence number >= 4")
    e = extendable_defects(graph)
    s = [coefficient(poly, alpha - d) for d in range(5)]
    b = [s[d] - e[d] for d in range(5)]
    require(min(b) >= 0 and b[0] == 0, "invalid extendable/blocked split")
    blocked_turan = b[3] ** 2 - b[2] * b[4]
    cross = 2 * e[3] * b[3] - e[2] * b[4] - e[4] * b[2]
    correction = blocked_turan + cross
    reserve = e[3] ** 2 - e[2] * e[4]
    full_margin = s[3] ** 2 - s[2] * s[4]
    require(full_margin == reserve + correction, "decomposition mismatch")
    pascal_slack = 27 * (alpha - 3) * e[3] ** 2 - 32 * (alpha - 2) * e[2] * e[4]
    shadow_slack = 6 * b[3] - (alpha - 4) * b[2]
    incidence_slack = (alpha - 3) * e[3] - 4 * e[4]
    require(min(pascal_slack, shadow_slack, incidence_slack) >= 0,
            "supplied forest inequality or counting implementation failed")
    n = len(graph)
    in_window = 33 <= n <= 38 and 17 <= alpha <= 19 and 2 * alpha <= n + 5
    low_density = 3 * (alpha - 3) * b[4] <= (alpha - 7) * e[4]
    sufficient_margin = (5 * alpha + 17) * e[2] * e[4] + 27 * (alpha - 3) * correction
    endpoint_margin = (5 * alpha + 17) * e[2] * e[4] - 27 * (alpha - 3) * (
        e[2] * b[4] + e[4] * b[2] + b[2] * b[4]
    )
    return {
        "n": n,
        "alpha": alpha,
        "deficiency": 2 * alpha - n,
        "graph6": nx.to_graph6_bytes(graph, header=False).decode().strip(),
        "independence_polynomial": poly,
        "e_0_to_4": e,
        "b_0_to_4": b,
        "s_0_to_4": s,
        "blocked_turan": blocked_turan,
        "cross_term": cross,
        "combined_correction": correction,
        "extendable_turan": reserve,
        "full_margin": full_margin,
        "scaled_sufficient_joint_margin_S": sufficient_margin,
        "scaled_endpoint_margin_E": endpoint_margin,
        "pascal_reserve_slack": pascal_slack,
        "blocked_shadow_slack": shadow_slack,
        "extendable_incidence_slack": incidence_slack,
        "in_target_window": in_window,
        "low_density": low_density,
        "in_residual_regime_R": in_window and b[1] > 0 and not low_density and correction < 0,
        "target_W_holds_for_this_graph": full_margin >= 0,
        "unimodal": indpoly.is_unimodal(poly),
    }


def selected_witnesses() -> list[tuple[str, dict]]:
    """Select exactly the three witnesses documented in the brief."""
    signs = json.loads((ROOT / "evidence/b1_positive_sign_obstructions_20260905.json").read_text())
    lifts = json.loads((ROOT / "evidence/cross_reserve_witness_lifts_20260904.json").read_text())
    return [
        ("pair_only_defect_one_adverse", signs["seed"]),
        ("unary_only_defect_one_adverse", signs["unary_b1_positive_counterexample"]),
        ("blocked_profile_not_log_concave", lifts["first_witnesses"]["blocked_lc_failure"]),
    ]


def replay() -> dict:
    rows = []
    for label, expected in selected_witnesses():
        graph = nx.from_graph6_bytes(expected["graph6"].encode("ascii"))
        row = analyze(graph)
        for key in ("n", "alpha", "deficiency", "e_0_to_4", "b_0_to_4",
                    "combined_correction", "unimodal"):
            require(row[key] == expected[key], f"{label}: mismatch in {key}")
        expected_margin = expected.get("full_depth3_margin", expected.get("full_margin"))
        require(row["full_margin"] == expected_margin, f"{label}: full margin mismatch")
        for key in ("independence_polynomial", "s_0_to_4", "blocked_turan",
                    "cross_term", "extendable_turan"):
            if key in expected:
                require(row[key] == expected[key], f"{label}: mismatch in {key}")
        if "unary_b_0_to_4" in expected:
            allowed, forced = roles(graph, row["alpha"])
            allowed_poly = polynomial(graph.subgraph(allowed))
            unary = [row["s_0_to_4"][d] - coefficient(allowed_poly, row["alpha"] - d)
                     for d in range(5)]
            pair = [row["b_0_to_4"][d] - unary[d] for d in range(5)]
            require(unary == expected["unary_b_0_to_4"], f"{label}: unary split mismatch")
            require(pair == expected["pair_b_0_to_4"], f"{label}: pair split mismatch")
            require(sorted(forced) == expected["forced_vertices"], f"{label}: forced mismatch")
            require(sorted(set(graph) - allowed) == expected["forbidden_vertices"],
                    f"{label}: forbidden mismatch")
            require(row["scaled_endpoint_margin_E"] == expected["joint_margin_numerator"],
                    f"{label}: historical endpoint margin mismatch")
            row.update(unary_b_0_to_4=unary, pair_b_0_to_4=pair)
        if label == "blocked_profile_not_log_concave":
            require(row["blocked_turan"] == -3960, "blocked-profile counterexample lost")
            require(row["scaled_endpoint_margin_E"] ==
                    expected["pascal_combined_endpoint_margin_numerator"],
                    "lift: historical endpoint margin mismatch")
        require(row["combined_correction"] < 0, f"{label}: correction should be adverse")
        require(row["full_margin"] > 0 and row["unimodal"],
                f"{label}: historical target/unimodality claim changed")
        row["label"] = label
        rows.append(row)
    return {
        "status": "PASS",
        "fixed_witnesses_recomputed": len(rows),
        "arithmetic": "Python arbitrary-precision integers; NumPy disabled",
        "warning": "These refute auxiliary shortcuts only; no universal theorem is proved.",
        "witnesses": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph6", help="analyze one additional graph instead of the saved witnesses")
    args = parser.parse_args()
    if args.graph6:
        result = analyze(nx.from_graph6_bytes(args.graph6.encode("ascii")))
    else:
        result = replay()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
