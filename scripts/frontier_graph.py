#!/usr/bin/env python3
"""Validate and inspect the Erdős #993 proof-frontier graph."""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]

CLAIM_STATUSES = {
    "open",
    "candidate",
    "speculative",
    "verified_bounded",
    "proved",
    "formalized",
    "cited",
    "refuted",
}
GROUNDED_STATUSES = {"verified_bounded", "proved", "formalized", "cited", "refuted"}
RESOLVED_STATUSES = {"proved", "formalized", "cited", "refuted"}
EDGE_TYPES = {
    "depends_on",
    "reduces_to",
    "formalizes",
    "cites",
    "refutes",
    "tests",
}
DEPENDENCY_EDGE_TYPES = {"depends_on", "reduces_to"}
LOCAL_EVIDENCE_TYPES = {"note", "script", "result", "lean", "source"}
BRIDGE_STATUSES = {"candidate", "speculative", "adapted", "refuted"}
LEVEL_ORDER = {"high": 0, "medium": 1, "low": 2}
EFFORT_ORDER = {"low": 0, "medium": 1, "high": 2}
TRANSLATION_ORDER = {"complete": 0, "partial": 1, "analogy": 2}
BRIDGE_STATUS_ORDER = {"adapted": 0, "candidate": 1, "speculative": 2, "refuted": 3}


class GraphValidationError(ValueError):
    """Raised when the frontier graph violates its audit contract."""


def load_graph(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        graph = json.load(handle)
    if not isinstance(graph, dict):
        raise GraphValidationError("graph root must be a JSON object")
    return graph


def _require_text(item: dict[str, Any], key: str, context: str) -> None:
    if not isinstance(item.get(key), str) or not item[key].strip():
        raise GraphValidationError(f"{context}: missing nonempty {key!r}")


def _validate_evidence(claim: dict[str, Any], context: str) -> None:
    evidence = claim.get("evidence")
    if not isinstance(evidence, list):
        raise GraphValidationError(f"{context}: evidence must be a list")
    if claim["status"] in GROUNDED_STATUSES and not evidence:
        raise GraphValidationError(
            f"{context}: status {claim['status']!r} requires grounded evidence"
        )
    for index, entry in enumerate(evidence):
        entry_context = f"{context}.evidence[{index}]"
        if not isinstance(entry, dict):
            raise GraphValidationError(f"{entry_context}: must be an object")
        _require_text(entry, "type", entry_context)
        _require_text(entry, "ref", entry_context)
        if entry["type"] in LOCAL_EVIDENCE_TYPES:
            evidence_path = REPO_ROOT / entry["ref"]
            if not evidence_path.exists():
                raise GraphValidationError(
                    f"{entry_context}: local evidence does not exist: {entry['ref']}"
                )


def _dependency_adjacency(
    claim_ids: set[str], edges: list[dict[str, Any]]
) -> dict[str, set[str]]:
    adjacency = {claim_id: set() for claim_id in claim_ids}
    for edge in edges:
        if edge["type"] in DEPENDENCY_EDGE_TYPES:
            adjacency[edge["from"]].add(edge["to"])
    return adjacency


def _validate_acyclic(adjacency: dict[str, set[str]]) -> None:
    state = {node: 0 for node in adjacency}
    trail: list[str] = []

    def visit(node: str) -> None:
        if state[node] == 1:
            cycle_start = trail.index(node)
            cycle = trail[cycle_start:] + [node]
            raise GraphValidationError(
                "dependency cycle: " + " -> ".join(cycle)
            )
        if state[node] == 2:
            return
        state[node] = 1
        trail.append(node)
        for dependency in sorted(adjacency[node]):
            visit(dependency)
        trail.pop()
        state[node] = 2

    for node in sorted(adjacency):
        visit(node)


def validate_graph(graph: dict[str, Any]) -> list[str]:
    """Validate the graph and return non-fatal warnings."""
    for key in ("schema_version", "graph_id", "title", "root_claim", "claims", "edges", "bridges"):
        if key not in graph:
            raise GraphValidationError(f"graph: missing {key!r}")
    if graph["schema_version"] != 1:
        raise GraphValidationError("graph: unsupported schema_version")

    claims = graph["claims"]
    if not isinstance(claims, list) or not claims:
        raise GraphValidationError("graph: claims must be a nonempty list")
    claim_ids: set[str] = set()
    for index, claim in enumerate(claims):
        context = f"claims[{index}]"
        if not isinstance(claim, dict):
            raise GraphValidationError(f"{context}: must be an object")
        for key in ("id", "kind", "slogan", "statement", "scope", "falsification"):
            _require_text(claim, key, context)
        if claim["id"] in claim_ids:
            raise GraphValidationError(f"{context}: duplicate id {claim['id']!r}")
        claim_ids.add(claim["id"])
        if claim.get("status") not in CLAIM_STATUSES:
            raise GraphValidationError(f"{context}: invalid status {claim.get('status')!r}")
        priority = claim.get("priority")
        if not isinstance(priority, dict):
            raise GraphValidationError(f"{context}: priority must be an object")
        if priority.get("falsifiability") not in LEVEL_ORDER:
            raise GraphValidationError(f"{context}: invalid falsifiability")
        if priority.get("verification_cost") not in EFFORT_ORDER:
            raise GraphValidationError(f"{context}: invalid verification_cost")
        search_profiles = claim.get("search_profiles", {})
        if not isinstance(search_profiles, dict) or any(
            not isinstance(values, list) or not all(isinstance(value, str) for value in values)
            for values in search_profiles.values()
        ):
            raise GraphValidationError(f"{context}: search_profiles must map to string lists")
        _validate_evidence(claim, context)

    if graph["root_claim"] not in claim_ids:
        raise GraphValidationError("graph: root_claim is not a claim id")

    edges = graph["edges"]
    if not isinstance(edges, list):
        raise GraphValidationError("graph: edges must be a list")
    seen_edges: set[tuple[str, str, str]] = set()
    for index, edge in enumerate(edges):
        context = f"edges[{index}]"
        if not isinstance(edge, dict):
            raise GraphValidationError(f"{context}: must be an object")
        for key in ("from", "to", "type"):
            _require_text(edge, key, context)
        if edge["from"] not in claim_ids or edge["to"] not in claim_ids:
            raise GraphValidationError(f"{context}: references an unknown claim")
        if edge["from"] == edge["to"]:
            raise GraphValidationError(f"{context}: self-edge")
        if edge["type"] not in EDGE_TYPES:
            raise GraphValidationError(f"{context}: invalid type {edge['type']!r}")
        signature = (edge["from"], edge["to"], edge["type"])
        if signature in seen_edges:
            raise GraphValidationError(f"{context}: duplicate edge {signature}")
        seen_edges.add(signature)

    _validate_acyclic(_dependency_adjacency(claim_ids, edges))

    bridges = graph["bridges"]
    if not isinstance(bridges, list):
        raise GraphValidationError("graph: bridges must be a list")
    bridge_ids: set[str] = set()
    for index, bridge in enumerate(bridges):
        context = f"bridges[{index}]"
        if not isinstance(bridge, dict):
            raise GraphValidationError(f"{context}: must be an object")
        for key in (
            "id",
            "title",
            "source_theory",
            "shared_object",
            "translation",
            "expected_payoff",
            "kill_test",
            "falsification",
        ):
            _require_text(bridge, key, context)
        if bridge["id"] in bridge_ids:
            raise GraphValidationError(f"{context}: duplicate id {bridge['id']!r}")
        bridge_ids.add(bridge["id"])
        if bridge.get("status") not in BRIDGE_STATUSES:
            raise GraphValidationError(f"{context}: invalid bridge status")
        targets = bridge.get("target_claims")
        if not isinstance(targets, list) or not targets:
            raise GraphValidationError(f"{context}: target_claims must be nonempty")
        if any(target not in claim_ids for target in targets):
            raise GraphValidationError(f"{context}: target_claims contains unknown id")
        missing = bridge.get("missing_hypotheses")
        if not isinstance(missing, list) or not all(
            isinstance(item, str) and item.strip() for item in missing
        ):
            raise GraphValidationError(f"{context}: missing_hypotheses must be a string list")
        priority = bridge.get("priority")
        if not isinstance(priority, dict):
            raise GraphValidationError(f"{context}: priority must be an object")
        if priority.get("falsifiability") not in LEVEL_ORDER:
            raise GraphValidationError(f"{context}: invalid falsifiability")
        if priority.get("verification_cost") not in EFFORT_ORDER:
            raise GraphValidationError(f"{context}: invalid verification_cost")
        if priority.get("translation_completeness") not in TRANSLATION_ORDER:
            raise GraphValidationError(f"{context}: invalid translation_completeness")

    warnings: list[str] = []
    for claim in claims:
        if claim["status"] not in GROUNDED_STATUSES and claim["evidence"]:
            warnings.append(
                f"{claim['id']}: evidence is candidate/supporting only because status={claim['status']}"
            )
    return warnings


def _claim_maps(graph: dict[str, Any]) -> tuple[dict[str, dict[str, Any]], dict[str, set[str]], dict[str, set[str]]]:
    claims = {claim["id"]: claim for claim in graph["claims"]}
    dependencies = _dependency_adjacency(set(claims), graph["edges"])
    dependents = {claim_id: set() for claim_id in claims}
    for claim_id, deps in dependencies.items():
        for dependency in deps:
            dependents[dependency].add(claim_id)
    return claims, dependencies, dependents


def _transitive_reach(start: str, adjacency: dict[str, set[str]]) -> set[str]:
    reached: set[str] = set()
    stack = list(adjacency[start])
    while stack:
        node = stack.pop()
        if node in reached:
            continue
        reached.add(node)
        stack.extend(adjacency[node] - reached)
    return reached


def graph_report(graph: dict[str, Any]) -> str:
    claims, dependencies, dependents = _claim_maps(graph)
    unresolved = {
        claim_id for claim_id, claim in claims.items() if claim["status"] not in RESOLVED_STATUSES
    }
    frontier = [
        claim_id
        for claim_id in unresolved
        if not any(dependency in unresolved for dependency in dependencies[claim_id])
    ]
    ranked_claims = sorted(
        unresolved,
        key=lambda claim_id: (
            -len(_transitive_reach(claim_id, dependents)),
            LEVEL_ORDER[claims[claim_id]["priority"]["falsifiability"]],
            EFFORT_ORDER[claims[claim_id]["priority"]["verification_cost"]],
            claim_id,
        ),
    )

    bridge_rows = []
    for bridge in graph["bridges"]:
        cut_value = max(
            len(_transitive_reach(target, dependents)) for target in bridge["target_claims"]
        )
        bridge_rows.append((bridge, cut_value))
    bridge_rows.sort(
        key=lambda row: (
            -row[1],
            LEVEL_ORDER[row[0]["priority"]["falsifiability"]],
            EFFORT_ORDER[row[0]["priority"]["verification_cost"]],
            TRANSLATION_ORDER[row[0]["priority"]["translation_completeness"]],
            BRIDGE_STATUS_ORDER[row[0]["status"]],
            row[0]["id"],
        )
    )

    lines = [
        f"Graph: {graph['title']}",
        f"Root: {graph['root_claim']}",
        f"Claims: {len(claims)} ({len(unresolved)} unresolved)",
        "",
        "Ready frontier nodes:",
    ]
    for claim_id in sorted(
        frontier,
        key=lambda item: (
            -len(_transitive_reach(item, dependents)),
            item,
        ),
    ):
        claim = claims[claim_id]
        reach = len(_transitive_reach(claim_id, dependents))
        lines.append(f"- {claim_id} [{claim['status']}; downstream={reach}]: {claim['slogan']}")
    lines.extend(["", "Unresolved nodes by proof leverage:"])
    for claim_id in ranked_claims:
        claim = claims[claim_id]
        reach = len(_transitive_reach(claim_id, dependents))
        lines.append(f"- {claim_id} [{claim['status']}; downstream={reach}]")
    lines.extend(["", "Bridge candidates:"])
    for bridge, cut_value in bridge_rows:
        lines.append(
            f"- {bridge['id']} [{bridge['status']}; cut={cut_value}; "
            f"falsify={bridge['priority']['falsifiability']}; "
            f"cost={bridge['priority']['verification_cost']}; "
            f"map={bridge['priority']['translation_completeness']}]: {bridge['title']}"
        )
    return "\n".join(lines)


def search_packet(graph: dict[str, Any], claim_id: str) -> str:
    claims, dependencies, _ = _claim_maps(graph)
    if claim_id not in claims:
        raise GraphValidationError(f"unknown claim id: {claim_id}")
    claim = claims[claim_id]
    lines = [
        f"# Search packet: {claim_id}",
        "",
        f"Status: {claim['status']}",
        f"Scope: {claim['scope']}",
        "",
        "## Exact target",
        "",
        claim["statement"],
        "",
        "## Refutation condition",
        "",
        claim["falsification"],
        "",
        "## Dependencies",
        "",
    ]
    for dependency in sorted(dependencies[claim_id]):
        lines.append(f"- {dependency} [{claims[dependency]['status']}]: {claims[dependency]['slogan']}")
    lines.extend(["", "## Search formulations", ""])
    for profile, queries in claim.get("search_profiles", {}).items():
        lines.append(f"### {profile.replace('_', ' ').title()}")
        lines.append("")
        lines.extend(f"- {query}" for query in queries)
        lines.append("")
    relevant_bridges = [
        bridge for bridge in graph["bridges"] if claim_id in bridge["target_claims"]
    ]
    if relevant_bridges:
        lines.extend(["## Candidate bridges", ""])
        for bridge in relevant_bridges:
            lines.extend(
                [
                    f"### {bridge['title']}",
                    "",
                    f"- Shared object: {bridge['shared_object']}",
                    f"- Translation: {bridge['translation']}",
                    f"- Missing hypotheses: {'; '.join(bridge['missing_hypotheses']) or 'none identified'}",
                    f"- Kill-test: {bridge['kill_test']}",
                    f"- Falsification: {bridge['falsification']}",
                    "",
                ]
            )
    return "\n".join(lines).rstrip() + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for name in ("validate", "report"):
        command = subparsers.add_parser(name)
        command.add_argument("graph", type=Path)
    packet = subparsers.add_parser("packet")
    packet.add_argument("graph", type=Path)
    packet.add_argument("claim_id")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        graph = load_graph(args.graph)
        warnings = validate_graph(graph)
        if args.command == "validate":
            print(
                f"PASS: {graph['graph_id']} ({len(graph['claims'])} claims, "
                f"{len(graph['edges'])} edges, {len(graph['bridges'])} bridges)"
            )
            for warning in warnings:
                print(f"WARNING: {warning}")
        elif args.command == "report":
            print(graph_report(graph))
        else:
            print(search_packet(graph, args.claim_id), end="")
    except (OSError, json.JSONDecodeError, GraphValidationError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
