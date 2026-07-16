#!/usr/bin/env python3
"""Compare Will Blair's published rooted-bush witnesses with local hard corpora.

This is an ingestion and collision-checking tool, not a new tree search.  It
reconstructs the 25 witnesses serialized in Blair's public ``results.json``,
checks their orders and exact independence polynomials, and compares them with
the local family grids, retained graph6 certificates, and explicit-edge
literature root-stress cases that can be rebuilt from committed artifacts.

The full set of 4,445 distinct non-log-concave polynomials claimed by the
upstream search is not serialized publicly.  This script therefore does not
claim to compare that full set and does not invoke the upstream V2 sweep.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import platform
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

import networkx as nx

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from graph6 import parse_graph6  # noqa: E402
from indpoly import independence_poly, is_log_concave, is_unimodal  # noqa: E402
from scripts.analyze_prufer_corpus import (  # noqa: E402
    canonical_graph6_batch,
    graph6_from_adj,
    iter_family_records,
    make_galvin_tree,
    make_li_tree,
)
from scripts.valley_search import bouquet_adj, bouquet_poly  # noqa: E402


PUBLIC_URL = "https://github.com/willblair0708/verified-combinatorics/tree/main/erdos-993"
LABEL_RE = re.compile(r"^bush\((?P<body>[^)]*)\)\|n=(?P<n>\d+)$")
TOKEN_RE = re.compile(r"^(?P<kind>[us])(?P<c>\d+)p(?P<length>\d+)$")
MATCH_LIMIT = 25

FAMILY_ARTIFACTS = (
    "results/galvin_family_grid_m1-50_t1-50_nle500_summary.json",
    "results/bautista_ramos_family_grid_m1-8_t1-8_summary.json",
    "results/li_family_grid_m1-50_n1-50_summary.json",
    "results/li_star_family_grid_m1-50_n1-50_summary.json",
)

EXACT_RETAINED_ARTIFACTS = (
    "results/analysis_n26.json",
    "results/analysis_n27_modal_lc_nm.json",
    "results/analysis_n28_modal_lc_nm.json",
)

EXPLICIT_EDGE_ARTIFACTS = (
    "results/literature_root_stress_20260716.json",
)

EVOLUTIONARY_TREE_GLOB = "best_evolutionary_tree*.json"
LC_BREAKER_GLOB = "lc_breaker_evo_n26*.json"


@dataclass(frozen=True)
class BushToken:
    kind: str
    child_count: int
    pendant_length: int

    @property
    def legs(self) -> tuple[int, ...]:
        """Translate Blair's token to the local center-to-leaf leg lengths."""
        legs = [self.pendant_length + 1] * self.child_count
        if self.kind == "s":
            legs[0] += 1
        return tuple(sorted(legs))

    def payload(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "child_count": self.child_count,
            "pendant_length": self.pendant_length,
            "local_spider_legs": list(self.legs),
        }


@dataclass
class LocalCandidate:
    n: int
    adj: list[list[int]]
    graph6: str
    sources: list[dict[str, Any]]


class LocalIndex:
    """Streaming exact indexes for local polynomials, signatures, and trees."""

    def __init__(self, remote_orders: set[int]) -> None:
        self.remote_orders = remote_orders
        self.polynomials: dict[tuple[int, ...], list[dict[str, Any]]] = defaultdict(list)
        self.signatures: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
        self.tree_candidates: list[LocalCandidate] = []
        self.source_counts: Counter[str] = Counter()
        self.record_count = 0
        self.source_reference_count = 0
        self.max_near_miss = Fraction(0, 1)
        self.max_near_miss_sources: list[dict[str, Any]] = []
        self.max_near_miss_by_order: dict[int, Fraction] = {}
        self.max_near_miss_sources_by_order: dict[int, list[dict[str, Any]]] = {}

    def add(
        self,
        *,
        adj: list[list[int]],
        poly: list[int],
        sources: list[dict[str, Any]],
        source_group: str,
        graph6: str | None = None,
    ) -> None:
        n = len(adj)
        validate_tree(adj)
        if len(poly) < 2 or poly[0] != 1 or poly[1] != n:
            raise ValueError(f"invalid independence polynomial header for local order {n}")

        metrics = exact_metrics(poly)
        poly_key = tuple(poly)
        signature_key = exact_signature(n, metrics)
        self.polynomials[poly_key].extend(sources)
        self.signatures[signature_key].extend(sources)
        self.record_count += 1
        self.source_reference_count += len(sources)
        self.source_counts[source_group] += 1

        near = fraction_from_payload(metrics["post_descent_max_ratio"])
        if near > self.max_near_miss:
            self.max_near_miss = near
            self.max_near_miss_sources = list(sources)
        elif near == self.max_near_miss:
            self.max_near_miss_sources.extend(sources)

        order_best = self.max_near_miss_by_order.get(n, Fraction(-1, 1))
        if near > order_best:
            self.max_near_miss_by_order[n] = near
            self.max_near_miss_sources_by_order[n] = list(sources)
        elif near == order_best:
            self.max_near_miss_sources_by_order[n].extend(sources)

        if n in self.remote_orders:
            effective_graph6 = graph6 if graph6 is not None else graph6_from_adj(adj)
            self.tree_candidates.append(
                LocalCandidate(
                    n=n,
                    adj=adj,
                    graph6=effective_graph6,
                    sources=list(sources),
                )
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--blair-repo",
        type=Path,
        required=True,
        help="Path to verified-combinatorics or its erdos-993 subdirectory",
    )
    parser.add_argument("--output", type=Path, required=True, help="Comparison JSON to write")
    return parser.parse_args()


def resolve_blair_dir(path: Path) -> Path:
    path = path.expanduser().resolve()
    if (path / "results.json").is_file():
        return path
    nested = path / "erdos-993"
    if (nested / "results.json").is_file():
        return nested
    raise FileNotFoundError(f"could not find results.json below {path}")


def git_revision(path: Path) -> str | None:
    proc = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        check=False,
    )
    return proc.stdout.strip() if proc.returncode == 0 else None


def parse_bush_label(label: str) -> tuple[list[BushToken], int]:
    match = LABEL_RE.fullmatch(label)
    if match is None:
        raise ValueError(f"unrecognized Blair bush label: {label!r}")
    tokens: list[BushToken] = []
    for raw in match.group("body").split(","):
        token_match = TOKEN_RE.fullmatch(raw)
        if token_match is None:
            raise ValueError(f"unrecognized Blair branch token {raw!r} in {label!r}")
        token = BushToken(
            kind=token_match.group("kind"),
            child_count=int(token_match.group("c")),
            pendant_length=int(token_match.group("length")),
        )
        if token.child_count < 1 or token.pendant_length < 0:
            raise ValueError(f"invalid Blair branch token {raw!r}")
        tokens.append(token)
    if not tokens:
        raise ValueError(f"empty Blair bush label: {label!r}")
    return tokens, int(match.group("n"))


def reconstruct_bush(tokens: list[BushToken]) -> tuple[list[list[int]], list[int]]:
    gadgets = [token.legs for token in tokens]
    n, adj = bouquet_adj(tuple(gadgets))
    fast_poly = bouquet_poly(gadgets)
    generic_poly = independence_poly(n, adj)
    if fast_poly != generic_poly:
        raise AssertionError("local bouquet formula disagrees with generic tree DP")
    return adj, generic_poly


def validate_tree(adj: list[list[int]]) -> None:
    n = len(adj)
    if n == 0:
        raise ValueError("empty graph is not a tree witness")
    if sum(len(row) for row in adj) != 2 * (n - 1):
        raise ValueError(f"order-{n} graph does not have n-1 edges")
    seen = {0}
    stack = [0]
    while stack:
        vertex = stack.pop()
        for neighbor in adj[vertex]:
            if neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    if len(seen) != n:
        raise ValueError(f"order-{n} graph is disconnected")


def rational_payload(value: Fraction) -> dict[str, Any]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def fraction_from_payload(payload: dict[str, Any]) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def polynomial_sha256(poly: Iterable[int]) -> str:
    raw = json.dumps(list(poly), separators=(",", ":")).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def exact_metrics(poly: list[int]) -> dict[str, Any]:
    peak = max(poly)
    modes = [idx for idx, value in enumerate(poly) if value == peak]
    first_descent = next(
        (idx for idx in range(1, len(poly)) if poly[idx] < poly[idx - 1]),
        None,
    )

    near_ratio = Fraction(0, 1)
    near_pos: int | None = None
    if first_descent is not None:
        for idx in range(first_descent, len(poly) - 1):
            if poly[idx] == 0:
                continue
            ratio = Fraction(poly[idx + 1], poly[idx])
            if ratio > near_ratio:
                near_ratio = ratio
                near_pos = idx

    lc_defects = []
    max_lc_ratio = Fraction(0, 1)
    max_lc_pos: int | None = None
    severity = Fraction(0, 1)
    severity_pos: int | None = None
    for idx in range(1, len(poly) - 1):
        center_square = poly[idx] * poly[idx]
        flanks = poly[idx - 1] * poly[idx + 1]
        if center_square == 0:
            continue
        lc_ratio = Fraction(flanks, center_square)
        signed_severity = Fraction(center_square - flanks, center_square)
        if lc_ratio > max_lc_ratio:
            max_lc_ratio = lc_ratio
            max_lc_pos = idx
        if signed_severity < severity:
            severity = signed_severity
            severity_pos = idx
        if flanks > center_square:
            lc_defects.append(
                {
                    "index": idx,
                    "absolute_defect": flanks - center_square,
                    "flank_product_over_center_square": rational_payload(lc_ratio),
                }
            )

    return {
        "alpha": len(poly) - 1,
        "mode_first": modes[0],
        "mode_last": modes[-1],
        "first_descent": first_descent,
        "unimodal": is_unimodal(poly),
        "log_concave": is_log_concave(poly),
        "lc_defect_count": len(lc_defects),
        "lc_defects": lc_defects,
        "max_lc_ratio": rational_payload(max_lc_ratio),
        "max_lc_ratio_index": max_lc_pos,
        "normalized_lc_severity": rational_payload(severity),
        "normalized_lc_severity_index": severity_pos,
        "post_descent_max_ratio": rational_payload(near_ratio),
        "post_descent_max_ratio_index": near_pos,
    }


def exact_signature(n: int, metrics: dict[str, Any]) -> tuple[Any, ...]:
    near = metrics["post_descent_max_ratio"]
    return (
        n,
        metrics["alpha"],
        metrics["mode_first"],
        metrics["mode_last"],
        metrics["first_descent"],
        tuple(defect["index"] for defect in metrics["lc_defects"]),
        metrics["post_descent_max_ratio_index"],
        near["numerator"],
        near["denominator"],
    )


def signature_payload(n: int, metrics: dict[str, Any]) -> dict[str, Any]:
    return {
        "n": n,
        "alpha": metrics["alpha"],
        "mode_interval": [metrics["mode_first"], metrics["mode_last"]],
        "first_descent": metrics["first_descent"],
        "lc_defect_indices": [defect["index"] for defect in metrics["lc_defects"]],
        "post_descent_max_ratio_index": metrics["post_descent_max_ratio_index"],
        "post_descent_max_ratio": metrics["post_descent_max_ratio"],
    }


def adjacency_to_networkx(adj: list[list[int]]) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(len(adj)))
    graph.add_edges_from((u, v) for u, row in enumerate(adj) for v in row if u < v)
    return graph


def source_id(prefix: str, **fields: Any) -> str:
    tail = ",".join(f"{key}={fields[key]}" for key in sorted(fields))
    return f"{prefix}:{tail}" if tail else prefix


def compact_matches(matches: list[dict[str, Any]]) -> dict[str, Any]:
    unique = {item["id"]: item for item in matches}
    ordered = [unique[key] for key in sorted(unique)]
    return {
        "count": len(ordered),
        "matches": ordered[:MATCH_LIMIT],
        "truncated": len(ordered) > MATCH_LIMIT,
    }


def reconstruct_family_grids(index: LocalIndex) -> list[dict[str, Any]]:
    reports = []
    for relative in FAMILY_ARTIFACTS:
        path = REPO / relative
        summary = json.loads(path.read_text(encoding="utf-8"))
        source = summary["source"]
        if source.get("kind") != "generated_family":
            raise ValueError(f"{relative} is not a generated-family artifact")
        records = iter_family_records(
            family=source["family"],
            m_values=[int(value) for value in source["m_values"]],
            t_values=[int(value) for value in source["t_values"]],
            max_n=source.get("max_n"),
        )
        if len(records) != summary["processed"]:
            raise AssertionError(
                f"{relative}: reconstructed {len(records)} rows, expected {summary['processed']}"
            )

        non_lc = 0
        non_unimodal = 0
        for record in records:
            poly = independence_poly(len(record.adj), record.adj)
            non_lc += not is_log_concave(poly)
            non_unimodal += not is_unimodal(poly)
            params = dict(record.params)
            ref = {
                "id": source_id(f"grid:{record.family}", **params),
                "kind": "generated_family_grid",
                "artifact": relative,
                "family": record.family,
                "params": params,
                "rank": record.rank,
            }
            index.add(
                adj=record.adj,
                poly=poly,
                sources=[ref],
                source_group=f"grid:{record.family}",
            )

        expected_counts = summary["counts"]
        if non_lc != expected_counts["non_log_concave"]:
            raise AssertionError(f"{relative}: non-log-concave count changed")
        if non_unimodal != expected_counts["non_unimodal"]:
            raise AssertionError(f"{relative}: non-unimodal count changed")
        reports.append(
            {
                "artifact": relative,
                "family": source["family"],
                "reconstructed": len(records),
                "non_log_concave": non_lc,
                "non_unimodal": non_unimodal,
                "summary_counts_verified": True,
            }
        )
    return reports


def json_pointer_part(value: Any) -> str:
    return str(value).replace("~", "~0").replace("/", "~1")


def walk_graph6_rows(value: Any, pointer: str = "") -> Iterable[tuple[str, str, dict[str, Any]]]:
    if isinstance(value, dict):
        graph6 = value.get("graph6")
        if isinstance(graph6, str):
            yield pointer or "/", graph6, value
        for key, child in value.items():
            yield from walk_graph6_rows(child, f"{pointer}/{json_pointer_part(key)}")
    elif isinstance(value, list):
        for idx, child in enumerate(value):
            yield from walk_graph6_rows(child, f"{pointer}/{idx}")


def retained_artifact_paths() -> list[Path]:
    paths = [REPO / relative for relative in EXACT_RETAINED_ARTIFACTS]
    paths.extend(sorted((REPO / "results").glob("ramos_sun_60_epoch*_summary.json")))
    return paths


def ingest_retained_graph6(index: LocalIndex) -> list[dict[str, Any]]:
    by_graph6: dict[str, dict[str, Any]] = {}
    artifact_reports = []
    for path in retained_artifact_paths():
        payload = json.loads(path.read_text(encoding="utf-8"))
        relative = str(path.relative_to(REPO))
        found = 0
        for pointer, graph6, row in walk_graph6_rows(payload):
            found += 1
            ref = {
                "id": f"retained:{relative}#{pointer}",
                "kind": "retained_graph6",
                "artifact": relative,
                "json_pointer": pointer,
            }
            entry = by_graph6.setdefault(
                graph6,
                {"sources": [], "stored_polynomials": [], "stored_orders": []},
            )
            entry["sources"].append(ref)
            for key in ("poly", "independence_polynomial"):
                if isinstance(row.get(key), list):
                    entry["stored_polynomials"].append((ref["id"], [int(x) for x in row[key]]))
            if isinstance(row.get("n"), int):
                entry["stored_orders"].append((ref["id"], int(row["n"])))
        artifact_reports.append({"artifact": relative, "graph6_references": found})

    verified_polynomials = 0
    for graph6 in sorted(by_graph6):
        entry = by_graph6[graph6]
        n, adj = parse_graph6(graph6.encode("ascii"))
        poly = independence_poly(n, adj)
        for ref_id, stored_n in entry["stored_orders"]:
            if stored_n != n:
                raise AssertionError(f"{ref_id}: stored order {stored_n} != graph6 order {n}")
        for ref_id, stored_poly in entry["stored_polynomials"]:
            if stored_poly != poly:
                raise AssertionError(f"{ref_id}: stored polynomial disagrees with exact DP")
            verified_polynomials += 1
        index.add(
            adj=adj,
            poly=poly,
            sources=entry["sources"],
            source_group="retained_graph6",
            graph6=graph6,
        )

    return [
        *artifact_reports,
        {
            "unique_labelled_graph6": len(by_graph6),
            "stored_polynomials_verified": verified_polynomials,
        },
    ]


def adjacency_from_edge_list(n: int, edges: Iterable[Iterable[Any]]) -> list[list[int]]:
    adj: list[list[int]] = [[] for _ in range(n)]
    seen_edges: set[tuple[int, int]] = set()
    for raw_edge in edges:
        edge = list(raw_edge)
        if len(edge) != 2:
            raise ValueError(f"expected a two-vertex edge, got {edge!r}")
        u, v = int(edge[0]), int(edge[1])
        if not (0 <= u < n and 0 <= v < n) or u == v:
            raise ValueError(f"invalid order-{n} edge {(u, v)!r}")
        normalized = (min(u, v), max(u, v))
        if normalized in seen_edges:
            raise ValueError(f"duplicate edge {normalized!r}")
        seen_edges.add(normalized)
        adj[u].append(v)
        adj[v].append(u)
    validate_tree(adj)
    return adj


def ingest_explicit_edge_artifacts(index: LocalIndex) -> list[dict[str, Any]]:
    """Replay committed exact-polynomial cases whose trees are stored as edges."""
    reports = []
    for relative in EXPLICIT_EDGE_ARTIFACTS:
        path = REPO / relative
        payload = json.loads(path.read_text(encoding="utf-8"))
        cases = payload.get("cases")
        if not isinstance(cases, list):
            raise ValueError(f"{relative}: expected a cases list")

        lane_counts: Counter[str] = Counter()
        for idx, row in enumerate(cases):
            pointer = f"/cases/{idx}"
            label = str(row["label"])
            lane = str(row["lane"])
            tree = row["tree"]
            stored_polynomial = row["polynomial"]["coefficients"]
            n = int(tree["order"])
            edges = tree["edge_list"]
            if int(tree["edges"]) != len(edges):
                raise AssertionError(f"{relative}#{pointer}: stored edge count changed")
            adj = adjacency_from_edge_list(n, edges)
            poly = independence_poly(n, adj)
            if poly != [int(value) for value in stored_polynomial]:
                raise AssertionError(
                    f"{relative}#{pointer}: stored polynomial disagrees with exact DP"
                )
            ref = {
                "id": f"explicit-edge:{relative}#{pointer}",
                "kind": "literature_root_stress_case",
                "artifact": relative,
                "json_pointer": pointer,
                "label": label,
                "lane": lane,
            }
            index.add(
                adj=adj,
                poly=poly,
                sources=[ref],
                source_group="literature_root_stress",
            )
            lane_counts[lane] += 1

        reports.append(
            {
                "artifact": relative,
                "cases_replayed": len(cases),
                "stored_polynomials_verified": len(cases),
                "lane_counts": dict(sorted(lane_counts.items())),
            }
        )
    return reports


def ingest_evolutionary_artifacts(index: LocalIndex) -> dict[str, Any]:
    """Replay saved evolutionary trees and resolve polynomial-only breaker rows."""
    tree_reports = []
    by_graph6: dict[str, dict[str, Any]] = {}
    for path in sorted((REPO / "results").glob(EVOLUTIONARY_TREE_GLOB)):
        relative = str(path.relative_to(REPO))
        payload = json.loads(path.read_text(encoding="utf-8"))
        row = payload if isinstance(payload.get("adj"), list) else payload.get("best")
        if not isinstance(row, dict) or not isinstance(row.get("adj"), list):
            raise ValueError(f"{relative}: no saved best adjacency list")
        adj = [[int(v) for v in neighbors] for neighbors in row["adj"]]
        n = int(row.get("n", len(adj)))
        if n != len(adj):
            raise AssertionError(f"{relative}: stored order disagrees with adjacency")
        validate_tree(adj)
        poly = independence_poly(n, adj)
        if poly != [int(value) for value in row["poly"]]:
            raise AssertionError(f"{relative}: stored polynomial disagrees with exact DP")
        graph6 = graph6_from_adj(adj)
        ref = {
            "id": f"evolutionary:{relative}",
            "kind": "saved_evolutionary_tree",
            "artifact": relative,
            "run_class": "smoke" if "smoke" in path.name else "production",
        }
        entry = by_graph6.setdefault(graph6, {"adj": adj, "poly": poly, "sources": []})
        if entry["poly"] != poly:
            raise AssertionError(f"{relative}: graph identity has inconsistent polynomial")
        entry["sources"].append(ref)
        tree_reports.append(
            {
                "artifact": relative,
                "order": n,
                "polynomial_verified": True,
                "graph6": graph6,
            }
        )

    for graph6 in sorted(by_graph6):
        entry = by_graph6[graph6]
        index.add(
            adj=entry["adj"],
            poly=entry["poly"],
            sources=entry["sources"],
            source_group="evolutionary_tree",
            graph6=graph6,
        )

    breaker_reports = []
    for path in sorted((REPO / "results").glob(LC_BREAKER_GLOB)):
        relative = str(path.relative_to(REPO))
        payload = json.loads(path.read_text(encoding="utf-8"))
        poly = [int(value) for value in payload["best"]["poly"]]
        if poly[0] != 1 or poly[1] != int(payload["n"]):
            raise AssertionError(f"{relative}: invalid polynomial header")
        matches = index.polynomials.get(tuple(poly), [])
        if not matches:
            raise AssertionError(
                f"{relative}: polynomial-only breaker is absent from its retained analysis"
            )
        breaker_reports.append(
            {
                "artifact": relative,
                "order": int(payload["n"]),
                "polynomial_sha256": polynomial_sha256(poly),
                "resolved_to_index": compact_matches(matches),
            }
        )

    return {
        "saved_tree_artifacts": tree_reports,
        "unique_saved_trees_indexed": len(by_graph6),
        "polynomial_only_breakers": breaker_reports,
    }


def exact_polynomial_quotient(
    dividend: list[int], divisor: list[int]
) -> list[int] | None:
    """Return the exact integer quotient when divisor (constant 1) divides dividend."""
    if not divisor or divisor[0] != 1 or len(divisor) > len(dividend):
        return None
    quotient_degree = len(dividend) - len(divisor)
    quotient: list[int] = []
    for degree in range(quotient_degree + 1):
        correction = sum(
            divisor[idx] * quotient[degree - idx]
            for idx in range(1, min(degree, len(divisor) - 1) + 1)
        )
        quotient.append(dividend[degree] - correction)
    product = [0] * (len(divisor) + len(quotient) - 1)
    for left_degree, left_value in enumerate(divisor):
        for right_degree, right_value in enumerate(quotient):
            product[left_degree + right_degree] += left_value * right_value
    return quotient if product == dividend else None


def path_polynomial(order: int) -> list[int]:
    adj: list[list[int]] = [[] for _ in range(order)]
    for vertex in range(order - 1):
        adj[vertex].append(vertex + 1)
        adj[vertex + 1].append(vertex)
    return independence_poly(order, adj)


def multiply_polynomials(left: list[int], right: list[int]) -> list[int]:
    product = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            product[left_degree + right_degree] += left_value * right_value
    return product


def closure_analysis(poly: list[int], index: LocalIndex) -> dict[str, Any]:
    """Check the V2 pair/triple, power, and one-to-three-path closure surface."""
    n = poly[1]
    indexed = []
    for factor_key, factor_sources in index.polynomials.items():
        factor = list(factor_key)
        factor_order = factor[1]
        if not (0 < factor_order < n):
            continue
        indexed.append(
            {
                "poly": factor,
                "order": factor_order,
                "sha256": polynomial_sha256(factor),
                "sources": factor_sources,
            }
        )

    pair_matches: dict[tuple[str, str], dict[str, Any]] = {}
    triple_matches: dict[tuple[str, str, str], dict[str, Any]] = {}
    for factor_row in indexed:
        factor = factor_row["poly"]
        quotient = exact_polynomial_quotient(poly, factor)
        if quotient is None:
            continue
        quotient_sources = index.polynomials.get(tuple(quotient))
        if quotient_sources:
            quotient_hash = polynomial_sha256(quotient)
            factors = sorted(
                [
                    {
                        "order": factor_row["order"],
                        "polynomial_sha256": factor_row["sha256"],
                        "sources": compact_matches(factor_row["sources"]),
                    },
                    {
                        "order": quotient[1],
                        "polynomial_sha256": quotient_hash,
                        "sources": compact_matches(quotient_sources),
                    },
                ],
                key=lambda item: item["polynomial_sha256"],
            )
            key = tuple(item["polynomial_sha256"] for item in factors)
            pair_matches[key] = {"factor_count": 2, "factors": factors}

        quotient_order = quotient[1]
        for second_row in indexed:
            if second_row["order"] >= quotient_order:
                continue
            residual = exact_polynomial_quotient(quotient, second_row["poly"])
            if residual is None:
                continue
            residual_sources = index.polynomials.get(tuple(residual))
            if not residual_sources:
                continue
            factors = sorted(
                [
                    {
                        "order": factor_row["order"],
                        "polynomial_sha256": factor_row["sha256"],
                        "sources": compact_matches(factor_row["sources"]),
                    },
                    {
                        "order": second_row["order"],
                        "polynomial_sha256": second_row["sha256"],
                        "sources": compact_matches(second_row["sources"]),
                    },
                    {
                        "order": residual[1],
                        "polynomial_sha256": polynomial_sha256(residual),
                        "sources": compact_matches(residual_sources),
                    },
                ],
                key=lambda item: item["polynomial_sha256"],
            )
            key = tuple(item["polynomial_sha256"] for item in factors)
            triple_matches[key] = {"factor_count": 3, "factors": factors}

    power_matches: dict[tuple[str, int], dict[str, Any]] = {}
    for factor_row in indexed:
        factor = factor_row["poly"]
        power = factor
        for exponent in range(2, 21):
            if factor_row["order"] * exponent > n:
                break
            power = multiply_polynomials(power, factor)
            if power == poly:
                key = (factor_row["sha256"], exponent)
                power_matches[key] = {
                    "base_order": factor_row["order"],
                    "base_polynomial_sha256": factor_row["sha256"],
                    "base_sources": compact_matches(factor_row["sources"]),
                    "exponent": exponent,
                }

    paths = [(order, path_polynomial(order)) for order in range(1, 17)]
    single_path_factors = []
    for order, path_poly in paths:
        quotient = exact_polynomial_quotient(poly, path_poly)
        if quotient is None:
            continue
        single_path_factors.append(
            {
                "path_order": order,
                "residual_order": quotient[1] if len(quotient) > 1 else 0,
                "residual_polynomial_sha256": polynomial_sha256(quotient),
                "residual_in_local_index": tuple(quotient) in index.polynomials,
                "residual_sources": compact_matches(index.polynomials.get(tuple(quotient), [])),
            }
        )

    path_product_matches: dict[tuple[int, ...], dict[str, Any]] = {}
    for path_count in range(1, 4):
        for combo in itertools.combinations_with_replacement(paths, path_count):
            path_orders = tuple(item[0] for item in combo)
            path_product = [1]
            for _, path_poly in combo:
                path_product = multiply_polynomials(path_product, path_poly)
            quotient = exact_polynomial_quotient(poly, path_product)
            if quotient is None:
                continue
            residual_sources = index.polynomials.get(tuple(quotient), [])
            if residual_sources:
                path_product_matches[path_orders] = {
                    "path_orders": list(path_orders),
                    "residual_order": quotient[1],
                    "residual_polynomial_sha256": polynomial_sha256(quotient),
                    "residual_sources": compact_matches(residual_sources),
                }

    ordered_products = [pair_matches[key] for key in sorted(pair_matches)]
    ordered_products.extend(triple_matches[key] for key in sorted(triple_matches))
    ordered_powers = [power_matches[key] for key in sorted(power_matches)]
    ordered_path_products = [
        path_product_matches[key] for key in sorted(path_product_matches)
    ]
    return {
        "known_indexed_product": {
            "count": len(ordered_products),
            "matches": ordered_products[:MATCH_LIMIT],
            "truncated": len(ordered_products) > MATCH_LIMIT,
        },
        "known_indexed_power": {
            "count": len(ordered_powers),
            "matches": ordered_powers[:MATCH_LIMIT],
            "truncated": len(ordered_powers) > MATCH_LIMIT,
        },
        "known_indexed_component_times_up_to_three_paths": {
            "count": len(ordered_path_products),
            "matches": ordered_path_products[:MATCH_LIMIT],
            "truncated": len(ordered_path_products) > MATCH_LIMIT,
        },
        "single_path_factors": single_path_factors,
        "known_closure": bool(
            ordered_products or ordered_powers or ordered_path_products
        ),
    }


def pressure_reference(
    index: LocalIndex, *, order: int, at_most: bool
) -> tuple[Fraction | None, list[dict[str, Any]], list[int]]:
    eligible_orders = [
        candidate_order
        for candidate_order in index.max_near_miss_by_order
        if (
            candidate_order <= order
            if at_most
            else candidate_order == order
        )
    ]
    if not eligible_orders:
        return None, [], []
    best = max(index.max_near_miss_by_order[item] for item in eligible_orders)
    best_orders = [
        item for item in eligible_orders if index.max_near_miss_by_order[item] == best
    ]
    sources = [
        source
        for item in best_orders
        for source in index.max_near_miss_sources_by_order[item]
    ]
    return best, sources, sorted(best_orders)


def pressure_comparison_payload(
    candidate: Fraction,
    reference: Fraction | None,
    sources: list[dict[str, Any]],
    orders: list[int],
) -> dict[str, Any]:
    if reference is None:
        return {
            "available": False,
            "candidate": rational_payload(candidate),
            "reference": None,
            "reference_orders": [],
            "reference_sources": compact_matches([]),
            "beats_reference": None,
            "ties_reference": None,
            "candidate_minus_reference": None,
        }
    return {
        "available": True,
        "candidate": rational_payload(candidate),
        "reference": rational_payload(reference),
        "reference_orders": orders,
        "reference_sources": compact_matches(sources),
        "beats_reference": candidate > reference,
        "ties_reference": candidate == reference,
        "candidate_minus_reference": rational_payload(candidate - reference),
    }


def verify_upstream_summary_arithmetic(payload: dict[str, Any]) -> dict[str, Any]:
    core = 80
    paths = 16
    seed_count = int(payload["distinct_non_log_concave_seeds"])
    pair_count = math.comb(core + 1, 2)
    triple_count = math.comb(core + 2, 3)
    power_count = seed_count * 19
    path_count = core * (math.comb(paths, 1) + math.comb(paths + 1, 2) + math.comb(paths + 2, 3))
    computed = pair_count + triple_count + power_count + path_count
    claimed = int(payload["forest_objects_tested"])
    return {
        "assumptions": {
            "core_seeds": core,
            "path_orders": [1, paths],
            "powers": [2, 20],
            "path_factor_counts": [1, 3],
        },
        "pair_products": pair_count,
        "triple_products": triple_count,
        "proper_powers": power_count,
        "seed_times_path_products": path_count,
        "computed_total": computed,
        "claimed_total": claimed,
        "verified": computed == claimed,
    }


def run_upstream_bounded_verifier(blair_dir: Path) -> dict[str, Any]:
    command = [sys.executable, "verify_993_result.py"]
    proc = subprocess.run(
        command,
        cwd=blair_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    return {
        "command": " ".join(command),
        "returncode": proc.returncode,
        "passed": proc.returncode == 0,
        "stdout": proc.stdout,
        "scope": "bounded verifier; does not rerun the 112,916-tree V2 enumeration",
    }


def verify_seed_files(blair_dir: Path) -> list[dict[str, Any]]:
    reports = []
    expected_families = {
        "seed_T_3_4_4.json": {
            "family": "li",
            "params": {"m": 4, "n": 4},
            "adj": make_li_tree(4, 4, starred=False),
        },
        "seed_Tstar_3_3_4.json": {
            "family": "li-star",
            "params": {"m": 3, "n": 4},
            "adj": make_li_tree(3, 4, starred=True),
        },
    }
    for name in ("seed_T_3_4_4.json", "seed_Tstar_3_3_4.json"):
        path = blair_dir / name
        payload = json.loads(path.read_text(encoding="utf-8"))
        n = int(payload["order"])
        adj: list[list[int]] = [[] for _ in range(n)]
        for raw_u, raw_v in payload["edges"]:
            u, v = int(raw_u), int(raw_v)
            adj[u].append(v)
            adj[v].append(u)
        validate_tree(adj)
        poly = independence_poly(n, adj)
        if poly != [int(value) for value in payload["independence_polynomial"]]:
            raise AssertionError(f"{name}: stored polynomial disagrees with exact DP")
        expected = expected_families[name]
        known_family_match = nx.is_isomorphic(
            adjacency_to_networkx(adj),
            adjacency_to_networkx(expected["adj"]),
        )
        if not known_family_match:
            raise AssertionError(f"{name}: expected local known-family match failed")
        reports.append(
            {
                "artifact": name,
                "name": payload["name"],
                "order": n,
                "edge_count": len(payload["edges"]),
                "polynomial_sha256": polynomial_sha256(poly),
                "verified": True,
                "known_family_collision": {
                    "family": expected["family"],
                    "params": expected["params"],
                    "exact_unrooted_tree": known_family_match,
                    "exact_independence_polynomial": (
                        poly == independence_poly(len(expected["adj"]), expected["adj"])
                    ),
                },
            }
        )
    return reports


def classify_ancestry(tokens: list[BushToken], adj: list[list[int]]) -> dict[str, Any]:
    ancestries: list[dict[str, Any]] = []
    if (
        all(token.kind == "u" and token.pendant_length == 1 for token in tokens)
        and len({token.child_count for token in tokens}) == 1
    ):
        m = len(tokens)
        t = tokens[0].child_count
        expected = make_galvin_tree(m, t)
        if not nx.is_isomorphic(adjacency_to_networkx(adj), adjacency_to_networkx(expected)):
            raise AssertionError("syntactic Galvin classification failed isomorphism check")
        ancestries.append({"family": "galvin", "params": {"m": m, "t": t}})

    if len(tokens) == 3 and all(
        token.kind == "u" and token.pendant_length == 1 for token in tokens
    ):
        counts = [token.child_count for token in tokens]
        if 3 in counts:
            remainder = list(counts)
            remainder.remove(3)
            m, n = sorted(remainder)
            expected = make_li_tree(m, n, starred=False)
            if nx.is_isomorphic(adjacency_to_networkx(adj), adjacency_to_networkx(expected)):
                ancestries.append({"family": "li", "params": {"m": m, "n": n}})

    has_single_extension = any(token.kind == "s" for token in tokens)
    if ancestries:
        structural_class = "known_family_member"
    elif has_single_extension:
        structural_class = "mixed_spider_bouquet_with_single_leg_extension"
    else:
        structural_class = "mixed_uniform_spider_bouquet"
    return {
        "local_grammar": "scripts.valley_search spider-bouquet",
        "local_grammar_relation": "exact_special_case",
        "structural_grammar_novel": False,
        "structural_class": structural_class,
        "known_family_ancestry": ancestries,
        "contains_blair_single_extension_token": has_single_extension,
        "single_extension_is_actual_li_starred_family": False,
    }


def canonical_tree_matches(
    remote_rows: list[dict[str, Any]],
    candidates: list[LocalCandidate],
) -> tuple[str, str | None, dict[int, list[dict[str, Any]]]]:
    labelg = shutil.which("labelg")
    matches: dict[int, list[dict[str, Any]]] = defaultdict(list)
    if labelg is not None:
        all_graph6 = [row["graph6"] for row in remote_rows]
        all_graph6.extend(candidate.graph6 for candidate in candidates)
        unique_graph6 = sorted(set(all_graph6))
        canonical = canonical_graph6_batch(unique_graph6)
        canonical_by_raw = dict(zip(unique_graph6, canonical, strict=True))
        local_by_canonical: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for candidate in candidates:
            local_by_canonical[canonical_by_raw[candidate.graph6]].extend(candidate.sources)
        for idx, row in enumerate(remote_rows):
            row["canonical_graph6"] = canonical_by_raw[row["graph6"]]
            matches[idx].extend(local_by_canonical[row["canonical_graph6"]])
        return "nauty_labelg_canonical_graph6", labelg, matches

    local_by_order_hash: dict[tuple[int, str], list[LocalCandidate]] = defaultdict(list)
    for candidate in candidates:
        graph = adjacency_to_networkx(candidate.adj)
        wl_hash = nx.weisfeiler_lehman_graph_hash(graph)
        local_by_order_hash[(candidate.n, wl_hash)].append(candidate)
    for idx, row in enumerate(remote_rows):
        graph = adjacency_to_networkx(row["adj"])
        wl_hash = nx.weisfeiler_lehman_graph_hash(graph)
        for candidate in local_by_order_hash[(row["n"], wl_hash)]:
            if nx.is_isomorphic(graph, adjacency_to_networkx(candidate.adj)):
                matches[idx].extend(candidate.sources)
        row["canonical_graph6"] = None
    return "wl_hash_prefilter_plus_exact_networkx_isomorphism", None, matches


def build_remote_rows(
    blair_payload: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    public_rows = blair_payload.get("top25_most_severe_non_log_concave")
    if not isinstance(public_rows, list) or len(public_rows) != 25:
        raise ValueError("expected exactly 25 public top witnesses in Blair results.json")

    internal_rows = []
    validation_reports = []
    for rank, stored in enumerate(public_rows, start=1):
        label = str(stored["label"])
        tokens, label_n = parse_bush_label(label)
        adj, poly = reconstruct_bush(tokens)
        n = len(adj)
        stored_poly = [int(value) for value in stored["independence_polynomial"]]
        if n != label_n or n != int(stored["order"]):
            raise AssertionError(f"{label}: computed/label/stored orders disagree")
        if poly != stored_poly:
            raise AssertionError(f"{label}: reconstructed polynomial disagrees with results.json")
        metrics = exact_metrics(poly)
        stored_severity = float(stored["lc_defect"])
        exact_severity = fraction_from_payload(metrics["normalized_lc_severity"])
        if round(float(exact_severity), 4) != stored_severity:
            raise AssertionError(f"{label}: stored rounded LC severity disagrees with exact value")
        if metrics["unimodal"] is not True or metrics["log_concave"] is not False:
            raise AssertionError(f"{label}: public top witness classification changed")
        graph6 = graph6_from_adj(adj)
        internal_rows.append(
            {
                "rank": rank,
                "label": label,
                "n": n,
                "tokens": tokens,
                "adj": adj,
                "poly": poly,
                "metrics": metrics,
                "graph6": graph6,
                "stored_lc_defect_rounded": stored_severity,
            }
        )
        validation_reports.append(
            {
                "rank": rank,
                "label": label,
                "order_verified": True,
                "polynomial_verified": True,
                "rounded_lc_severity_verified": True,
            }
        )
    return internal_rows, validation_reports


def main() -> int:
    args = parse_args()
    blair_dir = resolve_blair_dir(args.blair_repo)
    blair_payload = json.loads((blair_dir / "results.json").read_text(encoding="utf-8"))
    remote_rows, remote_validation = build_remote_rows(blair_payload)
    remote_orders = {row["n"] for row in remote_rows}

    local_index = LocalIndex(remote_orders)
    family_reports = reconstruct_family_grids(local_index)
    retained_reports = ingest_retained_graph6(local_index)
    explicit_edge_reports = ingest_explicit_edge_artifacts(local_index)
    evolutionary_report = ingest_evolutionary_artifacts(local_index)
    canonical_method, labelg_path, tree_matches = canonical_tree_matches(
        remote_rows,
        local_index.tree_candidates,
    )

    witnesses = []
    for idx, row in enumerate(remote_rows):
        poly_matches = local_index.polynomials.get(tuple(row["poly"]), [])
        sig_key = exact_signature(row["n"], row["metrics"])
        signature_matches = local_index.signatures.get(sig_key, [])
        exact_tree = compact_matches(tree_matches[idx])
        exact_poly = compact_matches(poly_matches)
        exact_sig = compact_matches(signature_matches)
        near = fraction_from_payload(row["metrics"]["post_descent_max_ratio"])
        global_comparison = pressure_comparison_payload(
            near,
            local_index.max_near_miss,
            local_index.max_near_miss_sources,
            sorted(
                order
                for order, value in local_index.max_near_miss_by_order.items()
                if value == local_index.max_near_miss
            ),
        )
        same_best, same_sources, same_orders = pressure_reference(
            local_index, order=row["n"], at_most=False
        )
        same_order_comparison = pressure_comparison_payload(
            near, same_best, same_sources, same_orders
        )
        pareto_best, pareto_sources, pareto_orders = pressure_reference(
            local_index, order=row["n"], at_most=True
        )
        pareto_comparison = pressure_comparison_payload(
            near, pareto_best, pareto_sources, pareto_orders
        )
        closures = closure_analysis(row["poly"], local_index)
        classification = classify_ancestry(row["tokens"], row["adj"])
        if exact_tree["count"] > 0:
            comparison_class = "exact duplicate"
        elif exact_poly["count"] > 0:
            comparison_class = "polynomial duplicate"
        elif closures["known_closure"]:
            comparison_class = "known closure"
        else:
            comparison_class = "known grammar/new parameter"

        if exact_tree["count"] > 0 or exact_poly["count"] > 0:
            scientific_reason = (
                "The exact tree or polynomial is already present in the bounded local index; "
                "the grammar is also a local spider-bouquet special case."
            )
        elif closures["known_closure"]:
            scientific_reason = (
                "The coefficient vector is an exact product, power, or bounded path-product "
                "closure of indexed components; the tree grammar is also already local."
            )
        elif pareto_comparison["beats_reference"]:
            scientific_reason = (
                "The witness is a new parameter point in the existing spider-bouquet grammar "
                "and extends the bounded order-versus-pressure Pareto frontier; this is not a "
                "new construction mechanism."
            )
        elif same_order_comparison["beats_reference"]:
            scientific_reason = (
                "The witness improves the available same-order pressure comparator but not the "
                "order-at-most Pareto frontier; its spider-bouquet grammar is already local."
            )
        else:
            scientific_reason = (
                "The exact tree and polynomial are absent from the bounded local index, but "
                "the grammar is already a local spider-bouquet special case and the witness "
                "does not extend its global, same-order, or order-at-most pressure frontiers."
            )
        classification.update(
            {
                "comparison_class": comparison_class,
                "exact_tree_novel_relative_to_index": exact_tree["count"] == 0,
                "exact_polynomial_novel_relative_to_index": exact_poly["count"] == 0,
                "exact_signature_novel_relative_to_index": exact_sig["count"] == 0,
                "candidate_exact_data_novel_relative_to_index": (
                    exact_tree["count"] == 0 and exact_poly["count"] == 0
                ),
                "beats_global_local_post_descent_ratio": global_comparison["beats_reference"],
                "beats_same_order_local_post_descent_ratio": same_order_comparison[
                    "beats_reference"
                ],
                "extends_order_pressure_pareto_frontier": pareto_comparison[
                    "beats_reference"
                ],
                "scientifically_new_hard_family": False,
                "scientific_novelty_reason": scientific_reason,
            }
        )
        witnesses.append(
            {
                "rank": row["rank"],
                "label": row["label"],
                "order": row["n"],
                "branch_tokens": [token.payload() for token in row["tokens"]],
                "normalized_local_bouquet_spec": [list(token.legs) for token in row["tokens"]],
                "graph6": row["graph6"],
                "canonical_graph6": row["canonical_graph6"],
                "independence_polynomial": row["poly"],
                "polynomial_sha256": polynomial_sha256(row["poly"]),
                "stored_lc_defect_rounded": row["stored_lc_defect_rounded"],
                "exact_metrics": row["metrics"],
                "exact_signature": signature_payload(row["n"], row["metrics"]),
                "matches": {
                    "exact_unrooted_tree": exact_tree,
                    "exact_independence_polynomial": exact_poly,
                    "exact_metric_signature": exact_sig,
                    "polynomial_closures": closures,
                },
                "pressure_comparisons": {
                    "global_size_unrestricted": global_comparison,
                    "same_order": same_order_comparison,
                    "order_at_most_pareto": pareto_comparison,
                },
                "classification": classification,
            }
        )

    exact_tree_collisions = sum(
        witness["matches"]["exact_unrooted_tree"]["count"] > 0 for witness in witnesses
    )
    polynomial_collisions = sum(
        witness["matches"]["exact_independence_polynomial"]["count"] > 0
        for witness in witnesses
    )
    signature_collisions = sum(
        witness["matches"]["exact_metric_signature"]["count"] > 0 for witness in witnesses
    )
    known_family_members = sum(
        bool(witness["classification"]["known_family_ancestry"]) for witness in witnesses
    )
    known_closure_witnesses = sum(
        witness["matches"]["polynomial_closures"]["known_closure"]
        for witness in witnesses
    )
    single_path_factor_witnesses = sum(
        bool(witness["matches"]["polynomial_closures"]["single_path_factors"])
        for witness in witnesses
    )
    global_improvements = sum(
        witness["classification"]["beats_global_local_post_descent_ratio"]
        for witness in witnesses
    )
    same_order_improvements = sum(
        witness["classification"]["beats_same_order_local_post_descent_ratio"]
        for witness in witnesses
    )
    pareto_improvements = sum(
        witness["classification"]["extends_order_pressure_pareto_frontier"]
        for witness in witnesses
    )
    max_remote = max(
        (
            (
                fraction_from_payload(witness["exact_metrics"]["post_descent_max_ratio"]),
                witness,
            )
            for witness in witnesses
        ),
        key=lambda item: item[0],
    )

    report = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "purpose": "Bounded exact comparison; no new bouquet enumeration.",
        "reproduction_command": (
            "python3 scripts/compare_blair_rooted_bush.py "
            f"--blair-repo {args.blair_repo} --output {args.output}"
        ),
        "upstream": {
            "public_url": PUBLIC_URL,
            "local_path": str(blair_dir),
            "git_revision": git_revision(blair_dir),
            "results_sha256": hashlib.sha256((blair_dir / "results.json").read_bytes()).hexdigest(),
            "claimed_single_trees_scanned": blair_payload.get("single_trees_scanned"),
            "claimed_distinct_non_log_concave_seeds": blair_payload.get(
                "distinct_non_log_concave_seeds"
            ),
            "claimed_forest_objects_tested": blair_payload.get("forest_objects_tested"),
            "public_witnesses_materialized": len(remote_rows),
            "claimed_library_fully_materialized": False,
            "upstream_v2_executed_by_this_comparison": False,
            "forest_inventory_arithmetic": verify_upstream_summary_arithmetic(blair_payload),
            "bounded_verifier": run_upstream_bounded_verifier(blair_dir),
            "seed_files": verify_seed_files(blair_dir),
        },
        "method": {
            "tree_identity": canonical_method,
            "labelg_path": labelg_path,
            "polynomial_identity": "exact integer coefficient tuple",
            "signature_identity": (
                "exact tuple of order, alpha, mode interval, first descent, LC indices, "
                "and reduced post-descent maximum ratio"
            ),
            "remote_reconstruction": (
                "Blair u(c,pL) -> c local legs L+1; Blair s(c,pL) -> one leg L+2 "
                "and c-1 legs L+1; generic tree DP cross-checked against local bouquet formula"
            ),
            "tree_fallback_limitation": (
                None
                if labelg_path is not None
                else "WL hashes are only prefilters; exact NetworkX isomorphism decides matches, "
                "but no portable canonical graph6 certificate is emitted."
            ),
            "match_lists_capped_at": MATCH_LIMIT,
            "closure_scope": (
                "all pairs and triples of indexed components, direct indexed-component "
                "powers 2..20, and one to three path factors P1..P16 against an indexed residual"
            ),
        },
        "validation": {
            "remote_top25": remote_validation,
            "family_grids": family_reports,
            "retained_graph6": retained_reports,
            "explicit_edge_artifacts": explicit_edge_reports,
            "evolutionary_artifacts": evolutionary_report,
        },
        "local_index": {
            "principal_family_grid_records": sum(item["reconstructed"] for item in family_reports),
            "total_tree_records_indexed": local_index.record_count,
            "source_references_indexed": local_index.source_reference_count,
            "source_record_counts": dict(sorted(local_index.source_counts.items())),
            "unique_exact_polynomials": len(local_index.polynomials),
            "unique_exact_metric_signatures": len(local_index.signatures),
            "tree_candidates_at_public_witness_orders": len(local_index.tree_candidates),
            "public_witness_orders": sorted(remote_orders),
            "max_post_descent_ratio": rational_payload(local_index.max_near_miss),
            "max_post_descent_ratio_sources": compact_matches(
                local_index.max_near_miss_sources
            ),
        },
        "summary": {
            "public_witness_count": len(witnesses),
            "exact_tree_collision_witnesses": exact_tree_collisions,
            "exact_polynomial_collision_witnesses": polynomial_collisions,
            "exact_signature_collision_witnesses": signature_collisions,
            "known_family_members": known_family_members,
            "known_product_power_or_path_closure_witnesses": known_closure_witnesses,
            "witnesses_with_single_path_polynomial_factor": single_path_factor_witnesses,
            "structurally_new_grammars": 0,
            "witnesses_beating_global_size_unrestricted_pressure": global_improvements,
            "witnesses_beating_same_order_pressure": same_order_improvements,
            "witnesses_extending_order_pressure_pareto_frontier": pareto_improvements,
            "best_public_post_descent_ratio": rational_payload(max_remote[0]),
            "best_public_post_descent_ratio_label": max_remote[1]["label"],
            "verdict": (
                "The public top 25 add exact parameter points but no new grammar or known "
                f"polynomial closure. They include {same_order_improvements} same-order and "
                f"{pareto_improvements} order-pressure Pareto improvements, so they should be "
                "ingested as bounded frontier data rather than treated as a new mechanism. "
                "Full 4,445-way collision claims require materializing the unpublished "
                "remainder of the V2 library."
            ),
        },
        "known_upstream_scope_caveats": [
            (
                "search_993_v2.py's s-token lengthens one leg by one edge, whereas the actual "
                "T* builder in search_993.py and the local Li-star generator lengthen it by two; "
                "the V2 grammar therefore does not literally contain the advertised T* family."
            ),
            (
                "search_993_v3_wide.py widens the parameter ranges but has no matching committed "
                "results artifact; the 112,916/4,445 claims belong to V2."
            ),
        ],
        "limitations": [
            "Only 25 of the claimed 4,445 distinct non-log-concave polynomials are public records.",
            (
                "The 253,695 forest-object total is arithmetically reproduced, but the "
                "112,916-tree scan and 4,445-library cardinality were not regenerated."
            ),
            (
                "The raw 50,000-row Ramos-Sun Prüfer corpora are not retained locally; this index "
                "uses every committed ranked graph6 row from all eleven summaries."
            ),
            (
                "Exact-tree comparisons are unrooted. Normalized bouquet specs separately record "
                "rooted structural ancestry."
            ),
        ],
        "witnesses": witnesses,
        "environment": {
            "python": platform.python_version(),
            "networkx": nx.__version__,
            "platform": platform.platform(),
        },
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")
    print(
        f"top25: tree collisions={exact_tree_collisions}, "
        f"polynomial collisions={polynomial_collisions}, "
        f"known-family members={known_family_members}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
