#!/usr/bin/env python3
"""Analyze adversarial tree corpora and structured non-LC families.

The first use case is the Ramos--Sun PatternBoost output, whose
``search_output_*.txt`` files contain one 1-based Prüfer code per line.
The script also generates the Galvin and Bautista--Ramos non-log-concave
tree families.  In all modes it recomputes independence polynomials with
this repository's exact DP and writes a compact JSON summary.
"""

from __future__ import annotations

import argparse
import ast
import json
import math
import os
import statistics
import subprocess
import sys
import time
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from indpoly import (  # noqa: E402
    independence_poly,
    is_log_concave,
    is_unimodal,
    log_concavity_ratio,
    near_miss_ratio,
)


@dataclass(frozen=True)
class PruferRecord:
    line_no: int
    code: tuple[int, ...]


@dataclass(frozen=True)
class GeneratedRecord:
    rank: int
    family: str
    params: dict[str, int]
    adj: list[list[int]]
    default_target_k: int | None
    expected_lc_breaks: tuple[int, ...]


def parse_prufer_line(line: str) -> tuple[int, ...]:
    """Parse a 1-based Prüfer code from Ramos--Sun text formats."""
    line = line.strip()
    if not line:
        raise ValueError("empty line")
    if line.startswith("["):
        values = ast.literal_eval(line)
        return tuple(int(x) for x in values)
    return tuple(int(token.strip().lstrip("V")) for token in line.split(",") if token.strip())


def iter_prufer_records(path: Path, limit: int | None) -> list[PruferRecord]:
    records: list[PruferRecord] = []
    with path.open("r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, start=1):
            if not line.strip():
                continue
            records.append(PruferRecord(line_no=line_no, code=parse_prufer_line(line)))
            if limit is not None and len(records) >= limit:
                break
    return records


def parse_score_distribution(path: Path, limit: int | None) -> list[float]:
    """Expand ``Score: x, Count: y`` rows into rank-aligned score values."""
    scores: list[float] = []
    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            raw = raw.strip()
            if not raw:
                continue
            # Format: Score: 9.919602e8, Count: 1
            parts = raw.replace("Score:", "").replace("Count:", "").split(",")
            if len(parts) != 2:
                raise ValueError(f"unexpected score row: {raw!r}")
            score = float(parts[0].strip())
            count = int(parts[1].strip())
            for _ in range(count):
                scores.append(score)
                if limit is not None and len(scores) >= limit:
                    return scores
    return scores


def prufer_to_adj(code_1_based: tuple[int, ...]) -> list[list[int]]:
    """Convert a 1-based Prüfer code to a 0-based adjacency list."""
    code = [x - 1 for x in code_1_based]
    n = len(code) + 2
    if any(v < 0 or v >= n for v in code):
        raise ValueError(f"Prüfer code has value outside 1..{n}")

    degree = [1] * n
    for v in code:
        degree[v] += 1

    # n is small in the intended corpora; a sorted leaf list keeps this simple.
    leaves = sorted(i for i, d in enumerate(degree) if d == 1)
    adj: list[list[int]] = [[] for _ in range(n)]

    for v in code:
        leaf = leaves.pop(0)
        adj[leaf].append(v)
        adj[v].append(leaf)
        degree[leaf] -= 1
        degree[v] -= 1
        if degree[v] == 1:
            insert_at = 0
            while insert_at < len(leaves) and leaves[insert_at] < v:
                insert_at += 1
            leaves.insert(insert_at, v)

    if len(leaves) != 2:
        raise ValueError("invalid Prüfer code did not leave two leaves")
    u, v = leaves
    adj[u].append(v)
    adj[v].append(u)
    for neighbors in adj:
        neighbors.sort()
    return adj


def _add_edge(adj: list[list[int]], u: int, v: int) -> None:
    adj[u].append(v)
    adj[v].append(u)


def make_galvin_tree(m: int, t: int) -> list[list[int]]:
    """Return Galvin's rooted tree T_{m,t,1}.

    The root has m children; each child has t pendant P2 paths.  The tree
    has 1 + m + 2mt vertices, independence number m(t + 1), and its
    asymptotic LC break is at coefficient index mt + 2.
    """
    if m < 0 or t < 0:
        raise ValueError("m and t must be nonnegative")
    n = 1 + m + 2 * m * t
    adj: list[list[int]] = [[] for _ in range(n)]
    for i in range(m):
        w = 1 + i
        _add_edge(adj, 0, w)
        for j in range(t):
            x = 1 + m + i * t + j
            y = 1 + m + m * t + i * t + j
            _add_edge(adj, w, x)
            _add_edge(adj, x, y)
    return adj


def make_bautista_ramos_tree(m: int, t: int) -> list[list[int]]:
    """Return Bautista--Ramos' tree TG_{m,t}.

    TG_{m,t} is formed from m disjoint copies of T_{3,t}, joining each
    copy root to a new root v0, and adding one extra leaf adjacent to v0.
    """
    if m < 0 or t < 0:
        raise ValueError("m and t must be nonnegative")
    adj: list[list[int]] = [[], []]
    _add_edge(adj, 0, 1)
    for _ in range(m):
        copy = make_galvin_tree(3, t)
        offset = len(adj)
        adj.extend([] for _ in range(len(copy)))
        for u, neighbors in enumerate(copy):
            for v in neighbors:
                if u < v:
                    _add_edge(adj, offset + u, offset + v)
        _add_edge(adj, 0, offset)
    return adj


def make_li_tree(m: int, n: int, *, starred: bool = False) -> list[list[int]]:
    """Return Li's Kadrawi--Levit family tree T_{3,m,n} or T^*_{3,m,n}.

    The root has three children.  Their numbers of pendant P2 paths are
    3, m, and n.  In the starred family, the edge from the third P2 path
    in the first branch to its leaf is replaced by a path on four vertices.
    """
    if m < 0 or n < 0:
        raise ValueError("m and n must be nonnegative")
    branch_sizes = [3, m, n]
    adj: list[list[int]] = [[]]
    branch_roots: list[int] = []
    for size in branch_sizes:
        branch_root = len(adj)
        adj.append([])
        _add_edge(adj, 0, branch_root)
        branch_roots.append(branch_root)
        for child_idx in range(size):
            child = len(adj)
            leaf = child + 1
            adj.extend(([], []))
            _add_edge(adj, branch_root, child)
            if starred and branch_root == branch_roots[0] and child_idx == 2:
                mid = len(adj)
                final_leaf = mid + 1
                adj.extend(([], []))
                _add_edge(adj, child, leaf)
                _add_edge(adj, leaf, mid)
                _add_edge(adj, mid, final_leaf)
            else:
                _add_edge(adj, child, leaf)
    return adj


def parse_int_values(spec: str) -> list[int]:
    """Parse comma-separated integers and inclusive ranges such as 1,3-5."""
    values: set[int] = set()
    for raw_part in spec.split(","):
        part = raw_part.strip()
        if not part:
            continue
        if "-" in part:
            start_s, end_s = part.split("-", 1)
            start = int(start_s)
            end = int(end_s)
            if start > end:
                raise ValueError(f"empty descending range {part!r}")
            values.update(range(start, end + 1))
        else:
            values.add(int(part))
    if not values:
        raise ValueError("no integer values parsed")
    if any(value < 0 for value in values):
        raise ValueError("family parameters must be nonnegative")
    return sorted(values)


def iter_family_records(
    *,
    family: str,
    m_values: list[int],
    t_values: list[int],
    max_n: int | None,
) -> list[GeneratedRecord]:
    families = ["galvin", "bautista-ramos", "li", "li-star"] if family == "all" else [family]
    records: list[GeneratedRecord] = []
    rank = 1
    for family_name in families:
        for m in m_values:
            for t in t_values:
                if family_name == "galvin":
                    adj = make_galvin_tree(m, t)
                    expected = (m * t + 2,) if m > 0 else ()
                    default_target_k = expected[0] if expected else None
                elif family_name == "bautista-ramos":
                    adj = make_bautista_ramos_tree(m, t)
                    alpha = 3 * (t + 1) * m + 1
                    expected = tuple(alpha - (2 * j + 1) for j in range(m))
                    default_target_k = expected[0] if expected else None
                elif family_name == "li":
                    adj = make_li_tree(m, t, starred=False)
                    expected = ()
                    default_target_k = None
                elif family_name == "li-star":
                    adj = make_li_tree(m, t, starred=True)
                    expected = ()
                    default_target_k = None
                else:
                    raise ValueError(f"unknown family {family_name!r}")
                if max_n is not None and len(adj) > max_n:
                    continue
                records.append(
                    GeneratedRecord(
                        rank=rank,
                        family=family_name,
                        params={"m": m, "n" if family_name in {"li", "li-star"} else "t": t},
                        adj=adj,
                        default_target_k=default_target_k,
                        expected_lc_breaks=expected,
                    )
                )
                rank += 1
    return records


def _graph6_order_header(n: int) -> str:
    """Return the standard graph6 order header."""
    if n < 0:
        raise ValueError("graph order must be nonnegative")
    if n <= 62:
        return chr(n + 63)
    if n <= 258047:
        return "~" + "".join(chr(((n >> shift) & 0x3F) + 63) for shift in (12, 6, 0))
    if n <= 68719476735:
        return "~~" + "".join(
            chr(((n >> shift) & 0x3F) + 63) for shift in (30, 24, 18, 12, 6, 0)
        )
    raise ValueError("graph6 only supports graph orders below 2^36")


def graph6_from_adj(adj: list[list[int]]) -> str:
    """Encode a simple graph in graph6 format."""
    n = len(adj)
    neighbor_sets = [set(neighbors) for neighbors in adj]
    chars = [_graph6_order_header(n)]
    value = 0
    used = 0
    for j in range(1, n):
        for i in range(j):
            value = (value << 1) | (1 if j in neighbor_sets[i] else 0)
            used += 1
            if used == 6:
                chars.append(chr(value + 63))
                value = 0
                used = 0
    if used:
        chars.append(chr((value << (6 - used)) + 63))
    return "".join(chars)


def canonical_graph6_batch(graph6_values: list[str]) -> list[str]:
    """Canonicalize graph6 strings with nauty labelg, preserving order."""
    if not graph6_values:
        return []
    proc = subprocess.run(
        ["labelg", "-q", "-g"],
        input=("\n".join(graph6_values) + "\n").encode("ascii"),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    out = [line.decode("ascii").strip() for line in proc.stdout.splitlines() if line.strip()]
    if len(out) != len(graph6_values):
        raise RuntimeError(
            f"labelg returned {len(out)} graph6 rows for {len(graph6_values)} inputs; "
            f"stderr={proc.stderr.decode('utf-8', 'replace')}"
        )
    return out


def lc_defects(poly: list[int]) -> list[dict[str, Any]]:
    defects = []
    for k in range(1, len(poly) - 1):
        defect = poly[k - 1] * poly[k + 1] - poly[k] * poly[k]
        if defect > 0:
            defects.append(
                {
                    "k": k,
                    "defect": int(defect),
                    "ratio": float(Fraction(poly[k - 1] * poly[k + 1], poly[k] * poly[k])),
                }
            )
    return defects


def mode_interval(poly: list[int]) -> tuple[int, int]:
    peak = max(poly)
    modes = [idx for idx, value in enumerate(poly) if value == peak]
    return modes[0], modes[-1]


def first_descent(poly: list[int]) -> int | None:
    for k in range(1, len(poly)):
        if poly[k] < poly[k - 1]:
            return k
    return None


def mean_independent_set_size(poly: list[int]) -> float:
    total = sum(poly)
    return sum(k * value for k, value in enumerate(poly)) / total


def structural_metrics(adj: list[list[int]]) -> dict[str, Any]:
    degrees = [len(neighbors) for neighbors in adj]
    leaves = [v for v, degree in enumerate(degrees) if degree == 1]
    leaf_set = set(leaves)
    leaf_children = [sum(1 for u in adj[v] if u in leaf_set) for v in range(len(adj))]
    degree2_adjacent_leaf = sum(
        1 for v, degree in enumerate(degrees) if degree == 2 and any(u in leaf_set for u in adj[v])
    )
    return {
        "leaves": len(leaves),
        "max_degree": max(degrees) if degrees else 0,
        "degree2": sum(1 for degree in degrees if degree == 2),
        "degree2_adjacent_leaf": degree2_adjacent_leaf,
        "max_d_leaf": max(leaf_children) if leaf_children else 0,
        "multi_leaf_hubs": sum(1 for value in leaf_children if value >= 2),
        "d_leaf_le1": all(value <= 1 for value in leaf_children),
    }


def evaluate_tree(
    adj: list[list[int]],
    *,
    rank: int,
    target_k: int | None,
    distribution_score: float | None,
    canonical_g6: str | None,
    metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    n = len(adj)
    poly = independence_poly(n, adj)
    if target_k is None:
        target_k = n // 2
    if target_k <= 0 or target_k >= len(poly) - 1:
        exact_score = None
    else:
        exact_score = poly[target_k - 1] * poly[target_k + 1] - poly[target_k] * poly[target_k]

    lc_ratio, lc_pos = log_concavity_ratio(poly)
    nm_ratio, nm_pos = near_miss_ratio(poly)
    mode_first, mode_last = mode_interval(poly)
    mu = mean_independent_set_size(poly)
    defects = lc_defects(poly)
    struct = structural_metrics(adj)
    graph6 = graph6_from_adj(adj)
    alpha = len(poly) - 1
    descent = first_descent(poly)
    expected_lc_breaks = tuple(metadata.get("expected_lc_breaks", ())) if metadata else ()
    actual_lc_breaks = {defect["k"] for defect in defects}
    low_mode_threshold = n // 3 + 1
    mean_gap_n_over_3 = n / 3 - mu
    defect_distances_from_mode = [defect["k"] - mode_first for defect in defects]
    defect_distances_from_near_miss = [
        defect["k"] - nm_pos for defect in defects if nm_pos != -1
    ]

    out: dict[str, Any] = {
        "rank": rank,
        "n": n,
        "alpha": alpha,
        "mode_first": mode_first,
        "mode_last": mode_last,
        "mu": mu,
        "ceil_mu": math.ceil(mu),
        "mode_minus_mu": mode_first - mu,
        "ceil_mu_minus_mode": math.ceil(mu) - mode_first,
        "low_mode_threshold": low_mode_threshold,
        "low_mode_surplus": low_mode_threshold - mode_first,
        "mean_gap_n_over_3": mean_gap_n_over_3,
        "mean_ratio_n_over_3": mu / (n / 3) if n else None,
        "first_descent": descent,
        "unimodal": is_unimodal(poly),
        "log_concave": is_log_concave(poly),
        "lc_ratio": lc_ratio,
        "lc_pos": lc_pos,
        "lc_defects": defects,
        "lc_defect_count": len(defects),
        "lc_defect_min_distance_from_mode": min(defect_distances_from_mode)
        if defect_distances_from_mode
        else None,
        "lc_defect_min_distance_from_near_miss": min(defect_distances_from_near_miss)
        if defect_distances_from_near_miss
        else None,
        "near_miss_ratio": nm_ratio,
        "near_miss_pos": nm_pos,
        "near_miss_reserve": 1.0 - nm_ratio if nm_pos != -1 else None,
        "target_k": target_k,
        "target_score": int(exact_score) if exact_score is not None else None,
        "distribution_score": distribution_score,
        "graph6": graph6,
        "canonical_graph6": canonical_g6,
        "expected_lc_breaks": list(expected_lc_breaks),
        "expected_lc_breaks_missing": [
            break_k for break_k in expected_lc_breaks if break_k not in actual_lc_breaks
        ],
        **struct,
    }
    if metadata:
        out.update(metadata)
    return out


def evaluate_prufer_record(
    record: PruferRecord,
    *,
    target_k: int | None,
    distribution_score: float | None,
    canonical_g6: str | None,
) -> dict[str, Any]:
    return evaluate_tree(
        prufer_to_adj(record.code),
        rank=record.line_no,
        target_k=target_k,
        distribution_score=distribution_score,
        canonical_g6=canonical_g6,
        metadata={"source_kind": "prufer", "line_no": record.line_no},
    )


def evaluate_generated_record(
    record: GeneratedRecord,
    *,
    target_k: int | None,
    canonical_g6: str | None,
) -> dict[str, Any]:
    effective_target_k = target_k if target_k is not None else record.default_target_k
    return evaluate_tree(
        record.adj,
        rank=record.rank,
        target_k=effective_target_k,
        distribution_score=None,
        canonical_g6=canonical_g6,
        metadata={
            "source_kind": "generated_family",
            "family": record.family,
            "params": record.params,
            "expected_lc_breaks": list(record.expected_lc_breaks),
        },
    )


def summarize_numeric(values: list[float]) -> dict[str, float | int | None]:
    if not values:
        return {"min": None, "max": None, "mean": None, "median": None}
    return {
        "min": min(values),
        "max": max(values),
        "mean": statistics.fmean(values),
        "median": statistics.median(values),
    }


def _histogram_sort_key(value: Any) -> tuple[int, Any]:
    if value is None:
        return (0, "")
    if isinstance(value, (int, float)):
        return (1, value)
    return (2, str(value))


def histogram(values: list[Any]) -> dict[str, int]:
    return {
        str(key): count
        for key, count in sorted(Counter(values).items(), key=lambda item: _histogram_sort_key(item[0]))
    }


def top_rows(rows: list[dict[str, Any]], key: str, limit: int) -> list[dict[str, Any]]:
    keep = sorted(
        rows,
        key=lambda row: float("-inf") if row[key] is None else row[key],
        reverse=True,
    )[:limit]
    return compact_rows(keep)


def bottom_rows(rows: list[dict[str, Any]], key: str, limit: int) -> list[dict[str, Any]]:
    keep = sorted(
        rows,
        key=lambda row: float("inf") if row[key] is None else row[key],
    )[:limit]
    return compact_rows(keep)


def compact_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    compact = []
    for row in rows:
        compact.append(
            {
                "rank": row["rank"],
                "source_kind": row.get("source_kind"),
                "family": row.get("family"),
                "params": row.get("params"),
                "n": row["n"],
                "alpha": row["alpha"],
                "mode_first": row["mode_first"],
                "mu": row["mu"],
                "ceil_mu_minus_mode": row["ceil_mu_minus_mode"],
                "low_mode_surplus": row["low_mode_surplus"],
                "mean_gap_n_over_3": row["mean_gap_n_over_3"],
                "mean_ratio_n_over_3": row["mean_ratio_n_over_3"],
                "lc_ratio": row["lc_ratio"],
                "lc_pos": row["lc_pos"],
                "lc_defect_count": row["lc_defect_count"],
                "lc_defect_min_distance_from_mode": row["lc_defect_min_distance_from_mode"],
                "lc_defect_min_distance_from_near_miss": row["lc_defect_min_distance_from_near_miss"],
                "near_miss_ratio": row["near_miss_ratio"],
                "near_miss_pos": row["near_miss_pos"],
                "near_miss_reserve": row["near_miss_reserve"],
                "target_score": row["target_score"],
                "expected_lc_breaks": row.get("expected_lc_breaks", []),
                "expected_lc_breaks_missing": row.get("expected_lc_breaks_missing", []),
                "leaves": row["leaves"],
                "degree2": row["degree2"],
                "degree2_adjacent_leaf": row["degree2_adjacent_leaf"],
                "max_d_leaf": row["max_d_leaf"],
                "graph6": row["graph6"],
                "canonical_graph6": row["canonical_graph6"],
            }
        )
    return compact


def build_summary(
    rows: list[dict[str, Any]],
    *,
    source: dict[str, Any],
    limit: int | None,
    canonicalize: bool,
    elapsed_s: float,
    top_limit: int,
) -> dict[str, Any]:
    lc_positions = []
    for row in rows:
        lc_positions.extend(defect["k"] for defect in row["lc_defects"])
    canon_values = [row["canonical_graph6"] for row in rows if row["canonical_graph6"]]
    graph6_values = {row["graph6"] for row in rows}
    expected_rows = [row for row in rows if row.get("expected_lc_breaks")]
    return {
        "source": {
            **source,
            "limit": limit,
            "canonicalize": canonicalize,
        },
        "elapsed_s": elapsed_s,
        "processed": len(rows),
        "n_hist": histogram([row["n"] for row in rows]),
        "alpha_hist": histogram([row["alpha"] for row in rows]),
        "family_hist": histogram([row.get("family") for row in rows if row.get("family")]),
        "mode_first_hist": histogram([row["mode_first"] for row in rows]),
        "first_descent_hist": histogram([row["first_descent"] for row in rows]),
        "lc_position_hist": histogram(lc_positions),
        "lc_defect_count_hist": histogram([row["lc_defect_count"] for row in rows]),
        "counts": {
            "non_unimodal": sum(1 for row in rows if not row["unimodal"]),
            "non_log_concave": sum(1 for row in rows if not row["log_concave"]),
            "low_mode": sum(1 for row in rows if row["low_mode_surplus"] >= 0),
            "mean_lt_n_over_3": sum(1 for row in rows if row["mean_gap_n_over_3"] > 0),
            "mode_le_ceil_mu": sum(1 for row in rows if row["ceil_mu_minus_mode"] >= 0),
            "lc_defects_all_after_mode": sum(
                1
                for row in rows
                if row["lc_defects"] and all(defect["k"] > row["mode_first"] for defect in row["lc_defects"])
            ),
            "lc_defects_all_after_near_miss": sum(
                1
                for row in rows
                if row["lc_defects"]
                and row["near_miss_pos"] != -1
                and all(defect["k"] > row["near_miss_pos"] for defect in row["lc_defects"])
            ),
            "positive_target_score": sum(
                1 for row in rows if row["target_score"] is not None and row["target_score"] > 0
            ),
            "expected_lc_breaks_all_hit": sum(
                1 for row in expected_rows if not row["expected_lc_breaks_missing"]
            ),
            "expected_lc_breaks_any_missing": sum(
                1 for row in expected_rows if row["expected_lc_breaks_missing"]
            ),
            "d_leaf_le1": sum(1 for row in rows if row["d_leaf_le1"]),
            "multi_leaf_hub": sum(1 for row in rows if row["multi_leaf_hubs"] > 0),
            "unique_graph6": len(graph6_values),
            "unique_labelled_prufer": len(graph6_values) if source.get("kind") == "prufer" else None,
            "unique_canonical_graph6": len(set(canon_values)) if canonicalize else None,
        },
        "numeric": {
            "lc_ratio": summarize_numeric([row["lc_ratio"] for row in rows]),
            "near_miss_ratio": summarize_numeric([row["near_miss_ratio"] for row in rows]),
            "near_miss_reserve": summarize_numeric(
                [row["near_miss_reserve"] for row in rows if row["near_miss_reserve"] is not None]
            ),
            "mu": summarize_numeric([row["mu"] for row in rows]),
            "mode_minus_mu": summarize_numeric([row["mode_minus_mu"] for row in rows]),
            "ceil_mu_minus_mode": summarize_numeric([row["ceil_mu_minus_mode"] for row in rows]),
            "low_mode_surplus": summarize_numeric([row["low_mode_surplus"] for row in rows]),
            "mean_gap_n_over_3": summarize_numeric([row["mean_gap_n_over_3"] for row in rows]),
            "mean_ratio_n_over_3": summarize_numeric(
                [row["mean_ratio_n_over_3"] for row in rows if row["mean_ratio_n_over_3"] is not None]
            ),
            "lc_defect_min_distance_from_mode": summarize_numeric(
                [
                    row["lc_defect_min_distance_from_mode"]
                    for row in rows
                    if row["lc_defect_min_distance_from_mode"] is not None
                ]
            ),
            "lc_defect_min_distance_from_near_miss": summarize_numeric(
                [
                    row["lc_defect_min_distance_from_near_miss"]
                    for row in rows
                    if row["lc_defect_min_distance_from_near_miss"] is not None
                ]
            ),
            "leaves": summarize_numeric([row["leaves"] for row in rows]),
            "degree2": summarize_numeric([row["degree2"] for row in rows]),
            "degree2_adjacent_leaf": summarize_numeric([row["degree2_adjacent_leaf"] for row in rows]),
            "max_d_leaf": summarize_numeric([row["max_d_leaf"] for row in rows]),
        },
        "top_by_lc_ratio": top_rows(rows, "lc_ratio", top_limit),
        "top_by_near_miss_ratio": top_rows(rows, "near_miss_ratio", top_limit),
        "top_by_target_score": top_rows(rows, "target_score", top_limit),
        "top_by_mean_ratio_n_over_3": top_rows(rows, "mean_ratio_n_over_3", top_limit),
        "bottom_by_mean_gap_n_over_3": bottom_rows(rows, "mean_gap_n_over_3", top_limit),
        "bottom_by_low_mode_surplus": bottom_rows(rows, "low_mode_surplus", top_limit),
        "bottom_by_near_miss_reserve": bottom_rows(rows, "near_miss_reserve", top_limit),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("codes", nargs="?", type=Path, help="Path to a search_output_*.txt Prüfer-code file")
    parser.add_argument("--scores", type=Path, help="Optional score distribution aligned to sorted codes")
    parser.add_argument(
        "--family",
        choices=["galvin", "bautista-ramos", "li", "li-star", "all"],
        help="Generate a structured family instead of reading a Prüfer-code corpus",
    )
    parser.add_argument(
        "--m-values",
        default="1-8",
        help="Family m parameters, e.g. 1-8 or 2,4,6; default 1-8",
    )
    parser.add_argument(
        "--t-values",
        default="1-8",
        help="Family t parameters, or n parameters for Li families; e.g. 1-8 or 5,8-10; default 1-8",
    )
    parser.add_argument("--max-n", type=int, help="Skip generated family trees above this order")
    parser.add_argument("--limit", type=int, default=1000, help="Rows to analyze; use 0 for all rows")
    parser.add_argument("--target-k", type=int, help="Coefficient index for target LC score; default n//2")
    parser.add_argument("--canonicalize", action="store_true", help="Canonicalize graph6 rows with nauty labelg")
    parser.add_argument("--top", type=int, default=10, help="Number of top rows retained in summary")
    parser.add_argument("--out", type=Path, required=True, help="Output JSON summary path")
    args = parser.parse_args()

    if args.codes and args.family:
        parser.error("pass either a Prüfer corpus path or --family, not both")
    if not args.codes and not args.family:
        parser.error("pass a Prüfer corpus path or choose --family")
    if args.scores and args.family:
        parser.error("--scores is only valid for Prüfer corpus input")
    if args.max_n is not None and args.max_n < 1:
        parser.error("--max-n must be positive")

    limit = None if args.limit == 0 else args.limit
    start = time.time()

    if args.family:
        try:
            m_values = parse_int_values(args.m_values)
            t_values = parse_int_values(args.t_values)
        except ValueError as exc:
            parser.error(str(exc))

        generated_records = iter_family_records(
            family=args.family,
            m_values=m_values,
            t_values=t_values,
            max_n=args.max_n,
        )
        if limit is not None:
            generated_records = generated_records[:limit]

        canonical_values: list[str | None] = [None] * len(generated_records)
        if args.canonicalize:
            graph6_values = [graph6_from_adj(record.adj) for record in generated_records]
            canonical_values = canonical_graph6_batch(graph6_values)

        rows = [
            evaluate_generated_record(
                record,
                target_k=args.target_k,
                canonical_g6=canonical_values[idx],
            )
            for idx, record in enumerate(generated_records)
        ]
        source = {
            "kind": "generated_family",
            "family": args.family,
            "m_values": m_values,
            "t_values": t_values,
            "max_n": args.max_n,
        }
    else:
        prufer_records = iter_prufer_records(args.codes, limit)
        scores = parse_score_distribution(args.scores, limit) if args.scores else []
        if scores and len(scores) != len(prufer_records):
            raise SystemExit(f"score count {len(scores)} does not match record count {len(prufer_records)}")

        canonical_values = [None] * len(prufer_records)
        if args.canonicalize:
            graph6_values = [graph6_from_adj(prufer_to_adj(record.code)) for record in prufer_records]
            canonical_values = canonical_graph6_batch(graph6_values)

        rows = []
        for idx, record in enumerate(prufer_records):
            rows.append(
                evaluate_prufer_record(
                    record,
                    target_k=args.target_k,
                    distribution_score=scores[idx] if scores else None,
                    canonical_g6=canonical_values[idx],
                )
            )
        source = {
            "kind": "prufer",
            "codes_path": str(args.codes),
            "scores_path": str(args.scores) if args.scores else None,
        }

    if not rows:
        raise SystemExit("no rows to analyze")

    elapsed_s = time.time() - start
    summary = build_summary(
        rows,
        source=source,
        limit=limit,
        canonicalize=args.canonicalize,
        elapsed_s=elapsed_s,
        top_limit=args.top,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write("\n")
    print(
        f"Analyzed {summary['processed']} rows in {elapsed_s:.2f}s; "
        f"non-LC={summary['counts']['non_log_concave']}, "
        f"non-unimodal={summary['counts']['non_unimodal']}, "
        f"wrote {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
