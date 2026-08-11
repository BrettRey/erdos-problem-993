#!/usr/bin/env python3
"""Exact adversarial test of the high-degree vertex-cover LC conjecture.

For each input tree, find a minimum-cost vertex cover P, where selecting v
costs max(0, 4-deg(v)).  Attach exactly that many new leaves to each selected
vertex.  The resulting tree has a certified vertex cover P whose vertices all
have degree at least four.  We then recompute its independence polynomial and
all log-concavity minors using arbitrary-precision integer arithmetic.

This is evidence, not a proof of the conjecture.  Its purpose is to aim the
strongest known non-log-concave tree families directly at the proposed class.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

from graph6 import parse_graph6
from indpoly import independence_poly, is_log_concave, is_unimodal
from scripts.analyze_prufer_corpus import (
    make_bautista_ramos_tree,
    make_galvin_tree,
    make_li_tree,
)


def minimum_completion_cover(
    adj: list[list[int]], threshold: int = 4,
) -> tuple[list[int], int]:
    """Return a deterministic minimum-cost vertex cover and its cost."""
    n = len(adj)
    parent = [-2] * n
    children: list[list[int]] = [[] for _ in range(n)]
    roots: list[int] = []
    order: list[int] = []
    for root in range(n):
        if parent[root] != -2:
            continue
        roots.append(root)
        parent[root] = -1
        stack = [root]
        while stack:
            v = stack.pop()
            order.append(v)
            for u in adj[v]:
                if u == parent[v]:
                    continue
                if parent[u] != -2:
                    raise ValueError("input is not a forest")
                parent[u] = v
                children[v].append(u)
                stack.append(u)

    # Each state is (leaf-completion cost, cover cardinality).  The second
    # coordinate makes ties deterministic and mildly favors smaller covers.
    dp_out = [(0, 0)] * n
    dp_in = [(0, 0)] * n
    for v in reversed(order):
        in_state = (max(0, threshold - len(adj[v])), 1)
        out_state = (0, 0)
        for child in children[v]:
            in_child = min(dp_out[child], dp_in[child])
            in_state = (
                in_state[0] + in_child[0],
                in_state[1] + in_child[1],
            )
            # If v is excluded, every child must be in the cover.
            out_state = (
                out_state[0] + dp_in[child][0],
                out_state[1] + dp_in[child][1],
            )
        dp_in[v] = in_state
        dp_out[v] = out_state

    selected: list[int] = []
    stack: list[tuple[int, bool]] = []
    for root in roots:
        stack.append((root, dp_in[root] < dp_out[root]))
    while stack:
        v, take = stack.pop()
        if take:
            selected.append(v)
        for child in children[v]:
            child_take = True if not take else dp_in[child] < dp_out[child]
            stack.append((child, child_take))

    cost = sum(max(0, threshold - len(adj[v])) for v in selected)
    expected = sum(min(dp_out[root], dp_in[root])[0] for root in roots)
    if cost != expected:
        raise AssertionError((cost, expected))
    return sorted(selected), cost


def complete_at_cover(
    adj: list[list[int]], cover: Iterable[int], threshold: int = 4,
) -> tuple[list[list[int]], int]:
    """Clone a forest and attach the minimum required leaves at cover nodes."""
    out = [neighbors.copy() for neighbors in adj]
    added = 0
    for v in cover:
        for _ in range(max(0, threshold - len(adj[v]))):
            leaf = len(out)
            out.append([v])
            out[v].append(leaf)
            added += 1
    return out, added


def validate_completion(
    original: list[list[int]], completed: list[list[int]], cover: list[int],
    threshold: int = 4,
) -> None:
    chosen = set(cover)
    for u, neighbors in enumerate(original):
        for v in neighbors:
            if u < v and u not in chosen and v not in chosen:
                raise AssertionError(f"uncovered original edge {u}-{v}")
    for u, neighbors in enumerate(completed):
        for v in neighbors:
            if u < v and u not in chosen and v not in chosen:
                raise AssertionError(f"uncovered completed edge {u}-{v}")
    if any(len(completed[v]) < threshold for v in cover):
        raise AssertionError("cover degree below threshold")
    edges = sum(map(len, completed)) // 2
    if completed and edges != len(completed) - 1:
        raise AssertionError("completed graph is not a tree")


def worst_minor(poly: list[int]) -> tuple[int, int, int]:
    """Return k,lhs,rhs maximizing rhs/lhs in the LC inequalities."""
    worst = (1, poly[1] * poly[1], poly[0] * poly[2])
    for k in range(2, len(poly) - 1):
        lhs = poly[k] * poly[k]
        rhs = poly[k - 1] * poly[k + 1]
        if rhs * worst[1] > worst[2] * lhs:
            worst = (k, lhs, rhs)
    return worst


def analyze(label: str, params: dict[str, int], adj: list[list[int]]) -> dict:
    original_poly = independence_poly(len(adj), adj)
    cover, cost = minimum_completion_cover(adj)
    completed, added = complete_at_cover(adj, cover)
    if added != cost:
        raise AssertionError((added, cost))
    validate_completion(adj, completed, cover)
    poly = independence_poly(len(completed), completed)
    k, lhs, rhs = worst_minor(poly)
    return {
        "label": label,
        "params": params,
        "original_order": len(adj),
        "original_log_concave": is_log_concave(original_poly),
        "original_unimodal": is_unimodal(original_poly),
        "cover": cover,
        "cover_size": len(cover),
        "leaves_added": added,
        "completed_order": len(completed),
        "completed_log_concave": lhs >= rhs,
        "completed_unimodal": is_unimodal(poly),
        "worst_minor": {"k": k, "lhs": lhs, "rhs": rhs},
        "worst_ratio": rhs / lhs,
    }


def known_failures() -> list[tuple[str, dict[str, int], list[list[int]]]]:
    records: list[tuple[str, dict[str, int], list[list[int]]]] = []
    sources = (
        (Path("results/analysis_n26.json"), "lc_failures"),
        (Path("results/analysis_n28_modal_lc_nm.json"), "top_lc_failures"),
    )
    for path, key in sources:
        payload = json.loads(path.read_text())
        for rank, item in enumerate(payload[key], start=1):
            n, adj = parse_graph6(item["graph6"].encode("ascii"))
            if n != item["n"]:
                raise AssertionError((n, item["n"]))
            records.append((
                f"known_n{n}_lc_failure_{rank}",
                {"n": n, "rank": rank},
                adj,
            ))
    return records


def family_instances(max_parameter: int, bautista_max_m: int):
    for m in range(1, max_parameter + 1):
        for t in range(1, max_parameter + 1):
            yield "galvin", {"m": m, "t": t}, make_galvin_tree(m, t)
    for m in range(1, bautista_max_m + 1):
        for t in range(1, max_parameter + 1):
            yield (
                "bautista_ramos", {"m": m, "t": t},
                make_bautista_ramos_tree(m, t),
            )
    for m in range(1, max_parameter + 1):
        for n in range(1, max_parameter + 1):
            yield "li", {"m": m, "n": n}, make_li_tree(m, n)
            yield (
                "li_star", {"m": m, "n": n},
                make_li_tree(m, n, starred=True),
            )


def summarize(records: list[dict]) -> dict:
    worst = max(records, key=lambda item: item["worst_ratio"])
    return {
        "instances": len(records),
        "original_lc_failures": sum(
            not item["original_log_concave"] for item in records
        ),
        "completed_lc_failures": sum(
            not item["completed_log_concave"] for item in records
        ),
        "completed_unimodality_failures": sum(
            not item["completed_unimodal"] for item in records
        ),
        "worst_completed_case": worst,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-parameter", type=int, default=20)
    parser.add_argument("--bautista-max-m", type=int, default=8)
    parser.add_argument("--known-only", action="store_true")
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/high_degree_cover_completion_20260808.json"),
    )
    args = parser.parse_args()

    known_records = [analyze(*item) for item in known_failures()]
    family_records: list[dict] = []
    if not args.known_only:
        for index, item in enumerate(
            family_instances(args.max_parameter, args.bautista_max_m), start=1,
        ):
            record = analyze(*item)
            # The adversarial corpus consists of the original LC failures;
            # retaining only them keeps the certificate compact.
            if not record["original_log_concave"]:
                family_records.append(record)
            if index % 100 == 0:
                print(json.dumps({
                    "event": "progress", "family_instances_examined": index,
                    "original_lc_failures": len(family_records),
                }), flush=True)

    result = {
        "claim_scope": "exact finite computation; evidence only",
        "statement_tested": (
            "Completing a minimum-cost vertex cover to minimum degree four "
            "repairs the selected known non-log-concave tree instances."
        ),
        "threshold": 4,
        "known_failures": summarize(known_records),
        "family_parameter_bounds": None if args.known_only else {
            "galvin_m_t": [1, args.max_parameter],
            "bautista_m": [1, args.bautista_max_m],
            "bautista_t": [1, args.max_parameter],
            "li_m_n": [1, args.max_parameter],
        },
        "family_failures": None if args.known_only else summarize(family_records),
        "known_records": known_records,
        "family_records": family_records,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    completed_failures = (
        result["known_failures"]["completed_lc_failures"]
        + (0 if args.known_only else result["family_failures"]["completed_lc_failures"])
    )
    print(json.dumps({
        "event": "done",
        "known": result["known_failures"],
        "family": result["family_failures"],
        "out": str(args.out),
    }), flush=True)
    raise SystemExit(1 if completed_failures else 0)


if __name__ == "__main__":
    main()
