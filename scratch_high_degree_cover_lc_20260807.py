#!/usr/bin/env python3
"""Exact stress test for the high-degree vertex-cover LC conjecture.

Conjecture: if a forest has a vertex cover P and every vertex of P has
degree at least four, then its independence polynomial is log-concave.

This is a candidate graph-theoretic class suggested by the observed EADS
orientations: in the exact small-tree corpus, right-oriented bad vertices are
independent and left-oriented bad vertices have degree at least four.  Those
two local statements are not proved, so this script does not establish a
reduction from an EADS counterexample.  The generator builds a random
P/N-labelled tree with N independent, then adds N-leaves until every P-vertex
has degree at least four.  All arithmetic is exact.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from indpoly import independence_poly, is_unimodal


def random_labelled_core(
    p: int, r: int, rng: random.Random,
) -> tuple[list[list[int]], list[str]]:
    """Build a recursive P/N tree in which no two N vertices are adjacent."""
    adj: list[list[int]] = [[]]
    labels = ["P"]
    remaining = ["P"] * (p - 1) + ["N"] * r
    rng.shuffle(remaining)
    for label in remaining:
        allowed = [
            v for v, old_label in enumerate(labels)
            if label == "P" or old_label == "P"
        ]
        parent = rng.choice(allowed)
        vertex = len(adj)
        adj.append([parent])
        adj[parent].append(vertex)
        labels.append(label)
    return adj, labels


def complete_degrees(
    adj: list[list[int]], labels: list[str], extra: int, rng: random.Random,
) -> None:
    """Attach N-leaves to make every P degree >=4, then add extra leaves."""
    p_vertices = [v for v, label in enumerate(labels) if label == "P"]
    loads = {v: max(0, 4 - len(adj[v])) for v in p_vertices}
    for _ in range(extra):
        loads[rng.choice(p_vertices)] += 1
    for support, count in loads.items():
        for _ in range(count):
            leaf = len(adj)
            adj.append([support])
            adj[support].append(leaf)
            labels.append("N")


def worst_minor(poly: list[int]) -> tuple[int, int, int]:
    worst = (1, poly[1] * poly[1], poly[0] * poly[2])
    for k in range(2, len(poly) - 1):
        lhs = poly[k] * poly[k]
        rhs = poly[k - 1] * poly[k + 1]
        if rhs * worst[1] > worst[2] * lhs:
            worst = k, lhs, rhs
    return worst


def encode(
    trial: int, adj: list[list[int]], labels: list[str], poly: list[int],
) -> dict[str, object]:
    k, lhs, rhs = worst_minor(poly)
    return {
        "trial": trial,
        "n": len(adj),
        "labels": labels,
        "edges": [
            [u, v]
            for u, neighbors in enumerate(adj)
            for v in neighbors if u < v
        ],
        "polynomial": poly,
        "unimodal": is_unimodal(poly),
        "worst_minor": {"k": k, "lhs": lhs, "rhs": rhs},
        "worst_ratio": rhs / lhs,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=100_000)
    parser.add_argument("--max-p", type=int, default=40)
    parser.add_argument("--max-extra-per-p", type=int, default=8)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/high_degree_cover_lc_20260807.json"),
    )
    args = parser.parse_args()
    rng = random.Random(args.seed)
    first_lc_failure = None
    first_unimodality_failure = None
    worst = None
    worst_lhs, worst_rhs = 1, 0

    for trial in range(args.trials):
        p = rng.randint(1, args.max_p)
        r = rng.randint(0, max(0, p - 1))
        adj, labels = random_labelled_core(p, r, rng)
        complete_degrees(
            adj, labels, rng.randint(0, args.max_extra_per_p * p), rng,
        )
        poly = independence_poly(len(adj), adj)
        k, lhs, rhs = worst_minor(poly)
        if rhs * worst_lhs > worst_rhs * lhs:
            worst_lhs, worst_rhs = lhs, rhs
            worst = encode(trial, adj, labels, poly)
        if lhs < rhs and first_lc_failure is None:
            first_lc_failure = encode(trial, adj, labels, poly)
            print(json.dumps({
                "event": "lc_failure", "trial": trial,
                "n": len(adj), "ratio": rhs / lhs,
            }), flush=True)
        if not is_unimodal(poly):
            first_unimodality_failure = encode(trial, adj, labels, poly)
            break
        if trial and trial % 10_000 == 0:
            print(json.dumps({
                "event": "progress", "trial": trial,
                "worst_ratio": worst_rhs / worst_lhs,
                "lc_failure": first_lc_failure is not None,
            }), flush=True)

    result = {
        "claim_scope": "exact randomized computation; evidence only",
        "statement_tested": (
            "A tree with a vertex cover P such that deg(v)>=4 for every "
            "v in P has a log-concave independence polynomial."
        ),
        "seed": args.seed,
        "trials_completed": (
            first_unimodality_failure["trial"] + 1
            if first_unimodality_failure is not None else args.trials
        ),
        "first_lc_failure": first_lc_failure,
        "first_unimodality_failure": first_unimodality_failure,
        "worst": worst,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "trials": result["trials_completed"],
        "lc_failure": first_lc_failure is not None,
        "unimodality_failure": first_unimodality_failure is not None,
        "worst_ratio": worst_rhs / worst_lhs,
        "out": str(args.out),
    }), flush=True)
    raise SystemExit(1 if first_unimodality_failure is not None else 0)


if __name__ == "__main__":
    main()
