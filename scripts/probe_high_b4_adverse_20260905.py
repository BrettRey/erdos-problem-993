#!/usr/bin/env python3
"""Target high b4/e4 with D<0 in the remaining b1-positive window.

Search evidence is not a theorem. All graph counts and reported margins are
exact; randomness only chooses mutations and permits occasional downhill moves.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path
from time import perf_counter

import networkx as nx

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from indpoly import is_unimodal
from scripts.audit_b1_zero_forbidden_core_20260905 import (
    archived_rows,
    coefficient,
    extendable_defects,
    polynomial,
)
from scripts.replay_b1_positive_sign_obstructions_20260905 import MUTATED_GRAPH6

SOURCES = ["cross_reserve_witness_lifts_20260904.json",
           "cross_reserve_multicell_lifts_20260904.json",
           "negative_correction_seeds_n18_20260904.json"]


def profile(graph: nx.Graph) -> dict | None:
    p = polynomial(graph)
    a = len(p) - 1
    delta = 2 * a - len(graph)
    if a not in (17, 18, 19) or not 33 <= len(graph) <= 38 or delta > 5:
        return None
    e = extendable_defects(graph)
    b = [coefficient(p, a - d) - e[d] for d in range(5)]
    assert min(b) >= 0 and b[0] == 0
    assert 6 * b[3] >= (a - 4) * b[2]
    assert (a - 3) * e[3] >= 4 * e[4]
    d = b[3] ** 2 - b[2] * b[4] + 2 * e[3] * b[3] - e[2] * b[4] - e[4] * b[2]
    full = p[a - 3] ** 2 - p[a - 2] * p[a - 4]
    threshold = Fraction(a - 7, 3 * (a - 3))
    density = Fraction(b[4], e[4])
    if density <= threshold:
        assert full > 0
    return {"graph6": nx.to_graph6_bytes(graph, header=False).decode().strip(),
            "n": len(graph), "alpha": a, "deficiency": delta,
            "e_0_to_4": e, "b_0_to_4": b, "combined_correction": d,
            "full_depth3_margin": full, "unimodal": is_unimodal(p),
            "b4_over_e4": str(density), "proved_density_threshold": str(threshold),
            "density_over_threshold": str(density / threshold),
            "high_density_adverse": d < 0 and density > threshold}


def mutate(graph: nx.Graph, rng: random.Random) -> nx.Graph:
    candidate = graph.copy()
    if rng.random() < 0.75:
        leaf = rng.choice([v for v in graph if graph.degree(v) == 1])
        old = next(iter(graph[leaf]))
        new = rng.choice([v for v in graph if v not in (leaf, old)])
        candidate.remove_edge(leaf, old)
        candidate.add_edge(leaf, new)
    else:
        left, right = rng.choice(list(graph.edges()))
        candidate.remove_edge(left, right)
        component = nx.node_connected_component(candidate, left)
        u = rng.choice(sorted(component))
        v = rng.choice(sorted(set(graph) - component))
        candidate.add_edge(u, v)
    assert nx.is_tree(candidate)
    return candidate


def run(steps: int, seed: int) -> dict:
    start = perf_counter()
    archive = {MUTATED_GRAPH6}
    for name in SOURCES:
        for row in archived_rows(json.loads((REPO / "results" / name).read_text())):
            if row["b_0_to_4"][1]:
                archive.add(row["graph6"])
    seeds = []
    for graph6 in sorted(archive):
        graph = nx.from_graph6_bytes(graph6.encode())
        row = profile(graph)
        if row is not None and row["combined_correction"] < 0:
            seeds.append((graph, row))
    assert len(seeds) == 6
    best_graph, best = max(seeds, key=lambda pair: Fraction(pair[1]["density_over_threshold"]))
    rng = random.Random(seed)
    counters = Counter(seed_trees=len(seeds))
    current_graph, current = best_graph, best
    unique = set()
    found = None
    for index in range(steps):
        if index % 200 == 0:
            current_graph, current = (best_graph, best) if rng.random() < 0.5 else rng.choice(seeds)
        candidate = mutate(current_graph, rng)
        counters["mutation_attempts"] += 1
        row = profile(candidate)
        if row is None:
            counters["outside_window"] += 1
            continue
        counters["in_window_profile_checks"] += 1
        unique.add(row["graph6"])
        if row["full_depth3_margin"] < 0 or not row["unimodal"]:
            found = row
            counters["target_failures"] += 1
            break
        if row["b_0_to_4"][1] == 0 or row["combined_correction"] >= 0:
            continue
        counters["b1_positive_adverse_visits"] += 1
        score = Fraction(row["density_over_threshold"])
        if score > Fraction(best["density_over_threshold"]):
            best_graph, best = candidate, row
            counters["best_improvements"] += 1
        if row["high_density_adverse"]:
            found = row
            counters["high_density_adverse_found"] += 1
            break
        if score >= Fraction(current["density_over_threshold"]) or rng.random() < 0.08:
            current_graph, current = candidate, row
    # Fresh graph6 replay is independent of the mutation object's state.
    assert profile(nx.from_graph6_bytes(best["graph6"].encode())) == best
    counters["distinct_labeled_graph6_profiles"] = len(unique)
    return {"kind": "high_b4_adverse_targeted_probe", "certificate_date": "2026-09-05",
            "scope": {"n": [33, 38], "alpha": [17, 18, 19], "delta_max": 5,
                      "requested_mutations": steps, "seed": seed, "archive_sources": SOURCES},
            "claim_status": "bounded adaptive search, not exhaustive evidence or a theorem",
            "counters": dict(counters), "best_adverse": best, "first_target_found": found,
            "elapsed_seconds": round(perf_counter() - start, 3)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260905)
    parser.add_argument("--output", type=Path,
                        default=REPO / "results/high_b4_adverse_probe_20260905.json")
    args = parser.parse_args()
    result = run(args.steps, args.seed)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"counters": result["counters"], "best_adverse": result["best_adverse"],
                      "first_target_found": result["first_target_found"],
                      "elapsed_seconds": result["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
