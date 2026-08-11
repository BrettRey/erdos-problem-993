#!/usr/bin/env python3
"""Exact adversarial EADS search over caterpillars.

A caterpillar is encoded by the numbers of leaves attached to consecutive
vertices of its spine.  The search mutates those loads directly, which is much
more effective against EADS than generic leaf moves: the current low-good
example has loads ``(3, 0, 4, 1, 5, 6)`` and only three good vertices.

Every score uses exact independence-polynomial coefficients.  A zero-good
record is an exact counterexample to EADS, not necessarily to Erdos #993.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from scratch_eads_optimizer_20260807 import evaluate


def caterpillar(loads: tuple[int, ...]) -> list[list[int]]:
    if not loads:
        raise ValueError("the spine must be nonempty")
    n = len(loads) + sum(loads)
    adj: list[list[int]] = [[] for _ in range(n)]
    for u in range(len(loads) - 1):
        adj[u].append(u + 1)
        adj[u + 1].append(u)
    leaf = len(loads)
    for support, count in enumerate(loads):
        for _ in range(count):
            adj[support].append(leaf)
            adj[leaf].append(support)
            leaf += 1
    return adj


def score(ev: dict[str, object]) -> tuple[object, ...]:
    distances = ev["distances"]
    assert isinstance(distances, list)
    good = sorted(d for d in distances if d <= 1)
    # First eliminate good vertices.  At a tie, prefer fewer distance-zero
    # vertices, then push the remaining good distances and the whole profile
    # upward.  The final order term mildly favors smaller certificates.
    return (
        len(good),
        sum(d == 0 for d in good),
        tuple(-d for d in good),
        tuple(-d for d in sorted(distances)),
        len(distances),
    )


def mutate(
    loads: tuple[int, ...], rng: random.Random, max_load: int,
    min_spine: int,
) -> tuple[int, ...]:
    out = list(loads)
    move = rng.random()
    if move < 0.40:
        i = rng.randrange(len(out))
        delta = rng.choice((-2, -1, 1, 2))
        out[i] = min(max_load, max(0, out[i] + delta))
    elif move < 0.72 and sum(out) > 0:
        donors = [i for i, value in enumerate(out) if value]
        i = rng.choice(donors)
        j = rng.randrange(len(out))
        amount = 1 if out[i] == 1 else rng.choice((1, 1, 2))
        amount = min(amount, out[i], max_load - out[j])
        out[i] -= amount
        out[j] += amount
    elif move < 0.86 and len(out) < 18:
        i = rng.randrange(len(out) + 1)
        out.insert(i, rng.randint(0, min(6, max_load)))
    elif len(out) > min_spine:
        i = rng.randrange(len(out))
        value = out.pop(i)
        if out and value:
            j = min(i, len(out) - 1)
            out[j] = min(max_load, out[j] + value)
    return tuple(out)


def encode(loads: tuple[int, ...], ev: dict[str, object]) -> dict[str, object]:
    adj = caterpillar(loads)
    return {
        "loads": list(loads),
        "spine_order": len(loads),
        "n": len(adj),
        "edges": [
            [u, v]
            for u, neighbors in enumerate(adj)
            for v in neighbors
            if u < v
        ],
        **ev,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--restarts", type=int, default=100)
    parser.add_argument("--steps", type=int, default=500)
    parser.add_argument("--min-spine", type=int, default=3)
    parser.add_argument("--max-spine", type=int, default=12)
    parser.add_argument("--max-load", type=int, default=20)
    parser.add_argument("--min-n", type=int, default=20)
    parser.add_argument("--max-n", type=int, default=180)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/eads_caterpillar_search_20260807.json"),
    )
    args = parser.parse_args()
    rng = random.Random(args.seed)

    seed_loads = (3, 0, 4, 1, 5, 6)
    best_loads = seed_loads
    best_ev = evaluate(caterpillar(best_loads))
    evaluated = 1

    for restart in range(args.restarts):
        if restart == 0:
            loads = seed_loads
        else:
            while True:
                spine = rng.randint(args.min_spine, args.max_spine)
                loads = tuple(rng.randint(0, args.max_load) for _ in range(spine))
                order = len(loads) + sum(loads)
                if args.min_n <= order <= args.max_n:
                    break
        ev = evaluate(caterpillar(loads))
        evaluated += 1
        temperature = 0.06
        for _ in range(args.steps):
            candidate = mutate(loads, rng, args.max_load, args.min_spine)
            candidate_n = len(candidate) + sum(candidate)
            if candidate_n < args.min_n or candidate_n > args.max_n:
                continue
            candidate_ev = evaluate(caterpillar(candidate))
            evaluated += 1
            if score(candidate_ev) <= score(ev) or rng.random() < temperature:
                loads, ev = candidate, candidate_ev
            temperature *= 0.994
            if score(ev) < score(best_ev):
                best_loads, best_ev = loads, ev
                print(json.dumps({
                    "event": "improvement",
                    "restart": restart,
                    "evaluated": evaluated,
                    "loads": best_loads,
                    "n": len(caterpillar(best_loads)),
                    "good_count": best_ev["good_count"],
                    "good_vertices": best_ev["good_vertices"],
                }), flush=True)
            if best_ev["good_count"] == 0:
                break
        if restart % 10 == 0:
            print(json.dumps({
                "event": "progress",
                "restart": restart,
                "evaluated": evaluated,
                "best_good_count": best_ev["good_count"],
                "best_loads": best_loads,
            }), flush=True)
        if best_ev["good_count"] == 0:
            break

    result = {
        "claim_scope": "exact caterpillar search; evidence only",
        "seed": args.seed,
        "evaluated": evaluated,
        "counterexample_found": best_ev["good_count"] == 0,
        "best": encode(best_loads, best_ev),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "evaluated": evaluated,
        "counterexample_found": result["counterexample_found"],
        "best_good_count": best_ev["good_count"],
        "best_loads": best_loads,
        "out": str(args.out),
    }), flush=True)
    raise SystemExit(1 if result["counterexample_found"] else 0)


if __name__ == "__main__":
    main()
