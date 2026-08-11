#!/usr/bin/env python3
"""Stress-test an independent-majority sufficient condition.

Conjectural lemma under test: if a graph has an independent set N of size m
and its complement P has size p, with m >= 2p + 2, then its independence
polynomial is unimodal (and perhaps even log-concave).

The graph is represented without expanding N: ``p_adj`` records the graph on
P and ``n_masks`` records, for each vertex of N, its neighborhood in P.  For
an independent S subseteq P, exactly those vertices of N whose masks avoid S
remain available, so

    I_G(x) = sum_{S independent in G[P]}
                 x^|S| (1+x)^#{u in N : N(u) cap S = empty}.

Everything is evaluated with exact integer arithmetic.  This file is a
counterexample search, not a proof of the lemma.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path

from indpoly import is_unimodal


def independent_polynomial(
    p: int, m: int, p_adj: list[int], n_masks: list[int],
) -> list[int]:
    """Return the exact independence polynomial of the encoded graph."""
    if len(p_adj) != p or len(n_masks) != m:
        raise ValueError("representation has the wrong dimensions")

    size = 1 << p
    independent = [False] * size
    independent[0] = True
    for mask in range(1, size):
        bit = mask & -mask
        v = bit.bit_length() - 1
        rest = mask ^ bit
        independent[mask] = independent[rest] and not (p_adj[v] & rest)

    # Zeta transform the multiplicities of N-neighborhood masks.  Then
    # subset_sum[complement(S)] is the number of N vertices still available.
    subset_sum = [0] * size
    for mask in n_masks:
        subset_sum[mask] += 1
    for v in range(p):
        bit = 1 << v
        for mask in range(size):
            if mask & bit:
                subset_sum[mask] += subset_sum[mask ^ bit]

    full = size - 1
    poly = [0] * (p + m + 1)
    binomial_rows = [[math.comb(q, k) for k in range(q + 1)] for q in range(m + 1)]
    for mask in range(size):
        if not independent[mask]:
            continue
        r = mask.bit_count()
        available = subset_sum[full ^ mask]
        for k, value in enumerate(binomial_rows[available]):
            poly[r + k] += value

    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def worst_lc_minor(poly: list[int]) -> tuple[int, int, int]:
    """Return (index, square, adjacent product) for the worst LC ratio."""
    if len(poly) < 3:
        return (0, 1, 0)
    worst = (1, poly[1] * poly[1], poly[0] * poly[2])
    for k in range(2, len(poly) - 1):
        lhs = poly[k] * poly[k]
        rhs = poly[k - 1] * poly[k + 1]
        if rhs * worst[1] > worst[2] * lhs:
            worst = (k, lhs, rhs)
    return worst


def first_valley(poly: list[int]) -> int | None:
    """Return the first strict fall-then-rise index, allowing plateaux."""
    falling = False
    last_strict_fall = -1
    for k in range(1, len(poly)):
        if poly[k] < poly[k - 1]:
            falling = True
            last_strict_fall = k
        elif poly[k] > poly[k - 1] and falling:
            return last_strict_fall
    return None


def worst_rebound(poly: list[int]) -> tuple[int, int, int]:
    """Return (valley index, later maximum, valley value).

    Only indices reached by a strict one-step fall are eligible.  The ratio
    ``later maximum / valley value`` exceeds one exactly when that fall is
    followed somewhere by a strict rebound, hence exactly when ``poly`` is
    nonunimodal.  Unlike an LC objective, this directly targets the property
    at issue.
    """
    suffix_max = [0] * len(poly)
    running = 0
    for k in range(len(poly) - 1, -1, -1):
        running = max(running, poly[k])
        suffix_max[k] = running

    worst = (0, 0, 1)
    for k in range(1, len(poly) - 1):
        if poly[k] >= poly[k - 1]:
            continue
        later = suffix_max[k + 1]
        if later * worst[2] > worst[1] * poly[k]:
            worst = (k, later, poly[k])
    return worst


def p_graph(p: int, style: str, rng: random.Random) -> list[int]:
    adj = [0] * p
    if style == "complete":
        probability = 1.0
    elif style == "empty":
        probability = 0.0
    else:
        probability = rng.choice((0.1, 0.25, 0.5, 0.75, 0.9))
    for u in range(p):
        for v in range(u + 1, p):
            if rng.random() < probability:
                adj[u] |= 1 << v
                adj[v] |= 1 << u
    return adj


def cross_masks(p: int, m: int, style: str, rng: random.Random) -> list[int]:
    full = (1 << p) - 1
    if style == "complete":
        return [full] * m
    if style == "empty":
        return [0] * m
    if style == "singletons":
        return [1 << (j % p) for j in range(m)]
    if style == "nested":
        order = list(range(p))
        rng.shuffle(order)
        out = []
        for j in range(m):
            width = (j * (p + 1)) // m
            mask = 0
            for v in order[:width]:
                mask |= 1 << v
            out.append(mask)
        rng.shuffle(out)
        return out
    if style == "blocks":
        width = rng.randint(1, max(1, p // 2))
        return [
            sum(1 << ((j + d) % p) for d in range(width))
            for j in range(m)
        ]
    if style == "polarized":
        choices = (0, full) + tuple(1 << v for v in range(p))
        return [rng.choice(choices) for _ in range(m)]

    probability = rng.choice((0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.95))
    return [
        sum(1 << v for v in range(p) if rng.random() < probability)
        for _ in range(m)
    ]


def encode_edges(p_adj: list[int]) -> list[list[int]]:
    return [
        [u, v]
        for u in range(len(p_adj))
        for v in range(u + 1, len(p_adj))
        if p_adj[u] & (1 << v)
    ]


def make_certificate(
    trial: int,
    source: str,
    p: int,
    m: int,
    p_adj: list[int],
    n_masks: list[int],
    poly: list[int],
) -> dict[str, object]:
    k, lhs, rhs = worst_lc_minor(poly)
    return {
        "trial": trial,
        "source": source,
        "p": p,
        "m": m,
        "p_edges": encode_edges(p_adj),
        "n_neighborhood_masks": n_masks,
        "polynomial": poly,
        "unimodal": is_unimodal(poly),
        "first_valley": first_valley(poly),
        "worst_lc_minor": {"k": k, "lhs": lhs, "rhs": rhs},
        "worst_lc_ratio": rhs / lhs,
    }


def mutate(
    p: int, p_adj: list[int], n_masks: list[int], rng: random.Random,
) -> tuple[list[int], list[int]]:
    new_adj = p_adj.copy()
    new_masks = n_masks.copy()
    if p >= 2 and rng.random() < 0.25:
        u = rng.randrange(p)
        v = rng.randrange(p - 1)
        if v >= u:
            v += 1
        new_adj[u] ^= 1 << v
        new_adj[v] ^= 1 << u
    else:
        j = rng.randrange(len(new_masks))
        new_masks[j] ^= 1 << rng.randrange(p)
    return new_adj, new_masks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=30_000)
    parser.add_argument("--max-p", type=int, default=12)
    parser.add_argument("--hill-restarts", type=int, default=40)
    parser.add_argument("--hill-steps", type=int, default=500)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument(
        "--hill-objective",
        choices=("lc", "valley", "mixed"),
        default="lc",
        help="quantity optimized by the local search",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/independent_majority_search_20260807.json"),
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    rng = random.Random(args.seed)

    p_styles = ("empty", "complete", "random")
    cross_styles = (
        "empty", "complete", "singletons", "nested", "blocks",
        "polarized", "random",
    )
    first_lc_failure: dict[str, object] | None = None
    first_unimodality_failure: dict[str, object] | None = None
    worst: dict[str, object] | None = None
    worst_lhs, worst_rhs = 1, 0
    worst_valley: dict[str, object] | None = None
    worst_valley_den, worst_valley_num = 1, 0
    completed = 0

    def test(
        source: str, p: int, m: int, adj: list[int], masks: list[int], trial: int,
    ) -> tuple[int, int]:
        nonlocal completed, first_lc_failure, first_unimodality_failure
        nonlocal worst, worst_lhs, worst_rhs
        nonlocal worst_valley, worst_valley_den, worst_valley_num
        poly = independent_polynomial(p, m, adj, masks)
        k, lhs, rhs = worst_lc_minor(poly)
        completed += 1
        if rhs * worst_lhs > worst_rhs * lhs:
            worst_lhs, worst_rhs = lhs, rhs
            worst = make_certificate(trial, source, p, m, adj, masks, poly)
        valley_k, later, valley = worst_rebound(poly)
        if later * worst_valley_den > worst_valley_num * valley:
            worst_valley_den, worst_valley_num = valley, later
            worst_valley = make_certificate(
                trial, source, p, m, adj, masks, poly,
            )
            worst_valley["rebound_minor"] = {
                "valley_k": valley_k,
                "later_max": later,
                "valley_value": valley,
                "ratio": later / valley,
            }
        if rhs > lhs and first_lc_failure is None:
            first_lc_failure = make_certificate(
                trial, source, p, m, adj, masks, poly,
            )
        if not is_unimodal(poly) and first_unimodality_failure is None:
            first_unimodality_failure = make_certificate(
                trial, source, p, m, adj, masks, poly,
            )
        return lhs, rhs

    # Random and deliberately extreme constructions at and above the boundary.
    for trial in range(args.trials):
        p = rng.randint(1, args.max_p)
        m = 2 * p + 2 + rng.randint(0, p)
        ps = p_styles[trial % len(p_styles)] if trial < 3 * len(cross_styles) else rng.choice(p_styles)
        cs = cross_styles[trial % len(cross_styles)] if trial < 3 * len(cross_styles) else rng.choice(cross_styles)
        test(
            f"structured-random:{ps}:{cs}",
            p,
            m,
            p_graph(p, ps, rng),
            cross_masks(p, m, cs, rng),
            trial,
        )
        if first_unimodality_failure is not None:
            break

    # Greedy local search maximizes either the exact worst LC ratio or the
    # exact post-descent rebound ratio at m=2p+2.
    # Occasional downhill moves prevent immediate trapping on plateaux.
    if first_unimodality_failure is None:
        for restart in range(args.hill_restarts):
            p = rng.randint(max(2, args.max_p // 2), args.max_p)
            m = 2 * p + 2
            adj = p_graph(p, "random", rng)
            masks = cross_masks(p, m, "random", rng)
            poly = independent_polynomial(p, m, adj, masks)
            _, lhs, rhs = worst_lc_minor(poly)
            _, rebound_num, rebound_den = worst_rebound(poly)
            for step in range(args.hill_steps):
                candidate_adj, candidate_masks = mutate(p, adj, masks, rng)
                candidate_poly = independent_polynomial(
                    p, m, candidate_adj, candidate_masks,
                )
                _, candidate_lhs, candidate_rhs = worst_lc_minor(candidate_poly)
                _, candidate_rebound_num, candidate_rebound_den = worst_rebound(
                    candidate_poly,
                )
                lc_improves = candidate_rhs * lhs >= rhs * candidate_lhs
                valley_improves = (
                    candidate_rebound_num * rebound_den
                    >= rebound_num * candidate_rebound_den
                )
                improves = (
                    lc_improves
                    if args.hill_objective == "lc"
                    else valley_improves
                    if args.hill_objective == "valley"
                    else lc_improves or valley_improves
                )
                if improves or rng.random() < 0.01:
                    adj, masks = candidate_adj, candidate_masks
                    lhs, rhs = candidate_lhs, candidate_rhs
                    rebound_num = candidate_rebound_num
                    rebound_den = candidate_rebound_den
                test(
                    "hill-climb",
                    p,
                    m,
                    adj,
                    masks,
                    args.trials + restart * args.hill_steps + step,
                )
                if first_unimodality_failure is not None:
                    break
            if first_unimodality_failure is not None:
                break

    summary = {
        "claim_scope": "exact adversarial computation; evidence only",
        "statement_tested": (
            "If N is independent, |N|=m, |V\\N|=p, and m>=2p+2, "
            "then I_G is unimodal; log-concavity was also tested."
        ),
        "seed": args.seed,
        "random_trials_requested": args.trials,
        "hill_restarts_requested": args.hill_restarts,
        "hill_steps_requested": args.hill_steps,
        "hill_objective": args.hill_objective,
        "cases_completed": completed,
        "first_lc_failure": first_lc_failure,
        "first_unimodality_failure": first_unimodality_failure,
        "worst_lc_case": worst,
        "worst_rebound_case": worst_valley,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    compact = {
        "cases_completed": completed,
        "lc_failure": first_lc_failure is not None,
        "unimodality_failure": first_unimodality_failure is not None,
        "worst_lc_ratio": worst["worst_lc_ratio"] if worst else None,
        "worst_rebound_ratio": (
            worst_valley["rebound_minor"]["ratio"]
            if worst_valley else None
        ),
    }
    print(json.dumps(compact if args.quiet else summary, indent=2))
    raise SystemExit(1 if first_unimodality_failure is not None else 0)


if __name__ == "__main__":
    main()
