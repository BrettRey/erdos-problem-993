#!/usr/bin/env python3
"""Adversarial pendant-load search on two-facet independence complexes.

The independence complex of the complete bipartite graph ``K_{r,b-r}`` is the
union of two disjoint simplices.  This gives a compact exact formula even when
the face set is enormous, and supplies a sharp test of whether pendant-leaf
majority statements genuinely need the core to be a tree.
"""

from __future__ import annotations

import argparse
import json
import random
from math import comb
from pathlib import Path

from indpoly import _polymul, is_unimodal


def binomial_poly(power: int) -> list[int]:
    return [comb(power, k) for k in range(power + 1)]


def poly_add(a: list[int], b: list[int], scale: int = 1) -> list[int]:
    out = [0] * max(len(a), len(b))
    for k, value in enumerate(a):
        out[k] += value
    for k, value in enumerate(b):
        out[k] += scale * value
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def simplex_transform(indices: range, q: list[int], ell: int) -> list[int]:
    """Transform of the full simplex on ``indices`` under private loads q."""
    total_on_facet = sum(q[i] for i in indices)
    poly = binomial_poly(ell - total_on_facet)
    for i in indices:
        factor = binomial_poly(q[i])
        if len(factor) == 1:
            factor.append(0)
        factor[1] += 1
        poly = _polymul(poly, factor)
    return poly


def two_facet_polynomial(b: int, cut: int, q: list[int]) -> list[int]:
    """Pendant transform for the independence complex of K_{cut,b-cut}."""
    ell = sum(q)
    left = simplex_transform(range(0, cut), q, ell)
    right = simplex_transform(range(cut, b), q, ell)
    # The empty core face occurs in both simplex transforms.
    return poly_add(poly_add(left, right), binomial_poly(ell), -1)


def worst_lc_minor(poly: list[int]) -> tuple[int, int, int]:
    worst = (1, poly[1] * poly[1], poly[0] * poly[2])
    for k in range(2, len(poly) - 1):
        lhs = poly[k] * poly[k]
        rhs = poly[k - 1] * poly[k + 1]
        if rhs * worst[1] > worst[2] * lhs:
            worst = (k, lhs, rhs)
    return worst


def load_vector(
    b: int, ell: int, cut: int, style: int, rng: random.Random,
) -> list[int]:
    q = [0] * b
    if style == 0:
        q[rng.randrange(b)] = ell
    elif style == 1:
        for _ in range(ell):
            q[rng.randrange(b)] += 1
    elif style == 2:
        side = range(0, cut) if rng.random() < 0.5 else range(cut, b)
        choices = tuple(side)
        for _ in range(ell):
            q[rng.choice(choices)] += 1
    else:
        for i in range(ell):
            q[i % b] += 1
    return q


def certificate(
    trial: int, b: int, ell: int, cut: int, q: list[int], poly: list[int],
) -> dict[str, object]:
    k, lhs, rhs = worst_lc_minor(poly)
    return {
        "trial": trial,
        "core_graph": f"K_{{{cut},{b-cut}}}",
        "b": b,
        "ell": ell,
        "q": q,
        "polynomial": poly,
        "unimodal": is_unimodal(poly),
        "worst_lc_minor": {"k": k, "lhs": lhs, "rhs": rhs},
        "worst_lc_ratio": rhs / lhs,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=100_000)
    parser.add_argument("--min-b", type=int, default=6)
    parser.add_argument("--max-b", type=int, default=60)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/pendant_two_facet_search_20260807.json"),
    )
    args = parser.parse_args()

    rng = random.Random(args.seed)
    first_lc_failure: dict[str, object] | None = None
    first_unimodality_failure: dict[str, object] | None = None
    worst: dict[str, object] | None = None
    worst_lhs = 1
    worst_rhs = 0

    for trial in range(args.trials):
        b = rng.randint(args.min_b, args.max_b)
        ell = b + 2 + rng.randint(0, b)
        cut = rng.randint(1, b - 1)
        q = load_vector(b, ell, cut, rng.randrange(4), rng)
        poly = two_facet_polynomial(b, cut, q)
        k, lhs, rhs = worst_lc_minor(poly)
        if rhs * worst_lhs > worst_rhs * lhs:
            worst_lhs, worst_rhs = lhs, rhs
            worst = certificate(trial, b, ell, cut, q, poly)
        if lhs < rhs and first_lc_failure is None:
            first_lc_failure = certificate(trial, b, ell, cut, q, poly)
        if not is_unimodal(poly):
            first_unimodality_failure = certificate(
                trial, b, ell, cut, q, poly,
            )
            break

    summary = {
        "claim_scope": "exact adversarial computation; evidence only",
        "family": "private-leaf attachments to complete bipartite cores",
        "seed": args.seed,
        "trials_requested": args.trials,
        "trials_completed": (
            first_unimodality_failure["trial"] + 1
            if first_unimodality_failure is not None
            else args.trials
        ),
        "first_lc_failure": first_lc_failure,
        "first_unimodality_failure": first_unimodality_failure,
        "worst_lc_case": worst,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    if not args.quiet:
        print(json.dumps(summary, indent=2))
    else:
        print({
            "trials_completed": summary["trials_completed"],
            "lc_failure": first_lc_failure is not None,
            "unimodality_failure": first_unimodality_failure is not None,
            "worst_lc_ratio": worst["worst_lc_ratio"] if worst else None,
        })
    raise SystemExit(1 if first_unimodality_failure is not None else 0)


if __name__ == "__main__":
    main()
