#!/usr/bin/env python3
"""Exact fixed-alpha extremal probes for the D13 prefix-GSB ratio.

For a tree T with independence coefficients i_j, define

    R_r(T) = (r+2)i_r i_{r+2} /
             ((r+1)i_{r+1}^2 + i_r i_{r+1}).

Prefix GSB is R_r(T) <= 1 for r <= ceil((2 alpha(T)-1)/3)-2.
This scratch driver groups exact maxima by (alpha,r) and compares them with
the perfect-matching corona of a star,

    C_alpha = S(2^(alpha-1),1),
    I(C_alpha;x)=(1+x)(1+2x)^(alpha-1)+x(1+x)^(alpha-1).

It is a falsification tool, not a proof of extremality.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from math import ceil, comb

from indpoly import independence_poly
from trees import trees_geng_raw


Adj = list[list[int]]


@dataclass(frozen=True)
class Row:
    num: int
    den: int
    n: int
    alpha: int
    rank: int
    graph6: str
    degrees: tuple[int, ...]

    @property
    def ratio(self) -> Fraction:
        return Fraction(self.num, self.den)


def gsb_ratio(poly: list[int], rank: int) -> Fraction:
    num = (rank + 2) * poly[rank] * poly[rank + 2]
    den = ((rank + 1) * poly[rank + 1] * poly[rank + 1]
           + poly[rank] * poly[rank + 1])
    return Fraction(num, den)


def corona_star_poly(alpha: int) -> list[int]:
    """I(S(2^(alpha-1),1);x), by a closed coefficient formula."""
    k = alpha - 1
    out: list[int] = []
    for j in range(alpha + 1):
        first = (2**j * comb(k, j)) if j <= k else 0
        second = ((2 ** (j - 1) + 1) * comb(k, j - 1)) if j >= 1 else 0
        out.append(first + second)
    return out


def prefix_last(alpha: int) -> int:
    return ceil((2 * alpha - 1) / 3) - 2


def scan(max_n: int) -> dict:
    best: dict[tuple[int, int], Row] = {}
    best_boundary: dict[int, Row] = {}
    counts: dict[int, int] = {}

    for n in range(1, max_n + 1):
        count = 0
        for _, adj, raw in trees_geng_raw(n):
            count += 1
            poly = independence_poly(n, adj)
            alpha = len(poly) - 1
            last = prefix_last(alpha)
            if last < 0:
                continue
            degrees = tuple(sorted((len(ns) for ns in adj), reverse=True))
            for rank in range(last + 1):
                ratio = gsb_ratio(poly, rank)
                row = Row(
                    ratio.numerator, ratio.denominator, n, alpha, rank,
                    raw.decode(), degrees,
                )
                key = (alpha, rank)
                if key not in best or row.ratio > best[key].ratio:
                    best[key] = row
                if rank == last and (
                    alpha not in best_boundary
                    or row.ratio > best_boundary[alpha].ratio
                ):
                    best_boundary[alpha] = row
        counts[n] = count
        print(json.dumps({"n": n, "trees": count}), flush=True)

    comparisons = []
    for (alpha, rank), row in sorted(best.items()):
        corona = gsb_ratio(corona_star_poly(alpha), rank)
        comparisons.append({
            "alpha": alpha,
            "rank": rank,
            "exhaustive_best": [row.num, row.den],
            "exhaustive_best_float": float(row.ratio),
            "corona": [corona.numerator, corona.denominator],
            "corona_float": float(corona),
            "best_exceeds_corona": row.ratio > corona,
            "attained_by_corona_ratio": row.ratio == corona,
            "witness": {
                "n": row.n,
                "graph6": row.graph6,
                "degrees": row.degrees,
            },
        })

    return {
        "max_n": max_n,
        "counts": counts,
        "comparisons": comparisons,
        "boundary": [
            {
                "alpha": alpha,
                "rank": row.rank,
                "ratio": [row.num, row.den],
                "ratio_float": float(row.ratio),
                "n": row.n,
                "graph6": row.graph6,
                "degrees": row.degrees,
            }
            for alpha, row in sorted(best_boundary.items())
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=18)
    parser.add_argument("--out")
    args = parser.parse_args()
    result = scan(args.max_n)
    payload = json.dumps(result, indent=2)
    if args.out:
        # Output redirection is handled by the caller; avoid persistent output
        # unless an explicit path is requested.
        from pathlib import Path
        Path(args.out).write_text(payload + "\n")
    print(payload)


if __name__ == "__main__":
    main()
