#!/usr/bin/env python3
"""Falsify the oriented fork inequality in the abstract P=E+xR algebra."""

from __future__ import annotations

import argparse
import random

from indpoly import _polyadd, _polymul


def add_many(polys: list[list[int]]) -> list[int]:
    out = [0]
    for poly in polys:
        out = _polyadd(out, poly)
    return out


def product(polys: list[list[int]]) -> list[int]:
    out = [1]
    for poly in polys:
        out = _polymul(out, poly)
    return out


def c(poly: list[int], rank: int) -> int:
    return poly[rank] if rank < len(poly) else 0


def trial(rng: random.Random, degree: int, max_len: int, max_coeff: int):
    excluded = []
    reserve = []
    total = []
    for _ in range(degree):
        length = rng.randint(1, max_len)
        e = [rng.randint(0, max_coeff) for _ in range(length)]
        # In a rooted tree R_i counts a subfamily of the configurations
        # counted by E_i, grade by grade.  Enforce that load-bearing fact.
        r = [rng.randint(0, value) for value in e]
        # Keep constant terms nonzero, as for graph polynomials.
        e[0] = 1
        r[0] = 1
        p = _polyadd(e, [0] + r)
        excluded.append(e)
        reserve.append(r)
        total.append(p)

    f = product(total)
    g = product(excluded)
    whole = _polyadd(f, [0] + g)
    a_center = g
    a_neighbor = [
        product([reserve[i], *[total[j] for j in range(degree) if j != i]])
        for i in range(degree)
    ]
    joint = add_many(
        [
            product(
                [
                    reserve[i],
                    reserve[j],
                    *[total[k] for k in range(degree) if k not in (i, j)],
                ]
            )
            for i in range(degree)
            for j in range(i + 1, degree)
        ]
    )
    for rank, n_r in enumerate(whole):
        av = [c(poly, rank) for poly in a_neighbor]
        s = c(a_center, rank) + sum(av)
        parent = c(a_center, rank) + av[0]
        gap = s * s - parent * parent - 2 * n_r * c(joint, rank)
        if gap < 0:
            return {
                "degree": degree,
                "rank": rank,
                "gap": gap,
                "excluded": excluded,
                "reserve": reserve,
                "N": n_r,
                "a_center": c(a_center, rank),
                "a_neighbor": av,
                "joint": c(joint, rank),
            }
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=10000)
    parser.add_argument("--max-degree", type=int, default=6)
    parser.add_argument("--max-len", type=int, default=5)
    parser.add_argument("--max-coeff", type=int, default=20)
    parser.add_argument("--seed", type=int, default=252993)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    for index in range(args.trials):
        witness = trial(
            rng,
            rng.randint(2, args.max_degree),
            args.max_len,
            args.max_coeff,
        )
        if witness:
            print({"trial": index, "failure": witness}, flush=True)
            raise SystemExit(1)
    print({"trials": args.trials, "failure": None}, flush=True)


if __name__ == "__main__":
    main()
