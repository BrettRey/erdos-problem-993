#!/usr/bin/env python3
"""Exact falsification driver for the D25 aggregate nonedge covariance lemma.

For a tree T and rank r let N=i_r(T) and

    a_v = [x^r] I(T-N[v]).

The candidate lemma is

    N * C(r+2,2) * i_(r+2) <= sum_({u,v} nonedge) a_u a_v

in the early prefix r>=4, r<=ceil((2 alpha-1)/3)-2.  The difference of the
left and right sides is exactly q2+qfar in the D18 distance-class profile.
All comparisons in this driver use Python integers.
"""

from __future__ import annotations

import argparse
import random
from fractions import Fraction

from scratch_d18_covariance_search_20260711 import covariance_profile, prufer_tree
from targeted import make_T_m_t_1, make_T_m_t_d
from trees import trees


def candidate_rows(adj: list[list[int]]) -> tuple[dict, list[dict]]:
    profile = covariance_profile(adj)
    rows = [
        row
        for row in profile["rows"]
        if row["prefix"] and row["rank"] >= 4
    ]
    for row in rows:
        row["q_nonedge"] = row["q2"] + row["qfar"]
    return profile, rows


def update_best(best: tuple | None, profile: dict, row: dict, label: str):
    # Normalize by N*d1 only for reporting; exact cross multiplication selects.
    denominator = profile["rows"][row["rank"]]["N"] * row["d1"]
    if denominator == 0:
        return best
    item = (row["q_nonedge"], denominator, label, profile["n"], profile["alpha"], row)
    if best is None or item[0] * best[1] > best[0] * item[1]:
        return item
    return best


def report(best: tuple | None) -> dict | None:
    if best is None:
        return None
    num, den, label, n, alpha, row = best
    return {
        "label": label,
        "n": n,
        "alpha": alpha,
        "rank": row["rank"],
        "q_nonedge": num,
        "normalizer": den,
        "ratio": str(Fraction(num, den)),
        "ratio_float": float(Fraction(num, den)),
        "q2": row["q2"],
        "qfar": row["qfar"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exhaustive-max-n", type=int, default=0)
    parser.add_argument("--random", type=int, default=0)
    parser.add_argument("--min-random-n", type=int, default=20)
    parser.add_argument("--max-random-n", type=int, default=100)
    parser.add_argument("--galvin-max-n", type=int, default=0)
    parser.add_argument("--generalized-max-n", type=int, default=0)
    parser.add_argument("--seed", type=int, default=250993)
    args = parser.parse_args()

    best = None
    for n in range(1, args.exhaustive_max_n + 1):
        count = 0
        for graph6, adj in trees(n):
            profile, rows = candidate_rows(adj)
            for row in rows:
                best = update_best(best, profile, row, f"tree:{graph6}")
                if row["q_nonedge"] > 0:
                    print({"failure": report(best), "adj": adj}, flush=True)
                    raise SystemExit(1)
            count += 1
        print({"scan_n": n, "trees": count, "best": report(best)}, flush=True)

    rng = random.Random(args.seed)
    for trial in range(args.random):
        n = rng.randint(args.min_random_n, args.max_random_n)
        adj = prufer_tree(n, rng)
        profile, rows = candidate_rows(adj)
        for row in rows:
            best = update_best(best, profile, row, f"prufer:{trial}")
            if row["q_nonedge"] > 0:
                print({"failure": report(best), "adj": adj}, flush=True)
                raise SystemExit(1)
        if (trial + 1) % 100 == 0:
            print({"random": trial + 1, "best": report(best)}, flush=True)

    if args.galvin_max_n:
        count = 0
        for t in range(2, 30):
            for m in range(1, 80):
                n, adj = make_T_m_t_1(m, t)
                if n > args.galvin_max_n:
                    break
                profile, rows = candidate_rows(adj)
                for row in rows:
                    best = update_best(best, profile, row, f"T({m},{t},1)")
                    if row["q_nonedge"] > 0:
                        print({"failure": report(best), "adj": adj}, flush=True)
                        raise SystemExit(1)
                count += 1
        print({"galvin": count, "best": report(best)}, flush=True)

    if args.generalized_max_n:
        count = 0
        for d in range(2, 6):
            for t in range(2, 16):
                for m in range(1, 30):
                    n, adj = make_T_m_t_d(m, t, d)
                    if n > args.generalized_max_n:
                        break
                    profile, rows = candidate_rows(adj)
                    for row in rows:
                        best = update_best(best, profile, row, f"T({m},{t},{d})")
                        if row["q_nonedge"] > 0:
                            print({"failure": report(best), "adj": adj}, flush=True)
                            raise SystemExit(1)
                    count += 1
        print({"generalized": count, "best": report(best)}, flush=True)

    print({"certificate": "no obstruction", "best": report(best)}, flush=True)


if __name__ == "__main__":
    main()
