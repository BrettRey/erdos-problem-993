#!/usr/bin/env python3
"""Search synthetic (f,g) counterexamples to H_CW => H_R.

H_R: ratio g_k/f_k nonincreasing on k<=d(g)-2.
H_CW: covariance window condition Cov_{P_f,lambda}(k, g_k/f_k) <= 0 whenever
      mu_f(lambda) <= d(g)-1.

This script searches generic positive log-concave sequence pairs (f,g) to show
the implication can fail outside tree-realizable structure.
"""

from __future__ import annotations

import argparse
import json
import random
from typing import Any


def first_descent(seq: list[int]) -> int:
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            return i
    return len(seq)


def is_log_concave(seq: list[int]) -> bool:
    for k in range(1, len(seq) - 1):
        if seq[k] * seq[k] < seq[k - 1] * seq[k + 1]:
            return False
    return True


def mu_f(f: list[int], lam: float) -> float:
    vals = [c * (lam**k) for k, c in enumerate(f)]
    z = sum(vals)
    return sum(k * v for k, v in enumerate(vals)) / z


def cov_window_value(f: list[int], g: list[int], lam: float) -> tuple[float, float]:
    fk = [c * (lam**k) for k, c in enumerate(f)]
    zf = sum(fk)
    p = [v / zf for v in fk]
    r = [g[k] / f[k] for k in range(len(f))]
    mu = sum(k * p[k] for k in range(len(f)))
    er = sum(r[k] * p[k] for k in range(len(f)))
    ekr = sum(k * r[k] * p[k] for k in range(len(f)))
    cov = ekr - mu * er
    return mu, cov


def hr_failure_index(f: list[int], g: list[int]) -> int:
    dg = first_descent(g)
    for k in range(max(0, dg - 1)):  # k<=dg-2
        if g[k + 1] * f[k] > g[k] * f[k + 1]:
            return k
    return -1


def hcw_holds(f: list[int], g: list[int], d_g: int, lambdas: list[float], tol: float) -> bool:
    for lam in lambdas:
        mu, cov = cov_window_value(f, g, lam)
        if mu <= d_g - 1 + tol and cov > tol:
            return False
    return True


def random_logconcave_sequence(m: int, rng: random.Random) -> list[int]:
    # crude sampler: random walk + LC filter (fast enough for search scale)
    while True:
        seq = [1]
        last = rng.randint(2, 12)
        seq.append(last)
        for _ in range(2, m + 1):
            last = max(1, last + rng.randint(-5, 5))
            seq.append(last)
        if is_log_concave(seq):
            return seq


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-deg", type=int, default=4)
    ap.add_argument("--max-deg", type=int, default=8)
    ap.add_argument("--trials-per-deg", type=int, default=200000)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    lambdas = [10 ** (-6 + 12 * i / 80) for i in range(81)]

    result: dict[str, Any] = {
        "claim": "H_CW does not imply H_R for generic log-concave positive pairs",
        "seed": args.seed,
        "min_deg": args.min_deg,
        "max_deg": args.max_deg,
        "trials_per_deg": args.trials_per_deg,
        "found": False,
        "counterexample": None,
    }

    for m in range(args.min_deg, args.max_deg + 1):
        for t in range(args.trials_per_deg):
            f = random_logconcave_sequence(m, rng)
            g = [rng.randint(1, f[k]) for k in range(m + 1)]
            g[0] = 1
            if not is_log_concave(g):
                continue
            d_g = first_descent(g)
            if d_g <= 2:
                continue

            bad_k = hr_failure_index(f, g)
            if bad_k < 0:
                continue

            if hcw_holds(f, g, d_g, lambdas, args.tol):
                result["found"] = True
                result["counterexample"] = {
                    "deg": m,
                    "trial": t,
                    "f": f,
                    "g": g,
                    "d_g": d_g,
                    "hr_bad_k": bad_k,
                    "i1_f": f[1],
                    "i2_f": f[2] if len(f) > 2 else None,
                    "note": "Typically violates tree constraints (e.g., i2 must equal C(n,2)-(n-1) for trees).",
                }
                if args.out:
                    with open(args.out, "w", encoding="utf-8") as fobj:
                        json.dump(result, fobj, indent=2)
                print(json.dumps(result, indent=2))
                return

    if args.out:
        with open(args.out, "w", encoding="utf-8") as fobj:
            json.dump(result, fobj, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
