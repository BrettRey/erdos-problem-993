#!/usr/bin/env python3
"""Exact collision checker for the realizable 3-type canonical-root family.

Family ratios at fixed lambda:
  R_E = (1+lambda)/(1+2lambda)
  R_P = (1+2lambda)/(1+3lambda+lambda^2)
  R_F = (1+2lambda)^2/(1+5lambda+6lambda^2+lambda^3)

For parameters (e,p,f) >= 0:
  rho = lambda * R_E^e * R_P^p * R_F^f
  N   = 2e + 3p + 5f

This script searches for any same-(lambda,rho) collisions with different
parameter vectors and reports whether N can split.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from fractions import Fraction
from typing import Any


def parse_lambda_list(raw: str | None, max_den: int) -> list[Fraction]:
    vals: list[Fraction] = []
    seen: set[tuple[int, int]] = set()
    if raw:
        for tok in raw.split(","):
            tok = tok.strip()
            if not tok:
                continue
            if "/" not in tok:
                raise ValueError(f"bad lambda token: {tok}")
            a, b = tok.split("/", 1)
            f = Fraction(int(a), int(b))
            if not (0 < f < 1):
                raise ValueError(f"lambda must be in (0,1): {tok}")
            key = (f.numerator, f.denominator)
            if key not in seen:
                seen.add(key)
                vals.append(f)
    else:
        for b in range(2, max_den + 1):
            for a in range(1, b):
                if math.gcd(a, b) != 1:
                    continue
                f = Fraction(a, b)
                key = (a, b)
                if key not in seen:
                    seen.add(key)
                    vals.append(f)
    vals.sort(key=lambda z: (z.denominator, z.numerator))
    return vals


def ratios(lam: Fraction) -> tuple[Fraction, Fraction, Fraction]:
    one = Fraction(1, 1)
    re = (one + lam) / (one + 2 * lam)
    rp = (one + 2 * lam) / (one + 3 * lam + lam * lam)
    rf = (one + 2 * lam) ** 2 / (one + 5 * lam + 6 * lam * lam + lam * lam * lam)
    return re, rp, rf


def run(max_e: int, max_p: int, max_f: int, lambdas: list[Fraction]) -> dict[str, Any]:
    started = time.time()
    total_triples = (max_e + 1) * (max_p + 1) * (max_f + 1)

    lambda_summaries: list[dict[str, Any]] = []
    first_collision_any: dict[str, Any] | None = None

    for lam in lambdas:
        re, rp, rf = ratios(lam)
        seen: dict[tuple[int, int], tuple[int, int, int, int]] = {}
        collisions = 0
        split_collisions = 0
        first_collision = None
        first_split = None

        for e in range(max_e + 1):
            re_pow = re**e
            for p in range(max_p + 1):
                rep_pow = re_pow * (rp**p)
                for f in range(max_f + 1):
                    rho = lam * rep_pow * (rf**f)
                    key = (rho.numerator, rho.denominator)
                    n = 2 * e + 3 * p + 5 * f

                    prev = seen.get(key)
                    if prev is None:
                        seen[key] = (e, p, f, n)
                        continue

                    collisions += 1
                    e0, p0, f0, n0 = prev
                    if first_collision is None:
                        first_collision = {
                            "A": {"e": e0, "p": p0, "f": f0, "N": n0},
                            "B": {"e": e, "p": p, "f": f, "N": n},
                            "rho": [rho.numerator, rho.denominator],
                        }
                    if (e, p, f) != (e0, p0, f0):
                        if n != n0:
                            split_collisions += 1
                            if first_split is None:
                                first_split = {
                                    "A": {"e": e0, "p": p0, "f": f0, "N": n0},
                                    "B": {"e": e, "p": p, "f": f, "N": n},
                                    "rho": [rho.numerator, rho.denominator],
                                }
                                if first_collision_any is None:
                                    first_collision_any = {
                                        "lambda": [lam.numerator, lam.denominator],
                                        "witness": first_split,
                                    }

        lam_rec = {
            "lambda": [lam.numerator, lam.denominator],
            "triples": total_triples,
            "unique_rho": len(seen),
            "collisions": collisions,
            "split_collisions": split_collisions,
            "first_collision": first_collision,
            "first_split": first_split,
        }
        lambda_summaries.append(lam_rec)

    return {
        "scan": "three_type_family_no_go_check",
        "family": {
            "types": ["E(edge endpoint)", "P(path3 endpoint)", "F(fork depth2 root)"],
            "rho_formula": "rho=lambda*RE^e*RP^p*RF^f",
            "N_formula": "N=2e+3p+5f",
        },
        "params": {
            "max_e": max_e,
            "max_p": max_p,
            "max_f": max_f,
            "lambda_count": len(lambdas),
        },
        "totals": {
            "triples_per_lambda": total_triples,
            "all_split_collisions": sum(x["split_collisions"] for x in lambda_summaries),
            "elapsed_sec": time.time() - started,
        },
        "first_split_any": first_collision_any,
        "per_lambda": lambda_summaries,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-e", type=int, default=10)
    ap.add_argument("--max-p", type=int, default=10)
    ap.add_argument("--max-f", type=int, default=10)
    ap.add_argument(
        "--lambdas",
        default="",
        help="Comma-separated rationals a/b in (0,1). If empty, use all reduced a/b with b<=max-den.",
    )
    ap.add_argument("--max-den", type=int, default=20)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    lambdas = parse_lambda_list(args.lambdas or None, args.max_den)
    payload = run(args.max_e, args.max_p, args.max_f, lambdas)

    print(
        f"done lambdas={payload['params']['lambda_count']} triples/lambda={payload['totals']['triples_per_lambda']} "
        f"split_collisions={payload['totals']['all_split_collisions']} elapsed={payload['totals']['elapsed_sec']:.2f}s"
    )
    if payload["first_split_any"] is None:
        print("first_split_any=None")
    else:
        print("first_split_any=", payload["first_split_any"])

    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
