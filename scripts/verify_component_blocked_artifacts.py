#!/usr/bin/env python3
"""Verify component-level BLOCKED artifacts reported by 5.2 instances.

Verifies two fixed-lambda aggregate collisions:

Case A:
  lambda = 1
  same (d, mu1, mu2, rho, sigma) but different N
  using custom rooted component types C0..C8 and two explicit multisets.

Case C:
  lambda = 1/2
  same (d, mu1, mu2, rho, sigma) but different N
  using rooted-path component types and two explicit multisets.
"""

from __future__ import annotations

import argparse
import json
import os
from collections import defaultdict
from fractions import Fraction
from typing import Any


def frac_pair(x: Fraction) -> list[int]:
    return [x.numerator, x.denominator]


def factor_int(n: int) -> dict[int, int]:
    out: dict[int, int] = defaultdict(int)
    x = n
    p = 2
    while p * p <= x:
        while x % p == 0:
            out[p] += 1
            x //= p
        p += 1 if p == 2 else 2
    if x > 1:
        out[x] += 1
    return dict(out)


def add_factor_exp(exp: dict[int, int], val: Fraction, mult: int) -> None:
    fn = factor_int(val.numerator)
    fd = factor_int(val.denominator)
    for p, e in fn.items():
        exp[p] = exp.get(p, 0) + mult * e
    for p, e in fd.items():
        exp[p] = exp.get(p, 0) - mult * e


def norm_exp(exp: dict[int, int]) -> dict[int, int]:
    return {p: e for p, e in sorted(exp.items()) if e != 0}


def aggregate_from_types(
    lam: Fraction,
    types: dict[Any, dict[str, Any]],
    counts: dict[Any, int],
) -> dict[str, Any]:
    d = 0
    mu1 = Fraction(0, 1)
    csum = Fraction(0, 1)  # sum (b - a^2)
    gsum = Fraction(0, 1)
    n_total = 0

    rho_exp: dict[int, int] = {}
    add_factor_exp(rho_exp, lam, 1)  # rho includes lambda prefactor

    for t, cnt in counts.items():
        rec = types[t]
        a = rec["a"]
        b = rec.get("b")
        c = rec.get("c")
        r = rec["r"]
        g = rec["g"]
        d += cnt * rec["d"]
        mu1 += cnt * a
        if c is None:
            assert b is not None
            c = b - a * a
        csum += cnt * c
        gsum += cnt * g
        n_total += cnt * rec["n"]
        add_factor_exp(rho_exp, r, cnt)

    mu2 = mu1 * mu1 + csum
    sigma_over_rho = Fraction(1, 1) + gsum
    return {
        "d": d,
        "mu1": mu1,
        "mu2": mu2,
        "rho_exp": norm_exp(rho_exp),
        "gsum": gsum,
        "sigma_over_rho": sigma_over_rho,
        "N": n_total,
    }


def build_case_a() -> dict[str, Any]:
    lam = Fraction(1, 1)
    types = {
        "C0": {"n": 2, "d": 1, "a": Fraction(2, 3), "b": Fraction(0, 1), "r": Fraction(2, 3), "g": Fraction(1, 2)},
        "C1": {"n": 3, "d": 2, "a": Fraction(1, 1), "b": Fraction(2, 5), "r": Fraction(3, 5), "g": Fraction(2, 3)},
        "C2": {"n": 4, "d": 2, "a": Fraction(5, 4), "b": Fraction(3, 4), "r": Fraction(5, 8), "g": Fraction(1, 1)},
        "C3": {"n": 4, "d": 2, "a": Fraction(5, 4), "b": Fraction(3, 4), "r": Fraction(3, 4), "g": Fraction(7, 6)},
        "C4": {"n": 5, "d": 3, "a": Fraction(23, 14), "b": Fraction(12, 7), "r": Fraction(4, 7), "g": Fraction(5, 4)},
        "C5": {"n": 6, "d": 3, "a": Fraction(38, 21), "b": Fraction(44, 21), "r": Fraction(5, 7), "g": Fraction(5, 3)},
        "C6": {"n": 6, "d": 3, "a": Fraction(38, 21), "b": Fraction(44, 21), "r": Fraction(16, 21), "g": Fraction(7, 4)},
        "C7": {"n": 7, "d": 4, "a": Fraction(13, 6), "b": Fraction(10, 3), "r": Fraction(5, 6), "g": Fraction(13, 6)},
        "C8": {"n": 7, "d": 4, "a": Fraction(13, 6), "b": Fraction(10, 3), "r": Fraction(2, 3), "g": Fraction(23, 12)},
    }
    counts_a = {"C0": 152, "C3": 89, "C4": 29, "C5": 161, "C8": 79}
    counts_b = {"C1": 70, "C2": 209, "C6": 190, "C7": 22}
    return {"name": "A", "lambda": lam, "types": types, "A": counts_a, "B": counts_b}


def build_case_c() -> dict[str, Any]:
    lam = Fraction(1, 2)
    types = {
        "P2t1": {"n": 2, "d": 1, "a": Fraction(1, 2), "c": Fraction(-1, 4), "r": Fraction(3, 4), "g": Fraction(1, 3)},
        "P3t1": {"n": 3, "d": 2, "a": Fraction(8, 11), "c": Fraction(-42, 121), "r": Fraction(8, 11), "g": Fraction(1, 2)},
        "P3t2": {"n": 3, "d": 2, "a": Fraction(8, 11), "c": Fraction(-42, 121), "r": Fraction(9, 11), "g": Fraction(2, 3)},
        "P4t1": {"n": 4, "d": 2, "a": Fraction(14, 15), "c": Fraction(-106, 225), "r": Fraction(11, 15), "g": Fraction(8, 11)},
        "P4t2": {"n": 4, "d": 2, "a": Fraction(14, 15), "c": Fraction(-106, 225), "r": Fraction(4, 5), "g": Fraction(5, 6)},
        "P5t1": {"n": 5, "d": 3, "a": Fraction(47, 41), "c": Fraction(-979, 1681), "r": Fraction(30, 41), "g": Fraction(14, 15)},
        "P5t2": {"n": 5, "d": 3, "a": Fraction(47, 41), "c": Fraction(-979, 1681), "r": Fraction(33, 41), "g": Fraction(35, 33)},
        "P6t1": {"n": 6, "d": 3, "a": Fraction(19, 14), "c": Fraction(-137, 196), "r": Fraction(41, 56), "g": Fraction(47, 41)},
        "P6t2": {"n": 6, "d": 3, "a": Fraction(19, 14), "c": Fraction(-137, 196), "r": Fraction(45, 56), "g": Fraction(19, 15)},
        "P7t1": {"n": 7, "d": 4, "a": Fraction(80, 51), "c": Fraction(-2116, 2601), "r": Fraction(112, 153), "g": Fraction(19, 14)},
        "P7t2": {"n": 7, "d": 4, "a": Fraction(80, 51), "c": Fraction(-2116, 2601), "r": Fraction(41, 51), "g": Fraction(182, 123)},
    }
    counts_a = {
        "P3t2": 9489887,
        "P4t1": 18979774,
        "P4t2": 11368676,
        "P5t1": 41890196,
        "P7t2": 21031633,
    }
    counts_b = {
        "P2t1": 9489887,
        "P3t1": 20858563,
        "P5t2": 30348450,
        "P6t1": 9489887,
        "P6t2": 11541746,
        "P7t1": 21031633,
    }
    return {"name": "C", "lambda": lam, "types": types, "A": counts_a, "B": counts_b}


def verify_case(case: dict[str, Any]) -> dict[str, Any]:
    agg_a = aggregate_from_types(case["lambda"], case["types"], case["A"])
    agg_b = aggregate_from_types(case["lambda"], case["types"], case["B"])
    equal = {
        "d": agg_a["d"] == agg_b["d"],
        "mu1": agg_a["mu1"] == agg_b["mu1"],
        "mu2": agg_a["mu2"] == agg_b["mu2"],
        "rho_exp": agg_a["rho_exp"] == agg_b["rho_exp"],
        "sigma_over_rho": agg_a["sigma_over_rho"] == agg_b["sigma_over_rho"],
        "N": agg_a["N"] == agg_b["N"],
    }
    return {
        "case": case["name"],
        "lambda": frac_pair(case["lambda"]),
        "equal_checks": equal,
        "A": {
            "d": agg_a["d"],
            "mu1": frac_pair(agg_a["mu1"]),
            "mu2": frac_pair(agg_a["mu2"]),
            "rho_exp": agg_a["rho_exp"],
            "sigma_over_rho": frac_pair(agg_a["sigma_over_rho"]),
            "N": agg_a["N"],
        },
        "B": {
            "d": agg_b["d"],
            "mu1": frac_pair(agg_b["mu1"]),
            "mu2": frac_pair(agg_b["mu2"]),
            "rho_exp": agg_b["rho_exp"],
            "sigma_over_rho": frac_pair(agg_b["sigma_over_rho"]),
            "N": agg_b["N"],
        },
        "N_diff_A_minus_B": agg_a["N"] - agg_b["N"],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Verify reported component-level blocked artifacts.")
    ap.add_argument(
        "--case",
        choices=["A", "C", "both"],
        default="both",
        help="Which reported case to verify.",
    )
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    payload: dict[str, Any] = {"scan": "verify_component_blocked_artifacts", "cases": []}
    if args.case in ("A", "both"):
        payload["cases"].append(verify_case(build_case_a()))
    if args.case in ("C", "both"):
        payload["cases"].append(verify_case(build_case_c()))

    for c in payload["cases"]:
        eq = c["equal_checks"]
        print(
            f"case={c['case']} equal(d,mu1,mu2,rho,sigma)="
            f"({eq['d']},{eq['mu1']},{eq['mu2']},{eq['rho_exp']},{eq['sigma_over_rho']}) "
            f"N_equal={eq['N']} N_diff={c['N_diff_A_minus_B']}",
            flush=True,
        )

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
