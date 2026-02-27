#!/usr/bin/env python3
"""Symbolic no-go identity for T0 + (X/Y kernel) + (Fork2/P4 scaffolds).

Computes exact truncated degree-5 coefficient lock:
  E := 9*i4 - 7*i5
and defines
  S := -120*E
then emits:
  - identity text
  - JSON term list for S
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Any

import sympy as sp
from sympy import binomial


def poly_list(coeffs: list[int], kmax: int) -> list[sp.Expr]:
    out = [sp.Integer(0)] * (kmax + 1)
    for i, c in enumerate(coeffs):
        if i > kmax:
            break
        out[i] = sp.Integer(c)
    return out


def mul_trunc(p: list[sp.Expr], q: list[sp.Expr], kmax: int) -> list[sp.Expr]:
    out = [sp.Integer(0)] * (kmax + 1)
    for i, pi in enumerate(p):
        if pi == 0:
            continue
        for j, qj in enumerate(q):
            if qj == 0:
                continue
            ij = i + j
            if ij > kmax:
                continue
            out[ij] += pi * qj
    return out


def pow_binom_trunc(F: list[sp.Expr], n: sp.Symbol, kmax: int) -> list[sp.Expr]:
    """Compute (1+u)^n truncated to degree kmax where F=1+u."""
    u = [sp.Integer(0)] * (kmax + 1)
    for i in range(1, kmax + 1):
        if i < len(F):
            u[i] = F[i]
    out = [sp.Integer(0)] * (kmax + 1)
    u_pow = [sp.Integer(0)] * (kmax + 1)
    u_pow[0] = sp.Integer(1)
    for k in range(0, kmax + 1):
        ck = binomial(n, k)
        for i in range(0, kmax + 1):
            if u_pow[i] != 0:
                out[i] += ck * u_pow[i]
        if k < kmax:
            u_pow = mul_trunc(u_pow, u, kmax)
    return out


def prod_trunc(polys: list[list[sp.Expr]], kmax: int) -> list[sp.Expr]:
    out = [sp.Integer(0)] * (kmax + 1)
    out[0] = sp.Integer(1)
    for p in polys:
        out = mul_trunc(out, p, kmax)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out-txt",
        default="results/family_attack_t0_xy_forkp4_no_go_identity.txt",
        help="Path for identity text output.",
    )
    ap.add_argument(
        "--out-json",
        default="results/family_attack_t0_xy_forkp4_no_go_identity_terms.json",
        help="Path for JSON terms output.",
    )
    args = ap.parse_args()

    x = sp.Symbol("x")
    a, b, c, d = sp.symbols("a b c d", integer=True, nonnegative=True)

    kmax = 5

    P0 = poly_list([1, 13, 66, 169, 235, 177, 67, 10], kmax)
    Q0 = poly_list([0, 1, 12, 58, 146, 206, 162, 65, 10], kmax)

    FX = poly_list([1, 6, 10, 5], kmax)
    GX = poly_list([1, 5, 7, 3], kmax)
    FY = poly_list([1, 8, 21, 20, 5], kmax)
    GY = poly_list([1, 7, 16, 13, 3], kmax)
    FF = poly_list([1, 5, 6, 1], kmax)
    GF = poly_list([1, 4, 4], kmax)
    FP4 = poly_list([1, 4, 3], kmax)
    GP4 = poly_list([1, 3, 1], kmax)

    P = prod_trunc(
        [
            P0,
            pow_binom_trunc(FX, a, kmax),
            pow_binom_trunc(FY, b, kmax),
            pow_binom_trunc(FF, c, kmax),
            pow_binom_trunc(FP4, d, kmax),
        ],
        kmax,
    )
    Q = prod_trunc(
        [
            Q0,
            pow_binom_trunc(GX, a, kmax),
            pow_binom_trunc(GY, b, kmax),
            pow_binom_trunc(GF, c, kmax),
            pow_binom_trunc(GP4, d, kmax),
        ],
        kmax,
    )

    def i_k(k: int) -> sp.Expr:
        pk = P[k]
        pk1 = P[k - 1] if k - 1 >= 0 else sp.Integer(0)
        qk = Q[k]
        qk1 = Q[k - 1] if k - 1 >= 0 else sp.Integer(0)
        return sp.simplify(pk + 2 * pk1 + qk + qk1)

    i4 = i_k(4)
    i5 = i_k(5)
    E = sp.expand(9 * i4 - 7 * i5)
    S = sp.expand(-120 * E)

    poly = sp.Poly(S, a, b, c, d)
    terms = []
    for exp, coeff in poly.terms():
        term = {
            "exp": {"a": int(exp[0]), "b": int(exp[1]), "c": int(exp[2]), "d": int(exp[3])},
            "coeff": str(int(coeff)),
        }
        terms.append(term)

    all_positive = all(int(t["coeff"]) > 0 for t in terms)
    has_constant = any(
        t["exp"]["a"] == 0 and t["exp"]["b"] == 0 and t["exp"]["c"] == 0 and t["exp"]["d"] == 0
        for t in terms
    )
    min_coeff = min(int(t["coeff"]) for t in terms) if terms else 0

    text = []
    text.append("Family-level no-go identity for T0 + X/Y kernel + Fork2/P4 scaffolds.")
    text.append("")
    text.append("Symbols: a,b,c,d are nonnegative integers.")
    text.append("Define P,Q and I=(1+2x)P+(1+x)Q as in the family construction.")
    text.append("")
    text.append(f"9*i4 - 7*i5 = {sp.factor(E)}")
    text.append("")
    text.append(f"S := -120*(9*i4 - 7*i5) = {S}")
    text.append("")
    text.append(
        "Checks: "
        f"num_terms={len(terms)}, min_coeff={min_coeff}, "
        f"all_positive_nonconstant={all_positive and not has_constant}, "
        f"has_constant_term={has_constant}"
    )
    text_blob = "\n".join(text) + "\n"

    out_txt = args.out_txt
    out_json = args.out_json
    if out_txt:
        out_dir = os.path.dirname(out_txt)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(out_txt, "w", encoding="utf-8") as f:
            f.write(text_blob)
    if out_json:
        out_dir = os.path.dirname(out_json)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        payload: dict[str, Any] = {
            "identity": str(sp.factor(E)),
            "S_terms": terms,
            "summary": {
                "num_terms": len(terms),
                "min_coeff": min_coeff,
                "all_positive_nonconstant": bool(all_positive and not has_constant),
                "has_constant_term": bool(has_constant),
            },
        }
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)

    print(
        "done "
        f"num_terms={len(terms)} min_coeff={min_coeff} "
        f"all_positive_nonconstant={all_positive and not has_constant}",
        flush=True,
    )
    if out_txt:
        print(f"wrote {out_txt}", flush=True)
    if out_json:
        print(f"wrote {out_json}", flush=True)


if __name__ == "__main__":
    main()
