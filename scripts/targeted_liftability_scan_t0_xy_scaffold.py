#!/usr/bin/env python3
"""Targeted liftability scan for the T0 + (X/Y kernel) + (E/P3 scaffold) family.

Family (all attached at canonical root u of fixed core T0):
  - Kernel children count t, split as b copies of Y and (t-b) copies of X
  - Scaffold children cE copies of E and dP3 copies of P3

Polynomials:
  H1 = 1+x
  H2 = 1+3x+x^2
  F_red = 1+5x+5x^2
  G_red = 1+4x+3x^2

  Hb = H1^(t-b) * H2^b
  Pbase = P0 * F_red^t * (1+2x)^cE * (1+3x+x^2)^dP3
  Qbase = Q0 * G_red^t * (1+x)^cE * (1+2x)^dP3
  I = Hb * ((1+2x)Pbase + (1+x)Qbase)

Collision criterion (liftability trigger at fixed t,cE,dP3):
  two different b values with same derived (m,lambda), m>=m_min.
Then rho is automatically equal (b-independent), while N differs by 2*delta_b.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from fractions import Fraction
from typing import Any


def poly_add(a: list[int], b: list[int]) -> list[int]:
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        if i < len(a):
            out[i] += a[i]
        if i < len(b):
            out[i] += b[i]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0:
                continue
            out[i + j] += ai * bj
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_pow_table(base: list[int], max_exp: int) -> list[list[int]]:
    tbl: list[list[int]] = [[1]]
    for _ in range(max_exp):
        tbl.append(poly_mul(tbl[-1], base))
    return tbl


def leftmost_mode_and_lambda(coeffs: list[int]) -> tuple[int, Fraction]:
    max_c = max(coeffs)
    m = 0
    for i, c in enumerate(coeffs):
        if c == max_c:
            m = i
            break
    if m == 0:
        return m, Fraction(0, 1)
    return m, Fraction(coeffs[m - 1], coeffs[m])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-t", type=int, default=20)
    ap.add_argument("--max-e", type=int, default=20, help="max cE")
    ap.add_argument("--max-p3", type=int, default=20, help="max dP3")
    ap.add_argument("--m-min", type=int, default=4)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    # Fixed core from verified canonical record
    P0 = [1, 13, 66, 169, 235, 177, 67, 10]
    Q0 = [0, 1, 12, 58, 146, 206, 162, 65, 10]
    I0 = [1, 16, 105, 371, 777, 999, 789, 371, 95, 10]

    H1 = [1, 1]
    H2 = [1, 3, 1]
    F_RED = [1, 5, 5]
    G_RED = [1, 4, 3]

    E_F = [1, 2]
    E_G = [1, 1]
    P3_F = [1, 3, 1]
    P3_G = [1, 2]

    t0 = time.time()

    h1_pow = poly_pow_table(H1, args.max_t)
    h2_pow = poly_pow_table(H2, args.max_t)
    fred_pow = poly_pow_table(F_RED, args.max_t)
    gred_pow = poly_pow_table(G_RED, args.max_t)
    ef_pow = poly_pow_table(E_F, args.max_e)
    eg_pow = poly_pow_table(E_G, args.max_e)
    p3f_pow = poly_pow_table(P3_F, args.max_p3)
    p3g_pow = poly_pow_table(P3_G, args.max_p3)

    combos_scanned = 0
    combos_passing_m_gate = 0
    collision_count = 0
    first_split: dict[str, Any] | None = None

    for t in range(args.max_t + 1):
        pcore_t = poly_mul(P0, fred_pow[t])
        qcore_t = poly_mul(Q0, gred_pow[t])
        for c_e in range(args.max_e + 1):
            pcore_te = poly_mul(pcore_t, ef_pow[c_e])
            qcore_te = poly_mul(qcore_t, eg_pow[c_e])
            for d_p3 in range(args.max_p3 + 1):
                pbase = poly_mul(pcore_te, p3f_pow[d_p3])
                qbase = poly_mul(qcore_te, p3g_pow[d_p3])
                ibase = poly_add(poly_mul([1, 2], pbase), poly_mul([1, 1], qbase))

                # key=(m,lam_num,lam_den) -> representative b
                seen_b: dict[tuple[int, int, int], tuple[int, int]] = {}

                for b in range(t + 1):
                    combos_scanned += 1
                    hb = poly_mul(h1_pow[t - b], h2_pow[b])
                    i_poly = poly_mul(hb, ibase)
                    m, lam = leftmost_mode_and_lambda(i_poly)

                    if m < args.m_min:
                        continue
                    combos_passing_m_gate += 1
                    key = (m, lam.numerator, lam.denominator)

                    # N = [x] (Hb * Pbase)
                    p_poly = poly_mul(hb, pbase)
                    n_val = p_poly[1] if len(p_poly) > 1 else 0

                    prev = seen_b.get(key)
                    if prev is None:
                        seen_b[key] = (b, n_val)
                    else:
                        collision_count += 1
                        if first_split is None and prev[0] != b:
                            first_split = {
                                "fixed_params": {
                                    "t": t,
                                    "cE": c_e,
                                    "dP3": d_p3,
                                },
                                "key": {
                                    "m": m,
                                    "lambda": f"{lam.numerator}/{lam.denominator}",
                                },
                                "A": {"b": prev[0], "N": prev[1]},
                                "B": {"b": b, "N": n_val},
                                "delta_N": abs(n_val - prev[1]),
                            }

    payload = {
        "family": {
            "core_g6": "O??????_A?C?E?@_WG@j?",
            "core_triplet": {"leaf": 4, "support": 12, "root": 3},
            "core_P": P0,
            "core_Q": Q0,
            "core_I": I0,
            "core_m": 5,
            "core_lambda": "7/9",
            "core_N": 13,
            "kernel_reduced": {
                "F_red": F_RED,
                "G_red": G_RED,
                "H1": H1,
                "H2": H2,
            },
            "scaffold_types": {
                "E": {"F": E_F, "G": E_G},
                "P3": {"F": P3_F, "G": P3_G},
            },
        },
        "scan_params": {
            "max_t": args.max_t,
            "max_e": args.max_e,
            "max_p3": args.max_p3,
            "m_min": args.m_min,
        },
        "totals": {
            "combos_scanned": combos_scanned,
            "combos_passing_m_gate": combos_passing_m_gate,
            "collision_count": collision_count,
            "split_found": first_split is not None,
            "runtime_seconds": time.time() - t0,
        },
        "first_split": first_split,
    }

    print(
        "done "
        f"scanned={combos_scanned} pass_m={combos_passing_m_gate} "
        f"collisions={collision_count} split_found={first_split is not None} "
        f"time={time.time()-t0:.3f}s",
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
