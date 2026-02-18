#!/usr/bin/env python3
"""Verify Sub-claim B compensation gap on the j=1 mixed-spider branch.

For T = S(2^k,1), unit-leaf decomposition gives:
  gap_comp = beta*c2 - alpha*deficit,
where deficit=max(0,-c1), and c1,c2 are unit-leaf slacks at lambda_m(T).

In the c1<0 regime, gap_comp = margin(T), so positivity of margin proves
Sub-claim B for this branch.
"""

from __future__ import annotations

import argparse
from fractions import Fraction

from conjecture_a_mixed_spider_exact_margin import unit_leaf_decomp


def main() -> None:
    ap = argparse.ArgumentParser(description="Check Sub-claim B compensation on j=1 branch.")
    ap.add_argument("--k-max", type=int, default=300)
    args = ap.parse_args()

    if args.k_max < 1:
        raise ValueError("k-max must be >= 1")

    min_gap = None
    min_row = None
    n_neg_c1 = 0
    bad = 0

    for k in range(1, args.k_max + 1):
        d = unit_leaf_decomp(k, 1)
        if d is None:
            continue

        c1 = d["c1"]
        c2 = d["c2"]
        alpha = d["alpha"]
        beta = d["beta"]
        margin = d["margin"]
        deficit = -c1 if c1 < 0 else Fraction(0)

        gap_comp = beta * c2 - alpha * deficit

        if c1 < 0:
            n_neg_c1 += 1
            # In failing-cond1 regime, compensation gap equals margin exactly.
            if gap_comp != margin:
                raise AssertionError(
                    f"identity mismatch at k={k}: gap_comp={gap_comp}, margin={margin}"
                )

        if gap_comp <= 0:
            bad += 1

        if min_gap is None or gap_comp < min_gap:
            min_gap = gap_comp
            min_row = {
                "k": k,
                "m": d["m"],
                "lambda": d["lam"],
                "c1": c1,
                "c2": c2,
                "alpha": alpha,
                "beta": beta,
                "gap_comp": gap_comp,
                "margin": margin,
            }

    if min_gap is None or min_row is None:
        raise AssertionError("No rows checked")

    print("Sub-claim B compensation check on T=S(2^k,1)")
    print(f"k=1..{args.k_max}")
    print(f"c1<0 cases: {n_neg_c1}")
    print(f"gap_comp<=0 failures: {bad}")
    print("minimum gap witness:")
    print(
        f"  k={min_row['k']} m={min_row['m']} "
        f"lambda={float(min_row['lambda']):.12f}"
    )
    print(
        f"  c1={float(min_row['c1']):.12f} c2={float(min_row['c2']):.12f} "
        f"alpha={float(min_row['alpha']):.12f} beta={float(min_row['beta']):.12f}"
    )
    print(
        f"  gap_comp={float(min_row['gap_comp']):.12f} "
        f"margin={float(min_row['margin']):.12f}"
    )


if __name__ == "__main__":
    main()

