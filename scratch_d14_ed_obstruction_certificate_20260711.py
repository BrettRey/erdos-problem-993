#!/usr/bin/env python3
"""Replay the exact obstruction to the D14 aggregate energy bound.

The tree consists of two identical rooted spherical trees joined at their
roots.  Reading away from either central root, the child multiplicities are

    (2, 3, 2, 1, 2, 1, 1).

At prefix rank r=81 this script computes, with integer/rational arithmetic,
the two unnormalised sides D and E of

    Dir_r(e) <= E_r(e+h)/r  iff  D <= E.

The asserted exact values make the script a regression certificate rather
than a floating-point search.
"""

from __future__ import annotations

from fractions import Fraction

from indpoly import independence_poly
from scratch_spectral_energy_dp_20260711 import (
    energy_distribution,
    make_double_spherical_tree,
)


BRANCHING = (2, 3, 2, 1, 2, 1, 1)
RANK = 81

EXPECTED_E = 535256584637996471952494945973518820229192
EXPECTED_D = Fraction(
    245700013545993183538657579915009525438626090145082199405497542315007719,
    450170153991699740612130417300,
)
EXPECTED_DEFICIT = Fraction(
    -42691269729915413999161105330329696919890850296588501772067055418675071,
    50018905999077748956903379700,
)
EXPECTED_RATIO = Fraction(
    245700013545993183538657579915009525438626090145082199405497542315007719,
    240956539131558137538750790433861781336415995667683476986378980601821600,
)


def main() -> None:
    tree = make_double_spherical_tree(BRANCHING)
    polynomial = independence_poly(len(tree), tree)
    alpha = len(polynomial) - 1
    prefix_boundary = (2 * alpha + 1) // 3

    assert len(tree) == 210
    assert sum(map(len, tree)) == 2 * (len(tree) - 1)
    assert alpha == 124
    assert prefix_boundary == 83
    assert RANK == prefix_boundary - 2

    print(
        {
            "branching": BRANCHING,
            "n": len(tree),
            "alpha": alpha,
            "prefix_boundary": prefix_boundary,
            "rank": RANK,
        },
        flush=True,
    )

    distribution = energy_distribution(tree, max_rank=RANK)
    rank_counts = [0] * (RANK + 1)
    for (rank, _), moment in distribution.items():
        rank_counts[rank] += moment.count
    assert rank_counts == polynomial[: RANK + 1]

    # D = sum_{C in I_{r-1}} (S2_C - 4 h_C^2/q_C).
    d_total = Fraction(0)
    for (rank, q), moment in distribution.items():
        if rank == RANK - 1 and q:
            d_total += (
                moment.degree_square_sum
                - Fraction(4 * moment.edge_square_sum, q)
            )

    # E = sum_{A in I_r} (e(A)+h(A)).  The extension part is
    # (r+1)i_{r+1}; the marked DP supplies the residual-edge part.
    residual_edge_total = sum(
        moment.edge_sum
        for (rank, _), moment in distribution.items()
        if rank == RANK
    )
    extension_total = (RANK + 1) * polynomial[RANK + 1]
    e_total = extension_total + residual_edge_total
    deficit = RANK * (e_total - d_total)
    ratio = d_total / e_total

    assert polynomial[81] == 12150264177683772579409779766331299105518
    assert polynomial[82] == 4933589294165997470864362459890912420846
    assert polynomial[83] == 1914913869650769840173000753914910371834
    assert extension_total == 404554322121611792610877721711054818509372
    assert residual_edge_total == 130702262516384679341617224262464001719820
    assert e_total == EXPECTED_E
    assert d_total == EXPECTED_D
    assert deficit == EXPECTED_DEFICIT
    assert ratio == EXPECTED_RATIO
    assert deficit < 0
    assert ratio > 1

    # The actual prefix-GSB inequality is nevertheless safely true here.
    mu = Fraction((RANK + 1) * polynomial[RANK + 1], polynomial[RANK])
    eta = Fraction(residual_edge_total, polynomial[RANK])
    variance = (
        Fraction(
            (RANK + 2) * (RANK + 1) * polynomial[RANK + 2],
            polynomial[RANK],
        )
        + 2 * eta
        + mu
        - mu * mu
    )
    gsb_ratio = variance / (2 * (mu + eta))
    assert gsb_ratio == Fraction(
        4456090562654887266294141184427378259936953170060911002599884558253076453306371,
        20844579827552598787630694022550327446979294710701437563672195144185711388948338,
    )
    assert gsb_ratio < 1

    print(
        {
            "E": e_total,
            "D": str(d_total),
            "81(E-D)": str(deficit),
            "D/E": str(ratio),
            "D/E_float": float(ratio),
            "GSB_ratio": str(gsb_ratio),
            "GSB_ratio_float": float(gsb_ratio),
            "certificate": "passed",
        },
        flush=True,
    )


if __name__ == "__main__":
    main()
