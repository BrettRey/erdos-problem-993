"""Domination diagnostics for proving Route-2 on S(2^a,8).

The intended proof splits

    I_{a,8}(x) = F_a(x) + G_a(x)
             = P_8(x)(1+2x)^a + x P_7(x)(1+x)^a.

F_a is the hub-off term and G_a is exponentially smaller near the target mode.
This probe checks the quantitative conditions behind the prospective proof:

* the proposed mode matches the actual full mode;
* G_a has already peaked before the proposed mode;
* after the G_a peak, the ratio G_a/F_a is tiny enough that F_a's increasing
  margin dominates any negative contribution from G_a.

It is diagnostic evidence for the fixed-r proof, not the proof itself.
"""

from __future__ import annotations

import argparse
import json
from math import comb
from pathlib import Path

from route2_spider_lane_scan import (
    path_polys,
    route2_float_for_removed_arm,
    spider_poly_from_counts,
)


P8 = [1, 8, 21, 20, 5]
P7 = [1, 7, 15, 10, 1]
SHIFTS = {0: 9, 1: 7, 2: 8}


def proposed_mode(a: int) -> int:
    return (2 * a + SHIFTS[a % 3]) // 3


def coeff(poly: list[int], base_coeff: int, a: int, k: int) -> int:
    total = 0
    for j, pj in enumerate(poly):
        kk = k - j
        if 0 <= kk <= a:
            total += pj * (base_coeff**kk) * comb(a, kk)
    return total


def f_coeff(a: int, k: int) -> int:
    return coeff(P8, 2, a, k)


def g_coeff(a: int, k: int) -> int:
    total = 0
    for j, qj in enumerate(P7):
        kk = k - 1 - j
        if 0 <= kk <= a:
            total += qj * comb(a, kk)
    return total


def mode(poly: list[int]) -> int:
    mx = max(poly)
    return next(i for i, v in enumerate(poly) if v == mx)


def diagnostics(a: int, paths: list[list[int]]) -> dict:
    m = proposed_mode(a)
    full_poly = spider_poly_from_counts({2: a, 8: 1}, paths)
    full_mode = mode(full_poly)
    g_values = [g_coeff(a, k) for k in range(a + 6)]
    g_mode = mode(g_values)

    delta_left = 1.0 - f_coeff(a, m - 1) / f_coeff(a, m)
    tail_ratios = [
        (g_coeff(a, k) / f_coeff(a, k), k)
        for k in range(g_mode + 1, m + 1)
        if f_coeff(a, k) > 0
    ]
    max_tail_ratio, max_tail_k = max(tail_ratios) if tail_ratios else (0.0, None)
    route2 = route2_float_for_removed_arm({2: a, 8: 1}, 2, paths)

    return {
        "a": a,
        "proposed_mode": m,
        "full_mode": full_mode,
        "g_mode": g_mode,
        "delta_left": delta_left,
        "max_tail_ratio": max_tail_ratio,
        "max_tail_ratio_k": max_tail_k,
        "tail_domination_margin": delta_left - 2 * max_tail_ratio,
        "route2_slack": float(route2["route2_slack"]),
        "mode_ok": full_mode == m,
        "g_peaked_before_mode": g_mode < m,
        "tail_domination_ok": delta_left > 2 * max_tail_ratio,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="S(2^a,8) domination diagnostics.")
    ap.add_argument("--a-max", type=int, default=250)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    paths = path_polys(8)
    rows = [diagnostics(a, paths) for a in range(1, args.a_max + 1)]
    failures = [
        r
        for r in rows
        if not (r["mode_ok"] and r["g_peaked_before_mode"] and r["tail_domination_ok"])
    ]
    min_route2 = min(rows, key=lambda r: r["route2_slack"])
    min_tail = min(rows, key=lambda r: r["tail_domination_margin"])

    print(f"checked a=1..{args.a_max}")
    print(f"diagnostic failures: {len(failures)}")
    if failures:
        print("first failures:")
        for row in failures[:10]:
            print(row)
    print(
        "min Route-2 slack: "
        f"a={min_route2['a']} slack={min_route2['route2_slack']:.12g}"
    )
    print(
        "min tail-domination margin: "
        f"a={min_tail['a']} margin={min_tail['tail_domination_margin']:.12g}"
    )

    for a in [12, 30, 50, 81, 100, 200, args.a_max]:
        if 1 <= a <= args.a_max:
            row = rows[a - 1]
            print(
                f"a={a:4d} m={row['proposed_mode']:4d} "
                f"gmode={row['g_mode']:4d} "
                f"tailE={row['max_tail_ratio']:.3e} "
                f"deltaL={row['delta_left']:.3e} "
                f"margin={row['tail_domination_margin']:.3e} "
                f"route2={row['route2_slack']:.12g}"
            )

    if args.out:
        path = Path(args.out)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"rows": rows}, indent=2))
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
