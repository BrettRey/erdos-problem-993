"""Asymptotic bound certificate scaffold for S(2^a,8).

This is not a complete proof by itself.  It certifies the numerical thresholds
for the intended proof split:

* exact finite Route-2 check for 1 <= a <= 200;
* for a >= 200, crude exponential bounds are already far below the polynomial
  reserves left by the hub-off calculation.

The missing formal lemmas are listed in the accompanying note.
"""

from __future__ import annotations

from fractions import Fraction

from route2_spider_lane_scan import path_polys, route2_for_removed_arm


THRESHOLD = 200


def decreasing_ratio(power: int, a: int = THRESHOLD) -> float:
    return ((a + 1) / a) ** power / (2**0.5)


def main() -> None:
    paths = path_polys(8)
    min_slack = None
    failures = []
    for a in range(1, THRESHOLD + 1):
        rec = route2_for_removed_arm({2: a, 8: 1}, 2, paths)
        slack = rec["route2_slack"]
        if slack <= 0:
            failures.append((a, slack))
        if min_slack is None or slack < min_slack[1]:
            min_slack = (a, slack, rec["m"])

    # Mode localization: G_max/F_m <= 918 (2/3)^(a/3).  Comparing
    # C_m-C_k=(F_m-F_k)+(G_m-G_k) requires |G_m-G_k| <= 2G_max.
    # For a=200 we upper-bound (2/3)^(200/3) by (2/3)^66.
    single_tail_bound = Fraction(918) * Fraction(2, 3) ** 66
    mode_perturb_bound = 2 * single_tail_bound
    domination_reserve = Fraction(1, 10 * THRESHOLD)

    # Route-2 perturbation: full fugacity shift term plus the hub-on mixture
    # remaining in B.  These are deliberately loose.
    lambda_shift_bound = Fraction(10_000 * THRESHOLD**8, 2 ** (THRESHOLD // 2))
    hub_on_mixture_bound = Fraction(10_000 * THRESHOLD**3 * 3 ** (THRESHOLD - 1), 4 ** (THRESHOLD - 1))
    route2_reserve = Fraction(1, 4 * THRESHOLD)

    print(f"finite exact Route-2 range: a=1..{THRESHOLD}")
    print(f"finite failures: {len(failures)}")
    print(
        "minimum finite slack: "
        f"a={min_slack[0]} slack={float(min_slack[1]):.12g} m={min_slack[2]}"
    )
    print()
    print(f"mode localization crude bound at a={THRESHOLD}:")
    print(
        "  918(2/3)^66 = "
        f"{single_tail_bound} ~= {float(single_tail_bound):.12g}"
    )
    print(
        "  2*918(2/3)^66 = "
        f"{mode_perturb_bound} ~= {float(mode_perturb_bound):.12g}"
    )
    print(f"  comparison reserve 1/(10a) = {float(domination_reserve):.12g}")
    print(f"  ratio for next a of a(2/3)^(a/3): {((THRESHOLD + 1) / THRESHOLD) * (2 / 3) ** (1 / 3):.12g}")
    print()
    print(f"perturbation crude bound at a={THRESHOLD}:")
    print(
        "  lambda shift term 10000 a^8 / 2^(a/2) = "
        f"{lambda_shift_bound} ~= {float(lambda_shift_bound):.12g}"
    )
    print(
        "  hub-on mixture term 10000 a^3 (3/4)^(a-1) = "
        f"{hub_on_mixture_bound} ~= {float(hub_on_mixture_bound):.12g}"
    )
    perturb_bound = lambda_shift_bound + hub_on_mixture_bound
    print(f"  total perturbation bound ~= {float(perturb_bound):.12g}")
    print(f"  comparison reserve 1/(4a) = {float(route2_reserve):.12g}")
    print(f"  ratio for next a, power=8: {decreasing_ratio(8):.12g}")
    print(f"  ratio for next a, hub-on term: {((THRESHOLD + 1) / THRESHOLD) ** 3 * (3 / 4):.12g}")

    assert not failures
    assert mode_perturb_bound < domination_reserve
    assert perturb_bound < route2_reserve
    assert ((THRESHOLD + 1) / THRESHOLD) * (2 / 3) ** (1 / 3) < 1
    assert decreasing_ratio(8) < 1
    assert ((THRESHOLD + 1) / THRESHOLD) ** 3 * (3 / 4) < 1


if __name__ == "__main__":
    main()
