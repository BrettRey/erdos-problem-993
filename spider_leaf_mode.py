#!/usr/bin/env python3
"""Analyze mode vs ℓ for spider trees.

A spider S(a_1, ..., a_d) has a central hub of degree d with arms of
lengths a_1 >= ... >= a_d >= 1.  n = 1 + sum(a_i), ℓ = d.

I(spider; x) = prod I(P_{a_i}; x)  +  x * prod I(P_{a_i - 1}; x)
             = [hub excluded]       +  [hub included]

where I(P_0; x) = 1 (empty path = single vertex removed = nothing left,
actually P_0 means no vertices, so I = [1]).

Since tight Leaf-Mode cases are all spiders, understanding mode/ℓ for
spiders is key to proving the Leaf-Mode Inequality.
"""

import sys
import time
from itertools import combinations_with_replacement

sys.path.insert(0, ".")
from indpoly import _polymul, _polyadd


def path_indpoly(m):
    """Independence polynomial of path P_m (m vertices, m >= 0)."""
    if m == 0:
        return [1]
    if m == 1:
        return [1, 1]
    # DP: I(P_m) = I(P_{m-1}) + x * I(P_{m-2})
    prev2 = [1]       # P_0
    prev1 = [1, 1]    # P_1
    for _ in range(2, m + 1):
        # I(P_k) = I(P_{k-1}) + x * I(P_{k-2})
        shifted = [0] + prev2  # x * I(P_{k-2})
        curr = _polyadd(prev1, shifted)
        prev2 = prev1
        prev1 = curr
    return prev1


def spider_indpoly(arms):
    """Independence polynomial of spider with given arm lengths.

    arms: list of positive integers [a_1, ..., a_d]
    Returns coefficient list.
    """
    d = len(arms)
    if d == 0:
        return [1, 1]  # single vertex (hub only)

    # Hub excluded: prod I(P_{a_i})
    prod_excl = [1]
    for a in arms:
        prod_excl = _polymul(prod_excl, path_indpoly(a))

    # Hub included: x * prod I(P_{a_i - 1})
    prod_incl = [1]
    for a in arms:
        prod_incl = _polymul(prod_incl, path_indpoly(a - 1))
    prod_incl = [0] + prod_incl  # multiply by x

    return _polyadd(prod_excl, prod_incl)


def mode_of(poly):
    """Return the mode index (position of max coefficient)."""
    return max(range(len(poly)), key=lambda i: poly[i])


def main():
    print("SPIDER TREES: MODE vs LEAVES ANALYSIS", flush=True)
    print("=" * 70, flush=True)
    print(flush=True)

    # Part 1: Uniform spiders S(a, a, ..., a) with d arms
    print("PART 1: UNIFORM SPIDERS S(a^d)", flush=True)
    print("-" * 70, flush=True)
    print(f"{'d':>4} {'a':>4} {'n':>6} {'ℓ':>4} {'mode':>5} {'ℓ-m':>5} "
          f"{'m/ℓ':>6} {'threshold':>9} {'high?':>5}", flush=True)

    for a in range(1, 8):
        for d in range(3, min(60, 200 // a)):
            n = 1 + a * d
            arms = [a] * d
            poly = spider_indpoly(arms)
            m = mode_of(poly)
            threshold = n // 3 + 1
            high = m > threshold
            gap = d - m
            ratio = m / d if d > 0 else 0

            if high or gap <= 5:
                marker = "HIGH" if high else ""
                print(f"{d:4d} {a:4d} {n:6d} {d:4d} {m:5d} {gap:5d} "
                      f"{ratio:6.3f} {threshold:9d} {marker:>5}", flush=True)

    # Part 2: Mixed-arm spiders -- find tightest mode/ℓ ratio for high-mode
    print(flush=True)
    print("PART 2: TIGHTEST HIGH-MODE SPIDERS BY n", flush=True)
    print("-" * 70, flush=True)

    for n in range(8, 61):
        best_gap = None
        best_info = None

        # Generate arm configurations: d arms summing to n-1, each >= 1
        # For efficiency, enumerate d and representative arm combos
        for d in range(3, n):
            total_arm = n - 1
            if total_arm < d:
                break
            # Generate a few representative configs for this d:
            # all-1, all-equal, and some mixed
            configs = set()

            # All equal (if divisible)
            if total_arm % d == 0:
                a = total_arm // d
                configs.add(tuple([a] * d))

            # Nearly equal
            base = total_arm // d
            extra = total_arm - base * d
            arms = sorted([base + 1] * extra + [base] * (d - extra), reverse=True)
            configs.add(tuple(arms))

            # One long arm, rest length 1
            if d >= 2:
                arms = sorted([total_arm - (d - 1)] + [1] * (d - 1), reverse=True)
                if arms[-1] >= 1 and arms[0] >= 1:
                    configs.add(tuple(arms))

            # One long arm, rest length 2
            remaining = total_arm - 2 * (d - 1)
            if remaining >= 1 and d >= 2:
                arms = sorted([remaining] + [2] * (d - 1), reverse=True)
                configs.add(tuple(arms))

            # Two long arms, rest length 1
            if d >= 3:
                half = (total_arm - (d - 2)) // 2
                other = total_arm - (d - 2) - half
                arms = sorted([half, other] + [1] * (d - 2), reverse=True)
                if arms[-1] >= 1:
                    configs.add(tuple(arms))

            for arms_tuple in configs:
                arms = list(arms_tuple)
                if any(a < 1 for a in arms):
                    continue
                if sum(arms) != total_arm:
                    continue

                poly = spider_indpoly(arms)
                m = mode_of(poly)
                threshold = n // 3 + 1

                if m > threshold:  # high-mode
                    gap = d - m  # ℓ - mode
                    if best_gap is None or gap < best_gap:
                        best_gap = gap
                        best_info = (d, arms[:5], m, threshold)

        if best_info:
            d, arms_show, m, thr = best_info
            print(f"n={n:3d}: ℓ={d:3d}, mode={m:3d}, gap={best_gap:3d}, "
                  f"arms={arms_show}, thr={thr}", flush=True)

    # Part 3: For uniform spiders, compute mode/ℓ ratio as d -> infinity
    print(flush=True)
    print("PART 3: ASYMPTOTIC mode/ℓ RATIO FOR UNIFORM SPIDERS", flush=True)
    print("-" * 70, flush=True)

    for a in range(1, 6):
        print(f"\nArm length a={a}:", flush=True)
        for d in [10, 20, 50, 100, 200]:
            n = 1 + a * d
            arms = [a] * d
            poly = spider_indpoly(arms)
            m = mode_of(poly)
            threshold = n // 3 + 1
            ratio = m / d
            high = m > threshold
            if high:
                print(f"  d={d:4d}: n={n:5d}, mode={m:4d}, ℓ={d:4d}, "
                      f"mode/ℓ={ratio:.4f}, HIGH", flush=True)
            else:
                print(f"  d={d:4d}: n={n:5d}, mode={m:4d}, ℓ={d:4d}, "
                      f"mode/ℓ={ratio:.4f}", flush=True)

    # Part 4: Path independence polynomial modes (for reference)
    print(flush=True)
    print("PART 4: PATH P_m MODES (reference)", flush=True)
    print("-" * 70, flush=True)
    for m in range(1, 31):
        poly = path_indpoly(m)
        md = mode_of(poly)
        print(f"P_{m:2d}: mode={md}, mode/m={md/m:.3f}, poly_deg={len(poly)-1}",
              flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
