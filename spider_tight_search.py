#!/usr/bin/env python3
"""Thorough search for tightest high-mode spiders at larger n.

Focus on configurations near the boundary: one or two long arms + many
short arms.  These are the patterns that produce gap = 0 (mode = ℓ).
Search for any violation (gap < 0, meaning mode > ℓ in a high-mode tree).
"""

import sys
sys.path.insert(0, ".")
from indpoly import _polymul, _polyadd


def path_indpoly(m):
    if m == 0:
        return [1]
    if m == 1:
        return [1, 1]
    prev2 = [1]
    prev1 = [1, 1]
    for _ in range(2, m + 1):
        curr = _polyadd(prev1, [0] + prev2)
        prev2 = prev1
        prev1 = curr
    return prev1


# Cache path polynomials
_path_cache = {}
def cached_path_indpoly(m):
    if m not in _path_cache:
        _path_cache[m] = path_indpoly(m)
    return _path_cache[m]


def spider_indpoly(arms):
    d = len(arms)
    if d == 0:
        return [1, 1]
    prod_excl = [1]
    for a in arms:
        prod_excl = _polymul(prod_excl, cached_path_indpoly(a))
    prod_incl = [1]
    for a in arms:
        prod_incl = _polymul(prod_incl, cached_path_indpoly(a - 1))
    return _polyadd(prod_excl, [0] + prod_incl)


def mode_of(poly):
    return max(range(len(poly)), key=lambda i: poly[i])


def main():
    print("THOROUGH SPIDER SEARCH: LOOKING FOR mode > ℓ IN HIGH-MODE TREES",
          flush=True)
    print("=" * 70, flush=True)
    print(flush=True)

    worst_gap = 999
    violations = 0

    # Strategy: for each n, try many (d, arm-config) combos
    # Key insight from Part 2: tightest cases have 1-2 long arms + unit arms
    for n in range(20, 151):
        threshold = n // 3 + 1
        best_for_n = None  # (gap, d, arms_short)

        for d in range(4, n):
            total_arm = n - 1
            if total_arm < d:
                break

            # Config 1: one long arm + (d-1) unit arms
            L = total_arm - (d - 1)
            if L >= 1:
                arms = [L] + [1] * (d - 1)
                poly = spider_indpoly(arms)
                m = mode_of(poly)
                if m > threshold:
                    gap = d - m
                    if best_for_n is None or gap < best_for_n[0]:
                        best_for_n = (gap, d, f"S({L},1^{d-1})")
                    if gap < 0:
                        violations += 1
                        print(f"  *** VIOLATION n={n}: ℓ={d}, mode={m}, "
                              f"gap={gap}, S({L},1^{d-1})", flush=True)

            # Config 2: two equal long arms + (d-2) unit arms
            remaining = total_arm - (d - 2)
            if remaining >= 2 and d >= 3:
                half = remaining // 2
                other = remaining - half
                arms = sorted([half, other] + [1] * (d - 2), reverse=True)
                if arms[-1] >= 1:
                    poly = spider_indpoly(arms)
                    m = mode_of(poly)
                    if m > threshold:
                        gap = d - m
                        if best_for_n is None or gap < best_for_n[0]:
                            best_for_n = (gap, d,
                                         f"S({half},{other},1^{d-2})")
                        if gap < 0:
                            violations += 1
                            print(f"  *** VIOLATION n={n}: ℓ={d}, mode={m}, "
                                  f"gap={gap}, S({half},{other},1^{d-2})",
                                  flush=True)

            # Config 3: three long arms + (d-3) unit arms
            remaining = total_arm - (d - 3)
            if remaining >= 3 and d >= 4:
                third = remaining // 3
                rem = remaining - 3 * third
                a1 = third + (1 if rem >= 1 else 0)
                a2 = third + (1 if rem >= 2 else 0)
                a3 = third
                arms = sorted([a1, a2, a3] + [1] * (d - 3), reverse=True)
                if arms[-1] >= 1:
                    poly = spider_indpoly(arms)
                    m = mode_of(poly)
                    if m > threshold:
                        gap = d - m
                        if best_for_n is None or gap < best_for_n[0]:
                            best_for_n = (gap, d,
                                         f"S({a1},{a2},{a3},1^{d-3})")
                        if gap < 0:
                            violations += 1

            # Config 4: one long arm + (d-1) arms of length 2
            if total_arm >= 2 * (d - 1) + 1 and d >= 2:
                L = total_arm - 2 * (d - 1)
                if L >= 1:
                    arms = [L] + [2] * (d - 1)
                    poly = spider_indpoly(arms)
                    m = mode_of(poly)
                    if m > threshold:
                        gap = d - m
                        if best_for_n is None or gap < best_for_n[0]:
                            best_for_n = (gap, d, f"S({L},2^{d-1})")
                        if gap < 0:
                            violations += 1

            # Config 5: uniform arms
            if total_arm % d == 0:
                a = total_arm // d
                arms = [a] * d
                poly = spider_indpoly(arms)
                m = mode_of(poly)
                if m > threshold:
                    gap = d - m
                    if best_for_n is None or gap < best_for_n[0]:
                        best_for_n = (gap, d, f"S({a}^{d})")
                    if gap < 0:
                        violations += 1

        if best_for_n:
            gap, d, desc = best_for_n
            if gap < worst_gap:
                worst_gap = gap
            marker = " <<<" if gap <= 0 else ""
            if n <= 65 or gap <= 1:
                print(f"n={n:4d}: ℓ={d:3d}, mode={d-gap:3d}, gap={gap:3d}, "
                      f"{desc}{marker}", flush=True)

    print(flush=True)
    print(f"Worst gap across all n: {worst_gap}", flush=True)
    print(f"Total violations (gap < 0): {violations}", flush=True)
    print(flush=True)

    # Deep dive into gap=0 cases
    print("DEEP DIVE: gap=0 CASES", flush=True)
    print("-" * 70, flush=True)

    for n in range(40, 101):
        threshold = n // 3 + 1
        gap0_cases = []

        for d in range(4, n):
            total_arm = n - 1
            if total_arm < d:
                break

            # Try all "interesting" configs at this d
            configs = []

            # One long arm
            L = total_arm - (d - 1)
            if L >= 1:
                configs.append([L] + [1] * (d - 1))

            # Two equal/near-equal long arms
            remaining = total_arm - (d - 2)
            if remaining >= 2 and d >= 3:
                half = remaining // 2
                other = remaining - half
                configs.append(sorted([half, other] + [1] * (d - 2),
                                     reverse=True))

            for arms in configs:
                if any(a < 1 for a in arms):
                    continue
                poly = spider_indpoly(arms)
                m = mode_of(poly)
                if m > threshold and d - m == 0:
                    gap0_cases.append((d, arms[:4], m))

        for d, arms_show, m in gap0_cases:
            print(f"n={n:3d}: ℓ=mode={m:3d}, arms={arms_show}..., thr={threshold}",
                  flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
