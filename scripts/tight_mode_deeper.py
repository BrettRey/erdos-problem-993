#!/usr/bin/env python3
"""Deeper analysis of tight mode cases from tight_mode_cases.json.

Focuses on:
1. Even/odd n parity effects
2. deg(w) relative to n
3. Coefficient symmetry in I(T-w) (palindromic patterns)
4. Component structure of T-w
5. Relationship between d_f and mode_f
6. How close I(T) coefficients are at the peak/descent boundary
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict

sys.path.insert(0, "/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993")


def main():
    with open("/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/tight_mode_cases.json") as f:
        data = json.load(f)

    cases = data["tight_cases_detail"]  # n <= 12 only
    by_n = data["by_n"]

    print("=" * 70)
    print("DEEPER ANALYSIS OF TIGHT MODE CASES")
    print("=" * 70)

    # 1. Even/odd parity
    print("\n--- Parity Analysis ---")
    for rec in by_n:
        n = rec["n"]
        parity = "even" if n % 2 == 0 else "odd"
        print(f"  n={n:2d} ({parity}): {rec['tight_cases']:5d} tight / {rec['vertex_cases']:6d} "
              f"= {rec['tight_fraction']:.6f}")

    odd_tight = sum(r["tight_cases"] for r in by_n if r["n"] % 2 == 1)
    even_tight = sum(r["tight_cases"] for r in by_n if r["n"] % 2 == 0)
    odd_total = sum(r["vertex_cases"] for r in by_n if r["n"] % 2 == 1)
    even_total = sum(r["vertex_cases"] for r in by_n if r["n"] % 2 == 0)
    print(f"\n  Odd n:  {odd_tight} tight / {odd_total} total = {odd_tight/odd_total:.6f}")
    print(f"  Even n: {even_tight} tight / {even_total} total = {even_tight/even_total:.6f}")

    # 2. Detailed case analysis (n <= 12)
    print(f"\n--- Detailed Case Analysis (n <= 12, {len(cases)} cases) ---")

    # Check palindromic / nearly-palindromic I(T-w) coefficients
    print("\n  Palindromic I(T-w)?")
    palindrome_count = 0
    for c in cases:
        g = c["coeffs_g"]
        if g == g[::-1]:
            palindrome_count += 1
            print(f"    n={c['n']} g6={c['graph6']} w={c['w']}: I(T-w) = {g} [PALINDROMIC]")
    print(f"  Total palindromic: {palindrome_count} / {len(cases)}")

    # Check: what are f[d_f-1] and f[d_f]?  How close is the descent?
    print("\n  Descent margin in I(T) at d_f:")
    for c in cases:
        f = c["coeffs_f"]
        d = c["d_f"]
        if d >= 1 and d < len(f):
            margin = f[d-1] - f[d]
            ratio = f[d] / f[d-1] if f[d-1] > 0 else 0
            print(f"    n={c['n']} g6={c['graph6']} w={c['w']}: "
                  f"f[{d-1}]={f[d-1]}, f[{d}]={f[d]}, margin={margin}, "
                  f"ratio={ratio:.6f}")

    # Check: What are g[d_f] and g[d_f-1]?  The tight condition means
    # mode(g) = d_f, so g[d_f] >= g[k] for all k.
    print("\n  I(T-w) at mode = d_f:")
    for c in cases:
        g = c["coeffs_g"]
        d = c["d_f"]
        if d < len(g) and d >= 1:
            # g[d] should be the max
            g_at_d = g[d]
            g_at_dm1 = g[d-1] if d-1 >= 0 else 0
            g_at_dp1 = g[d+1] if d+1 < len(g) else 0
            print(f"    n={c['n']} g6={c['graph6']} w={c['w']}: "
                  f"g[{d-1}]={g_at_dm1}, g[{d}]={g_at_d}, g[{d+1}]={g_at_dp1}, "
                  f"g[d]/g[d-1]={g_at_d/g_at_dm1:.4f}" + (f", g[d]=g[d-1]" if g_at_d == g_at_dm1 else ""))

    # Check: Is g[d_f] == g[d_f-1] (i.e., flat top at the boundary)?
    print("\n  Flat-top count (g[d_f] == g[d_f-1]):")
    flat_top = 0
    for c in cases:
        g = c["coeffs_g"]
        d = c["d_f"]
        if d < len(g) and d >= 1:
            if g[d] == g[d-1]:
                flat_top += 1
    print(f"    {flat_top} / {len(cases)} cases have g[d_f] == g[d_f-1]")

    # 3. Component structure analysis
    print("\n  Component structure of T-w in tight cases:")
    comp_patterns: Counter = Counter()
    for c in cases:
        pattern = tuple(c["comp_sizes_tw"])
        comp_patterns[pattern] += 1
    for pattern, count in comp_patterns.most_common(20):
        print(f"    sizes={list(pattern)}: {count}")

    # 4. deg(w) / n ratio
    print("\n  deg(w) / n ratio in tight cases:")
    for c in cases:
        ratio = c["deg_w"] / c["n"]
        print(f"    n={c['n']} deg(w)={c['deg_w']}: ratio={ratio:.4f}")

    # 5. Does w always have the highest degree in T?
    print("\n  Is w the highest-degree vertex?")
    w_is_max_deg = 0
    for c in cases:
        max_deg = c["deg_seq"][0]  # sorted descending
        if c["deg_w"] == max_deg:
            w_is_max_deg += 1
    print(f"    {w_is_max_deg} / {len(cases)} cases")

    # 6. d_f = mode_f + 1 always?
    print("\n  d(f) - mode(f) values:")
    diffs: Counter = Counter()
    for c in cases:
        diffs[c["d_f"] - c["mode_f"]] += 1
    for diff, count in sorted(diffs.items()):
        print(f"    d_f - mode_f = {diff}: {count}")

    # 7. Is mode_g actually the last max or just one max?
    print("\n  Checking if g has a unique mode at d_f or a plateau:")
    plateau_count = 0
    for c in cases:
        g = c["coeffs_g"]
        d = c["d_f"]
        maxg = max(g)
        indices_at_max = [i for i, v in enumerate(g) if v == maxg]
        if len(indices_at_max) > 1:
            plateau_count += 1
            print(f"    n={c['n']} g6={c['graph6']} w={c['w']}: "
                  f"max at indices {indices_at_max}, g={g}")
    print(f"  Plateau (non-unique mode): {plateau_count} / {len(cases)}")


if __name__ == "__main__":
    main()
