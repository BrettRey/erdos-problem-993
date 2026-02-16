#!/usr/bin/env python3
"""Classify tight cases into subcategories for proof strategy.

Two types of tight cases:
  Type A (palindromic): I(T-w) is palindromic, mode at both d_f-1 and d_f (flat top)
  Type B (non-palindromic): mode(g) = d_f but g is not palindromic

For Type A, the mechanism is clear: removing a high-degree vertex from a
"balanced" tree can produce a palindromic forest polynomial.

For Type B, what pushes the mode up to d_f?
"""

from __future__ import annotations

import json
import sys
from collections import Counter

sys.path.insert(0, "/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993")


def main():
    with open("/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/tight_mode_cases.json") as f:
        data = json.load(f)

    cases = data["tight_cases_detail"]

    print("=" * 70)
    print("CLASSIFICATION OF TIGHT CASES (n <= 12)")
    print("=" * 70)

    type_a = []  # palindromic g
    type_b = []  # non-palindromic g, but mode at d_f

    for c in cases:
        g = c["coeffs_g"]
        if g == g[::-1]:
            type_a.append(c)
        else:
            type_b.append(c)

    print(f"\nType A (palindromic I(T-w)): {len(type_a)} / {len(cases)}")
    print(f"Type B (non-palindromic):    {len(type_b)} / {len(cases)}")

    # Analyze Type A
    print("\n--- Type A Details ---")
    print("In all Type A cases, I(T-w) is palindromic, so the mode is at")
    print("floor((len-1)/2) which happens to equal d_f.")
    for c in type_a:
        g = c["coeffs_g"]
        glen = len(g) - 1  # degree of polynomial
        mid = glen / 2
        print(f"  n={c['n']} deg(w)={c['deg_w']} "
              f"alpha(T-w)={glen} mid={mid} d_f={c['d_f']} "
              f"|g|={len(g)} comp_tw={c['comp_sizes_tw']}")

    # Check: for palindromic cases, is T-w always a disjoint union of paths?
    print("\n  Component sizes in Type A:")
    for c in type_a:
        sizes = c["comp_sizes_tw"]
        n_minus_1 = c["n"] - 1
        total_comp = sum(sizes)
        print(f"    n={c['n']} comp={sizes} total={total_comp} n-1={n_minus_1} deg(w)={c['deg_w']}")

    # Analyze Type B
    print("\n--- Type B Details ---")
    print("Non-palindromic cases where mode(I(T-w)) = d(I(T)):")
    for c in type_b:
        g = c["coeffs_g"]
        f = c["coeffs_f"]
        d = c["d_f"]
        maxg = max(g)
        max_indices = [i for i, v in enumerate(g) if v == maxg]

        # How much does g[d_f] exceed g[d_f - 1]?
        g_at_d = g[d] if d < len(g) else 0
        g_at_dm1 = g[d-1] if d >= 1 else 0
        excess = g_at_d - g_at_dm1

        print(f"  n={c['n']} g6={c['graph6']} deg(w)={c['deg_w']} "
              f"d_f={d} type={c['tree_type']}")
        print(f"    f = {f}")
        print(f"    g = {g}")
        print(f"    g[{d}]-g[{d-1}] = {excess}, max_indices={max_indices}")
        print(f"    comp_tw={c['comp_sizes_tw']}")
        print()

    # Summary: what structural features distinguish Type A from Type B?
    print("\n--- Structural Comparison ---")
    for label, group in [("Type A", type_a), ("Type B", type_b)]:
        if not group:
            continue
        deg_dist = Counter(c["deg_w"] for c in group)
        type_dist = Counter(c["tree_type"] for c in group)
        center_count = sum(1 for c in group if c["w_in_center"])
        max_deg_count = sum(1 for c in group if c["deg_w"] == c["deg_seq"][0])
        avg_comp = sum(len(c["comp_sizes_tw"]) for c in group) / len(group)

        print(f"\n  {label} ({len(group)} cases):")
        print(f"    deg(w) dist: {dict(sorted(deg_dist.items()))}")
        print(f"    tree types:  {dict(sorted(type_dist.items()))}")
        print(f"    w is center: {center_count}/{len(group)}")
        print(f"    w is max-deg: {max_deg_count}/{len(group)}")
        print(f"    avg #components of T-w: {avg_comp:.1f}")

    # Check: is there a relationship between deg(w) and alpha(T)?
    # alpha(T) = len(f) - 1 = independence number
    print("\n--- deg(w) vs alpha(T) and alpha(T-w) ---")
    for c in cases:
        f = c["coeffs_f"]
        g = c["coeffs_g"]
        alpha_T = len(f) - 1
        alpha_Tw = len(g) - 1
        # alpha(T) - alpha(T-w) should be 0 or 1 (removing a vertex from a max IS)
        diff_alpha = alpha_T - alpha_Tw
        print(f"  n={c['n']} deg(w)={c['deg_w']} alpha(T)={alpha_T} "
              f"alpha(T-w)={alpha_Tw} diff={diff_alpha} "
              f"palin={'Y' if c['coeffs_g'] == c['coeffs_g'][::-1] else 'N'}")


if __name__ == "__main__":
    main()
