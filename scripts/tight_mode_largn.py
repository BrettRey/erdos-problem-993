#!/usr/bin/env python3
"""Check key properties for tight cases at larger n (13-16).

Since we don't store individual cases for n > 12, re-enumerate and check:
1. Is g[d_f] == g[d_f - 1] (flat-top) always?
2. Is alpha(T) == alpha(T-w) always?
3. Is w always non-leaf?
4. Is w always the max-degree vertex?
5. What fraction are palindromic?
"""

from __future__ import annotations

import sys
from collections import Counter

sys.path.insert(0, "/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993")
from indpoly import independence_poly, _polymul, _polyadd
from trees import trees


def mode_index(seq: list[int]) -> int:
    if not seq:
        return -1
    maxv = max(seq)
    idx = -1
    for i, v in enumerate(seq):
        if v == maxv:
            idx = i
    return idx


def first_descent(seq: list[int]) -> int:
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            return i
    return -1


def forest_poly_from_adj(adj: list[list[int]], removed: set[int]) -> list[int]:
    n = len(adj)
    remaining = [v for v in range(n) if v not in removed]
    if not remaining:
        return [1]
    rem_set = set(remaining)
    seen = set()
    out = [1]
    for start in remaining:
        if start in seen:
            continue
        comp = []
        stack = [start]
        seen.add(start)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v in rem_set and v not in seen:
                    seen.add(v)
                    stack.append(v)
        mapping = {old: i for i, old in enumerate(comp)}
        cadj: list[list[int]] = [[] for _ in range(len(comp))]
        for old in comp:
            ni = mapping[old]
            for v in adj[old]:
                j = mapping.get(v)
                if j is not None:
                    cadj[ni].append(j)
        poly_c = independence_poly(len(comp), cadj)
        out = _polymul(out, poly_c)
    return out


def main():
    for n in range(13, 17):
        print(f"\n{'='*60}")
        print(f"n = {n}")
        print(f"{'='*60}")

        tree_count = 0
        tight_count = 0
        palindromic = 0
        flat_top = 0
        alpha_preserved = 0
        w_is_leaf_count = 0
        w_is_max_deg = 0
        deg_w_dist: Counter = Counter()
        mode_g_h_dist: Counter = Counter()

        for _, adj in trees(n, backend="auto"):
            tree_count += 1
            f = independence_poly(n, adj)
            d_f = first_descent(f)
            if d_f == -1:
                continue

            degs = [len(adj[v]) for v in range(n)]
            max_deg = max(degs)

            for w in range(n):
                g = forest_poly_from_adj(adj, {w})
                mode_g = mode_index(g)

                if mode_g == d_f:
                    tight_count += 1

                    # Check properties
                    if g == g[::-1]:
                        palindromic += 1

                    if d_f < len(g) and d_f >= 1:
                        if g[d_f] == g[d_f - 1]:
                            flat_top += 1

                    alpha_T = len(f) - 1
                    alpha_Tw = len(g) - 1
                    if alpha_T == alpha_Tw:
                        alpha_preserved += 1

                    deg_w = degs[w]
                    if deg_w == 1:
                        w_is_leaf_count += 1
                    if deg_w == max_deg:
                        w_is_max_deg += 1

                    deg_w_dist[deg_w] += 1

                    # Compute mode(I(T - N[w]))
                    closed_nb = {w} | set(adj[w])
                    h = forest_poly_from_adj(adj, closed_nb)
                    mode_h = mode_index(h)
                    mode_g_h_dist[mode_g - mode_h] += 1

        print(f"Trees: {tree_count}")
        print(f"Tight cases: {tight_count}")
        if tight_count > 0:
            print(f"  Palindromic:       {palindromic} ({palindromic/tight_count*100:.1f}%)")
            print(f"  Flat-top:          {flat_top} ({flat_top/tight_count*100:.1f}%)")
            print(f"  alpha preserved:   {alpha_preserved} ({alpha_preserved/tight_count*100:.1f}%)")
            print(f"  w is leaf:         {w_is_leaf_count}")
            print(f"  w is max-deg:      {w_is_max_deg} ({w_is_max_deg/tight_count*100:.1f}%)")
            print(f"  deg(w) dist:       {dict(sorted(deg_w_dist.items()))}")
            print(f"  mode(g)-mode(h):   {dict(sorted(mode_g_h_dist.items()))}")


if __name__ == "__main__":
    main()
