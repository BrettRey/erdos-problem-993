#!/usr/bin/env python3
"""Verify d_f - mode_f = 1 for all tight cases at n = 13..16."""

from __future__ import annotations
import sys
sys.path.insert(0, "/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993")
from indpoly import independence_poly, _polymul
from trees import trees
from collections import Counter


def mode_index(seq):
    if not seq:
        return -1
    maxv = max(seq)
    idx = -1
    for i, v in enumerate(seq):
        if v == maxv:
            idx = i
    return idx


def first_descent(seq):
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            return i
    return -1


def forest_poly_from_adj(adj, removed):
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
        cadj = [[] for _ in range(len(comp))]
        for old in comp:
            ni = mapping[old]
            for v in adj[old]:
                j = mapping.get(v)
                if j is not None:
                    cadj[ni].append(j)
        from indpoly import independence_poly
        poly_c = independence_poly(len(comp), cadj)
        out = _polymul(out, poly_c)
    return out


for n in range(13, 17):
    df_mf_dist = Counter()
    tight = 0
    for _, adj in trees(n, backend="auto"):
        f = independence_poly(n, adj)
        d_f = first_descent(f)
        if d_f == -1:
            continue
        mode_f = mode_index(f)
        for w in range(n):
            g = forest_poly_from_adj(adj, {w})
            mode_g = mode_index(g)
            if mode_g == d_f:
                tight += 1
                df_mf_dist[d_f - mode_f] += 1
    print(f"n={n}: {tight} tight cases, d_f - mode_f distribution: {dict(sorted(df_mf_dist.items()))}")
