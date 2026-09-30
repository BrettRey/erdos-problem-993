#!/usr/bin/env python3
"""Exact checks of A, A', B, C on a few larger trees (integer 2-D DP, Fractions).
Trees recorded as parent arrays (1-indexed parents, 0 = root) in the output."""
import sys, json, math
from fractions import Fraction as Fr
from math import comb
from exact_AC import joint_dp, joint_brute, indep_dp, window, check_tree_k

def adj_from_parent(par):
    n = len(par); adj = [[] for _ in range(n)]
    for v, p in enumerate(par):
        if p:
            adj[v].append(p - 1); adj[p - 1].append(v)
    return adj

def hubstar(m, t):            # hub joined to m vertices each carrying t leaves
    par = [0]
    for i in range(m):
        par.append(1); mid = len(par)
        for _ in range(t): par.append(mid)
    return par

def path(n): return [0] + list(range(1, n))

def double_star(a, b):        # two adjacent centres with a and b leaves
    par = [0, 1] + [1] * a + [2] * b
    return par

def spider(legs, length):
    par = [0]
    for _ in range(legs):
        prev = 1
        for _ in range(length):
            par.append(prev); prev = len(par)
    return par

TREES = {"hubstar_H(9,2)_n28": hubstar(9, 2), "path_n40": path(40), "double_star_14_14_n30": double_star(14, 14),
         "spider_5x6_n31": spider(5, 6)}

if __name__ == "__main__":
    res = {}
    for name, par in TREES.items():
        adj = adj_from_parent(par); n = len(adj)
        Z = indep_dp(adj); alpha = len(Z) - 1
        lo, hi = window(n, alpha)
        stats = dict(A_cases=0, C_cases=0, cert_old_pos=0, cert_sharp_pos=0, mid_cases=0, mid_old_pos=0, mid_sharp_pos=0)
        recs = []
        brute_checked = []
        for side in (0, 1):
            N = joint_dp(adj, side)
            nR = sum(1 for key in N for _ in [0]) and None
            # independent cross-check: brute force over J subset R when |R| <= 20
            from exact_AC import bfs
            _, depth, _ = bfs(adj)
            sizeR = sum(1 for v in range(n) if depth[v] % 2 != side)
            if sizeR <= 20:
                assert N == joint_brute(adj, side); brute_checked.append(side)
            # lambda-free identity vs standard DP
            for j in range(alpha + 2):
                assert sum(c * comb(B, j - s) for (s, B), c in N.items() if 0 <= j - s <= B) == (Z[j] if j <= alpha else 0)
            ks = sorted(set([lo, (lo + hi) // 2, hi]))
            for k in ks:
                lam = Fr(Z[k - 1], Z[k])
                Zl = sum(c * lam ** j for j, c in enumerate(Z))
                EB = sum(c * lam ** s * (1 + lam) ** B * B for (s, B), c in N.items()) / Zl
                B0s = sorted(set(max(1, math.floor((1 - Fr(eta, 10)) * EB)) for eta in (3, 5)))
                sub = []
                check_tree_k(N, Z, k, B0s, stats, record=sub)
                for r in sub: r["side"] = side
                recs += sub
        res[name] = dict(parent_array=par, n=n, alpha=alpha, window=[lo, hi], brute_checked_sides=brute_checked,
                         stats=stats, rows=recs)
        print(name, n, alpha, (lo, hi), stats, flush=True)
    json.dump(res, open("exact_large.json", "w"), indent=1)
