#!/usr/bin/env python3
"""Exact (Fraction) checks of A, A0, A', B, C on three more trees, n = 60-63, DP only
(brute force over J is infeasible here; the DP was validated against brute force on all
trees n <= 12 and on n <= 40 trees with |R| <= 20). Reuses exact_AC.check_tree_k."""
import json, math, random
from fractions import Fraction as Fr
from math import comb
from exact_AC import joint_dp, indep_dp, window, check_tree_k
from exact_large import adj_from_parent

def complete_binary(depth):
    par = [0]
    for v in range(2, 2 ** (depth + 1)):
        par.append(v // 2)
    return par

def caterpillar(spine, legs):
    par = [0] + list(range(1, spine))
    for s in range(1, spine + 1):
        par += [s] * legs
    return par

def random_tree(n, seed):
    rng = random.Random(seed)
    return [0] + [rng.randint(1, v) for v in range(1, n)]

TREES = {"complete_binary_depth5_n63": complete_binary(5), "caterpillar_20x2_n60": caterpillar(20, 2),
         "random_recursive_n60_seed993": random_tree(60, 993)}
res = {}
for name, par in TREES.items():
    adj = adj_from_parent(par); n = len(adj)
    Z = indep_dp(adj); alpha = len(Z) - 1
    lo, hi = window(n, alpha)
    stats = dict(A_cases=0, C_cases=0, cert_old_pos=0, cert_sharp_pos=0, mid_cases=0, mid_old_pos=0, mid_sharp_pos=0)
    rows = []
    for side in (0, 1):
        N = joint_dp(adj, side)
        for j in range(alpha + 2):
            assert sum(c * comb(B, j - s) for (s, B), c in N.items() if 0 <= j - s <= B) == (Z[j] if j <= alpha else 0)
        for k in sorted(set([lo, (lo + hi) // 2, hi])):
            lam = Fr(Z[k - 1], Z[k]); Zl = sum(c * lam ** j for j, c in enumerate(Z))
            EB = sum(c * lam ** s * (1 + lam) ** B * B for (s, B), c in N.items()) / Zl
            B0s = sorted(set(max(1, math.floor((1 - Fr(e, 10)) * EB)) for e in (3, 5)))
            sub = []; check_tree_k(N, Z, k, B0s, stats, record=sub)
            for r in sub: r["side"] = side
            rows += sub
    res[name] = dict(parent_array=par, n=n, alpha=alpha, window=[lo, hi], stats=stats, rows=rows)
    print(name, n, alpha, (lo, hi), stats, flush=True)
json.dump(res, open("rerun/exact_large2.json", "w"), indent=1)
