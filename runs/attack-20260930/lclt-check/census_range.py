#!/usr/bin/env python3
"""Exact census of A0, A, A', B, C (exact_AC.check_tree_k) on all trees with n in [n1, n2]."""
import sys, json, subprocess
from exact_AC import par_to_adj, indep_brute, indep_dp, joint_brute, joint_dp, window, check_tree_k
from math import comb
n1, n2 = int(sys.argv[1]), int(sys.argv[2])
stats = dict(trees=0, A0_cases=0, A_cases=0, C_cases=0, cert_old_pos=0, cert_sharp_pos=0, mid_cases=0, mid_old_pos=0, mid_sharp_pos=0, dp_vs_brute=0)
per_n = {}
for n in range(n1, n2 + 1):
    t0 = stats["trees"]
    out = subprocess.run(["gentreeg", "-p", "-q", str(n)], capture_output=True, text=True).stdout
    for line in out.splitlines():
        par = list(map(int, line.split()))
        if len(par) != n: continue
        adj = par_to_adj(par); Z = indep_dp(adj)
        if n <= 14: assert Z == indep_brute(adj)
        alpha = len(Z) - 1; lo, hi = window(n, alpha); stats["trees"] += 1
        for side in (0, 1):
            N = joint_dp(adj, side)
            assert N == joint_brute(adj, side); stats["dp_vs_brute"] += 1
            for j in range(alpha + 2):
                assert sum(c * comb(B, j - s) for (s, B), c in N.items() if 0 <= j - s <= B) == (Z[j] if j <= alpha else 0)
            stats["A0_cases"] += 1
            Bmax = max(B for (s, B) in N)
            for k in range(lo, hi + 1):
                check_tree_k(N, Z, k, list(range(1, Bmax + 1)), stats)
    per_n[n] = stats["trees"] - t0
    print(n, per_n[n], flush=True)
print(json.dumps(dict(range=[n1, n2], per_n=per_n, **stats)))
