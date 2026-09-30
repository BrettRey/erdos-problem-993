#!/usr/bin/env python3
"""Independent cross-check of central_margin.c for 10 <= n <= 14.

Enumerates trees with networkx.nonisomorphic_trees (not nauty), computes the
independence polynomial with the Python DP in treepoly.py, and recomputes in
exact Fractions: tree count, min m5, min m4, violation counts, per-k minima,
window-top minimum, and the minimum Newton ratio. Compares with
data/census_n{n}.json produced from gentreeg + the C program.
"""
import json
import os
import sys
from fractions import Fraction

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from treepoly import ipoly_from_adj, margins, q_of  # noqa: E402

res = {}
ok_all = True
for n in range(10, 15):
    cnt = 0
    m5 = m4 = None
    v5 = v4 = 0
    perk = {}
    topq = None
    newton = None
    for T in nx.nonisomorphic_trees(n):
        cnt += 1
        adj = [list(T.neighbors(v)) for v in range(n)]
        p = ipoly_from_adj(adj)
        alpha = len(p) - 1
        q = q_of(alpha)
        ms5 = margins(p, n, "n5")
        ms4 = margins(p, n, "n4")
        a5 = min(ms5.values())
        a4 = min(ms4.values())
        m5 = a5 if m5 is None else min(m5, a5)
        m4 = a4 if m4 is None else min(m4, a4)
        v5 += any(d < 0 for d in ms5.values())
        v4 += any(d < 0 for d in ms4.values())
        for k, d in ms5.items():
            perk[k] = d if k not in perk else min(perk[k], d)
            r = d * Fraction((k + 1) * (alpha - k + 1), alpha + 1)
            newton = r if newton is None else min(newton, r)
        topq = ms5[q] if topq is None else min(topq, ms5[q])
    c = json.load(open(os.path.join(HERE, "data", f"census_n{n}.json")))
    checks = {
        "count": cnt == c["trees"],
        "m5": m5 == Fraction(c["TOP5"][0]["margin"]),
        "m4": m4 == Fraction(c["TOP4"][0]["margin"]),
        "viol5": v5 == c["viol5"],
        "viol4": v4 == c["viol4"],
        "perk": all(perk[int(k)] == Fraction(r["value"]) for k, r in c["per_k_min_delta"].items())
                and set(map(int, c["per_k_min_delta"])) == set(perk),
        "topq": topq == Fraction(c["window_top_min_delta_q"]["value"]),
        "newton": newton == Fraction(c["newton_ratio_min_exact"]["value"]),
    }
    ok = all(checks.values())
    ok_all &= ok
    res[n] = dict(trees=cnt, m5=str(m5), m4=str(m4), viol5=v5, viol4=v4, checks=checks, ok=ok)
    print(n, cnt, m5, m4, v5, v4, checks, flush=True)

json.dump(dict(all_ok=ok_all, per_n=res), open(os.path.join(HERE, "data", "crosscheck_python_n10_14.json"), "w"), indent=1)
print("ALL_OK" if ok_all else "MISMATCH")
