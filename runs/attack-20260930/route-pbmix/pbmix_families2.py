#!/usr/bin/env python3
"""Second batch: (a) random trees with B = branch set (path components, enumeration)
versus B = non-leaves (binomial DP); (b) combs / sparse combs at scale.
Exact arithmetic for rho/cert/LC; *_float columns are diagnostics."""
import json
import random
import sys
from fractions import Fraction as Fr

import pbmix_killtest as P


def comb_tree(L, every=1):
    """Spine of L vertices; a pendant leaf on spine vertex i when i % every == 0."""
    adj = [[] for _ in range(L)]
    for i in range(L - 1):
        adj[i].append(i + 1)
        adj[i + 1].append(i)
    for i in range(L):
        if i % every == 0:
            v = len(adj)
            adj.append([i])
            adj[i].append(v)
    return adj, len(adj)


def branch_run(name, adj, n, cert=True):
    B = {v for v in range(n) if len(adj[v]) >= 3}
    comps = P.components_branch_set(adj, n, B)
    f = [0] * (n + 1)
    for mult, s, Q in comps:
        for i, c in enumerate(Q):
            f[s + i] += mult * c
    while len(f) > 1 and f[-1] == 0:
        f.pop()
    assert f == P.indpoly_tree(adj, n)
    alpha = P.alpha_of(f)
    lo, q, _, _ = P.windows(n, alpha)
    rows = []
    for k in range(lo, q + 1):
        st = P.mixture_stats(comps, k, want_cert=cert)
        lam = Fr(f[k], f[k + 1])
        bf, Ew, Bv = P.between_fraction(comps, lam)
        rows.append({"k": k, "LC_exact": st["delta_f"] > 0, "rho_float": P.fl(st["rho"]),
                     "cert_pb_pass": (st["cert_pb"] > 0) if cert else None,
                     "cert_best_pass": (st["cert_best"] > 0) if cert else None,
                     "between_frac_float": P.fl(bf)})
    rh = [x["rho_float"] for x in rows]
    bfr = [x["between_frac_float"] for x in rows]
    print(name, "B=branch |B|", len(B), "configs", sum(c[0] for c in comps), "window", [lo, q],
          "allLC", all(x["LC_exact"] for x in rows),
          "rho[min,max]", (round(min(rh), 4), round(max(rh), 4)),
          "bf[min,max]", (round(min(bfr), 4), round(max(bfr), 4)),
          "cert_pb_fail", sum(1 for x in rows if x["cert_pb_pass"] is False),
          "cert_best_fail", sum(1 for x in rows if x["cert_best_pass"] is False),
          "pairs", len(rows), flush=True)
    return {"family": name, "B": "branch", "n": n, "window": [lo, q], "rows": rows}


def nonleaf_run(name, adj, n, cert=False):
    Bint = {v for v in range(n) if len(adj[v]) >= 2}
    r = P.analyse_binomial_family(name, adj, n, Bint, cert=cert)
    rh = [x["rho_float"] for x in r["rows"]]
    bfr = [x["between_frac_float"] for x in r["rows"]]
    print(name, "B=nonleaf |B|", len(Bint), "window", r["window"],
          "allLC", all(x["LC_exact"] for x in r["rows"]),
          "rho[min,max]", (round(min(rh), 4), round(max(rh), 4)),
          "bf[min,max]", (round(min(bfr), 4), round(max(bfr), 4)),
          "cert_pb_fail", sum(1 for x in r["rows"] if x["cert_pb_pass"] is False),
          "pairs", len(r["rows"]), flush=True)
    return r


def main():
    out = []
    rng = random.Random(993)
    for n in (30, 40, 50):
        for rep in range(3):
            adj = P.random_tree(n, rng)
            nb = sum(1 for v in range(n) if len(adj[v]) >= 3)
            if nb > 14:
                continue
            out.append(branch_run(f"rand n={n} rep={rep}", adj, n))
            out.append(nonleaf_run(f"rand n={n} rep={rep}", adj, n, cert=True))
    for L in (20, 40, 80):
        adj, n = comb_tree(L, 1)
        out.append(nonleaf_run(f"comb L={L} every=1", adj, n))
        adj, n = comb_tree(L, 2)
        out.append(nonleaf_run(f"comb L={L} every=2", adj, n))
        adj, n = comb_tree(L, 3)
        out.append(nonleaf_run(f"comb L={L} every=3", adj, n))
    json.dump(out, open("results_families2.json", "w"), indent=1)


if __name__ == "__main__":
    main()
