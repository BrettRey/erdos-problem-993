#!/usr/bin/env python3
"""KT1 + KT2 exact census over all trees with gentreeg -p, n in [nmin, nmax].

For every tree and every k in the central window [ceil(n/4), min(q, alpha-1)],
q = ceil((2 alpha - 1)/3):
  delta_k = 1 - i_{k-1} i_{k+1} / i_k^2                      (exact Fraction)
  Gamma_k = sigma^2(lam) * delta_k at lam in {i_{k-1}/i_k, i_k/i_{k+1}}  (exact)
  rho_k   = max over bipartition sides of E Var(X | other side)/sigma^2 (exact)
Usage: kt12_census.py nmin nmax res mod  (shard res of mod; gentreeg res/mod)
Output: one JSON line per n with exact extremes (as strings) and witnesses.
"""
import sys, json, subprocess
from fractions import Fraction
from lclt_lib import parent_line_to_adj, analyse

def run(n, res, mod):
    cmd = ["gentreeg", "-p", "-q", str(n), f"{res}/{mod}"]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    st = dict(n=n, shard=f"{res}/{mod}", trees=0, pairs=0, delta_nonpos=0,
              min_Gamma=None, max_Gamma=None, min_rho=None, max_lam_hi=None,
              w_min_Gamma=None, w_max_Gamma=None, w_min_rho=None, w_max_lam=None)
    for line in p.stdout:
        par = list(map(int, line.split()))
        if len(par) != n:
            continue
        adj = parent_line_to_adj(par)
        Z, alpha, (lo, hi, q), rows = analyse(adj)
        st["trees"] += 1
        for r in rows:
            st["pairs"] += 1
            if not r["delta_pos"]:
                st["delta_nonpos"] += 1
            tag = dict(par=line.strip(), k=r["k"], alpha=alpha, q=q)
            if st["min_Gamma"] is None or r["Gamma_lo"] < st["min_Gamma"]:
                st["min_Gamma"] = r["Gamma_lo"]; st["w_min_Gamma"] = tag
            if st["max_Gamma"] is None or r["Gamma_hi"] > st["max_Gamma"]:
                st["max_Gamma"] = r["Gamma_hi"]; st["w_max_Gamma"] = tag
            if st["min_rho"] is None or r["rho"] < st["min_rho"]:
                st["min_rho"] = r["rho"]; st["w_min_rho"] = tag
            if st["max_lam_hi"] is None or r["lam_hi"] > st["max_lam_hi"]:
                st["max_lam_hi"] = r["lam_hi"]; st["w_max_lam"] = tag
    p.wait()
    for key in ("min_Gamma", "max_Gamma", "min_rho", "max_lam_hi"):
        v = st[key]
        if v is not None:
            st[key + "_float"] = float(v)
            st[key] = f"{v.numerator}/{v.denominator}" if len(str(v)) < 400 else f"~{float(v)!r} (exact omitted, >400 chars)"
    return st

if __name__ == "__main__":
    nmin, nmax, res, mod = map(int, sys.argv[1:5])
    for n in range(nmin, nmax + 1):
        print(json.dumps(run(n, res, mod)), flush=True)
