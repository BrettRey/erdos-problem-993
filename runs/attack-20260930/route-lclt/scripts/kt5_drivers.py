#!/usr/bin/env python3
"""Drivers for KT4/KT5 (float64 DIAGNOSTICS; see kt4_mixing_clt.py).
  kt5_drivers.py families NMAX         -> certificate + mixing-CLT stats per family member
  kt5_drivers.py census nmin nmax res mod -> certificate over all gentreeg trees
"""
import sys, json, subprocess
import kt12_families as F
import kt4_mixing_clt as W
from lclt_lib import parent_line_to_adj

def fam(nmax):
    for name, lst in F.FAMILIES.items():
        for gen, args in lst:
            adj = gen(*args)
            if len(adj) > nmax: continue
            c = W.certificate(adj)
            if c is None: continue
            m = W.analyse_W(adj)
            out = dict(family=name, n=c["n"], alpha=c["alpha"], min_cert_rel=c["min_cert_rel"], all_cert=c["all_cert"],
                       max_dK_Mprime=max(r["dK_Mprime"] for r in c["rows"]),
                       max_dK_M=m["max_dK"], min_rho=m["min_rho"],
                       max_sdB_over_EB=max(r["sdB_over_EB"] or 0 for r in m["rows"]),
                       rows=[{k: (float(v) if hasattr(v, "__float__") and not isinstance(v, bool) else v) for k, v in r.items()} for r in c["rows"]])
            print(json.dumps(out), flush=True)

def census(nmin, nmax, res, mod):
    for n in range(nmin, nmax + 1):
        p = subprocess.Popen(["gentreeg", "-p", "-q", str(n), f"{res}/{mod}"], stdout=subprocess.PIPE, text=True)
        st = dict(n=n, shard=f"{res}/{mod}", trees=0, pairs=0, cert_fail_pairs=0, trees_all_cert=0,
                  min_cert_rel=None, w_min=None, max_dK_Mprime=0.0)
        for line in p.stdout:
            par = list(map(int, line.split()))
            if len(par) != n: continue
            c = W.certificate(parent_line_to_adj(par))
            if c is None: continue
            st["trees"] += 1
            ok = True
            for r in c["rows"]:
                st["pairs"] += 1
                st["max_dK_Mprime"] = max(st["max_dK_Mprime"], r["dK_Mprime"])
                if r["cert"] <= 0:
                    st["cert_fail_pairs"] += 1; ok = False
            st["trees_all_cert"] += ok
            if st["min_cert_rel"] is None or c["min_cert_rel"] < st["min_cert_rel"]:
                st["min_cert_rel"] = c["min_cert_rel"]; st["w_min"] = line.strip()
        p.wait()
        print(json.dumps(st), flush=True)

if __name__ == "__main__":
    if sys.argv[1] == "families":
        fam(int(sys.argv[2]))
    else:
        census(*map(int, sys.argv[2:6]))
