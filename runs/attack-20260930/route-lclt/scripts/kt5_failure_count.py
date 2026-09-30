#!/usr/bin/env python3
"""For all trees at n (gentreeg), recount certificate failures (Cert<=0) and how many
of them have a positive Kolmogorov part EhG - TV*dK > 0 (i.e. fail only through the
crude tail term 2P(E0^c)).  float64 DIAGNOSTIC."""
import sys, subprocess, json
import kt4_mixing_clt as W
from lclt_lib import parent_line_to_adj
n = int(sys.argv[1])
out = subprocess.run(["gentreeg", "-p", "-q", str(n)], capture_output=True, text=True).stdout
fails = 0; tail_only = 0; pairs = 0
for line in out.splitlines():
    par = list(map(int, line.split()))
    if len(par) != n: continue
    c = W.certificate(parent_line_to_adj(par))
    if c is None: continue
    for r in c["rows"]:
        pairs += 1
        if r["cert"] <= 0:
            fails += 1
            if r["EhG"] - r["TV_h"] * r["dK_Mprime"] > 0: tail_only += 1
print(json.dumps(dict(n=n, pairs=pairs, cert_fail_pairs=fails, fail_with_positive_kolmogorov_part=tail_only)))
