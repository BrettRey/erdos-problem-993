#!/usr/bin/env python3
"""DIAGNOSTIC: architecture floor with idealised constants: a in {2/3, 0.82}, C_delta = C_gamma = 1,
and c_sigma = c_F = 1/4 (versus the paper's 2.1e-4 and 4.2e-4). Uses the chain functions of n0_chains.py."""
import json, mpmath as mp
import n0_chains as nc
out={}
for a in [mp.mpf(2)/3, mp.mpf("0.82")]:
    cst=dict(a=a, C_delta=mp.mpf(1), C_gamma=mp.mpf(1))
    for label,(cs,cf) in {"paper_c":(nc.C_SIGMA_PAPER,nc.C_F_PAPER),"ideal_c=1/4":(mp.mpf(1)/4,mp.mpf(1)/4)}.items():
        L=nc.lean_style_chain(cst,cs,cf); P=nc.paper_route_chain(cst,cs,cf)
        key=f"a={mp.nstr(a,4)},{label}"
        out[key]=dict(lean_style_log10_N0=mp.nstr(L["log10_N0"],6), R=mp.nstr(L["R"],6), eta=mp.nstr(L["eta"],6),
                      paper_route_log10_N0=mp.nstr(P["log10_N0"],6), paper_route_dominant=P["dominant"])
        print(key, out[key])
json.dump(out, open("data/floor_ideal.json","w"), indent=1)
