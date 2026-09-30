#!/usr/bin/env python3
"""Summaries of KT1/KT2 families and KT4/KT5 families + census (reads data/*.jsonl)."""
import json, os
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
def jl(f):
    p = os.path.join(D, f)
    return [json.loads(l) for l in open(p)] if os.path.exists(p) else []
print("== KT1/KT2 families (Gamma, rho are 60-digit mpmath diagnostics; LC sign exact) ==")
print("family n alpha minGamma n*(1-minGamma) minRho maxLam lc_all")
fam = jl("kt12_families.jsonl")
for d in fam:
    print(d["family"], d["n"], d["alpha"], "%.4f" % d["min_Gamma"], "%.2f" % (d["n"] * d["max_dev"]), "%.3f" % d["min_rho"], "%.3f" % d["max_lam"], d["lc_all"])
print("overall: min rho %.4f ; max n*(1-Gamma) at n>=100: %.2f ; max lam (n>=100) %.3f ; all LC %s" % (
    min(d["min_rho"] for d in fam), max(d["n"] * d["max_dev"] for d in fam if d["n"] >= 100),
    max(d["max_lam"] for d in fam if d["n"] >= 100), all(d["lc_all"] for d in fam)))
print()
print("== KT4/KT5 families (float64 diagnostics) ==")
print("family n cert_all minCertRel maxdK(M') sqrt(n)*maxdK(M') maxdK(M) minRho maxsd(B)/E(B)")
k5 = jl("kt5_families.jsonl")
for d in k5:
    print(d["family"], d["n"], d["all_cert"], "%.3f" % d["min_cert_rel"], "%.4f" % d["max_dK_Mprime"],
          "%.3f" % (d["n"] ** 0.5 * d["max_dK_Mprime"]), "%.3f" % d["max_dK_M"], "%.3f" % d["min_rho"], "%.3f" % d["max_sdB_over_EB"])
if k5:
    print("overall: members %d, cert fails %d (n list %s); max sqrt(n) dK(M') = %.3f" % (
        len(k5), sum(not d["all_cert"] for d in k5), [ (d["family"], d["n"]) for d in k5 if not d["all_cert"]],
        max(d["n"] ** 0.5 * d["max_dK_Mprime"] for d in k5)))
print()
print("== KT5 census (all trees; certificate per (tree,k) at 5 window points) ==")
for d in jl("kt5_census_n8_16.jsonl"):
    print("n=%d trees=%d pairs=%d cert_fail_pairs=%d trees_all_cert=%d min_cert_rel=%.3f max_dK(M')=%.4f worst=%s" % (
        d["n"], d["trees"], d["pairs"], d["cert_fail_pairs"], d["trees_all_cert"], d["min_cert_rel"], d["max_dK_Mprime"], d["w_min"]))
print()
print("== KT1/KT2 census totals (exact) ==")
import glob as _g
from collections import defaultdict as _dd
_rows = _dd(list)
for _f in ["kt12_census_n4_14.jsonl"] + sorted(os.path.basename(x) for x in _g.glob(os.path.join(D, "kt12_census_n15_20_s*.jsonl"))):
    for _d in jl(_f): _rows[_d["n"]].append(_d)
_t = sum(d["trees"] for v in _rows.values() for d in v)
_p = sum(d["pairs"] for v in _rows.values() for d in v)
_b = sum(d["delta_nonpos"] for v in _rows.values() for d in v)
print("n range %d..%d: trees %d, window (tree,k) pairs %d, delta<=0 pairs %d" % (min(_rows), max(_rows), _t, _p, _b))
print("min rho over census %.4f ; min Gamma %.4f ; max Gamma %.4f ; max lam %.3f" % (
    min(d["min_rho_float"] for v in _rows.values() for d in v), min(d["min_Gamma_float"] for v in _rows.values() for d in v),
    max(d["max_Gamma_float"] for v in _rows.values() for d in v), max(d["max_lam_hi_float"] for v in _rows.values() for d in v)))
