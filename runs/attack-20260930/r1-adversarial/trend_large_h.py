"""Trend of R1 on S(h,m,2) as h grows (fixed m): exact verdict via r1lib, float F/C at the worst (u,k).
Usage: python trend_large_h.py > data/trend_large_h.jsonl"""
import json, time
from families_r1 import star_of_hubs
from r1lib import evaluate
for m in (10, 9, 11, 12, 8):
    for h in (12, 14, 16, 20, 24):
        t0 = time.time()
        r = evaluate(star_of_hubs(h, m, 2))
        print(json.dumps({"h": h, "m": m, "n": r["n"], "window": r["window"], "n_viol": len(r["viol"]),
                          "viol_k": sorted({v["k"] for v in r["viol"]}), "viol_u_deg": sorted({v["deg_u"] for v in r["viol"]}),
                          "lc_fail": r["lc_fail"], "F": r["F"], "F_at": r["F_at"], "C": r["C"], "secs": round(time.time() - t0, 1)}), flush=True)
