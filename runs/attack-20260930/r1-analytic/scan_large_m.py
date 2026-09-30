"""Exact window scans of H(m,s) (s=1..4) and TH(m,2) (plus TH(m,s)) at large m.
Per family and vertex class u: worst (max over window k) neighbourhood-normalised load
L_u / sum_{v in N[u]} |D_v|/(deg v+1)  [float diagnostic], exact count of R1 failures,
exact count of window levels where the hub fails PV, and the hub PV ratio max."""
import sys, json, fam
ms = [int(a) for a in sys.argv[2].split(",")]
kind = sys.argv[1]
for s in ([1, 2, 3, 4] if kind == "H" else [1, 2, 3]):
    for m in ms:
        F = fam.hubstar(m, s) if kind == "H" else fam.twohub(m, s)
        rows = F.scan()
        rec = {"fam": F.label, "n": F.n, "alpha": F.alpha, "window": [F.lo, F.q],
               "r1_fail_levels": sum(1 for r in rows if r["r1_fail"]),
               "lc_fail_levels": sum(1 for r in rows if not r["lc_ok"]),
               "hub_pv_fail_levels": sum(1 for r in rows if "hub" in r["pv_fail"]),
               "other_pv_fail_levels": sum(1 for r in rows if set(r["pv_fail"]) - {"hub"}),
               "hub_pv_ratio_max": max((r["pv_ratio"]["hub"] or 0) for r in rows),
               "worst_nbhd_norm": {u: round(max(r["r1_nbhd_norm"][u] for r in rows), 6) for u in rows[0]["r1_nbhd_norm"]},
               "worst_mean_norm": {u: round(max(r["r1_mean_norm"][u] for r in rows), 6) for u in rows[0]["r1_mean_norm"]}}
        print(json.dumps(rec), flush=True)
