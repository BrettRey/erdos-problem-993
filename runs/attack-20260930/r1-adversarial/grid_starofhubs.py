"""Grid scan of star-of-hubs S(h,m,t): centre joined to h hubs, each hub carrying m t-cherries.
Exact R1 verdict per tree (r1lib), float diagnostics. Usage: python grid_starofhubs.py T HMAX MMAX NMAX"""
import json, sys
from families_r1 import star_of_hubs
from r1lib import evaluate
t, hmax, mmax, nmax = map(int, sys.argv[1:5])
for m in range(3, mmax + 1):
    for h in range(2, hmax + 1):
        n = 1 + h * (1 + m * (t + 1))
        if n > nmax:
            continue
        r = evaluate(star_of_hubs(h, m, t))
        print(json.dumps({"h": h, "m": m, "t": t, "n": r["n"], "alpha": r["alpha"], "window": r["window"],
                          "n_viol": len(r["viol"]), "viol_k": sorted({v["k"] for v in r["viol"]}),
                          "viol_u_deg": sorted({v["deg_u"] for v in r["viol"]}), "lc_fail": r["lc_fail"],
                          "F": r["F"], "F_at": r["F_at"], "C": r["C"], "E": r["E"], "A": r["A"], "B": r["B"]}), flush=True)
