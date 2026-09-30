"""Exact scan of SH(h,m,s) along h (fixed m,s): R1 failure levels at z and worst z-load."""
import sys, json, fam
s = int(sys.argv[1]); m = int(sys.argv[2]); hs = [int(a) for a in sys.argv[3].split(",")]
for h in hs:
    F = fam.starofhubs(h, m, s)
    rows = F.scan()
    fk = [r["k"] for r in rows if "z" in r["r1_fail"]]
    other = [(r["k"], r["r1_fail"]) for r in rows if set(r["r1_fail"]) - {"z"}]
    wk = max(rows, key=lambda r: r["r1_nbhd_norm"]["z"])
    print(json.dumps({"fam": F.label, "n": F.n, "window": [F.lo, F.q], "z_fail_k": [fk[0], fk[-1], len(fk)] if fk else [],
                      "other_fail": other, "lc_fail_levels": sum(1 for r in rows if not r["lc_ok"]),
                      "z_worst_nbhd_norm": round(wk["r1_nbhd_norm"]["z"], 5), "z_worst_k": wk["k"],
                      "z_worst_k_per_unit": round(wk["k"] / h, 4)}), flush=True)
