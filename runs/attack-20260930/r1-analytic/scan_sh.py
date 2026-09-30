import sys, json, fam
s = int(sys.argv[1]); hmax = int(sys.argv[2]); ms = [int(a) for a in sys.argv[3].split(",")]
for m in ms:
    for h in range(2, hmax + 1):
        F = fam.starofhubs(h, m, s)
        rows = F.scan()
        fails = [(r["k"], r["r1_fail"]) for r in rows if r["r1_fail"]]
        pvf = [r["k"] for r in rows if "hub" in r["pv_fail"]]
        worst = max(r["r1_nbhd_norm"]["z"] for r in rows)
        print(json.dumps({"fam": F.label, "n": F.n, "win": [F.lo, F.q], "r1_fail": fails, "hub_pv_fail_k": pvf[:3] + (["..."] if len(pvf) > 3 else []), "z_nbhd_norm_max": round(worst, 4)}), flush=True)
        if fails:
            break
