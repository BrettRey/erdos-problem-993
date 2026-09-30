"""For each (m,s): smallest h with an R1 failure in SH(h,m,s) (exact), h <= HMAX."""
import sys, json, fam
HMAX = int(sys.argv[1])
res = []
for s in (1, 2, 3, 4):
    for m in range(2, 26):
        found = None
        for h in range(2, HMAX + 1):
            F = fam.starofhubs(h, m, s)
            if F.n > 700: break
            rows = F.scan()
            fails = [(r["k"], r["r1_fail"], r["lc_ok"]) for r in rows if r["r1_fail"]]
            if fails:
                found = {"s": s, "m": m, "h": h, "n": F.n, "window": [F.lo, F.q], "fails": fails}
                break
        print(json.dumps(found if found else {"s": s, "m": m, "none_up_to_h": h, "n_last": F.n}), flush=True)
