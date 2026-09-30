"""Do R1 failures occur strictly before the mode (ascending part, i_{k+1} > i_k)? Exact."""
import sys, json, fam
for (h, m, s) in [(8, 15, 3), (20, 16, 3), (40, 16, 3), (30, 11, 2), (40, 20, 2)]:
    F = fam.starofhubs(h, m, s)
    I = F.I
    mode = max(range(F.alpha + 1), key=lambda k: int(I[k]))
    rows = F.scan()
    fk = [r["k"] for r in rows if "z" in r["r1_fail"]]
    pre = [k for k in fk if fam.co(I, k + 1) > fam.co(I, k)]
    print(json.dumps({"fam": F.label, "n": F.n, "window": [F.lo, F.q], "mode": mode,
                      "z_fail_range": [fk[0], fk[-1]] if fk else None, "n_fail": len(fk),
                      "fail_levels_strictly_ascending": [pre[0], pre[-1], len(pre)] if pre else []}), flush=True)
