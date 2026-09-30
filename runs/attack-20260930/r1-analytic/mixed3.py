"""Capped search: centre joined to hubs of up to three types (m in 5..14, s in {2,3}), n < NMAX.
Top 16 window levels scanned exactly; hits printed (then fully rescanned + generic recheck separately)."""
import sys, json, itertools, time
import mixed
NMAX = int(sys.argv[1]); TLIM = float(sys.argv[2])
t0 = time.time()
kinds = [(m, 2) for m in range(5, 15)] + [(m, 3) for m in range(5, 13)]
size = lambda m, s: 1 + m * (s + 1)
cands = set()
def rec(start, h_left, acc, n):
    if acc and n < NMAX and sum(c for c, _, _ in acc) >= 3:
        cands.add(tuple(sorted(acc)))
    if h_left == 0 or len(acc) == 3: return
    for i in range(start, len(kinds)):
        m, s = kinds[i]
        for c in range(1, h_left + 1):
            nn = n + c * size(m, s)
            if nn >= NMAX: break
            rec(i + 1, h_left - c, acc + [(c, m, s)], nn)
rec(0, 10, [], 1)
cands = sorted(cands, key=lambda T: 1 + sum(c * size(m, s) for c, m, s in T))
print(json.dumps({"n_candidates": len(cands)}), flush=True)
best = None; done = 0
for T in cands:
    if time.time() - t0 > TLIM: break
    F = mixed.mixed(list(T)); done += 1
    rows = F.scan(range(max(F.lo, F.q - 15), F.q + 1))
    fails = [(r["k"], r["r1_fail"]) for r in rows if r["r1_fail"]]
    if fails:
        print(json.dumps({"hit": list(T), "n": F.n, "fails": fails}), flush=True)
        if best is None or F.n < best: best = F.n
print(json.dumps({"evaluated": done, "of": len(cands), "smallest_hit_n": best, "seconds": round(time.time() - t0)}), flush=True)
