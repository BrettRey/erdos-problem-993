"""Trend scans via the orbit evaluator. mode H: H(m,s) upper end; mode MSH: MSH(h;m,s) both ends.
Usage: python scan_grid.py MODE OUT.jsonl NMAX 'list1' 'list2' ['list3']"""
import sys, json, time
from orbit import analyse_spec, H, MSH
mode, out, nmax = sys.argv[1], open(sys.argv[2], 'w'), int(sys.argv[3])
L = [[int(x) for x in a.split(',') if x] for a in sys.argv[4:]]
def emit(spec, fam, par):
    t0 = time.time(); r = analyse_spec(spec); r.update(family=fam, params=par, secs=round(time.time() - t0, 2))
    out.write(json.dumps(r) + '\n'); out.flush()
    print(f"{fam}{tuple(par.values())} n={r['n']} low={r['low_f']:.6f} at{r.get('low_at')} up={r['up_f']:.6f} at{r.get('up_at')} win={r['window']} maxL23={r['maxL23']:.3e} viol={r['n_viol23']} {r['secs']}s", flush=True)
if mode == 'H':
    for s in L[0]:
        for m in L[1]:
            if 1 + m * (s + 1) <= nmax: emit(H(m, s), 'H', dict(m=m, s=s))
else:
    for s in L[0]:
        for m in L[1]:
            for h in L[2]:
                if 1 + h * (1 + m * (s + 1)) <= nmax: emit(MSH(h, m, s), 'MSH', dict(h=h, m=m, s=s))
