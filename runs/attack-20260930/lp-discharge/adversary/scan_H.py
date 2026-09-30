"""Upper-end trend for single hub-stars H(m,s) via the orbit evaluator."""
import sys, json, time
from orbit import analyse_spec, H
out = open(sys.argv[1], 'w')
for s in [int(x) for x in sys.argv[2].split(',') if x]:
    for m in [int(x) for x in sys.argv[3].split(',') if x]:
        t0 = time.time(); r = analyse_spec(H(m, s)); r.update(family='H', params=dict(m=m, s=s), secs=round(time.time() - t0, 2))
        out.write(json.dumps(r) + '\n'); out.flush()
        print(f"H({m},{s}) n={r['n']} up={r['up_f']:.6f} at{r['up_at']} win={r['window']} maxL23={r['maxL23']:.3e} viol={r['n_viol23']} {r['secs']}s", flush=True)
