"""Simulated annealing against S(sigma) over trees NMIN <= n <= NMAX (mutation set reused
from ../../r1-adversarial/sa_r1.py). OBJ: 'L23' maximise max normalised L_u at sigma=2/3;
'up' minimise the interval's upper end; 'low' maximise its lower end. Floats steer only;
every violation / interval end reported is the exact integer verdict from slib.analyse.
Usage: python sa_const.py OBJ NMIN NMAX SECONDS SEED OUT.json"""
import sys, os, json, math, random, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'r1-adversarial'))
from sa_r1 import mutate, edges_of
from slib import analyse, adj_from_edges, B
def hstar(m, s, extra=0):
    b = B()
    for _ in range(m):
        x = b.add(0)
        for _ in range(s): b.add(x)
    for _ in range(extra): b.add(0)
    return b.adj()
def msh(h, m, s):
    b = B()
    for _ in range(h):
        hub = b.add(0)
        for _ in range(m):
            x = b.add(hub)
            for _ in range(s): b.add(x)
    return b.adj()
def val(r, obj):
    return r['maxL23'] if obj == 'L23' else (-r['up_f'] if obj == 'up' else r['low_f'])
obj, nmin, nmax, secs, seed, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
rng = random.Random(seed); t0 = time.time(); evals = 0; best = None; viols = []; restarts = 0
def seeds():
    S = []
    for s in (3, 4, 5, 6):
        for m in range(10, 80):
            if nmin <= 1 + m * (s + 1) <= nmax: S.append(hstar(m, s))
    for s in (1, 2, 3):
        for m in range(6, 17):
            for h in range(2, 30):
                if nmin <= 1 + h * (1 + m * (s + 1)) <= nmax: S.append(msh(h, m, s))
    return S
SEEDS = seeds()
# start from the best seeds for this objective
scored = sorted(((val(analyse(a), obj), i) for i, a in enumerate(SEEDS)), reverse=True)
top = [SEEDS[i] for _, i in scored[:8]]
while time.time() - t0 < secs:
    restarts += 1
    adj = top[restarts % len(top)]
    r = analyse(adj); cur = val(r, obj); T0 = 0.02 * (abs(cur) + 1e-4)
    for step in range(600):
        if time.time() - t0 >= secs: break
        T = T0 * (0.01) ** (step / 600)
        try: adj2 = mutate(adj, rng, nmin, nmax)
        except (IndexError, ValueError): continue
        if not (nmin <= len(adj2) <= nmax): continue
        r2 = analyse(adj2); evals += 1; v2 = val(r2, obj)
        if r2['n_viol23'] and len(viols) < 10: viols.append(dict(n=r2['n'], edges=edges_of(adj2), viol=r2['viol23'][:5], up=r2.get('up'), low=r2.get('low')))
        if v2 >= cur or rng.random() < math.exp((v2 - cur) / T):
            adj, cur, r = adj2, v2, r2
            if best is None or cur > best['val']:
                best = dict(val=cur, n=r['n'], low_f=r['low_f'], up_f=r['up_f'], low=r.get('low'), up=r.get('up'), low_at=r.get('low_at'), up_at=r.get('up_at'),
                            maxL23=r['maxL23'], edges=edges_of(adj), degseq=sorted((len(x) for x in adj), reverse=True)[:6])
    json.dump(dict(obj=obj, nmin=nmin, nmax=nmax, seed=seed, elapsed=time.time() - t0, evals=evals, restarts=restarts,
                   seed_best=scored[0][0], best=best, violations=viols), open(out, 'w'))
