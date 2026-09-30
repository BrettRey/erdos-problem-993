"""Hill-climb / annealing in orbit-spec space (multiplicity vectors), exact via orbit.py.
MODE hub: root hub with a_s supports carrying s leaves (s = 1..9), e hub leaves, p pendant
          paths of length 2 at the hub, g sub-hubs H(4,2) at the hub.  Objective: maximise
          maxL23 (S(2/3) violation from the upper side) or minimise up_f.
MODE centre: centre with c_j hubs of type H(m_j, s_j) from a menu, plus e centre leaves.
          Objective: maximise low_f or maxL23.
Constraint n <= NMAX. Usage: python spec_search.py MODE OBJ NMAX SECONDS SEED OUT.json"""
import sys, json, math, random, time
from fractions import Fraction as F
from orbit import analyse_spec
sys.set_int_max_str_digits(0)
mode, obj, nmax, secs, seed, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
rng = random.Random(seed)
MENU = [(m, s) for s in (1, 2, 3) for m in range(6, 17)]
def spec_of(x):
    if mode == 'hub':
        sp = [(x[s - 1], [(s, [])]) for s in range(1, 10) if x[s - 1] > 0]
        if x[9]: sp.append((x[9], []))
        if x[10]: sp.append((x[10], [(1, [(1, [])])]))
        if x[11]: sp.append((x[11], [(4, [(2, [])])]))
        return sp
    sp = [(x[j], [(m, [(s, [])])]) for j, (m, s) in enumerate(MENU) if x[j] > 0]
    if x[-1]: sp.append((x[-1], []))
    return sp
def nsize(x):
    if mode == 'hub':
        return 1 + sum(x[s - 1] * (s + 1) for s in range(1, 10)) + x[9] + 2 * x[10] + x[11] * 13
    return 1 + sum(x[j] * (1 + m * (s + 1)) for j, (m, s) in enumerate(MENU)) + x[-1]
def score(r):
    if obj == 'L23': return r['maxL23']
    if obj == 'up': return -r['up_f']
    return r['low_f']
def start():
    if mode == 'hub':
        x = [0] * 12; s = rng.choice([4, 5, 6]); x[s - 1] = (nmax - 1) // (s + 1); return x
    x = [0] * (len(MENU) + 1); j = MENU.index((rng.choice([10, 11, 12]), 2)); x[j] = (nmax - 1) // (1 + MENU[j][0] * 3); return x
t0 = time.time(); best = None; evals = 0; hist = []
while time.time() - t0 < secs:
    x = start()
    while nsize(x) > nmax: x[[i for i, v in enumerate(x) if v][0]] -= 1
    r = analyse_spec(spec_of(x)); cur = score(r); T0 = 0.05 * abs(cur) + 1e-6
    for step in range(400):
        if time.time() - t0 >= secs: break
        y = list(x)
        for _ in range(rng.choice([1, 1, 2])):
            i = rng.randrange(len(y)); y[i] = max(0, y[i] + rng.choice([-3, -2, -1, 1, 2, 3]))
        if nsize(y) > nmax or nsize(y) < 30 or sum(y) < 3: continue
        try: r2 = analyse_spec(spec_of(y))
        except AssertionError: continue
        evals += 1; s2 = score(r2); T = T0 * 0.01 ** (step / 400)
        if s2 >= cur or rng.random() < math.exp((s2 - cur) / T):
            x, cur, r = y, s2, r2
            if best is None or cur > best['score']:
                best = dict(score=cur, x=x, spec=spec_of(x), n=r['n'], low_f=r['low_f'], up_f=r['up_f'], low=r.get('low'), up=r.get('up'),
                            low_at=r.get('low_at'), up_at=r.get('up_at'), maxL23=r['maxL23'], n_viol23=r['n_viol23'])
    hist.append(best and (best['score'], best['n']))
    json.dump(dict(mode=mode, obj=obj, nmax=nmax, seed=seed, evals=evals, elapsed=time.time() - t0, best=best, menu=MENU if mode == 'centre' else None), open(out, 'w'))
