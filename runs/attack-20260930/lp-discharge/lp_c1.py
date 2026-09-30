"""Family 1: degree-based one-hop discharging. A degree-d vertex keeps sigma(d)
of its defect and gives (1 - sigma(d))/d to each neighbour. Certificate for tree
T at level k: for every u, L_u = sigma(d_u) D_u + sum_{v~u} (1-sigma(d_v))/d_v D_v <= 0
(sum_u L_u = sum_v D_v, so this gives LC). R1 is sigma(d) = 1/(d+1).
LP: maximize t s.t. L_u/(k i_{k-1} i_k) <= -t, 0 <= sigma <= 1. Cutting planes over
all trees n <= 16 and the families. t* < 0 means NO degree-based rule works on
the rows collected (the dual names the obstructing rows)."""
import sys, json, time
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
import dlib

def tree_rows(name, adj):
    deg = [len(x) for x in adj]
    return [(name, k, u, deg, adj, Dn, D) for (k, D, Dn) in dlib.rows_for(adj) for u in range(len(adj))]

def row_value(sig, r):
    name, k, u, deg, adj, Dn, D = r
    return sig[deg[u]] * Dn[u] + sum((1 - sig[deg[v]]) / deg[v] * Dn[v] for v in adj[u])

def solve(rows, dmax):
    nv = dmax + 2                          # sigma[0..dmax] (0 unused), t at index dmax+1
    A = lil_matrix((len(rows), nv)); b = np.zeros(len(rows))
    for i, (name, k, u, deg, adj, Dn, D) in enumerate(rows):
        A[i, deg[u]] += Dn[u]; c0 = 0.0
        for v in adj[u]:
            A[i, deg[v]] += -Dn[v] / deg[v]; c0 += Dn[v] / deg[v]
        A[i, nv - 1] = 1.0; b[i] = -c0
    cost = np.zeros(nv); cost[-1] = -1.0
    bounds = [(0, 1)] * (nv - 1) + [(None, 1)]
    res = linprog(cost, A_ub=A.tocsr(), b_ub=b, bounds=bounds, method='highs')
    return res

t0 = time.time()
train = []
for n in range(4, 13):
    for adj in dlib.gentreeg(n): train += tree_rows(f'n{n}', adj)
fam = dlib.families(); fam['R1killer'] = dlib.load_r1_killer()
famrows = {nm: tree_rows(nm, adj) for nm, adj in fam.items()}
for nm in ('R1killer', 'H(9,2)', 'SoS(3,3)'): train += famrows[nm]
test = {}
for n in range(13, 17):
    test[n] = [r for adj in dlib.gentreeg(n) for r in tree_rows(f'n{n}', adj)]
print(f'rows: train {len(train)}, test n13-16 {sum(len(v) for v in test.values())}, families {sum(len(v) for v in famrows.values())} ({time.time()-t0:.0f}s)', flush=True)
dmax = max(max(r[3]) for rs in [train] + list(test.values()) + list(famrows.values()) for r in rs)
log = []
for it in range(30):
    res = solve(train, dmax)
    sig = list(res.x[:dmax + 1]); t = res.x[-1]
    viol = []
    for rs in list(test.values()) + list(famrows.values()):
        for r in rs:
            if row_value(sig, r) > -max(t, 0) + 1e-12: viol.append(r)
    worst = max((row_value(sig, r) for r in viol), default=None)
    print(f'iter {it}: status={res.status} t*={t:.6g} violated rows={len(viol)} worst={worst}', flush=True)
    log.append(dict(iter=it, t=t, n_viol=len(viol), sigma=sig[1:dmax + 1]))
    if t < 0 or not viol: break
    viol.sort(key=lambda r: -row_value(sig, r)); train += viol[:2000]
print('sigma(d), d=1..%d:' % dmax, [round(s, 4) for s in sig[1:dmax + 1]])
print('R1 sigma for comparison:', [round(1 / (d + 1), 4) for d in range(1, dmax + 1)])
json.dump(dict(log=log, dmax=dmax, final_t=t), open('lp_c1_result.json', 'w'), indent=1)
if t < 0:
    # dual: which rows carry the obstruction
    res2 = solve(train, dmax)
    y = -res2.ineqlin.marginals
    idx = np.argsort(-y)[:15]
    print('top dual-weighted rows (tree, k, u, deg_u, weight):')
    for i in idx:
        if y[i] > 1e-9: print('  ', train[i][0], train[i][1], train[i][2], train[i][3][train[i][2]], round(float(y[i]), 5))
