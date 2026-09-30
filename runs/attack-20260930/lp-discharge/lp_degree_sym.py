"""Degree-based discharging LP on the H and MSH grids: exists sigma(d) in [0,1]
with every row <= -t? Maximize t. If t* < 0, report the dual support (the rows
that jointly rule out every degree-based rule)."""
import numpy as np, json
from scipy.optimize import linprog
from scipy import sparse
Z = np.load('symrows.npz'); R, C, V, c, meta = Z['R'], Z['C'], Z['V'], Z['const'], Z['meta']
grid = [l.split('\t') for l in open('symrows_grid.txt').read().split('\n')]
degs = sorted(set(C.tolist())); col = {d: i for i, d in enumerate(degs)}; nd = len(degs)
A = sparse.csr_matrix((V, (R, [col[d] for d in C])), shape=(len(c), nd))
Aub = sparse.hstack([A, sparse.csr_matrix(np.ones((len(c), 1)))]).tocsr()
cost = np.zeros(nd + 1); cost[-1] = -1
res = linprog(cost, A_ub=Aub, b_ub=-c, bounds=[(0, 1)] * nd + [(None, 1)], method='highs')
t = res.x[-1]; sig = dict(zip(degs, res.x[:nd]))
print(f'rows={len(c)} degrees={nd} status={res.status} t*={t:.6g}')
out = dict(t=t, degrees=degs, sigma=[sig[d] for d in degs])
if t < 0:
    y = -res.ineqlin.marginals; sup = np.where(y > 1e-12)[0]
    print(f'dual support: {len(sup)} rows')
    typ = {('H', 0): 'hub', ('H', 1): 'sup', ('H', 2): 'leaf', ('MSH', 0): 'cen', ('MSH', 1): 'hub', ('MSH', 2): 'sup', ('MSH', 3): 'leaf'}
    rows = []
    for i in sup[np.argsort(-y[sup])]:
        ti, k, tt = meta[i]; fam, args = grid[ti][1], grid[ti][2]
        rows.append(dict(row=int(i), y=float(y[i]), family=fam, args=args, k=int(k), type=typ[(fam, int(tt))]))
    for r in rows[:25]: print('  ', r)
    out['dual_support'] = rows
else:
    print('sigma(d) for small degrees:', {d: round(sig[d], 4) for d in degs[:14]})
json.dump(out, open('lp_degree_sym_result.json', 'w'), indent=1)
