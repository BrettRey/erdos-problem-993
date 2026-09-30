"""Degree-based LP on the union of the base grid and the extended MSH grid
(h up to 120), rows unit-normalized. Reports t*, feasibility of t >= 0, and the
dual support (rows jointly excluding every degree-based rule when t* < 0)."""
import numpy as np, json, sys
from scipy.optimize import linprog
from scipy import sparse
parts = []
for f, g in (('symrows.npz', 'symrows_grid.txt'), ('symrows_ext.npz', 'symrows_ext_grid.txt')):
    Z = np.load(f); grid = [l.split('\t') for l in open(g).read().split('\n')]
    parts.append((Z, grid, f))
degs = sorted(set().union(*[set(Z['C'].tolist()) for Z, _, _ in parts])); col = {d: i for i, d in enumerate(degs)}; nd = len(degs)
blocks, consts, info = [], [], []
for Z, grid, f in parts:
    A = sparse.csr_matrix((Z['V'], (Z['R'], [col[d] for d in Z['C']])), shape=(len(Z['const']), nd)).tocsr()
    blocks.append(A); consts.append(Z['const']); info += [(f, grid, m) for m in Z['meta']]
A = sparse.vstack(blocks).tocsr(); c = np.concatenate(consts)
s = np.asarray(abs(A).sum(axis=1)).ravel() + np.abs(c); An = (sparse.diags(1 / s) @ A).tocsr(); cn = c / s
Aub = sparse.hstack([An, sparse.csr_matrix(np.ones((len(c), 1)))]).tocsr()
cost = np.zeros(nd + 1); cost[-1] = -1
res = linprog(cost, A_ub=Aub, b_ub=-cn, bounds=[(0, 1)] * nd + [(None, 1)], method='highs',
              options=dict(primal_feasibility_tolerance=1e-10, dual_feasibility_tolerance=1e-10))
t = res.x[-1]; print(f'rows={len(c)} degrees={nd} status={res.status} t*={t:.4e}')
y = -res.ineqlin.marginals; sup = np.where(y > 1e-10)[0]
typ = {('H', 0): 'hub', ('H', 1): 'sup', ('H', 2): 'leaf', ('MSH', 0): 'cen', ('MSH', 1): 'hub', ('MSH', 2): 'sup', ('MSH', 3): 'leaf'}
rows = []
for i in sup[np.argsort(-y[sup])]:
    f, grid, (ti, k, tt) = info[i]; fam, args = grid[ti][1], grid[ti][2]
    rows.append(dict(y=float(y[i]), family=fam, args=args, k=int(k), type=typ[(fam, int(tt))], file=f, ti=int(ti)))
print('dual support:', len(rows))
for r in rows[:12]: print('  ', r)
json.dump(dict(t=float(t), degrees=degs, sigma=res.x[:nd].tolist(), dual=rows), open('lp_degree_combined_result.json', 'w'), indent=1)
