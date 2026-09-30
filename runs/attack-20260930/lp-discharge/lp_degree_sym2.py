"""Degree-based LP with every row normalized to unit l1 norm (coefficients +
constant), so solver tolerances are relative. Maximize the relative margin t."""
import numpy as np, json
from scipy.optimize import linprog
from scipy import sparse
Z = np.load('symrows.npz'); R, C, V, c, meta = Z['R'], Z['C'], Z['V'], Z['const'], Z['meta']
degs = sorted(set(C.tolist())); col = {d: i for i, d in enumerate(degs)}; nd = len(degs)
A = sparse.csr_matrix((V, (R, [col[d] for d in C])), shape=(len(c), nd)).tocsr()
s = np.asarray(abs(A).sum(axis=1)).ravel() + np.abs(c)
Dg = sparse.diags(1 / s); An = (Dg @ A).tocsr(); cn = c / s
Aub = sparse.hstack([An, sparse.csr_matrix(np.ones((len(c), 1)))]).tocsr()
cost = np.zeros(nd + 1); cost[-1] = -1
opts = dict(primal_feasibility_tolerance=1e-10, dual_feasibility_tolerance=1e-10)
res = linprog(cost, A_ub=Aub, b_ub=-cn, bounds=[(0, 1)] * nd + [(None, 1)], method='highs', options=opts)
t = res.x[-1]; sig = res.x[:nd]
viol = An @ sig + cn
print(f'status={res.status} t*={t:.3e}; with t dropped, max normalized row = {viol.max():.3e}, rows > 0: {(viol > 0).sum()}')
json.dump(dict(t=float(t), degrees=degs, sigma=sig.tolist()), open('lp_degree_sym2_result.json', 'w'))
y = -res.ineqlin.marginals; sup = np.where(y > 1e-10)[0]
print('dual support size:', len(sup))
grid = [l.split('\t') for l in open('symrows_grid.txt').read().split('\n')]
typ = {('H', 0): 'hub', ('H', 1): 'sup', ('H', 2): 'leaf', ('MSH', 0): 'cen', ('MSH', 1): 'hub', ('MSH', 2): 'sup', ('MSH', 3): 'leaf'}
for i in sup[np.argsort(-y[sup])][:12]:
    ti, k, tt = meta[i]; print('   y=%.4g' % y[i], grid[ti][1], grid[ti][2], 'k=%d' % k, typ[(grid[ti][1], int(tt))])
