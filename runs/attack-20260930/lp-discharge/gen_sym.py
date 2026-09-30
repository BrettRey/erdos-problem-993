"""Generate degree-based discharging rows for the H and MSH grids (closed forms,
python-flint). Row for vertex type u at level k (normalized by k i_{k-1} i_k):
  sigma(d_u) Dn_u + sum_{nbr types} cnt (1 - sigma(d_v))/d_v Dn_v <= -t.
Stored sparse: (row, degree, coef), const, and metadata. Asserts the exact identity
sum_v D_v(k) = k(k+1)(i_{k-1}i_{k+1} - i_k^2) with type multiplicities."""
import sys, time, numpy as np
from math import ceil
from symrows import H, MSH, coeffs
t0 = time.time()
import os
EXT = os.environ.get('EXT') == '1'
grid = []
for s in (range(1, 9) if not EXT else []):
    for m in list(range(3, 41)) + list(range(45, 151, 5)) + [175, 200]:
        grid.append(('H', (m, s)))
for s in (1, 2, 3):
    for m in (range(3, 15) if not EXT else range(9, 15)):
        for h in ((2, 3, 4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 38, 44, 50) if not EXT else (56, 64, 72, 80, 96, 120)):
            grid.append(('MSH', (h, m, s)))
R, C, V, const, meta = [], [], [], [], []
mult = lambda fam, a: ({'hub': 1, 'sup': a[0], 'leaf': a[0] * a[1]} if fam == 'H' else
                       {'cen': 1, 'hub': a[0], 'sup': a[0] * a[1], 'leaf': a[0] * a[1] * a[2]})
TYPES = {'H': ['hub', 'sup', 'leaf'], 'MSH': ['cen', 'hub', 'sup', 'leaf']}
row = 0
for ti, (fam, args) in enumerate(grid):
    I, J, n = (H if fam == 'H' else MSH)(*args)
    I = coeffs(I); a = len(I) - 1; q = ceil((2 * a - 1) / 3); lo = ceil(n / 4)
    Jc = {t: (coeffs(Jp), d, nb) for t, (Jp, d, nb) in J.items()}; mu = mult(fam, args)
    co = lambda p, k: p[k] if 0 <= k < len(p) else 0
    for k in range(lo, min(q, a - 1) + 1):
        D = {t: k * co(I, k - 1) * co(Jc[t][0], k) - (k + 1) * co(I, k) * co(Jc[t][0], k - 1) for t in Jc}
        assert sum(mu[t] * D[t] for t in D) == k * (k + 1) * (co(I, k - 1) * co(I, k + 1) - co(I, k) ** 2)
        nrm = k * co(I, k - 1) * co(I, k); Dn = {t: D[t] / nrm for t in D}
        for u in TYPES[fam]:
            _, du, nb = Jc[u]; coef = {du: Dn[u]}; cc = 0.0
            for (vt, cnt, dv) in nb:
                coef[dv] = coef.get(dv, 0.0) - cnt * Dn[vt] / dv; cc += cnt * Dn[vt] / dv
            for d, x in coef.items(): R.append(row); C.append(d); V.append(x)
            const.append(cc); meta.append((ti, k, TYPES[fam].index(u))); row += 1
    if ti % 100 == 0: print(f'{ti}/{len(grid)} trees, {row} rows ({time.time()-t0:.0f}s)', flush=True)
np.savez_compressed('symrows_ext.npz' if EXT else 'symrows.npz', R=np.array(R), C=np.array(C), V=np.array(V), const=np.array(const), meta=np.array(meta))
open('symrows_ext_grid.txt' if EXT else 'symrows_grid.txt', 'w').write('\n'.join(f'{i}\t{f}\t{a}' for i, (f, a) in enumerate(grid)))
print(f'done: {len(grid)} trees, {row} rows ({time.time()-t0:.0f}s)')
