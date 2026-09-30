"""Compact row store. Each (tree, k, u) row becomes a vector over degrees:
L_u/norm = sum_d sigma[d] * M[row, d] + c[row], with
M[row, d] = Dn_u [d_u == d] - A[row, d]/d,  c[row] = sum_d A[row, d]/d,
A[row, d] = sum of Dn_v over neighbours v of u with degree d."""
import sys, time, pickle, numpy as np
import dlib
DMAX = 40
def rows_matrix(adj):
    deg = [len(x) for x in adj]; M = []; c = []
    for (k, D, Dn) in dlib.rows_for(adj):
        for u in range(len(adj)):
            m = np.zeros(DMAX + 1); cc = 0.0
            m[deg[u]] += Dn[u]
            for v in adj[u]:
                m[deg[v]] -= Dn[v] / deg[v]; cc += Dn[v] / deg[v]
            M.append(m); c.append(cc)
    return M, c
if __name__ == '__main__':
    lo, hi = int(sys.argv[1]), int(sys.argv[2]); t0 = time.time()
    store = {}
    for n in range(lo, hi + 1):
        M, c = [], []
        for adj in dlib.gentreeg(n):
            m, cc = rows_matrix(adj); M += m; c += cc
        store[f'n{n}'] = (np.array(M, dtype=np.float64), np.array(c))
        print(f'n={n}: {len(c)} rows ({time.time()-t0:.0f}s)', flush=True)
    if lo == 4:
        fam = dlib.families(); fam['R1killer'] = dlib.load_r1_killer()
        for nm, adj in fam.items():
            m, cc = rows_matrix(adj); store['F:' + nm] = (np.array(m), np.array(cc))
        print(f'families done ({time.time()-t0:.0f}s)', flush=True)
    pickle.dump(store, open(f'rows_{lo}_{hi}.pkl', 'wb'))
