"""Maximise the certificate margin t over structured rule families, using every
row in the store (all trees n <= 16 + families incl. the n=237 R1 killer).
Row value (normalized L_u) = M @ sigma + c; certificate needs <= -t < 0."""
import pickle, numpy as np, json
from scipy.optimize import linprog
from scipy import sparse
S = pickle.load(open('rows_4_16.pkl', 'rb'))
M = sparse.csr_matrix(np.vstack([v[0] for v in S.values()])); c = np.concatenate([v[1] for v in S.values()])
names = [k for k, v in S.items() for _ in range(len(v[1]))]
D = M.shape[1] - 1; d = np.arange(D + 1); used = sorted(set(M.nonzero()[1])); print('rows', M.shape[0], 'degrees used', used)
def worst(sig):
    val = M @ sig + c; i = int(np.argmax(val)); return val[i], names[i]
def lp(B, extra_ub=None, bounds=None):
    """sigma = B @ theta; maximize t: M B theta + t <= -c."""
    MB = M @ B; k = B.shape[1]
    A = sparse.hstack([sparse.csr_matrix(MB), sparse.csr_matrix(np.ones((M.shape[0], 1)))]).tocsr()
    Aub, bub = [A], [-c]
    # 0 <= sigma(d) <= 1 on degrees 1..D
    Bd = B[1:]; Aub += [sparse.csr_matrix(np.hstack([Bd, np.zeros((D, 1))])), sparse.csr_matrix(np.hstack([-Bd, np.zeros((D, 1))]))]
    bub += [np.ones(D), np.zeros(D)]
    if extra_ub is not None: Aub.append(sparse.csr_matrix(extra_ub[0])); bub.append(extra_ub[1])
    cost = np.zeros(k + 1); cost[-1] = -1
    r = linprog(cost, A_ub=sparse.vstack(Aub).tocsr(), b_ub=np.concatenate(bub), bounds=bounds or [(None, None)] * k + [(None, 1)], method='highs')
    th = r.x[:k]; return r.x[-1], B @ th, th
out = {}
# R1 itself
sR1 = np.array([0.0] + [1 / (x + 1) for x in range(1, D + 1)]); print('R1 (sigma=1/(d+1)) worst row:', worst(sR1))
# free sigma
I = np.eye(D + 1); t, s, _ = lp(I); print(f'free:      t*={t:.5f}  worst={worst(s)}'); out['free'] = dict(t=t, sigma=s[1:].tolist())
# monotone nonincreasing on degrees 1..D
mono = np.zeros((D - 1, D + 2))
for i in range(1, D): mono[i - 1, i] = -1; mono[i - 1, i + 1] = 1   # sigma(i+1) - sigma(i) <= 0
t, s, _ = lp(I, extra_ub=(mono, np.zeros(D - 1))); print(f'monotone:  t*={t:.5f}  sigma[1..16]={np.round(s[1:17],4).tolist()}'); out['monotone'] = dict(t=t, sigma=s[1:].tolist())
# constant sigma
B = np.ones((D + 1, 1)); t, s, th = lp(B); print(f'constant:  t*={t:.5f}  sigma={th[0]:.4f}'); out['constant'] = dict(t=t, theta=th.tolist())
# c/(d+1)
w = np.array([0.0] + [1 / (x + 1) for x in range(1, D + 1)])[:, None]; t, s, th = lp(w); print(f'c/(d+1):   t*={t:.5f}  c={th[0]:.4f}'); out['c_over_d1'] = dict(t=t, theta=th.tolist())
# alpha + beta/(d+1)
B = np.hstack([np.ones((D + 1, 1)), w]); t, s, th = lp(B); print(f'a+b/(d+1): t*={t:.5f}  alpha={th[0]:.4f} beta={th[1]:.4f}  worst={worst(s)}'); out['a_b_over_d1'] = dict(t=t, theta=th.tolist())
json.dump(out, open('lp_forms_result.json', 'w'), indent=1)
