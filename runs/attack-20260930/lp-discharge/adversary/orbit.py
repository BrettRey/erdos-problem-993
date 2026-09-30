"""Orbit-compressed exact evaluator for S(sigma) on symmetric rooted trees.
A spec is a nested list of (multiplicity, childspec); [] is a leaf. Example:
  H(m,s) = [(m, [(s, [])])];  MSH(h;m,s) = [(h, [(m, [(s, [])])])].
All vertices in one orbit (same spec node) share I(T - N[v]), so one rerooting
pass over the spec nodes with polynomial powers gives every D_v exactly.
Same rerooting identities as r1lib.all_seqs (A = I(sub - root), B = x prod A_c,
Qu[c] = excl(A+B) * Pu[t], Pu[c] = Qu[c] + x excl(A) Qu[t], j = R * Qu).
Output format matches slib.analyse. Exact integers throughout; floats only in
the labelled ranking fields (low_f, up_f, maxL23)."""
from math import lcm
import sys as _s; _s.set_int_max_str_digits(0)
import flint
from slib import win, lt, B as Builder

X = flint.fmpz_poly([0, 1]); ONE = flint.fmpz_poly([1])

def flatten(spec):
    """types: list of dict(parent, mult, children=[(idx, mult)])"""
    T = [dict(parent=-1, mult=1, children=[])]
    stack = [(0, spec)]
    while stack:
        t, sp = stack.pop()
        for mu, ch in sp:
            i = len(T); T.append(dict(parent=t, mult=mu, children=[])); T[t]['children'].append((i, mu)); stack.append((i, ch))
    return T

def ppow(p, e):
    return p ** e if e else ONE

def coeffs(p):
    return [int(c) for c in p.coeffs()]

def seqs(spec):
    T = flatten(spec); nt = len(T)
    order = list(range(nt))  # parents precede children by construction
    A = [None] * nt; Bp = [None] * nt; R = [None] * nt
    for t in reversed(order):
        a, r = ONE, ONE
        for c, mu in T[t]['children']:
            a = a * ppow(A[c] + Bp[c], mu); r = r * ppow(A[c], mu)
        A[t] = a; R[t] = r; Bp[t] = X * r
    Pu = [None] * nt; Qu = [None] * nt; Pu[0] = ONE; Qu[0] = ONE
    for t in order:
        ch = T[t]['children']
        for i, (c, mu) in enumerate(ch):
            eAB, eA = ONE, ONE
            for j, (c2, mu2) in enumerate(ch):
                e = mu2 - (1 if j == i else 0)
                eAB = eAB * ppow(A[c2] + Bp[c2], e); eA = eA * ppow(A[c2], e)
            Qu[c] = eAB * Pu[t]; Pu[c] = Qu[c] + X * eA * Qu[t]
    I = coeffs(A[0] + Bp[0]); J = [coeffs(R[t] * Qu[t]) for t in range(nt)]
    cnt = [1] * nt; deg = [0] * nt
    for t in order:
        if T[t]['parent'] >= 0: cnt[t] = cnt[T[t]['parent']] * T[t]['mult']; deg[t] += 1
        deg[t] += sum(mu for _, mu in T[t]['children'])
    co = lambda s, k: s[k] if 0 <= k < len(s) else 0
    for k in range(1, len(I) + 1):
        assert k * co(I, k) == sum(cnt[t] * co(J[t], k - 1) for t in range(nt)), ('double count', k)
    return T, I, J, cnt, deg

def analyse_spec(spec):
    T, I, J, cnt, deg = seqs(spec); nt = len(T)
    n = sum(cnt); alpha = len(I) - 1; lo, hi = win(n, alpha)
    co = lambda s, k: s[k] if 0 <= k < len(s) else 0
    M = 1
    for d in set(deg): M = lcm(M, d)
    low = up = empty_row = None; viol = []; best23 = (-1e300, None); lcfail = []
    for k in range(lo, hi + 1):
        c1, c2 = k * co(I, k - 1), (k + 1) * co(I, k)
        D = [c1 * co(J[t], k) - c2 * co(J[t], k - 1) for t in range(nt)]
        sD = sum(cnt[t] * D[t] for t in range(nt))
        assert sD == k * (k + 1) * (co(I, k - 1) * co(I, k + 1) - co(I, k) ** 2)
        if sD > 0: lcfail.append(k)
        WD = [(M // deg[t]) * D[t] for t in range(nt)]
        nrm = 3 * M * k * co(I, k - 1) * co(I, k)
        for t in range(nt):
            b = WD[T[t]['parent']] if T[t]['parent'] >= 0 else 0
            for c, mu in T[t]['children']: b += mu * WD[c]
            a = M * D[t] - b
            if a > 0:
                cand = (-b, a)
                if up is None or lt(cand, up[:2]): up = (cand[0], cand[1], k, t)
            elif a < 0:
                cand = (b, -a)
                if low is None or lt(low[:2], cand): low = (cand[0], cand[1], k, t)
            elif b > 0 and empty_row is None: empty_row = (k, t)
            s23 = 2 * a + 3 * b
            if s23 > 0: viol.append((t, k, deg[t]))
            f = s23 / nrm
            if f > best23[0]: best23 = (f, (t, k, deg[t]))
    out = dict(n=n, alpha=alpha, window=[lo, hi], lcfail=lcfail, viol23=viol[:20], n_viol23=len(viol),
               maxL23=best23[0], maxL23_at=best23[1], empty_row=empty_row, ntypes=nt)
    if low: out.update(low=[str(low[0]), str(low[1])], low_f=low[0] / low[1], low_at=(low[3], low[2], deg[low[3]]))
    else: out['low_f'] = -1e300
    if up: out.update(up=[str(up[0]), str(up[1])], up_f=up[0] / up[1], up_at=(up[3], up[2], deg[up[3]]))
    else: out['up_f'] = 1e300
    return out

def expand(spec):
    """explicit adjacency (vertex 0 = root) for exact rechecks"""
    b = Builder()
    def rec(v, sp):
        for mu, ch in sp:
            for _ in range(mu):
                w = b.add(v); rec(w, ch)
    rec(0, spec)
    return b.adj()

def H(m, s): return [(m, [(s, [])])]
def MSH(h, m, s): return [(h, [(m, [(s, [])])])]
