"""Cross-validate slib.analyse against dlib.rows_for (pv_lib, independent DP):
exact interval ends via Fractions from dlib's D values, all trees 4<=n<=NMAX
plus hub families; also time slib on MSH(10^8) and report its interval.
Run from lp-discharge/:  ../../../venv/bin/python adversary/validate.py 12"""
import sys, time, json
from fractions import Fraction
sys.path.insert(0, 'adversary')
import dlib
from slib import analyse, msh

def interval_dlib(adj):
    lo, hi, empty = None, None, False
    deg = [len(x) for x in adj]
    viol = 0
    for (k, D, Dn) in dlib.rows_for(adj):
        for u in range(len(adj)):
            A = sum(Fraction(D[v], deg[v]) for v in adj[u])
            a = D[u] - A; b = A
            if a > 0:
                x = -b / a; hi = x if hi is None or x < hi else hi
            elif a < 0:
                x = b / (-a); lo = x if lo is None or x > lo else lo
            elif b > 0:
                empty = True
            if Fraction(2, 3) * D[u] + Fraction(1, 3) * A > 0: viol += 1
    return lo, hi, empty, viol

def frac(p):
    return Fraction(int(p[0]), int(p[1]))

NMAX = int(sys.argv[1])
cnt = 0; t0 = time.time(); glo, ghi = None, None
for n in range(4, NMAX + 1):
    for adj in dlib.gentreeg(n):
        r = analyse(adj); lo, hi, emp, v = interval_dlib(adj)
        assert (lo is None) == ('low' not in r) and (hi is None) == ('up' not in r), n
        if lo is not None: assert frac(r['low']) == lo, (n, lo, r['low'])
        if hi is not None: assert frac(r['up']) == hi, (n, hi, r['up'])
        assert emp == (r['empty_row'] is not None) and v == r['n_viol23']
        if lo is not None: glo = lo if glo is None or lo > glo else glo
        if hi is not None: ghi = hi if ghi is None or hi < ghi else ghi
        cnt += 1
print(f'all {cnt} trees 4..{NMAX} agree exactly ({time.time()-t0:.1f}s); intersection [{float(glo):.5f}, {float(ghi):.5f}]')
F = dlib.families(); F['R1killer'] = dlib.load_r1_killer()
for nm in ['H(9,2)', 'SoS(4,3)', 'MSH(9,9)', 'MSH(5,5,5,5,5,5,5,5)']:
    adj = F[nm]; r = analyse(adj); lo, hi, emp, v = interval_dlib(adj)
    assert frac(r['low']) == lo and frac(r['up']) == hi and v == r['n_viol23'], nm
    print(f'{nm}: agree, [{float(lo):.5f}, {float(hi):.5f}]')
t0 = time.time(); r = analyse(msh([10] * 8)); print(f"MSH(10^8) n={r['n']}: [{r['low_f']:.5f} at {r['low_at']}, {r['up_f']:.5f} at {r['up_at']}], maxL23={r['maxL23']:.4g}  ({time.time()-t0:.2f}s)")
t0 = time.time(); adj = F['R1killer']; r = analyse(adj); print(f"R1killer: [{r['low_f']:.5f}, {r['up_f']:.5f}] ({time.time()-t0:.2f}s)")
