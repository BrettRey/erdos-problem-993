"""Exact (Fraction) check of the LP's degree rule sigma(d) against every H/MSH row
in both grids. sigma values rationalized with limit_denominator(10**9)."""
import json
from math import ceil
from fractions import Fraction
from symrows import H, MSH, coeffs
r = json.load(open('lp_degree_combined_result.json'))
sig = {d: Fraction(s).limit_denominator(10**9) for d, s in zip(r['degrees'], r['sigma'])}
sig = {d: min(max(v, Fraction(0)), Fraction(1)) for d, v in sig.items()}
sig = {d: (Fraction(0) if v < Fraction(1, 10**6) else Fraction(1) if v > 1 - Fraction(1, 10**6) else v) for d, v in sig.items()}
print("snapped sigma:", {d: str(v) for d, v in sorted(sig.items()) if d <= 20})
TYPES = {'H': ['hub', 'sup', 'leaf'], 'MSH': ['cen', 'hub', 'sup', 'leaf']}
co = lambda p, k: p[k] if 0 <= k < len(p) else 0
viol = []; checked = 0
for g in ('symrows_grid.txt', 'symrows_ext_grid.txt'):
    for line in open(g).read().split('\n'):
        ti, fam, args = line.split('\t'); args = eval(args)
        I, J, n = (H if fam == 'H' else MSH)(*args); I = coeffs(I); a = len(I) - 1; q = ceil((2 * a - 1) / 3); lo = ceil(n / 4)
        Jc = {t: (coeffs(Jp), d, nb) for t, (Jp, d, nb) in J.items()}
        for k in range(lo, min(q, a - 1) + 1):
            D = {t: k * co(I, k - 1) * co(Jc[t][0], k) - (k + 1) * co(I, k) * co(Jc[t][0], k - 1) for t in Jc}
            for u in TYPES[fam]:
                _, du, nb = Jc[u]
                L = sig[du] * D[u] + sum(cnt * (1 - sig[dv]) * Fraction(D[vt], dv) for (vt, cnt, dv) in nb)
                checked += 1
                if L > 0: viol.append((fam, args, k, u, float(L / (k * co(I, k - 1) * co(I, k)))))
print(f'rows checked exactly: {checked}; violations: {len(viol)}')
for v in sorted(viol, key=lambda x: -x[4])[:10]: print('  ', v)
json.dump(dict(checked=checked, violations=viol[:500], n_viol=len(viol)), open('exact_check_degree.json', 'w'))
