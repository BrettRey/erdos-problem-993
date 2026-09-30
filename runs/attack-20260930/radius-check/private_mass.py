"""Exact sufficient condition for a tree-specific radius-1 certificate in layered
trees: every positive vertex at depth i can be paid by its OWN descendants within
distance 2 (children and grandchildren, not shared with any other positive vertex
at depth i): D_i <= b_i*|D_{i+1}|*[D_{i+1}<0] + b_i*b_{i+1}*|D_{i+2}|*[D_{i+2}<0].
Integer arithmetic only. Also reports the size of the raw defects (float overflow check)."""
from math import ceil
from layered import layered, coeffs
co = lambda p, k: p[k] if 0 <= k < len(p) else 0
def check(b):
    I, J, N, deg = layered(b); I = coeffs(I); Jc = [coeffs(x) for x in J]
    n = sum(N); a = len(I) - 1; q = ceil((2 * a - 1) / 3); lo = ceil(n / 4); L = len(b)
    fails = 0; tested = 0; maxdigits = 0; worst = None
    for k in range(lo, min(q, a - 1) + 1):
        D = [k * co(I, k - 1) * co(Jc[i], k) - (k + 1) * co(I, k) * co(Jc[i], k - 1) for i in range(L + 1)]
        maxdigits = max(maxdigits, max(len(str(abs(d))) for d in D))
        for i, d in enumerate(D):
            if d <= 0: continue
            tested += 1
            pay = 0
            if i + 1 <= L and D[i + 1] < 0: pay += b[i] * -D[i + 1]
            if i + 2 <= L and D[i + 2] < 0: pay += b[i] * b[i + 1] * -D[i + 2]
            if pay < d: fails += 1
            r = pay // d if pay >= d else 0
            worst = r if worst is None else min(worst, r)
    return n, tested, fails, worst, maxdigits
for name, b in [('H(9,2)', [9, 2]), ('H(75,5)', [75, 5]), ('MSH(8;10,2)', [8, 10, 2]), ('MSH(38;11,2)', [38, 11, 2]),
                ('MSH(96;11,2)', [96, 11, 2]), ('MSH(120;14,3)', [120, 14, 3]), ('depth-6 worst', [2, 1, 4, 3, 6, 2])]:
    n, t, f, w, dg = check(b)
    print(f'{name}: n={n} positive (level,k) cases={t} not privately payable={f} min floor(private negative/positive)={w} raw defect digits up to {dg}')
