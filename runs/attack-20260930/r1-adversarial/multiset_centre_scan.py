"""Centre-only exact scan (closed form, flint) of stars of hubs with arbitrary hub sizes:
centre joined to hubs with m_1..m_h cherries (vertex + 2 leaves). Exact integer sign test of L_centre.
Usage: python multiset_centre_scan.py NMAX > data/multiset_centre_scan.jsonl"""
import itertools, json, sys
import flint
X = flint.fmpz_poly([0, 1]); O = flint.fmpz_poly([1, 1]); Ic = O * O + X
def co(p, k): return int(p[k]) if 0 <= k <= p.degree() else 0
def centre_check(ms):
    h = len(ms)
    A = [Ic ** m + X * O ** (2 * m) for m in ms]; Bs = [Ic ** m for m in ms]
    pA = flint.fmpz_poly([1]); pB = flint.fmpz_poly([1])
    for a in A: pA *= a
    for b in Bs: pB *= b
    I = pA + X * pB; Jc = pB
    Jh = []
    for i, m in enumerate(ms):
        p = O ** (2 * m)
        for j, a in enumerate(A):
            if j != i: p *= a
        Jh.append(p)
    n = 1 + sum(3 * m + 1 for m in ms); alpha = I.degree()
    lo, q = -((-n) // 4), -((-(2 * alpha - 1)) // 3)
    viol = []; bestF = -1e300
    from math import lcm
    M = h + 1
    for m in ms: M = lcm(M, m + 2)
    for k in range(lo, q + 1):
        a, b = k * co(I, k - 1), (k + 1) * co(I, k)
        LM = (a * co(Jc, k) - b * co(Jc, k - 1)) * (M // (h + 1))
        G2 = b * co(Jc, k - 1) * (M // (h + 1))
        for i, m in enumerate(ms):
            LM += (a * co(Jh[i], k) - b * co(Jh[i], k - 1)) * (M // (m + 2))
            G2 += b * co(Jh[i], k - 1) * (M // (m + 2))
        if LM > 0: viol.append(k)
        bestF = max(bestF, (k + 1) * LM / G2)
    return n, alpha, [lo, q], viol, bestF
if __name__ == '__main__':
    nmax = int(sys.argv[1])
    for h in range(5, 11):
        for ms in itertools.combinations_with_replacement(range(7, 14), h):
            n = 1 + sum(3 * m + 1 for m in ms)
            if n > nmax: continue
            n, alpha, w, viol, F = centre_check(list(ms))
            print(json.dumps({"ms": list(ms), "n": n, "alpha": alpha, "window": w, "viol_k": viol, "F_centre": F}), flush=True)
