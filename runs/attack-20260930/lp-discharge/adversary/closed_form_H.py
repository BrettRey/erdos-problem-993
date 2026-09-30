"""Third, closed-form check for hub-stars H(m,s) (no DP, no shared code): by hand,
I = P^m + x(1+x)^{sm}, P = (1+x)^s + x;  j_hub = (1+x)^{sm};  j_sup = P^{m-1};
j_leaf = (1+x)^{s-1} (P^{m-1} + x (1+x)^{s(m-1)}).  Coefficients by binomial sums.
Hub row at sigma=2/3: 3 L_hub = 2 D_hub + m D_sup/(s+1) (hub has m neighbours, all supports of degree s+1)."""
import sys
from math import comb
from fractions import Fraction
m, s = int(sys.argv[1]), int(sys.argv[2])
def c(N, k): return comb(N, k) if 0 <= k <= N else 0
def Ppow(e, k): return sum(c(e, a) * c(s * (e - a), k - a) for a in range(0, min(e, k) + 1))
def i(k): return Ppow(m, k) + c(s * m, k - 1)
def jh(k): return c(s * m, k)
def js(k): return Ppow(m - 1, k)
def jl(k): return sum(c(s - 1, t) * (Ppow(m - 1, k - t) + c(s * (m - 1), k - t - 1)) for t in range(0, s))
n = 1 + m * (s + 1); alpha = s * m + 1
assert i(alpha) > 0 and i(alpha + 1) == 0
lo = -(-n // 4); hi = min(-(-(2 * alpha - 1) // 3), alpha - 1)
tot = 0; up = None
for k in range(lo, hi + 1):
    D = lambda j: k * i(k - 1) * j(k) - (k + 1) * i(k) * j(k - 1)
    Dh, Ds, Dl = D(jh), D(js), D(jl)
    assert Dh + m * Ds + m * s * Dl == k * (k + 1) * (i(k - 1) * i(k + 1) - i(k) ** 2)
    L3 = 2 * Dh + Fraction(m * Ds, s + 1)
    if L3 > 0: tot += 1; print(f'  hub violation k={k}: 3*L_hub/(k i_(k-1) i_k) = {float(L3 / (k * i(k-1) * i(k))):.4e}')
    a = Dh - Fraction(m * Ds, s + 1)
    if a > 0:
        x = Fraction(-m * Ds, s + 1) / a
        up = (x, k) if up is None or x < up[0] else up
print(f'H({m},{s}) n={n} alpha={alpha} window=[{lo},{hi}] hub violations at sigma=2/3: {tot}; hub upper end {float(up[0]):.6f} at k={up[1]}')
