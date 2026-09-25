"""Slope conjecture tester (polynomial level).

mu_k := (k+1) i_{k+1} / i_k  (mean number of free vertices over uniform k-independent sets).
Slope conjecture:  mu_{k+1} - mu_k <= 1 for all 0 <= k <= alpha-2.
Equivalently:      (k+2) i_{k+2} i_k <= (k+1) i_{k+1}^2 + i_k i_{k+1}.
It implies unimodality: mu_k < k+1  =>  mu_{k+1} < k+2.

Also reports the 'Newton-type' statistic max slope (mu decreasing <=> slope <= 0).
All arithmetic exact (Fractions / integers).
"""
from fractions import Fraction

def mu_seq(poly):
    a = len(poly) - 1
    return [Fraction((k + 1) * poly[k + 1], poly[k]) for k in range(a)]

def slopes(poly):
    m = mu_seq(poly)
    return [m[k + 1] - m[k] for k in range(len(m) - 1)]

def max_slope(poly):
    s = slopes(poly)
    if not s:
        return None, None
    v = max(s)
    return v, s.index(v)

def is_unimodal(a):
    i = 0
    while i + 1 < len(a) and a[i] <= a[i + 1]:
        i += 1
    while i + 1 < len(a) and a[i] >= a[i + 1]:
        i += 1
    return i + 1 >= len(a)

def lc_breaks(a):
    return [k for k in range(1, len(a) - 1) if a[k] * a[k] < a[k - 1] * a[k + 1]]

def summarize(poly, label=""):
    v, k = max_slope(poly)
    alpha = len(poly) - 1
    return dict(label=label, n=None, alpha=alpha, max_slope=float(v) if v is not None else None,
                at=k, lc_breaks=lc_breaks(poly), unimodal=is_unimodal(poly))
