"""Refined slope diagnostics relative to the Levit-Mandrescu threshold.

thr = ceil((2 alpha - 1)/3).  Rises i_j < i_{j+1} can only occur for j <= thr-1,
so a valley needs slope_k := mu_{k+1}-mu_k > 1 at some k <= thr-2 (with mu_k<k+1).

For a polynomial we report:
  first_pos : min (k - thr) over k with slope_k > 0      (None if never)
  first_gt1 : min (k - thr) over k with slope_k > 1      (None if never)
  max_slope_prefix : max slope_k over k <= thr-2
"""
from fractions import Fraction
from mason_slope_20260902 import mu_seq

def thr_of(alpha):
    return -(-(2 * alpha - 1) // 3)

def diag(poly):
    alpha = len(poly) - 1
    thr = thr_of(alpha)
    mu = mu_seq(poly)
    sl = [mu[k + 1] - mu[k] for k in range(len(mu) - 1)]
    first_pos = min((k - thr for k, s in enumerate(sl) if s > 0), default=None)
    first_gt1 = min((k - thr for k, s in enumerate(sl) if s > 1), default=None)
    pref = [s for k, s in enumerate(sl) if k <= thr - 2]
    mx = max(pref) if pref else None
    kmx = pref.index(mx) if pref else None
    return dict(alpha=alpha, thr=thr, first_pos=first_pos, first_gt1=first_gt1,
                max_slope_prefix=(float(mx) if mx is not None else None), k_at=kmx,
                max_slope_all=float(max(sl)) if sl else None)
