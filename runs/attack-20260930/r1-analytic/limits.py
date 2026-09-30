"""Grand-canonical (linear-regime) limit criteria, evaluated EXACTLY (Fractions) at rational fugacities.

Linear regime (all pi_v = O(1), n -> infinity, k/n -> kappa):
  n*kappa*D_v(k)/(k i_{k-1} i_k) -> pi_v (eps_v - 1),
  eps_v = (E[|S| | v free] - E|S|) / (Var|S| / E|S|).
Star-of-hubs SH(h,m,s), m,s fixed, h -> infinity: units U = A^m + x B^m are i.i.d.,
centre z is free with probability -> 0 exponentially in h, so
  sign L_z -> sign(eps_hub - 1),  eps_hub - 1 > 0  <=>  m s lam/(1+lam) - t(lam) > Var_U/t(lam),
with t = per-unit mean size, Var_U = per-unit variance.
Window (h -> inf): per-unit t in [(1+m(s+1))/4, 2(sm+1)/3].
"""
import sys, json
from fractions import Fraction as Fr
from flint import fmpz_poly
X = fmpz_poly([0, 1])

def moments(P, lam):
    c = [int(P[j]) for j in range(P.degree() + 1)]
    Z = M1 = M2 = Fr(0); p = Fr(1)
    for j, a in enumerate(c):
        Z += a * p; M1 += j * a * p; M2 += j * j * a * p; p *= lam
    t = M1 / Z
    return t, M2 / Z - t * t

def sh_hub(m, s, lam):
    A = (1 + X) ** s + X; B = (1 + X) ** s
    U = A ** m + X * B ** m
    t, var = moments(U, lam)
    delta = Fr(m * s) * lam / (1 + lam) - t
    return t, var, delta, delta - var / t   # last > 0  <=>  eps_hub > 1

def grid():
    return [Fr(i, 20) for i in range(1, 400)]  # lam in (0, 20)

if __name__ == "__main__":
    out = []
    for s in (1, 2, 3, 4):
        for m in range(2, 31):
            tlo = Fr(1 + m * (s + 1), 4); thi = Fr(2 * (s * m + 1), 3)
            win_pos = []; best = None
            for lam in grid():
                t, var, delta, crit = sh_hub(m, s, lam)
                if tlo <= t <= thi:
                    eps = delta / (var / t)
                    if best is None or eps > best[1]: best = (lam, eps)
                    if crit > 0: win_pos.append(lam)
            rec = {"s": s, "m": m, "max_eps_hub_in_window_float": float(best[1]) if best else None,
                   "lam_at_max": str(best[0]) if best else None,
                   "eps_gt_1_lam_range_on_grid": [str(win_pos[0]), str(win_pos[-1])] if win_pos else None,
                   "predicts_R1_failure_for_large_h": bool(win_pos)}
            out.append(rec); print(json.dumps(rec), flush=True)
