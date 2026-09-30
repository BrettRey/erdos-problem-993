"""m -> infinity limit of R1 loads on H(m,s) and TH(m,s): the hub is free with probability
<= (B/A)^m-ish (exponentially small), so the tree is a forest of m stars K_{1,s} plus a blocked hub.
Per star (Fractions, exact at rational lam): A = (1+lam)^s + lam, t = mean, var = variance,
pi_mid = 1/A, pi_leaf = (1+lam)^(s-1)/A, Delta_mid = -t, Delta_leaf = (s-1)lam/(1+lam) - t,
eps_v = Delta_v/(var/t);  n*kappa*D_v/(k i_{k-1} i_k) -> pi_v (eps_v - 1).
Window: t in [(s+1)/4, 2s/3]. Reports sup over window of: eps_mid, eps_leaf, and the
mean|D|-normalised loads of mid and leaf (to compare with the exact large-m scans)."""
import json
from fractions import Fraction as Fr

def star(s, lam):
    A = (1 + lam) ** s + lam
    Ap = s * (1 + lam) ** (s - 1) + 1
    App = s * (s - 1) * (1 + lam) ** (s - 2) if s >= 2 else 0
    t = lam * Ap / A
    # var = lam d t / d lam
    var = t + lam * lam * App / A - t * t
    return A, t, var

for s in (1, 2, 3, 4):
    tlo, thi = Fr(s + 1, 4), Fr(2 * s, 3)
    best = {"eps_mid": None, "eps_leaf": None, "mid_mean_norm": None, "leaf_mean_norm": None}
    for i in range(1, 2000):
        lam = Fr(i, 200)
        A, t, var = star(s, lam)
        if not (tlo <= t <= thi):
            continue
        r = var / t
        pc, pl = 1 / A, (1 + lam) ** (s - 1) / A
        ec, el = -t / r, (Fr(s - 1) * lam / (1 + lam) - t) / r
        dc, dl = pc * (ec - 1), pl * (el - 1)
        mean = (abs(dc) + s * abs(dl)) / (s + 1)
        vals = {"eps_mid": ec, "eps_leaf": el, "mid_mean_norm": (dc / (s + 2) + s * dl / 2) / mean,
                "leaf_mean_norm": (dl / 2 + dc / (s + 2)) / mean}
        for k2, v in vals.items():
            if best[k2] is None or v > best[k2][0]:
                best[k2] = (v, lam)
    print(json.dumps({"s": s, **{k2: [round(float(v[0]), 6), str(v[1])] for k2, v in best.items()}}))
