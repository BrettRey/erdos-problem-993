"""End-to-end exact-arithmetic audit.

Part 1: Proposition 3.1 tested directly from the definitions (K, L_r, R_r,
A(delta)) with Fractions, at dense random and adversarial delta in (0, 1/4),
including points straddling every cell boundary H = m and triangular-cell
boundaries. Independent of the paper's cell decomposition.

Part 2: Theorem 1.1 stress-tested on random Poisson-binomial laws with exact
rational pmfs: check V*delta_D >= 1/4 whenever V >= 1, track the minimum.
Also check the log-concavity-based tail corollary on each sampled law.
"""
from fractions import Fraction as F
import random

random.seed(993)

def audit_scalar(delta):
    assert F(0) < delta < F(1, 4)
    # K = max{r>=1 : (r+1) delta < 1}
    K = 1
    while (K + 2) * delta < 1:
        K += 1
    assert (K + 1) * delta < 1 <= (K + 2) * delta
    a = (1 - 2 * delta) / (1 - delta)
    # incremental L_r, R_r (O(K) exact rational work)
    # R_r = a^r prod_{j=1}^{r-1}(1-j d);  L_r = (1-d)^{-r} prod_{j=2}^{r+1}(1-j d)
    S0 = S1 = S2 = F(0)   # sum w, sum i w, sum i^2 w over i = -K..K
    Rr, Lr = F(1), F(1)
    S0 += 1
    for r in range(1, K + 1):
        Rr *= a * (1 - (r - 1) * delta)
        Lr *= (1 - (r + 1) * delta) / (1 - delta)
        assert Rr > 0 and Lr > 0
        S0 += Rr + Lr
        S1 += r * Rr - r * Lr
        S2 += r * r * (Rr + Lr)
    A = S0 * S2 - S1 * S1   # = (1/2) sum_{i,j} w_i w_j (i-j)^2
    target = (3 + delta) / (4 * delta ** 2)
    return A >= target, A - target

fails = 0
tested = 0
margin_min = None
# adversarial: straddle K-transition points delta = 1/(m+1) (H = m) and 1/4, plus
# triangular boundaries delta = 1/(T_J+1)
special = []
for m in range(4, 60):
    for eps_num in (1, 7, 123):
        eps = F(eps_num, 10 ** 9)
        special += [F(1, m + 1) - eps, F(1, m + 1) + eps]
special += [F(1, 4) - F(1, 10 ** 12), F(1, 2001), F(1, 703)]
for J in range(5, 12):
    T = J * (J + 1) // 2
    special += [F(1, T + 2), F(1, T + 1), F(2, 2 * T + 1)]
rand = [F(random.randrange(1, 10 ** 9), 4 * 10 ** 9 + random.randrange(1, 10 ** 9))
        for _ in range(400)]
for d in special + rand:
    if not (F(1, 3000) < d < F(1, 4)):
        continue  # depth cap: smaller delta = larger K; that regime is the
                  # symbolically verified J>=6 family, spot-checked at 1/2001
    tested += 1
    passed, margin = audit_scalar(d)
    if not passed:
        fails += 1
        print("SCALAR FAIL at delta =", d)
    if margin_min is None or margin < margin_min:
        margin_min, argmin = margin, d
print(f"Part 1: Prop 3.1 direct exact test: {tested} deltas, {fails} failures")
print(f"        smallest margin A - (3+d)/(4d^2) = {float(margin_min):.6g} at delta = {argmin} (~{float(argmin):.6g})")

# ---------- Part 2: the theorem itself ----------
def pmf_exact(ps):
    f = [F(1)]
    for p in ps:
        g = [F(0)] * (len(f) + 1)
        for k, fk in enumerate(f):
            g[k] += fk * (1 - p)
            g[k + 1] += fk * p
        f = g
    return f

def theorem_check(ps):
    f = pmf_exact(ps)
    n = len(ps)
    V = sum(p * (1 - p) for p in ps)
    D = None
    for k in range(1, n + 1):
        if f[k] < f[k - 1]:
            D = k
            break
    assert D is not None, "no descent despite V>=1"
    fR = f[D + 1] if D + 1 <= n else F(0)
    delta = 1 - (f[D - 1] * fR) / (f[D] ** 2)
    score = V * delta
    # corollary checks: ratio drop and geometric tail (exact, using rho=1-1/(4V))
    ratio_ok = V * (1 - fR / f[D]) >= F(1, 4)
    rho = 1 - 1 / (4 * V)
    tail_ok = all(f[D + r] * f[D] ** 0 <= (rho ** r) * f[D] for r in range(1, n - D + 1))
    return score, ratio_ok, tail_ok

def rand_law():
    kind = random.choice(["uniform", "lopsided", "nearhalf", "mixture", "tiny"])
    n = random.randrange(2, 26)
    ps = []
    for _ in range(n):
        if kind == "uniform":
            ps.append(F(random.randrange(1, 999), 1000))
        elif kind == "lopsided":
            ps.append(random.choice([F(random.randrange(1, 60), 1000),
                                     F(random.randrange(940, 999), 1000)]))
        elif kind == "nearhalf":
            ps.append(F(500 + random.randrange(-40, 41), 1000))
        elif kind == "mixture":
            ps.append(random.choice([F(1, 50), F(49, 50), F(1, 2),
                                     F(random.randrange(1, 999), 1000)]))
        else:
            ps.append(F(random.randrange(1, 200), 1000))
    return ps

worst = None
checked = 0
viol = 0
for _ in range(600):
    ps = rand_law()
    V = sum(p * (1 - p) for p in ps)
    if V < 1:
        continue
    checked += 1
    score, ratio_ok, tail_ok = theorem_check(ps)
    if score < F(1, 4) or not ratio_ok or not tail_ok:
        viol += 1
        print("THEOREM FAIL:", [str(p) for p in ps], float(score), ratio_ok, tail_ok)
    if worst is None or score < worst[0]:
        worst = (score, len(ps))
# the paper's variance-one binomial family, exact-rational surrogate:
# choose rational p near (1-sqrt(1-4/n))/2 keeping V >= 1
import math
fam = []
for n in range(5, 90, 7):
    p_star = (1 - math.sqrt(1 - 4 / n)) / 2
    p = F(math.ceil(p_star * 10 ** 7), 10 ** 7)  # round UP so V >= 1
    ps = [p] * n
    V = sum(q * (1 - q) for q in ps)
    assert V >= 1
    score, ratio_ok, tail_ok = theorem_check(ps)
    assert score >= F(1, 4) and ratio_ok and tail_ok
    fam.append((n, float(score)))
    checked += 1
print(f"Part 2: theorem exact stress test: {checked} laws with V>=1, {viol} violations")
print(f"        worst random score {float(worst[0]):.6f} (n={worst[1]})")
print("        variance-~1 binomial family scores:", ", ".join(f"n={n}:{s:.5f}" for n, s in fam))
print("        limit (n+1)/(3(n-1)) -> 1/3 = 0.33333")
print("ALL END-TO-END CHECKS COMPLETE")
