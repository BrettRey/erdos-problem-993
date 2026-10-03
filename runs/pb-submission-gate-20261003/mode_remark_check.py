"""Exact test of the proposed mode remark: V*delta_c >= 1/5, c = D-1.

Derivation under test: the right inequality of (2.4) at k = D-1 (needs D >= 2)
gives delta_{D-1} >= delta_D/(1+delta_D); with Theorem 1.1 this gives
V delta_{D-1} >= 1/5. For D = 1, delta_0 = 1 and V delta_0 = V >= 1.
"""
from fractions import Fraction as F
import random

random.seed(20261003)


def pmf(ps):
    f = [F(1)]
    for p in ps:
        g = [F(0)] * (len(f) + 1)
        for k, x in enumerate(f):
            g[k] += x * (1 - p)
            g[k + 1] += x * p
        f = g
    return f


def deltas(f):
    n = len(f) - 1
    ext = [F(0)] + f + [F(0)]
    return [1 - ext[k] * ext[k + 2] / ext[k + 1] ** 2 for k in range(n + 1)]


worst_c = None
worst_D = None
viol_lemma = 0
tested = 0
for trial in range(4000):
    n = random.randint(2, 40)
    mode = random.random()
    ps = []
    for _ in range(n):
        if mode < 0.3:
            d = random.choice([7, 13, 29, 61, 127])
            x = F(random.randint(1, d - 1), d)
        elif mode < 0.6:
            x = F(random.choice([1, 2, 3]), random.choice([50, 100, 400]))
            if random.random() < 0.5:
                x = 1 - x
        else:
            x = F(random.randint(1, 999), 1000)
        ps.append(x)
    V = sum(p * (1 - p) for p in ps)
    if V < 1:
        continue
    f = pmf(ps)
    dl = deltas(f)
    D = next(k for k in range(1, len(f)) if f[k] < f[k - 1])
    c = D - 1
    tested += 1
    vD = V * dl[D]
    vc = V * dl[c]
    if D >= 2 and dl[D - 1] < dl[D] / (1 + dl[D]):
        viol_lemma += 1
    if worst_D is None or vD < worst_D[0]:
        worst_D = (vD, n, D)
    if worst_c is None or vc < worst_c[0]:
        worst_c = (vc, n, D)

print("tested", tested)
print("lemma violations", viol_lemma)
print("min V*delta_D", float(worst_D[0]), worst_D[1:])
print("min V*delta_c", float(worst_c[0]), worst_c[1:])
assert viol_lemma == 0 and worst_D[0] >= F(1, 4) and worst_c[0] >= F(1, 5)
print("OK")
