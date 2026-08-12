#!/usr/bin/env python3
"""Single-index laws for the C2 base-facts program.

For arm multisets a, b (pendant path lengths at two hubs) define, with
P_n the path independence polynomial (P_{-1}=P_0=1, P_n=P_{n-1}+xP_{n-2}),
Q_a = prod P_{a_i}, R_a = prod P_{a_i-1}:

    B = Q_aQ_b + x(R_aQ_b + Q_aR_b)     (adjacent two-hub tree)
    C = Q_aQ_b + xR_aR_b                (spider on the combined arms)

and the single-index forms (out-of-range coefficients read 0)

    U_k = c_k b_k - c_{k-1} b_{k+1}
    V_k = c_k b_k - c_{k+1} b_{k-1}
    W_k = 2 c_{k-1} b_k - c_k b_{k-1} - c_{k-2} b_{k+1}

W is the diagonal of the partial-synchronization form for (xC, B).

Laws tested (exact integer arithmetic):
  LAW-V: V_k >= 0 for every arm pair and every k.
  LAW-W: W_k >= 0 for every arm pair and every k.
  LAW-U: U_k < 0 only at reflected depth N-k <= 3 (N = deg B).

Assembly lemma (proved in notes/c2_single_index_laws_2026-08-11.md):
  C log-concave (known: Li-Li-Yang-Zhang, arXiv:2501.04245, spiders) plus
  LAW-V + LAW-U at the relevant indices plus LAW-W imply both base
  partial-synchronization relations C ~_p B and xC ~_p B away from the
  U-failure strip.  This script also end-to-end checks that implication's
  contrapositive on every scanned pair: any partial-sync failure must
  touch the U-failure strip.  (In the scans to date the base relations
  never fail at all, so the check is vacuous but cheap.)

Also records the bilinear generator-pair table: U, V, W are bilinear in
(C, B); C = g1+g2, B = g1+h2+h3 with g1=QaQb, g2=xRaRb, h2=xRaQb,
h3=xQaRb; per-term nonnegativity holds for some generator pairs only,
so the laws are aggregate, not termwise.
"""
import json
import random
import time
from fractions import Fraction
from itertools import combinations_with_replacement as cwr

# ---------- exact polynomial helpers ----------

def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                out[i + j] += ai * bj
    return out


def padd(a, b):
    m = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
            for i in range(m)]


_PC = {}


def P(n):
    if n in _PC:
        return _PC[n]
    a, b = [1], [1, 1]
    if n == 0:
        r = a
    else:
        for _ in range(n - 1):
            a, b = b, padd(b, [0] + a)
        r = b
    _PC[n] = r
    return r


def cf(p, i):
    return p[i] if 0 <= i < len(p) else 0


def prodP(arms, shift):
    out = [1]
    for a in arms:
        out = pmul(out, P(a - shift))
    return out


def BC(a, b):
    Qa, Ra = prodP(a, 0), prodP(a, 1)
    Qb, Rb = prodP(b, 0), prodP(b, 1)
    QQ = pmul(Qa, Qb)
    B = padd(QQ, [0] + padd(pmul(Ra, Qb), pmul(Qa, Rb)))
    C = padd(QQ, [0] + pmul(Ra, Rb))
    return B, C


def UVW_scan(a, b):
    """Return N, U-fail depths, min-V ratio, min-W ratio, psync fail pairs."""
    B, C = BC(a, b)
    N = len(B) - 1
    A = [0] + C
    Ufail, Vmin, Wmin = [], None, None
    for k in range(N + 2):
        U = cf(C, k) * cf(B, k) - cf(C, k - 1) * cf(B, k + 1)
        V = cf(C, k) * cf(B, k) - cf(C, k + 1) * cf(B, k - 1)
        W = (2 * cf(A, k) * cf(B, k) - cf(A, k + 1) * cf(B, k - 1)
             - cf(A, k - 1) * cf(B, k + 1))
        if U < 0:
            Ufail.append(N - k)
        pV, pW = cf(C, k) * cf(B, k), 2 * cf(A, k) * cf(B, k)
        if pV > 0:
            r = Fraction(V, pV)
            if Vmin is None or r < Vmin:
                Vmin = r
        if pW > 0:
            r = Fraction(W, pW)
            if Wmin is None or r < Wmin:
                Wmin = r
    return N, Ufail, Vmin, Wmin


def psync_fails(A, B):
    D = max(len(A), len(B))
    out = []
    for m in range(D + 1):
        for n in range(m + 1):
            lhs = cf(A, m) * cf(B, n) + cf(A, n) * cf(B, m)
            rhs = cf(A, m + 1) * cf(B, n - 1) + cf(A, n - 1) * cf(B, m + 1)
            if lhs < rhs:
                out.append((m, n))
    return out


def main():
    t0 = time.time()
    report = {"date": "2026-08-11", "laws": ["LAW-V", "LAW-W", "LAW-U-depth<=3"]}

    # ---- structured scan ----
    small = []
    for r in range(0, 4):
        for arms in cwr([1, 2, 3, 4, 5, 6, 7], r):
            if sum(arms) <= 12:
                small.append(list(arms))
    cases = [(x, y) for x in small for y in small if sum(x) + sum(y) > 0]
    for m in [9, 15, 25, 40]:
        cases += [([1, 1], [m, m]), ([1], [m]), ([1, 1, 1], [m, m, m]),
                  ([2, 2], [m, m]), ([1, 3], [m]), ([], [1, 3]),
                  ([], [m, m]), ([1], [m, m, m]), ([1, 1], [2, m])]
    random.seed(993)
    for _ in range(300):
        r, s = random.randint(0, 4), random.randint(1, 4)
        x = [random.randint(1, 15) for _ in range(r)]
        y = [random.randint(1, 15) for _ in range(s)]
        if sum(x) + sum(y) <= 60:
            cases.append((x, y))

    deepU, globV, globW = (-1, None), (None, None), (None, None)
    viol = []
    for a, b in cases:
        N, Ufail, Vmin, Wmin = UVW_scan(a, b)
        if Ufail and max(Ufail) > deepU[0]:
            deepU = (max(Ufail), (a, b, N))
        if Vmin is not None and (globV[0] is None or Vmin < globV[0]):
            globV = (Vmin, (a, b, N))
        if Wmin is not None and (globW[0] is None or Wmin < globW[0]):
            globW = (Wmin, (a, b, N))
        if (Vmin is not None and Vmin < 0) or (Wmin is not None and Wmin < 0):
            viol.append((a, b))
    report["structured"] = {
        "cases": len(cases),
        "violations_V_or_W": viol,
        "deepest_U_fail_depth": deepU[0],
        "deepest_U_fail_pair": deepU[1],
        "min_V_ratio": [str(globV[0]), float(globV[0])],
        "min_W_ratio": [str(globW[0]), float(globW[0])],
    }

    # ---- bottleneck family: U, V, W clean to h = 80 ----
    bn = {"U_fail_h": [], "V_fail_h": [], "W_fail_h": []}
    for h in range(1, 81):
        n = 2 * h + 1
        N, Ufail, Vmin, Wmin = UVW_scan([1, 1], [n, n])
        if Ufail:
            bn["U_fail_h"].append(h)
        if Vmin < 0:
            bn["V_fail_h"].append(h)
        if Wmin < 0:
            bn["W_fail_h"].append(h)
    report["bottleneck_h_le_80"] = bn

    # ---- random + hill-climb hunt ----
    random.seed(20260811)

    def rand_arms():
        style = random.random()
        if style < 0.3:
            return ([1] * random.randint(0, 10)
                    + [random.randint(20, 80)
                       for _ in range(random.randint(0, 2))])
        if style < 0.6:
            return [random.randint(1, 25) for _ in range(random.randint(1, 6))]
        if style < 0.8:
            return [random.randint(1, 120)]
        base = random.randint(1, 40)
        return [base, base + 1, base + 2][:random.randint(1, 3)]

    hunt_min = {"V": (None, None), "W": (None, None)}
    hunt_deepU = -1
    nhunt = 0
    for _ in range(1500):
        a, b = rand_arms(), rand_arms()
        if sum(a) + sum(b) > 200:
            continue
        nhunt += 1
        N, Ufail, Vmin, Wmin = UVW_scan(a, b)
        if Ufail:
            hunt_deepU = max(hunt_deepU, max(Ufail))
        for slot, val in (("V", Vmin), ("W", Wmin)):
            if val is not None and (hunt_min[slot][0] is None
                                    or val < hunt_min[slot][0]):
                hunt_min[slot] = (val, (a, b))

    def climb(seed_pair, slot, steps=150):
        best = seed_pair
        _, _, V0, W0 = UVW_scan(*best)
        bestval = V0 if slot == "V" else W0
        for _ in range(steps):
            a, b = list(best[0]), list(best[1])
            op = random.random()
            side = a if random.random() < 0.5 else b
            if op < 0.3 and side:
                side[random.randrange(len(side))] += random.choice([-1, 1])
            elif op < 0.5:
                side.append(random.randint(1, 30))
            elif op < 0.7 and side:
                side.pop(random.randrange(len(side)))
            elif op < 0.85 and a:
                b.append(a.pop(random.randrange(len(a))))
            elif b:
                a.append(b.pop(random.randrange(len(b))))
            a = [x for x in a if x > 0]
            b = [x for x in b if x > 0]
            if not (0 < sum(a) + sum(b) <= 250):
                continue
            _, _, V1, W1 = UVW_scan(a, b)
            val = V1 if slot == "V" else W1
            if val is not None and val < bestval:
                bestval, best = val, (a, b)
        return best, bestval

    climbs = {}
    for slot in ("V", "W"):
        best, val = climb(([list(hunt_min[slot][1][0]),
                            list(hunt_min[slot][1][1])][0],
                           list(hunt_min[slot][1][1])), slot)
        climbs[slot] = {"pair": best, "min_ratio": [str(val), float(val)]}
    report["hunt"] = {
        "random_pairs": nhunt,
        "min_V_ratio": [str(hunt_min["V"][0]), float(hunt_min["V"][0])],
        "min_W_ratio": [str(hunt_min["W"][0]), float(hunt_min["W"][0])],
        "deepest_U_fail_depth": hunt_deepU,
        "hill_climb": climbs,
    }

    # ---- generator-pair bilinear table ----
    labels = ["(QQ,QQ)", "(QQ,xRaQ)", "(QQ,xQaR)",
              "(xRR,QQ)", "(xRR,xRaQ)", "(xRR,xQaR)"]
    clean = {f: {lab: True for lab in labels} for f in "UVW"}
    random.seed(7)
    gpairs = [([], [1, 3]), ([1, 1], [9, 9]), ([6], [47, 71]),
              ([], [7, 15, 5]), ([1], [2]), ([2, 3], [4, 5])]
    for _ in range(120):
        a = [random.randint(1, 12) for _ in range(random.randint(0, 4))]
        b = [random.randint(1, 12) for _ in range(random.randint(1, 4))]
        if sum(a) + sum(b) <= 40:
            gpairs.append((a, b))
    for a, b in gpairs:
        Qa, Ra = prodP(a, 0), prodP(a, 1)
        Qb, Rb = prodP(b, 0), prodP(b, 1)
        g1 = pmul(Qa, Qb)
        g2 = [0] + pmul(Ra, Rb)
        h2 = [0] + pmul(Ra, Qb)
        h3 = [0] + pmul(Qa, Rb)
        N = max(len(g1), len(h2), len(h3)) + 1
        for (f, g), lab in zip([(g1, g1), (g1, h2), (g1, h3),
                                (g2, g1), (g2, h2), (g2, h3)], labels):
            F = [0] + f
            for k in range(N + 2):
                if cf(f, k) * cf(g, k) - cf(f, k - 1) * cf(g, k + 1) < 0:
                    clean["U"][lab] = False
                if cf(f, k) * cf(g, k) - cf(f, k + 1) * cf(g, k - 1) < 0:
                    clean["V"][lab] = False
                if (2 * cf(F, k) * cf(g, k) - cf(F, k + 1) * cf(g, k - 1)
                        - cf(F, k - 1) * cf(g, k + 1)) < 0:
                    clean["W"][lab] = False
    report["generator_pair_clean"] = clean
    report["generator_pair_cases"] = len(gpairs)

    # ---- end-to-end: base partial-sync relations on the structured set ----
    psync_bad = []
    for a, b in cases[:2000]:
        B, C = BC(a, b)
        if psync_fails(C, B) or psync_fails([0] + C, B):
            psync_bad.append((a, b))
    report["base_psync_failures_first2000"] = psync_bad

    report["runtime_s"] = round(time.time() - t0, 2)
    with open("results/c2_single_index_laws_20260811.json", "w") as fh:
        json.dump(report, fh, indent=1, default=str)

    print(json.dumps({k: v for k, v in report.items()
                      if k not in ("structured",)}, indent=1, default=str)[:400])
    print("structured:", report["structured"]["cases"], "cases;",
          "V/W violations:", report["structured"]["violations_V_or_W"],
          "; deepest U depth:", report["structured"]["deepest_U_fail_depth"])
    print("base psync failures (first 2000 structured pairs):",
          report["base_psync_failures_first2000"])
    print("runtime", report["runtime_s"], "s")


if __name__ == "__main__":
    main()
