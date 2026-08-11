#!/usr/bin/env python3
"""Independent audit of partial_sync_obstruction_20260811.

Re-derives everything claim-critical with freshly written code:
  - brute-force independence polynomials on small instances (bitmask
    enumeration, no DP), including an independently constructed double broom;
  - own psync-violation, Turan, and LC-break implementations;
  - full sweeps at the note's claimed parameter ranges;
  - per-finding checks for F1-F5 and the mechanism identity.
Exact integer arithmetic throughout.
"""
import importlib.util
import json
import sys

import os
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location(
    "vps", REPO + "/verify_partial_sync_obstruction_20260811.py")
vps = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vps)

fails = []


def check(label, ok):
    print(("PASS " if ok else "FAIL ") + label, flush=True)
    if not ok:
        fails.append(label)


# ---------- independent primitives ----------

def cf(P, i):
    return P[i] if 0 <= i < len(P) else 0


def my_psync_viol(A, B):
    D = max(len(A), len(B))
    return sorted(
        (m, n)
        for m in range(D + 1)
        for n in range(m + 1)
        if cf(A, m) * cf(B, n) + cf(A, n) * cf(B, m)
        < cf(A, m + 1) * cf(B, n - 1) + cf(A, n - 1) * cf(B, m + 1))


def my_turan(P, j):
    return cf(P, j) ** 2 - cf(P, j - 1) * cf(P, j + 1)


def my_lc_breaks(P):
    return [j for j in range(1, len(P) - 1)
            if P[j] * P[j] < P[j - 1] * P[j + 1]]


def brute_ipoly(n, adj):
    nb = [0] * n
    for v in range(n):
        for u in adj[v]:
            nb[v] |= 1 << u
    cnt = [0] * (n + 1)
    for S in range(1 << n):
        T = S
        ok = True
        while T:
            v = (T & -T).bit_length() - 1
            if nb[v] & S & ~(1 << v):
                ok = False
                break
            T &= T - 1
        if ok:
            cnt[bin(S).count("1")] += 1
    while cnt and cnt[-1] == 0:
        cnt.pop()
    return cnt


# ---------- 1. brute-force DP checks ----------

# my own double broom T_{5,5,2}: hubs a,b; 2 internal path vertices; 5+5 leaves
edges = []
a, i1, i2, b = 0, 1, 2, 3
edges += [(a, i1), (i1, i2), (i2, b)]
v = 4
for _ in range(5):
    edges.append((a, v)); v += 1
for _ in range(5):
    edges.append((b, v)); v += 1
nbb = [[] for _ in range(v)]
for x, y in edges:
    nbb[x].append(y); nbb[y].append(x)
_, F_broom, _, _ = vps.double_broom_instance(5, 5, 2)
check("brute T_{5,5,2} (independent construction) == their F",
      brute_ipoly(v, nbb) == F_broom)

T = vps.pattern_tree(lambda: vps.T1_ell(3), lambda: vps.path_tree(2), [2, 2])
check("brute KL k=2 tree (order %d) == their DP" % T.n,
      brute_ipoly(T.n, T.adj) == T.ipoly())

S = vps.spider2(4)
check("brute S_{2,4} == their DP", brute_ipoly(S.n, S.adj) == S.ipoly())
G = vps.T1_ell(2)
check("brute T_{1,2} == their DP", brute_ipoly(G.n, G.adj) == G.ipoly())

# ---------- 2. full sweeps with fresh implementations ----------

RANGES = [
    ("double_broom_ctrl", [("t", t) for t in range(2, 26)],
     lambda t: vps.double_broom_instance(5, 5, t)),
    ("KL_3kk", [("k", k) for k in range(2, 25)],
     lambda k: vps.kl_3kk_instance(k)),
    ("T12_S24_2", [("m", m) for m in range(2, 61)],
     lambda m: vps.m_family_instance(2, 4, 2, m)),
    ("T13_S24_2", [("m", m) for m in range(2, 31)],
     lambda m: vps.m_family_instance(3, 4, 2, m)),
    ("T17_S25_2", [("m", m) for m in range(2, 25)],
     lambda m: vps.m_family_instance(7, 5, 2, m)),
]

rows = {}
for fam, params, builder in RANGES:
    for pname, pval in params:
        order, F, A, B = builder(pval)   # asserts decomposition == DP inside
        alpha = len(F) - 1
        viol = my_psync_viol(A, B)
        breaks = my_lc_breaks(F)
        rows[(fam, pval)] = {
            "family": fam, pname: pval, "order": order, "alpha": alpha,
            "viol": viol,
            "count": len(viol),
            "width": max((m - n for m, n in viol), default=0),
            "depth": (alpha - min(n for _, n in viol)) if viol else 0,
            "head_diag": (alpha - 1, alpha - 1) in viol,
            "top_diag": (alpha, alpha) in viol,
            "breaks_refl": sorted(alpha - j for j in breaks),
            "unimodal": vps.is_unimodal(F),
        }
    print("swept " + fam, flush=True)

R = rows

# no instance violates at (alpha, alpha); all trees unimodal
check("(alpha,alpha) never violated (all %d instances)" % len(R),
      not any(r["top_diag"] for r in R.values()))
check("all instances unimodal", all(r["unimodal"] for r in R.values()))
check("double-broom control: zero violations at every t",
      all(r["count"] == 0 for r in R.values()
          if r["family"] == "double_broom_ctrl"))

# F1: KL violation set is exactly {(alpha-1, alpha-1)} for every k
check("F1: KL viol set == {(alpha-1,alpha-1)} for k=2..24",
      all(R[("KL_3kk", k)]["viol"]
          == [(R[("KL_3kk", k)]["alpha"] - 1, R[("KL_3kk", k)]["alpha"] - 1)]
          for k in range(2, 25)))

# F2: first sync failure vs first LC break, per family
def first(fam, rng, pred):
    return next((p for p in rng if pred(R[(fam, p)])), None)

f_sync_KL = first("KL_3kk", range(2, 25), lambda r: r["count"] > 0)
f_lc_KL = first("KL_3kk", range(2, 25), lambda r: r["breaks_refl"])
check("F2: KL sync fails from k=2, LC from k=4",
      f_sync_KL == 2 and f_lc_KL == 4)
f_sync_12 = first("T12_S24_2", range(2, 61), lambda r: r["count"] > 0)
f_lc_12 = first("T12_S24_2", range(2, 61), lambda r: r["breaks_refl"])
check("F2: T12 sync fails from m=3, LC from m=9",
      f_sync_12 == 3 and f_lc_12 == 9)
f_sync_13 = first("T13_S24_2", range(2, 31), lambda r: r["count"] > 0)
f_2brk_13 = first("T13_S24_2", range(2, 31), lambda r: len(r["breaks_refl"]) >= 2)
f_sync_17 = first("T17_S25_2", range(2, 25), lambda r: r["count"] > 0)
f_3brk_17 = first("T17_S25_2", range(2, 25), lambda r: len(r["breaks_refl"]) >= 3)
print(f"  T13: first sync fail m={f_sync_13}, first 2-break m={f_2brk_13}")
print(f"  T17: first sync fail m={f_sync_17}, first 3-break m={f_3brk_17}")
check("F2: T13/T17 sync fails by m=4 vs 2-/3-break thresholds 16 and 13",
      f_sync_13 <= 4 and f_sync_17 <= 4
      and f_2brk_13 == 16 and f_3brk_17 == 13)

# F3: head diagonal pair violated in every non-clean instance
check("F3: (alpha-1,alpha-1) violated in every non-clean instance",
      all(r["head_diag"] for r in R.values() if r["count"] > 0))
w60 = R[("T12_S24_2", 60)]
check("F3: T12 m=60 alpha=603, band width 26",
      w60["alpha"] == 603 and w60["width"] == 26)
print("  T12 widths m=14,30,40,50,60:",
      [R[("T12_S24_2", m)]["width"] for m in (14, 30, 40, 50, 60)])

# F4: depth/alpha ratios and marginal growth
ratios = {m: R[("T12_S24_2", m)]["depth"] / R[("T12_S24_2", m)]["alpha"]
          for m in (14, 30, 40, 50, 60)}
print("  T12 depth/alpha:", {m: round(v, 3) for m, v in ratios.items()})
check("F4: ratios match note (0.25, 0.42, 0.47, 0.51, 0.54)",
      [round(ratios[m], 2) for m in (14, 30, 40, 50, 60)]
      == [0.25, 0.42, 0.47, 0.51, 0.54])
d50, d60 = R[("T12_S24_2", 50)], R[("T12_S24_2", 60)]
marg = (d60["depth"] - d50["depth"]) / (d60["alpha"] - d50["alpha"])
print(f"  marginal depth ratio (m=50->60): {marg:.3f}")
check("F4: marginal ratio ~0.7", 0.6 <= marg <= 0.8)
# Levit-Mandrescu decreasing zone: i_k strictly decreasing for
# k >= ceil((2 alpha - 1)/3) (convention of paper/main.tex), i.e. the zone
# holds the reflected indices j <= alpha - ceil((2 alpha - 1)/3).
def lm_bound(alpha):
    return alpha - (-(-(2 * alpha - 1) // 3))

f_lm = first("T12_S24_2", range(2, 61),
             lambda r: r["depth"] > lm_bound(r["alpha"]))
r20 = R[("T12_S24_2", 20)]
print(f"  T12 first m with violations outside LM zone: {f_lm}; "
      f"m=20 depth {r20['depth']} vs boundary {lm_bound(r20['alpha'])}")
check("F4: band touches LM boundary at m=20, escapes from m=21 on",
      f_lm == 21 and r20["depth"] == lm_bound(r20["alpha"])
      and all(R[("T12_S24_2", m)]["depth"]
              > lm_bound(R[("T12_S24_2", m)]["alpha"])
              for m in range(21, 61)))

# F5: LC break sets; find non-contiguous-from-1 exceptions
exceptions = [(r["family"], p) for (f_, p), r in R.items()
              if r["breaks_refl"]
              and r["breaks_refl"] != list(range(1, len(r["breaks_refl"]) + 1))]
print("  F5 top-block exceptions:", exceptions)
for f_, p in exceptions:
    print("   ", f_, p, "breaks_refl:", R[(f_, p)]["breaks_refl"])

# ---------- 3. mechanism identity at T12 m=9 ----------

_, F9, A9, B9 = vps.m_family_instance(2, 4, 2, 9)
alpha9 = len(F9) - 1
ok = all(
    my_turan(F9, j) == my_turan(A9, j) + my_turan(B9, j)
    + 2 * cf(A9, j) * cf(B9, j) - cf(A9, j - 1) * cf(B9, j + 1)
    - cf(A9, j + 1) * cf(B9, j - 1)
    for j in range(alpha9 + 1))
check("mechanism: Turan(F)=Turan(A)+Turan(B)+S at all indices (T12 m=9)", ok)
j = alpha9 - 1
tA, tB = my_turan(A9, j), my_turan(B9, j)
Sj = (2 * cf(A9, j) * cf(B9, j) - cf(A9, j - 1) * cf(B9, j + 1)
      - cf(A9, j + 1) * cf(B9, j - 1))
print(f"  at break index {j}: Turan(A)={tA} Turan(B)={tB} S={Sj}")
check("mechanism: values 3945235 + 1089228 - 5272677 = -238214",
      sorted((tA, tB)) == [1089228, 3945235] and Sj == -5272677
      and tA + tB + Sj == -238214)
m = 9
check("Cor 4.12 closed form: Turan(F) at break == (456m+10)^2-5(47008m^2-41668m+6)",
      my_turan(F9, j)
      == (456 * m + 10) ** 2 - 5 * (47008 * m * m - 41668 * m + 6))

print()
if fails:
    print("AUDIT FAILED:", len(fails), "checks:", fails)
    sys.exit(1)
print("AUDIT PASSED: all checks")
