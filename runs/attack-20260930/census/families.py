#!/usr/bin/env python3
"""Central-window margins for structured tree families beyond the exhaustive range.

Exact integer independence polynomials (treepoly.ipoly_from_adj) and exact
Fraction margins delta_k = 1 - i_{k-1} i_{k+1}/i_k^2 over
W5 = [ceil(n/5), q], q = ceil((2 alpha - 1)/3) (and W4 = [ceil(n/4), q]).
The tilted quantities (lambda_k, V_k) are FLOAT DIAGNOSTICS (bisection on
log lambda in treepoly.tilt).

Family definitions (sourced):
  Galvin T_{m,t}: root v, m children, each child carries t pendant P2 arms
      (Bautista-Ramos arXiv:2511.00334 Def. 1(iii); Levit-Kadrawi
      arXiv:2603.17114 Def. 3.1(3)).
  Kadrawi-Levit T_{3,m,n'}: root v0 with children v1,v2,v3 carrying 3, m, n'
      P2 arms (Levit-Kadrawi arXiv:2603.17114 Def. 3.1(1)).
  TG_{m,t}: m disjoint copies of T_{3,t}, roots joined to a new root v0, plus
      one pendant leaf at v0 (Bautista-Ramos arXiv:2511.00334 Def. 1(iv)).
Other families are standard (star, path, spiders, brooms, double stars,
caterpillars, complete b-ary trees, edge-subdivided complete binary trees).

Usage: python3 families.py [--maxn 301] [--workers 5]
"""
import argparse
import json
import os
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from treepoly import ipoly_from_adj, margins, q_of, tilt, tilted_Q, window  # noqa: E402


class B:
    """tiny tree builder"""
    def __init__(self):
        self.adj = []

    def v(self):
        self.adj.append([])
        return len(self.adj) - 1

    def e(self, a, b):
        self.adj[a].append(b)
        self.adj[b].append(a)

    def path_from(self, a, length):
        """attach a path of `length` new vertices hanging from a"""
        prev = a
        for _ in range(length):
            w = self.v()
            self.e(prev, w)
            prev = w
        return prev


def star(n):
    t = B(); c = t.v()
    for _ in range(n - 1):
        t.e(c, t.v())
    return t.adj


def path(n):
    t = B(); c = t.v(); t.path_from(c, n - 1)
    return t.adj


def spider(legs):
    t = B(); c = t.v()
    for L in legs:
        t.path_from(c, L)
    return t.adj


def broom(p, s):
    """path on p vertices, s extra leaves at one end"""
    t = B(); a = t.v(); end = t.path_from(a, p - 1)
    for _ in range(s):
        t.e(end, t.v())
    return t.adj


def double_star(a, b):
    t = B(); x = t.v(); y = t.v(); t.e(x, y)
    for _ in range(a):
        t.e(x, t.v())
    for _ in range(b):
        t.e(y, t.v())
    return t.adj


def caterpillar(spine, leaves_each):
    t = B(); prev = None
    for i in range(spine):
        s = t.v()
        if prev is not None:
            t.e(prev, s)
        prev = s
        for _ in range(leaves_each):
            t.e(s, t.v())
    return t.adj


def complete_bary(b, h):
    t = B(); root = t.v(); level = [root]
    for _ in range(h):
        nxt = []
        for u in level:
            for _ in range(b):
                w = t.v(); t.e(u, w); nxt.append(w)
        level = nxt
    return t.adj


def subdivided_binary(h, s):
    """complete binary tree of depth h with every edge subdivided s times"""
    t = B(); root = t.v(); level = [root]
    for _ in range(h):
        nxt = []
        for u in level:
            for _ in range(2):
                end = t.path_from(u, s + 1)
                nxt.append(end)
        level = nxt
    return t.adj


def star_of_stars(m, s):
    """root joined to the centres of m copies of K_{1,s} (the exhaustive
    Newton-ratio minimiser shape at n=20,25)"""
    t = B(); r = t.v()
    for _ in range(m):
        c = t.v(); t.e(r, c)
        for _ in range(s):
            t.e(c, t.v())
    return t.adj


def bistar_mid(a, b):
    """two stars K_{1,a}, K_{1,b} whose centres are joined through one middle
    vertex (the exhaustive runner-up shape for min delta)"""
    t = B(); x = t.v(); m = t.v(); y = t.v(); t.e(x, m); t.e(m, y)
    for _ in range(a):
        t.e(x, t.v())
    for _ in range(b):
        t.e(y, t.v())
    return t.adj


def bundle_child(t, parent, arms):
    c = t.v(); t.e(parent, c)
    for _ in range(arms):
        t.path_from(c, 2)
    return c


def galvin_T(m, tt):
    t = B(); r = t.v()
    for _ in range(m):
        bundle_child(t, r, tt)
    return t.adj


def kl_T3(m, n2):
    t = B(); r = t.v()
    for arms in (3, m, n2):
        bundle_child(t, r, arms)
    return t.adj


def TG(m, tt):
    t = B(); r0 = t.v()
    for _ in range(m):
        v = t.v(); t.e(r0, v)
        for _ in range(3):
            bundle_child(t, v, tt)
    t.e(r0, t.v())
    return t.adj


def family_specs(maxn):
    S = []
    add = lambda fam, params, fn: S.append((fam, params, fn))
    for n in list(range(11, maxn + 1, 10)):
        add("star", {"n": n}, ("star", (n,)))
        add("path", {"n": n}, ("path", (n,)))
    for m in range(5, (maxn - 1) // 2 + 1, 10):
        add("spider_2^m", {"m": m}, ("spider", ([2] * m,)))
    for m in range(4, (maxn - 1) // 3 + 1, 8):
        add("spider_3^m", {"m": m}, ("spider", ([3] * m,)))
    for a in range(3, (maxn - 1) // 3 + 1, 8):
        add("spider_1^a2^a", {"a": a}, ("spider", ([1] * a + [2] * a,)))
    for m in range(2, (maxn - 1) // 5 + 1, 6):
        add("spider_5^m(star_of_P5)", {"m": m}, ("spider", ([5] * m,)))
    for m in range(2, (maxn - 1) // 10 + 1, 3):
        add("spider_10^m(star_of_P10)", {"m": m}, ("spider", ([10] * m,)))
    for L in range(5, (maxn - 1) // 3 + 1, 12):
        add("spider_3legs_L(long)", {"L": L}, ("spider", ([L] * 3,)))
    for p in range(6, maxn // 2 + 1, 12):
        add("broom_p=s", {"p": p, "s": p}, ("broom", (p, p)))
    for p in range(6, maxn // 4 + 1, 6):
        add("broom_s=3p", {"p": p, "s": 3 * p}, ("broom", (p, 3 * p)))
    for a in range(4, (maxn - 2) // 2 + 1, 12):
        add("double_star_a,a", {"a": a}, ("double_star", (a, a)))
    for a in range(3, (maxn - 2) // 4 + 1, 6):
        add("double_star_a,3a", {"a": a}, ("double_star", (a, 3 * a)))
    for s in range(5, maxn // 2 + 1, 12):
        add("caterpillar_1leaf(comb)", {"spine": s}, ("caterpillar", (s, 1)))
    for s in range(4, maxn // 3 + 1, 8):
        add("caterpillar_2leaves", {"spine": s}, ("caterpillar", (s, 2)))
    for s in range(3, maxn // 5 + 1, 5):
        add("caterpillar_4leaves", {"spine": s}, ("caterpillar", (s, 4)))
    for h in range(2, 10):
        if 2 ** (h + 1) - 1 <= maxn:
            add("complete_binary", {"h": h}, ("complete_bary", (2, h)))
    for h in range(2, 7):
        if (3 ** (h + 1) - 1) // 2 <= maxn:
            add("complete_ternary", {"h": h}, ("complete_bary", (3, h)))
    for h in range(2, 9):
        if (2 ** (h + 1) - 1) + (2 ** (h + 1) - 2) <= maxn:
            add("binary_subdiv1", {"h": h}, ("subdivided_binary", (h, 1)))
        if (2 ** (h + 1) - 1) + 2 * (2 ** (h + 1) - 2) <= maxn:
            add("binary_subdiv2", {"h": h}, ("subdivided_binary", (h, 2)))
    for s in (2, 3, 4, 6):
        for m in range(2, 200):
            if 1 + m * (1 + s) > maxn:
                break
            if m <= 8 or m % 4 == 0:
                add(f"star_of_stars_K1,{s}", {"m": m, "s": s}, ("star_of_stars", (m, s)))
    for a in range(4, maxn // 2, 12):
        add("bistar_mid_a,a", {"a": a}, ("bistar_mid", (a, a)))
    for a in range(4, maxn // 5, 6):
        add("bistar_mid_4a,a", {"a": a}, ("bistar_mid", (4 * a, a)))
    for tt in (2, 3, 4, 6):
        for m in range(2, 200):
            if 1 + m * (1 + 2 * tt) > maxn:
                break
            if m <= 6 or m % 3 == 0:
                add(f"galvin_T_m,{tt}", {"m": m, "t": tt}, ("galvin_T", (m, tt)))
    for tt in range(2, 200, 4):
        if 1 + 3 * (1 + 2 * tt) > maxn:
            break
        add("galvin_T_3,t", {"m": 3, "t": tt}, ("galvin_T", (3, tt)))
    for mm in range(3, 200, 4):
        if 1 + 7 + 2 * (1 + 2 * mm) > maxn:
            break
        add("kadrawi_levit_T_3,m,m", {"m": mm}, ("kl_T3", (mm, mm)))
    for mm in range(3, 200, 6):
        if 1 + 7 + (1 + 2 * mm) + (1 + 2 * (mm + 1)) > maxn:
            break
        add("kadrawi_levit_T_3,m,m+1", {"m": mm}, ("kl_T3", (mm, mm + 1)))
    for tt in (3, 5):
        for m in range(1, 50):
            if 2 + m * (1 + 3 * (1 + 2 * tt)) > maxn:
                break
            add(f"TG_m,{tt}", {"m": m, "t": tt}, ("TG", (m, tt)))
    return S


BUILDERS = dict(star=star, path=path, spider=spider, broom=broom, double_star=double_star,
                caterpillar=caterpillar, complete_bary=complete_bary,
                subdivided_binary=subdivided_binary, galvin_T=galvin_T,
                star_of_stars=star_of_stars, bistar_mid=bistar_mid, kl_T3=kl_T3, TG=TG)


def analyse(spec):
    fam, params, (bname, args) = spec
    adj = BUILDERS[bname](*args)
    n = len(adj)
    p = ipoly_from_adj(adj)
    alpha = len(p) - 1
    q = q_of(alpha)
    lo5, _ = window(n, alpha, "n5")
    lo4, _ = window(n, alpha, "n4")
    ms = margins(p, n, "n5")
    k5 = min(ms, key=lambda k: ms[k])
    d5 = ms[k5]
    ms4 = {k: v for k, v in ms.items() if k >= lo4}
    k4 = min(ms4, key=lambda k: ms4[k])
    newton = {k: d * Fraction((k + 1) * (alpha - k + 1), alpha + 1) for k, d in ms.items()}
    kn = min(newton, key=lambda k: newton[k])
    # FLOAT diagnostics
    vd = []
    for k, d in ms.items():
        lam, V = tilt(p, k, iters=90)
        vd.append((k, lam, V, V * float(d), tilted_Q(p, k, lam, V)))
    kq = min(vd, key=lambda r: r[4])
    vbins = {}
    for r in vd:
        b = 0 if r[2] < 2 else 2 if r[2] < 4 else 4 if r[2] < 8 else 8 if r[2] < 16 else 16
        vbins[b] = min(vbins.get(b, 9e9), r[3])
    kvd = min(vd, key=lambda r: r[3])
    at_argmin = next(r for r in vd if r[0] == k5)
    at_q = next(r for r in vd if r[0] == q)
    breaks = [k for k in range(1, alpha) if p[k - 1] * p[k + 1] > p[k] * p[k]]
    uni = True
    rising = True
    for i in range(1, len(p)):
        if rising:
            if p[i] < p[i - 1]:
                rising = False
        elif p[i] > p[i - 1]:
            uni = False
    return dict(
        family=fam, params=params, n=n, alpha=alpha, q=q, W5=[lo5, q], W4=[lo4, q],
        min_delta_W5=f"{d5.numerator}/{d5.denominator}" if n <= 60 else None,
        min_delta_W5_float=float(d5), argmin_k_W5=k5, argmin_rel_pos=(k5 - lo5) / max(1, q - lo5),
        n_min_delta_W5=float(n * d5),
        min_delta_W4_float=float(ms4[k4]), n_min_delta_W4=float(n * ms4[k4]), argmin_k_W4=k4,
        delta_q_float=float(ms[q]), n_delta_q=float(n * ms[q]),
        newton_ratio_min=float(newton[kn]), newton_argmin_k=kn, newton_q_minus_k=q - kn,
        V_at_argmin=at_argmin[2], lambda_at_argmin=at_argmin[1], Vdelta_at_argmin=at_argmin[3],
        V_at_q=at_q[2], lambda_at_q=at_q[1], Vdelta_at_q=at_q[3],
        Vdelta_min=kvd[3], Vdelta_min_k=kvd[0], Vdelta_min_V=kvd[2], Vdelta_min_lambda=kvd[1],
        Vdelta_min_by_Vbin={str(b): v for b, v in sorted(vbins.items())},
        Q_min=kq[4], Q_min_k=kq[0], Q_min_V=kq[2], Q_min_lambda=kq[1], Q_at_argmin=at_argmin[4], Q_at_q=at_q[4],
        central_violation=any(d < 0 for d in ms.values()),
        lc_breaks_minus_q=[k - q for k in breaks], unimodal=uni,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxn", type=int, default=301)
    ap.add_argument("--workers", type=int, default=5)
    a = ap.parse_args()
    specs = family_specs(a.maxn)
    t0 = time.time()
    with Pool(a.workers) as pool:
        rows = pool.map(analyse, specs, chunksize=1)
    out = os.path.join(HERE, "data", "families.json")
    json.dump(dict(maxn=a.maxn, elapsed_s=round(time.time() - t0, 1), rows=rows), open(out, "w"), indent=1)
    fams = {}
    for r in rows:
        fams.setdefault(r["family"], []).append(r)
    for fam, rs in fams.items():
        rs.sort(key=lambda r: r["n"])
        worst = min(rs, key=lambda r: r["n_min_delta_W5"])
        last = rs[-1]
        print(f"{fam:28s} members={len(rs):3d} maxn={last['n']:4d} "
              f"n*minD(last)={last['n_min_delta_W5']:.4f} worst n*minD={worst['n_min_delta_W5']:.4f}@n={worst['n']} "
              f"newton_min={min(r['newton_ratio_min'] for r in rs):.4f} "
              f"Vd_min={min(r['Vdelta_min'] for r in rs):.4f} Vd@argmin(last)={last['Vdelta_at_argmin']:.4f} "
              f"n*d_q(last)={last['n_delta_q']:.3f} Qmin={min(r['Q_min'] for r in rs):.4f} Q@argmin(last)={last['Q_at_argmin']:.4f} "
              f"viol={sum(r['central_violation'] for r in rs)} "
              f"nonuni={sum(not r['unimodal'] for r in rs)}", flush=True)
    print(f"elapsed {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
