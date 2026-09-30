#!/usr/bin/env python3
"""Route pbmix kill-test: trees as mixtures of shifted Poisson-binomial laws.

Exact arithmetic throughout (Python ints / Fractions).  Floats appear only in
the printed/JSON diagnostic columns, which are labelled *_float.

Setup.  Fix a tree T and a set B with T - B a linear forest.  For each
independent set sigma of T[B], the vertices outside B u N(sigma) induce a
linear forest, so the configurations extending sigma have generating
polynomial a^sigma(x) = x^{|sigma|} Q_sigma(x), Q_sigma a product of path
polynomials (real-rooted, hence Poisson-binomial after normalisation).
I(T;x) = f(x) = sum_sigma a^sigma(x).

At an index k, with posterior pi(sigma) = a^sigma_k / f_k (over sigma with
a^sigma_k > 0), u_sigma = a_{k+1}/a_k, v_sigma = a_{k-1}/a_k:
    f_{k+1}/f_k = E_pi[u] + beta_+ ,  f_{k-1}/f_k = E_pi[v] + beta_- ,
where beta_+- collect components with a_k = 0 (support boundary).
    delta(f) = 1 - (E u + beta_+)(E v + beta_-).
Each component has delta^sigma = 1 - u_sigma v_sigma.
Diagnostic mixture loss:  rho = 1 - delta(f) / dbar,  dbar = E_pi[delta^sigma].
LC at k  <=>  rho < 1.  A route that only knows delta^sigma >= d_sigma
(certified lower bound) certifies LC iff
    cert = 1 - (E u + beta_+)(E[(1-d_sigma)/u] + beta_-) > 0     (bulk comps)
(edge components, with u = 0 or v = 0, are kept exact).

Certified component bounds d_sigma:
  ULC (Newton):  delta_j >= (m+1)/((j+1)(m-j+1)), m = deg Q_sigma  [exact for binomials]
  PB  (thm:main + tilt invariance):  delta_j >= 1/(4 V(lambda)) for any
      lambda in [a_{j-2}/a_{j-1}, a_{j-1}/a_j) with V(lambda) >= 1
      (evaluated at the two closure endpoints; conservative).
"""
from __future__ import annotations

import itertools
import json
import math
import random
import subprocess
import sys
import time
from fractions import Fraction as Fr
from math import comb

# ---------------------------------------------------------------- polynomials

def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return r


def padd(a, b):
    if len(a) < len(b):
        a, b = b, a
    r = list(a)
    for i, y in enumerate(b):
        r[i] += y
    return r


PATH_CACHE = {0: [1], 1: [1, 1]}


def path_poly(l):
    """Independence polynomial of the path on l vertices."""
    if l in PATH_CACHE:
        return PATH_CACHE[l]
    a, b = [1], [1, 1]
    for i in range(2, l + 1):
        c = padd(b, [0] + a)
        a, b = b, c
        PATH_CACHE[i] = c
    return PATH_CACHE[l]


def coef(a, j):
    return a[j] if 0 <= j < len(a) else 0


def indpoly_tree(adj, n):
    """Standard tree DP (sanity)."""
    root = 0
    order, parent = [], [-1] * n
    seen = [False] * n
    stack = [root]
    seen[root] = True
    while stack:
        v = stack.pop()
        order.append(v)
        for w in adj[v]:
            if not seen[w]:
                seen[w] = True
                parent[w] = v
                stack.append(w)
    d0, d1 = [None] * n, [None] * n
    for v in reversed(order):
        p0, p1 = [1], [0, 1]
        for w in adj[v]:
            if w == parent[v]:
                continue
            p0 = pmul(p0, padd(d0[w], d1[w]))
            p1 = pmul(p1, d0[w])
        d0[v], d1[v] = p0, p1
    r = padd(d0[root], d1[root])
    while len(r) > 1 and r[-1] == 0:
        r.pop()
    return r


def alpha_of(poly):
    return len(poly) - 1


# ---------------------------------------------------------------- PB variance

def pb_variance(Q, lam):
    """Variance of the size law with pgf Q(lam x)/Q(lam): exact Fraction."""
    Z = sum(c * lam ** i for i, c in enumerate(Q))
    M1 = sum(i * c * lam ** i for i, c in enumerate(Q))
    M2 = sum(i * i * c * lam ** i for i, c in enumerate(Q))
    return Fr(M2) / Z - (Fr(M1) / Z) ** 2


def d_ulc(m, j):
    return Fr(m + 1, (j + 1) * (m - j + 1))


def d_pb(Q, j):
    """Certified lower bound for delta_j(Q) from thm:main via tilting, or None."""
    m = len(Q) - 1
    if not (1 <= j <= m - 1):
        return None
    cands = []
    # closed left endpoint lambda = Q_{j-2}/Q_{j-1} (j>=2) puts D=j exactly
    if j >= 2:
        lam = Fr(Q[j - 2], Q[j - 1])
        V = pb_variance(Q, lam)
        if V >= 1:
            cands.append(V)
    # open right endpoint lambda = Q_{j-1}/Q_j: valid by continuity if V > 1
    lam = Fr(Q[j - 1], Q[j])
    V = pb_variance(Q, lam)
    if V > 1:
        cands.append(V)
    if not cands:
        return None
    return 1 / (4 * min(cands))


# ---------------------------------------------------------------- mixture stats

def mixture_stats(comps, k, want_cert=True, binomial=False):
    """comps: list of (mult, shift, Q) ; component poly = mult * x^shift * Q.
    If binomial=True, Q is given as the integer m meaning (1+x)^m.
    Returns dict of exact Fractions."""
    def qc(Q, j):
        if binomial:
            return comb(Q, j) if 0 <= j <= Q else 0
        return coef(Q, j)

    fkm, fk, fkp = 0, 0, 0
    for mult, s, Q in comps:
        fkm += mult * qc(Q, k - 1 - s)
        fk += mult * qc(Q, k - s)
        fkp += mult * qc(Q, k + 1 - s)
    if fk == 0:
        return None
    delta_f = 1 - Fr(fkm * fkp, fk * fk)
    # posterior-weighted component curvature
    dbar_num = Fr(0)
    Eu_num, Ev_num, Evcert_num, Evcert2_num = 0, Fr(0), Fr(0), Fr(0)
    bplus, bminus = 0, 0
    n_bulk_uncert = 0
    for mult, s, Q in comps:
        j = k - s
        ak = qc(Q, j)
        am = qc(Q, j - 1)
        ap = qc(Q, j + 1)
        if ak == 0:
            bplus += mult * ap
            bminus += mult * am
            continue
        w = mult * ak
        dsig = 1 - Fr(am * ap, ak * ak)
        dbar_num += w * dsig
        Eu_num += mult * ap
        Ev_num += Fr(mult * am)
        if not want_cert:
            continue
        if am == 0 or ap == 0:
            # edge component: keep exact
            Evcert_num += Fr(mult * am)
            Evcert2_num += Fr(mult * am)
            continue
        m = Q if binomial else len(Q) - 1
        cands = [d_ulc(m, j)]
        Qpoly = [comb(Q, i) for i in range(Q + 1)] if binomial else Q
        dp = d_pb(Qpoly, j)
        # cert_PB: PB bound only (ULC only as fallback when PB inapplicable)
        d1 = dp if dp is not None else d_ulc(m, j)
        # cert_best: max(PB, ULC)
        d2 = max(cands + ([dp] if dp is not None else []))
        if dp is None:
            n_bulk_uncert += 1
        # v = (1-delta)/u ; certified v_cert = (1-d)/u = (1-d)*ak/ap
        Evcert_num += mult * (1 - d1) * Fr(ak * ak, ap)
        Evcert2_num += mult * (1 - d2) * Fr(ak * ak, ap)
    dbar = dbar_num / fk
    out = {
        "delta_f": delta_f,
        "dbar": dbar,
        "rho": (1 - delta_f / dbar) if dbar > 0 else None,
    }
    if want_cert:
        Eu = Fr(Eu_num + bplus, fk)
        out["cert_pb"] = 1 - Eu * (Evcert_num / fk + Fr(bminus, fk))
        out["cert_best"] = 1 - Eu * (Evcert2_num / fk + Fr(bminus, fk))
        out["n_bulk_no_pb"] = n_bulk_uncert
    return out


def between_fraction(comps, lam, binomial=False):
    """Var(E[X|sigma]) / Var(X) at fugacity lam (exact)."""
    Zs, means, vars_ = [], [], []
    for mult, s, Q in comps:
        if binomial:
            m = Q
            Z = mult * (1 + lam) ** m
            p = lam / (1 + lam)
            mu = s + m * p
            var = m * p * (1 - p)
        else:
            Zq = sum(c * lam ** i for i, c in enumerate(Q))
            M1 = sum(i * c * lam ** i for i, c in enumerate(Q))
            Z = mult * lam ** s * Zq
            mu = s + Fr(M1) / Zq
            var = pb_variance(Q, lam)
        if binomial:
            Z = Z * lam ** s
        Zs.append(Z)
        means.append(mu)
        vars_.append(var)
    Ztot = sum(Zs)
    Ew = sum(z * v for z, v in zip(Zs, vars_)) / Ztot
    mu_tot = sum(z * mu for z, mu in zip(Zs, means)) / Ztot
    Bv = sum(z * (mu - mu_tot) ** 2 for z, mu in zip(Zs, means)) / Ztot
    return Bv / (Bv + Ew), Ew, Bv


# ------------------------------------------------ decomposition over B (census)

def components_branch_set(adj, n, B):
    """Enumerate independent sets sigma of T[B]; return list of (1, |sigma|, Q_sigma)."""
    Bl = sorted(B)
    Bset = set(Bl)
    comps = {}
    for r in range(len(Bl) + 1):
        for sig in itertools.combinations(Bl, r):
            ok = True
            ss = set(sig)
            for v in sig:
                if any(w in ss for w in adj[v]):
                    ok = False
                    break
            if not ok:
                continue
            blocked = set(Bset)
            for v in sig:
                blocked.update(adj[v])
            rest = [v for v in range(n) if v not in blocked]
            rs = set(rest)
            # path components of T[rest]
            seen = set()
            lens = []
            for v in rest:
                if v in seen:
                    continue
                stack, cnt = [v], 0
                seen.add(v)
                while stack:
                    x = stack.pop()
                    cnt += 1
                    for y in adj[x]:
                        if y in rs and y not in seen:
                            seen.add(y)
                            stack.append(y)
                lens.append(cnt)
            key = (len(sig), tuple(sorted(lens)))
            comps[key] = comps.get(key, 0) + 1
    out = []
    for (s, lens), mult in comps.items():
        Q = [1]
        for l in lens:
            Q = pmul(Q, path_poly(l))
        out.append((mult, s, Q))
    return out


def degree_ok_linear_forest(adj, n, B):
    for v in range(n):
        if v in B:
            continue
        if sum(1 for w in adj[v] if w not in B) > 2:
            return False
    return True


# ---------------------------------- bivariate DP when T - B is an independent set

def bmul(A, Bd):
    R = {}
    for (s1, m1), c1 in A.items():
        for (s2, m2), c2 in Bd.items():
            key = (s1 + s2, m1 + m2)
            R[key] = R.get(key, 0) + c1 * c2
    return R


def badd(A, Bd):
    R = dict(A)
    for k_, c in Bd.items():
        R[k_] = R.get(k_, 0) + c
    return R


def bsub(A, Bd):
    R = dict(A)
    for k_, c in Bd.items():
        R[k_] = R.get(k_, 0) - c
    return {k_: c for k_, c in R.items() if c != 0}


def bshift(A, ds, dm):
    return {(s + ds, m + dm): c for (s, m), c in A.items()}


def binomial_mixture(adj, n, B, root=0):
    """Return {(s,m): count}: I(T;x) = sum count * x^s (1+x)^m.
    Requires V(T)-B independent."""
    order, parent = [], [-1] * n
    seen = [False] * n
    stack = [root]
    seen[root] = True
    while stack:
        v = stack.pop()
        order.append(v)
        for w in adj[v]:
            if not seen[w]:
                seen[w] = True
                parent[w] = v
                stack.append(w)
    G0, G1, HB, HF = {}, {}, {}, {}
    ONE = {(0, 0): 1}
    for v in reversed(order):
        ch = [w for w in adj[v] if w != parent[v]]
        if v in B:
            g0, g1 = ONE, {(1, 0): 1}
            for c in ch:
                if c in B:
                    g0 = bmul(g0, badd(G0[c], G1[c]))
                    g1 = bmul(g1, G0[c])
                else:
                    g0 = bmul(g0, badd(HB[c], bshift(HF[c], 0, 1)))
                    g1 = bmul(g1, badd(HB[c], HF[c]))
            G0[v], G1[v] = g0, g1
        else:
            allp, free = ONE, ONE
            for c in ch:
                assert c in B
                allp = bmul(allp, badd(G0[c], G1[c]))
                free = bmul(free, G0[c])
            HB[v] = bsub(allp, free)
            HF[v] = free
    if root in B:
        return badd(G0[root], G1[root])
    return badd(HB[root], bshift(HF[root], 0, 1))


def poly_from_binomial_mixture(mix):
    deg = max(s + m for (s, m) in mix)
    f = [0] * (deg + 1)
    for (s, m), c in mix.items():
        for j in range(m + 1):
            f[s + j] += c * comb(m, j)
    while len(f) > 1 and f[-1] == 0:
        f.pop()
    return f


# ---------------------------------------------------------------- tree families

def complete_binary(depth):
    n = 2 ** (depth + 1) - 1
    adj = [[] for _ in range(n)]
    for v in range(1, n):
        p = (v - 1) // 2
        adj[v].append(p)
        adj[p].append(v)
    return adj, n


def random_tree(n, rng):
    if n <= 2:
        adj = [[] for _ in range(n)]
        if n == 2:
            adj[0].append(1)
            adj[1].append(0)
        return adj
    prufer = [rng.randrange(n) for _ in range(n - 2)]
    deg = [1] * n
    for x in prufer:
        deg[x] += 1
    adj = [[] for _ in range(n)]
    import heapq
    leaves = [i for i in range(n) if deg[i] == 1]
    heapq.heapify(leaves)
    for x in prufer:
        leaf = heapq.heappop(leaves)
        adj[leaf].append(x)
        adj[x].append(leaf)
        deg[x] -= 1
        if deg[x] == 1:
            heapq.heappush(leaves, x)
    u = heapq.heappop(leaves)
    w = heapq.heappop(leaves)
    adj[u].append(w)
    adj[w].append(u)
    return adj


def spider(legs):
    n = 1 + sum(legs)
    adj = [[] for _ in range(n)]
    nxt = 1
    for l in legs:
        prev = 0
        for _ in range(l):
            adj[prev].append(nxt)
            adj[nxt].append(prev)
            prev = nxt
            nxt += 1
    return adj, n


def windows(n, alpha):
    q = -(-(2 * alpha - 1) // 3)
    lo = -(-n // 4)
    lo5 = -(-n // 5)
    hi17 = (17 * alpha) // 25
    return lo, q, lo5, hi17


def fl(x):
    return None if x is None else float(x)


# ---------------------------------------------------------------- experiments

def analyse_binomial_family(name, adj, n, B, kset=None, cert=True):
    t0 = time.time()
    mix = binomial_mixture(adj, n, B)
    f = poly_from_binomial_mixture(mix)
    g = indpoly_tree(adj, n)
    assert f == g, f"{name}: mixture does not reproduce I(T)"
    alpha = alpha_of(f)
    lo, q, lo5, hi17 = windows(n, alpha)
    comps = [(c, s, m) for (s, m), c in mix.items()]
    ks = kset if kset is not None else list(range(lo, q + 1))
    rows = []
    for k in ks:
        st = mixture_stats(comps, k, want_cert=cert, binomial=True)
        lam = Fr(f[k], f[k + 1])  # tilt with tilted f_k = f_{k+1}
        bf, Ew, Bv = between_fraction(comps, lam, binomial=True)
        rows.append({
            "k": k,
            "LC_exact": st["delta_f"] > 0,
            "delta_f_float": fl(st["delta_f"]),
            "dbar_float": fl(st["dbar"]),
            "rho_float": fl(st["rho"]),
            "cert_pb_pass": (st["cert_pb"] > 0) if cert else None,
            "cert_pb_float": fl(st.get("cert_pb")),
            "between_frac_float": fl(bf),
            "within_var_float": fl(Ew),
            "between_var_float": fl(Bv),
        })
    return {
        "family": name, "n": n, "alpha": alpha, "B_size": len(B),
        "groups": len(mix), "window": [lo, q], "rows": rows,
        "secs": round(time.time() - t0, 2),
    }


def census(nmax, nmin=8):
    """All trees n<=nmax (geng), B = vertices of degree >= 3 (components = paths)."""
    summary = []
    worst = []
    for n in range(nmin, nmax + 1):
        out = subprocess.run(["/opt/homebrew/bin/geng", "-c", "-q", str(n), f"{n-1}:{n-1}"],
                             capture_output=True, text=True).stdout.split()
        cnt = dict(trees=0, pairs=0, lc_fail=0, cert_pb_fail=0, cert_best_fail=0,
                   rho_max=0.0, rho_ge_quarter=0)
        for g6 in out:
            adj = decode_g6(g6)
            nn = len(adj)
            B = {v for v in range(nn) if len(adj[v]) >= 3}
            assert degree_ok_linear_forest(adj, nn, B)
            comps = components_branch_set(adj, nn, B)
            f = [0] * (nn + 1)
            for mult, s, Q in comps:
                for i, c in enumerate(Q):
                    f[s + i] += mult * c
            while len(f) > 1 and f[-1] == 0:
                f.pop()
            assert f == indpoly_tree(adj, nn)
            alpha = alpha_of(f)
            lo, q, _, _ = windows(nn, alpha)
            cnt["trees"] += 1
            for k in range(lo, q + 1):
                st = mixture_stats(comps, k, want_cert=True, binomial=False)
                cnt["pairs"] += 1
                if not st["delta_f"] > 0:
                    cnt["lc_fail"] += 1
                if not st["cert_pb"] > 0:
                    cnt["cert_pb_fail"] += 1
                if not st["cert_best"] > 0:
                    cnt["cert_best_fail"] += 1
                r = float(st["rho"]) if st["rho"] is not None else 0.0
                if r >= 0.25:
                    cnt["rho_ge_quarter"] += 1
                if r > cnt["rho_max"]:
                    cnt["rho_max"] = r
                    cnt["rho_max_at"] = [g6, k, len(B)]
        summary.append({"n": n, **cnt})
        print(json.dumps(summary[-1]), flush=True)
    return summary


def decode_g6(s):
    data = [ord(c) - 63 for c in s]
    n = data[0]
    bits = []
    for x in data[1:]:
        for i in range(5, -1, -1):
            bits.append((x >> i) & 1)
    adj = [[] for _ in range(n)]
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if bits[idx]:
                adj[i].append(j)
                adj[j].append(i)
            idx += 1
    return adj


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    results = {}
    if which in ("all", "cbt"):
        res = []
        for depth in range(2, 8):
            adj, n = complete_binary(depth)
            B3 = {v for v in range(n) if len(adj[v]) == 3}
            Bint = {v for v in range(n) if len(adj[v]) >= 2}
            for bname, B in (("deg3", B3), ("nonleaf", Bint)):
                r = analyse_binomial_family(f"CBT d={depth} B={bname}", adj, n, B,
                                            cert=(depth <= 6))
                res.append(r)
                rh = [x["rho_float"] for x in r["rows"]]
                bfr = [x["between_frac_float"] for x in r["rows"]]
                print(r["family"], "n", n, "|B|", len(B), "groups", r["groups"],
                      "window", r["window"],
                      "allLC", all(x["LC_exact"] for x in r["rows"]),
                      "rho[min,max]", (round(min(rh), 4), round(max(rh), 4)),
                      "betweenfrac[min,max]", (round(min(bfr), 4), round(max(bfr), 4)),
                      "cert_pb_fail", sum(1 for x in r["rows"] if x["cert_pb_pass"] is False),
                      "secs", r["secs"], flush=True)
        results["cbt"] = res
        json.dump(res, open("results_cbt.json", "w"), indent=1)
    if which in ("all", "rand"):
        rng = random.Random(20260930)
        res = []
        for n in (40, 80, 120, 160):
            for rep in range(3):
                adj = random_tree(n, rng)
                Bint = {v for v in range(n) if len(adj[v]) >= 2}
                r = analyse_binomial_family(f"rand n={n} rep={rep} B=nonleaf", adj, n, Bint,
                                            cert=(n <= 80))
                res.append(r)
                rh = [x["rho_float"] for x in r["rows"]]
                bfr = [x["between_frac_float"] for x in r["rows"]]
                print(r["family"], "|B|", len(Bint), "groups", r["groups"], "window", r["window"],
                      "allLC", all(x["LC_exact"] for x in r["rows"]),
                      "rho[min,max]", (round(min(rh), 4), round(max(rh), 4)),
                      "betweenfrac[min,max]", (round(min(bfr), 4), round(max(bfr), 4)),
                      "cert_pb_fail", sum(1 for x in r["rows"] if x["cert_pb_pass"] is False),
                      "secs", r["secs"], flush=True)
        results["rand"] = res
        json.dump(res, open("results_rand.json", "w"), indent=1)
    if which in ("all", "spider"):
        res = []
        for legs in ([2] * 10, [2] * 30, [3] * 15, [5] * 12, [8] * 10, [20] * 6, [40] * 5):
            adj, n = spider(legs)
            B = {0}
            comps = components_branch_set(adj, n, B)
            f = [0] * (n + 1)
            for mult, s, Q in comps:
                for i, c in enumerate(Q):
                    f[s + i] += mult * c
            while len(f) > 1 and f[-1] == 0:
                f.pop()
            assert f == indpoly_tree(adj, n)
            alpha = alpha_of(f)
            lo, q, _, _ = windows(n, alpha)
            rows = []
            for k in range(lo, q + 1):
                st = mixture_stats(comps, k, want_cert=True)
                lam = Fr(f[k], f[k + 1])
                bf, Ew, Bv = between_fraction(comps, lam)
                rows.append({"k": k, "LC_exact": st["delta_f"] > 0,
                             "rho_float": fl(st["rho"]), "cert_pb_pass": st["cert_pb"] > 0,
                             "cert_best_pass": st["cert_best"] > 0,
                             "between_frac_float": fl(bf)})
            rh = [x["rho_float"] for x in rows]
            print("spider", legs[0], "x", len(legs), "n", n, "window", [lo, q],
                  "rho[min,max]", (round(min(rh), 4), round(max(rh), 4)),
                  "cert_pb_fail", sum(1 for x in rows if not x["cert_pb_pass"]),
                  "cert_best_fail", sum(1 for x in rows if not x["cert_best_pass"]),
                  "pairs", len(rows), flush=True)
            res.append({"legs": legs, "n": n, "window": [lo, q], "rows": rows})
        results["spider"] = res
        json.dump(res, open("results_spider.json", "w"), indent=1, default=str)
    if which in ("all", "census"):
        nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 14
        s = census(nmax)
        json.dump(s, open("results_census.json", "w"), indent=1)


if __name__ == "__main__":
    main()
