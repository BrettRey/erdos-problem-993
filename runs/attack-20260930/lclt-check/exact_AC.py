#!/usr/bin/env python3
"""Independent exact (Fraction / integer) check of route-lclt claims A, B, C.

Written from scratch; does NOT import route-lclt's lclt_lib.

For a tree T with bipartition (L, R) and J subset R (every such J is independent):
  B_J = #{v in L : N(v) cap J = empty}.
(A0) combinatorial, lambda-free:  i_j(T) = sum_{J subset R} C(B_J, j - |J|).
(A)  at fugacity lam, p = lam/(1+lam): P(X=j) = sum_J w_J b_{B_J,p}(j-|J|),
     w_J = lam^|J| (1+lam)^{B_J} / Z(lam).
(A') fixed kernel B0: on E0 = {B_J >= B0}, P(X=j, E0) = sum_m P(M'=m, E0) b_{B0,p}(j-m),
     M' = |J| + Bin(B_J - B0, p).
(B)  Abel: |E h(k-M) - E h(k-G)| <= TV(h) sup_m |F_M(m) - F_G(m)| for lattice laws M, G.
(C)  D_k := 2p_k - p_{k-1} - p_{k+1}
        >= P(E0)[E h(k-G) - TV(h) dK(M'|E0, G)] - 2 P(E0^c)          (route-lclt form)
        >= ... - (a_{k-1} + a_{k+1}),  a_j = P(X=j, E0^c)              (sharp form, this check)
     G = any rational lattice law; here a rounded discretised normal with the
     conditional mean and variance of M' given E0 (rational after rounding).

Two data sources for the joint law of (|J|, B_J):
  brute force over all J subset R  (small trees), and
  an exact integer 2-D tree DP  N[s][B] = #{J subset R: |J| = s, B_J = B}  (any size),
cross-checked against each other and against brute-force / standard-DP i_k.
"""
import sys, json, math, subprocess
from fractions import Fraction as Fr
from math import comb

# ------------------------------------------------------------------ trees
def par_to_adj(par):
    n = len(par); adj = [[] for _ in range(n)]
    for v in range(1, n + 1):
        p = par[v - 1]
        if p:
            adj[v - 1].append(p - 1); adj[p - 1].append(v - 1)
    return adj

def bfs(adj, root=0):
    n = len(adj); parent = [-1] * n; depth = [0] * n; order = [root]; seen = [False] * n; seen[root] = True
    for u in order:
        for w in adj[u]:
            if not seen[w]:
                seen[w] = True; parent[w] = u; depth[w] = depth[u] + 1; order.append(w)
    return parent, depth, order

def indep_brute(adj):
    n = len(adj); nb = [0] * n
    for u in range(n):
        for w in adj[u]: nb[u] |= 1 << w
    cnt = [0] * (n + 1)
    for S in range(1 << n):
        ok = True; T = S
        while T:
            u = (T & -T).bit_length() - 1; T &= T - 1
            if nb[u] & S: ok = False; break
        if ok: cnt[bin(S).count("1")] += 1
    while cnt and cnt[-1] == 0: cnt.pop()
    return cnt

def indep_dp(adj):
    """standard two-state DP (independent of the J-mixture)."""
    parent, depth, order = bfs(adj)
    def mul(a, b):
        out = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b): out[i + j] += x * y
        return out
    def add(a, b):
        m = max(len(a), len(b)); return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(m)]
    ex = [None] * len(adj); inc = [None] * len(adj)
    for u in reversed(order):
        e = [1]; i = [0, 1]
        for c in adj[u]:
            if c == parent[u]: continue
            e = mul(e, add(ex[c], inc[c])); i = mul(i, ex[c])
        ex[u], inc[u] = e, i
    r = add(ex[order[0]], inc[order[0]])
    while r[-1] == 0: r.pop()
    return r

# ------------------------------------------------------------------ (|J|, B_J) joint counts
def joint_brute(adj, Lside):
    parent, depth, order = bfs(adj)
    L = [v for v in range(len(adj)) if depth[v] % 2 == Lside]
    R = [v for v in range(len(adj)) if depth[v] % 2 != Lside]
    idx = {v: i for i, v in enumerate(R)}
    masks = []
    for v in L:
        m = 0
        for w in adj[v]: m |= 1 << idx[w]
        masks.append(m)
    N = {}
    for J in range(1 << len(R)):
        B = sum(1 for m in masks if m & J == 0)
        key = (bin(J).count("1"), B); N[key] = N.get(key, 0) + 1
    return N

def joint_dp(adj, Lside):
    """exact integer 2-D DP; polys as dict {(s,B): count}."""
    parent, depth, order = bfs(adj)
    isL = [depth[v] % 2 == Lside for v in range(len(adj))]
    def mul(a, b):
        out = {}
        for (s1, b1), x in a.items():
            for (s2, b2), y in b.items():
                k = (s1 + s2, b1 + b2); out[k] = out.get(k, 0) + x * y
        return out
    def add(a, b, sign=1):
        out = dict(a)
        for k, y in b.items(): out[k] = out.get(k, 0) + sign * y
        return {k: v for k, v in out.items() if v}
    one = {(0, 0): 1}
    A = [None] * len(adj); Bv = [None] * len(adj)   # R vertex: A = not in J, Bv = in J ; L vertex: A = parent not in J, Bv = parent in J
    for u in reversed(order):
        ch = [c for c in adj[u] if c != parent[u]]
        if isL[u]:
            allp = one; outp = one
            for c in ch:
                allp = mul(allp, add(A[c], Bv[c])); outp = mul(outp, A[c])
            F1 = allp
            F0 = add(add(allp, outp, -1), {(s, b + 1): v for (s, b), v in outp.items()})
            A[u], Bv[u] = F0, F1
        else:
            p0 = one; p1 = one
            for c in ch:
                p0 = mul(p0, A[c]); p1 = mul(p1, Bv[c])
            A[u] = p0; Bv[u] = {(s + 1, b): v for (s, b), v in p1.items()}
    r = order[0]
    return A[r] if isL[r] else add(A[r], Bv[r])

# ------------------------------------------------------------------ exact probability pieces
def binpmf(B, p, j):
    if j < 0 or j > B: return Fr(0)
    return comb(B, j) * p ** j * (1 - p) ** (B - j)

def window(n, alpha):
    lo = -(-n // 4); q = -(-(2 * alpha - 1) // 3)
    return lo, min(q, alpha - 1)

def rational_normal(mean, var, lo, hi, den=10 ** 15):
    """rounded discretised normal on integers lo..hi (a rational lattice law; any lattice law is valid for B)."""
    s = math.sqrt(max(float(var), 1e-12)); m = float(mean)
    Phi = lambda x: 0.5 * math.erfc(-x / math.sqrt(2))
    w = []
    for j in range(lo, hi + 1):
        a = Phi((j - 0.5 - m) / s) if j > lo else 0.0
        b = Phi((j + 0.5 - m) / s) if j < hi else 1.0
        w.append(Fr(round((b - a) * den), den))
    tot = sum(w)
    return {j: w[j - lo] / tot for j in range(lo, hi + 1) if w[j - lo]}

def cdf_sup(P, Q):
    keys = sorted(set(P) | set(Q)); F = G = Fr(0); d = Fr(0)
    for m in keys:
        F += P.get(m, 0); G += Q.get(m, 0); d = max(d, abs(F - G))
    return d

def check_tree_k(N, Z, k, B0s, stats, record=None):
    """N: joint counts {(s,B):count}; Z: i_j list. All checks exact."""
    lam = Fr(Z[k - 1], Z[k]); p = lam / (1 + lam)
    Zl = sum(c * lam ** j for j, c in enumerate(Z))
    alpha = len(Z) - 1
    pk = lambda j: Fr(Z[j]) * lam ** j / Zl if 0 <= j <= alpha else Fr(0)
    w = {key: c * lam ** key[0] * (1 + lam) ** key[1] / Zl for key, c in N.items()}
    assert sum(w.values()) == 1
    # (A)
    for j in range(alpha + 2):
        assert pk(j) == sum(wt * binpmf(B, p, j - s) for (s, B), wt in w.items())
    D = 2 * pk(k) - pk(k - 1) - pk(k + 1)
    rhs = sum(wt * (2 * binpmf(B, p, k - s) - binpmf(B, p, k - 1 - s) - binpmf(B, p, k + 1 - s)) for (s, B), wt in w.items())
    assert D == rhs
    stats["A_cases"] += 1
    EB = sum(wt * B for (s, B), wt in w.items())
    for B0 in B0s:
        PE0 = sum(wt for (s, B), wt in w.items() if B >= B0)
        if PE0 == 0: continue
        lawM = {}
        for (s, B), wt in w.items():
            if B < B0: continue
            for t in range(B - B0 + 1):
                lawM[s + t] = lawM.get(s + t, 0) + wt * binpmf(B - B0, p, t)
        # (A') fixed-kernel split, all j
        for j in range(alpha + 2):
            lhs = sum(wt * binpmf(B, p, j - s) for (s, B), wt in w.items() if B >= B0)
            assert lhs == sum(pm * binpmf(B0, p, j - m) for m, pm in lawM.items())
        a = lambda j: sum(wt * binpmf(B, p, j - s) for (s, B), wt in w.items() if B < B0)
        h = lambda j: 2 * binpmf(B0, p, j) - binpmf(B0, p, j - 1) - binpmf(B0, p, j + 1)
        TV = sum(abs(h(j + 1) - h(j)) for j in range(-3, B0 + 3))
        condM = {m: v / PE0 for m, v in lawM.items()}
        mbar = sum(m * v for m, v in condM.items()); var = sum((m - mbar) ** 2 * v for m, v in condM.items())
        G = rational_normal(mbar, var, min(condM) - 4, max(condM) + 4)
        dK = cdf_sup(condM, G)
        EhM = sum(v * h(k - m) for m, v in condM.items()); EhG = sum(v * h(k - m) for m, v in G.items())
        # (B) Abel bound, exact
        assert abs(EhM - EhG) <= TV * dK
        # exact decomposition of D
        tailpart = 2 * a(k) - a(k - 1) - a(k + 1)
        assert D == PE0 * EhM + tailpart
        # tail bounds
        assert tailpart >= -(a(k - 1) + a(k + 1)) >= -(1 - PE0)
        cert_old = PE0 * (EhG - TV * dK) - 2 * (1 - PE0)
        cert_sharp = PE0 * (EhG - TV * dK) - (a(k - 1) + a(k + 1))
        assert D >= cert_sharp >= cert_old
        stats["C_cases"] += 1
        stats["cert_old_pos"] += cert_old > 0; stats["cert_sharp_pos"] += cert_sharp > 0
        if B0 == max(1, math.floor(Fr(1, 2) * EB)):
            stats["mid_cases"] += 1
            stats["mid_old_pos"] += cert_old > 0; stats["mid_sharp_pos"] += cert_sharp > 0
        if record is not None:
            record.append(dict(k=k, B0=B0, EB=float(EB), PE0=float(PE0), D=float(D), EhG=float(EhG), TVdK=float(TV * dK),
                               dK=float(dK), cert_old=float(cert_old), cert_sharp=float(cert_sharp),
                               TVv32=float(TV) * (float(B0 * p * (1 - p))) ** 1.5))
    return D

def run_census(nmax):
    stats = dict(trees=0, A0_cases=0, A_cases=0, C_cases=0, cert_old_pos=0, cert_sharp_pos=0, mid_cases=0, mid_old_pos=0, mid_sharp_pos=0, dp_vs_brute=0)
    for n in range(3, nmax + 1):
        out = subprocess.run(["gentreeg", "-p", "-q", str(n)], capture_output=True, text=True).stdout
        for line in out.splitlines():
            par = list(map(int, line.split()))
            if len(par) != n: continue
            adj = par_to_adj(par)
            Z = indep_brute(adj)
            assert Z == indep_dp(adj)
            alpha = len(Z) - 1
            lo, hi = window(n, alpha)
            stats["trees"] += 1
            for side in (0, 1):
                N = joint_brute(adj, side)
                assert N == joint_dp(adj, side); stats["dp_vs_brute"] += 1
                # (A0) lambda-free combinatorial identity
                for j in range(alpha + 2):
                    assert sum(c * comb(B, j - s) for (s, B), c in N.items() if 0 <= j - s <= B) == (Z[j] if j <= alpha else 0)
                stats["A0_cases"] += 1
                Bmax = max(B for (s, B) in N)
                for k in range(lo, hi + 1):
                    EBf = None
                    B0s = sorted(set([1, Bmax] + [max(1, b) for b in range(1, Bmax + 1)]))
                    check_tree_k(N, Z, k, B0s, stats)
    return stats

if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    st = run_census(nmax)
    print(json.dumps(st))
