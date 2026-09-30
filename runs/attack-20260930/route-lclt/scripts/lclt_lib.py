"""Exact tree DP helpers for the route-lclt kill tests (integer / Fraction only).

For a tree given by adjacency lists, root it, 2-colour it (L = even depth,
R = odd depth), and compute in one DP pass
    Z(x)   = sum_S x^|S|                (independence polynomial, i_k)
    ZL(x)  = sum_S |S cap L| x^|S|
Then, at a rational fugacity lam = a/b,
    mu     = lam Z'(lam)/Z(lam),  sigma^2 = lam d/dlam mu
    muL    = ZL(lam)/Z(lam),      muR = mu - muL
Smoothing identity (Fang Lemma 6.1 structure): conditional on I cap R, the
available L-vertices are iid Bernoulli(p), p = lam/(1+lam), so
    E[B_L] = muL / p,   E Var(X | I cap R) = p(1-p) E[B_L] = muL/(1+lam),
and the variance fraction carried by the binomial kernel is
    rho_L = muL / ((1+lam) sigma^2)       (similarly rho_R).
"""
from fractions import Fraction


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return out


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def shift(a):
    return [0] + a


def tree_polys(adj, root=0):
    """Return (Z, ZL, colour) for tree adjacency list adj (0-indexed)."""
    n = len(adj)
    parent = [-1] * n
    depth = [0] * n
    order = [root]
    seen = [False] * n
    seen[root] = True
    i = 0
    while i < len(order):
        u = order[i]; i += 1
        for w in adj[u]:
            if not seen[w]:
                seen[w] = True
                parent[w] = u
                depth[w] = depth[u] + 1
                order.append(w)
    # per vertex: (E0, E1, I0, I1): excluded/included polys, and their L-weighted versions
    E0 = [None] * n; E1 = [None] * n; I0 = [None] * n; I1 = [None] * n
    for u in reversed(order):
        # product over children of (E+I) with weighted derivative-style bookkeeping
        P, Pw = [1], [0]
        Q, Qw = [1], [0]
        for c in adj[u]:
            if c == parent[u]:
                continue
            T0 = padd(E0[c], I0[c]); T1 = padd(E1[c], I1[c])
            P, Pw = pmul(P, T0), padd(pmul(Pw, T0), pmul(P, T1))
            Q, Qw = pmul(Q, E0[c]), padd(pmul(Qw, E0[c]), pmul(Q, E1[c]))
        E0[u], E1[u] = P, Pw
        I0[u] = shift(Q)
        if depth[u] % 2 == 0:
            I1[u] = shift(padd(Qw, Q))
        else:
            I1[u] = shift(Qw)
    Z = padd(E0[root], I0[root])
    ZL = padd(E1[root], I1[root])
    while len(Z) > 1 and Z[-1] == 0:
        Z.pop()
    ZL = ZL[:len(Z)]
    colour = [d % 2 for d in depth]
    return Z, ZL, colour


def eval_scaled(coeffs, a, b, m=0):
    """Return sum_j j^m c_j a^j b^(D-j) (integer), D = len(coeffs)-1.
    Equals b^D * sum_j j^m c_j (a/b)^j."""
    D = len(coeffs) - 1
    tot = 0
    apow = 1
    bpows = [1] * (D + 1)
    for j in range(1, D + 1):
        bpows[j] = bpows[j - 1] * b
    for j, c in enumerate(coeffs):
        if c:
            tot += (j ** m) * c * apow * bpows[D - j]
        apow *= a
    return tot


def moments_at(Z, ZL, lam):
    """Exact (mu, sigma2, muL) at rational lam, as Fractions."""
    lam = Fraction(lam)
    a, b = lam.numerator, lam.denominator
    S0 = eval_scaled(Z, a, b, 0)
    S1 = eval_scaled(Z, a, b, 1)
    S2 = eval_scaled(Z, a, b, 2)
    SL = eval_scaled(ZL + [0] * (len(Z) - len(ZL)), a, b, 0)
    mu = Fraction(S1, S0)
    sig2 = Fraction(S2, S0) - mu * mu
    muL = Fraction(SL, S0)
    return mu, sig2, muL


def parent_line_to_adj(par):
    n = len(par)
    adj = [[] for _ in range(n)]
    for v in range(1, n + 1):
        p = par[v - 1]
        if p != 0:
            adj[v - 1].append(p - 1)
            adj[p - 1].append(v - 1)
    return adj


def window(n, alpha):
    lo = max(1, -(-n // 4))
    q = -(-(2 * alpha - 1) // 3)
    hi = min(q, alpha - 1)
    return lo, hi, q


def analyse(adj, ks=None):
    """Per-k exact diagnostics on the central window.  Returns list of dicts."""
    n = len(adj)
    Z, ZL, _ = tree_polys(adj)
    alpha = len(Z) - 1
    lo, hi, q = window(n, alpha)
    rows = []
    rng = range(lo, hi + 1) if ks is None else [k for k in ks if lo <= k <= hi]
    for k in rng:
        im, i0, ip = Z[k - 1], Z[k], Z[k + 1]
        delta = 1 - Fraction(im * ip, i0 * i0)   # exact, tilt invariant
        res = {"k": k, "delta_pos": delta > 0}
        g = []; rho = []; lams = []
        for lam in (Fraction(im, i0), Fraction(i0, ip)):
            mu, s2, muL = moments_at(Z, ZL, lam)
            muR = mu - muL
            rL = muL / ((1 + lam) * s2)
            rR = muR / ((1 + lam) * s2)
            g.append(s2 * delta)
            rho.append(max(rL, rR))
            lams.append(lam)
        res.update(Gamma_lo=min(g), Gamma_hi=max(g), rho=min(rho), lam_hi=lams[1], lam_lo=lams[0])
        rows.append(res)
    return Z, alpha, (lo, hi, q), rows
