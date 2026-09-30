"""Exact helpers for the central-margin census (2026-09-30).

Independence polynomial of a tree by the standard two-state DP in Python
integers; exact central-window margins as Fractions; an independent
(float, DIAGNOSTIC ONLY) solver for the hard-core fugacity with a given mean
and the variance there.
"""
import math
from fractions import Fraction


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def padd(a, b):
    if len(a) < len(b):
        a, b = b, a
    out = list(a)
    for i, y in enumerate(b):
        out[i] += y
    return out


def ipoly_from_adj(adj, root=0):
    """Exact independence polynomial of a tree (adjacency lists)."""
    n = len(adj)
    parent = [-1] * n
    order = []
    stack = [root]
    seen = [False] * n
    seen[root] = True
    while stack:
        v = stack.pop()
        order.append(v)
        for w in adj[v]:
            if not seen[w]:
                seen[w] = True
                parent[w] = v
                stack.append(w)
    assert len(order) == n, "not connected"
    F = [None] * n   # v excluded
    G = [None] * n   # v included
    for v in reversed(order):
        f = [1]
        g = [0, 1]
        for w in adj[v]:
            if w == parent[v]:
                continue
            f = pmul(f, padd(F[w], G[w]))
            g = pmul(g, F[w])
        F[v], G[v] = f, g
    p = padd(F[root], G[root])
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def adj_from_parent(par):
    """par: list, 1-indexed parent array as printed by gentreeg -p (par[0] of root = 0)."""
    n = len(par)
    adj = [[] for _ in range(n)]
    for i, p in enumerate(par):
        if p:
            adj[i].append(p - 1)
            adj[p - 1].append(i)
    return adj


def adj_from_edges(n, edges):
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    return adj


def q_of(alpha):
    return (2 * alpha + 1) // 3          # ceil((2 alpha - 1)/3)


def window(n, alpha, lower="n5"):
    lo = -(-n // 5) if lower == "n5" else -(-n // 4)
    lo = max(lo, 1)
    return lo, q_of(alpha)


def delta(p, k):
    return Fraction(p[k] * p[k] - p[k - 1] * p[k + 1], p[k] * p[k])


def margins(p, n, lower="n5"):
    alpha = len(p) - 1
    lo, q = window(n, alpha, lower)
    return {k: delta(p, k) for k in range(lo, q + 1)}


def min_margin(p, n, lower="n5"):
    ms = margins(p, n, lower)
    k = min(ms, key=lambda kk: ms[kk])
    return ms[k], k


# ---------------- FLOAT DIAGNOSTIC: tilted variance ----------------

def _moments(logc, t):
    lw = [lc + j * t for j, lc in enumerate(logc)]
    m = max(lw)
    w = [math.exp(x - m) for x in lw]
    s0 = sum(w)
    s1 = sum(j * x for j, x in enumerate(w))
    s2 = sum(j * j * x for j, x in enumerate(w))
    mu = s1 / s0
    return mu, s2 / s0 - mu * mu


def tilt(p, k, iters=200):
    """Return (lambda, V) with hard-core mean k, by bisection on log lambda.
    Independent of the C implementation (which uses safeguarded Newton)."""
    logc = [math.log(c) for c in p]
    lo, hi = -80.0, 80.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        mu, _ = _moments(logc, mid)
        if mu < k:
            lo = mid
        else:
            hi = mid
    t = 0.5 * (lo + hi)
    mu, V = _moments(logc, t)
    return math.exp(t), V


def tilted_Q(p, k, lam, V):
    """FLOAT DIAGNOSTIC: Q_k = V^{3/2} (2 p_k - p_{k-1} - p_{k+1}) under the law
    p_j proportional to i_j lam^j (arXiv:2609.20961, Eq. (6.4) quantity)."""
    t = math.log(lam)
    lw = [math.log(c) + j * t for j, c in enumerate(p)]
    m = max(lw)
    s0 = sum(math.exp(x - m) for x in lw)
    pr = lambda j: math.exp(lw[j] - m) / s0
    return V ** 1.5 * (2 * pr(k) - pr(k - 1) - pr(k + 1))


def vd_profile(p, n, lower="n5"):
    """list of (k, delta_float, lambda, V, V*delta) over the window."""
    out = []
    for k, d in margins(p, n, lower).items():
        lam, V = tilt(p, k)
        out.append((k, float(d), lam, V, V * float(d)))
    return out
