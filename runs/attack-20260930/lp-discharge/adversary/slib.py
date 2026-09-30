"""Constant-sigma discharging S(sigma): exact per-tree feasible interval and the
S(2/3) verdict.

Rule: L_u(k) = sigma D_u(k) + (1 - sigma) sum_{v~u} D_v(k)/deg v,
D_v(k) = k i_{k-1} j^v_k - (k+1) i_k j^v_{k-1},  j^v = I(T - N[v]),
window W(T) = [ceil(n/4), min(q, alpha-1)], q = ceil((2 alpha - 1)/3).

Per row (k, u), with M = lcm(degrees), AM_u = sum_{v~u} D_v * M/deg v (integer):
  M L_u = a sigma + b,  a = M D_u - AM_u,  b = AM_u.
  a > 0: sigma <= -b/a (upper end);  a < 0: sigma >= b/(-a) (lower end);
  a = 0: infeasible for every sigma iff b > 0.
Interval ends are kept as exact integer pairs (num, den>0) and compared by
cross-multiplication; floats appear only in labelled ranking fields.
S(2/3) violation at (k,u): 3 M L_u = 2a + 3b > 0 (exact integer test).
I and J come from the validated rerooting DP in ../../r1-adversarial/r1lib.py
(flint fmpz_poly, exact; asserts the double-count identity)."""
import os, sys
from math import lcm
import sys as _s; _s.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'r1-adversarial'))
from r1lib import all_seqs, adj_from_edges, B  # noqa: E402


def co(s, k):
    return s[k] if 0 <= k < len(s) else 0


def win(n, alpha):
    lo = -((-n) // 4)
    q = -((-(2 * alpha - 1)) // 3)
    return lo, min(q, alpha - 1)


def lt(p, q):
    """exact p < q for pairs (num, den>0)"""
    return p[0] * q[1] < q[0] * p[1]


def analyse(adj, keep_rows=False):
    n = len(adj)
    I, J = all_seqs(adj)
    alpha = len(I) - 1
    lo, hi = win(n, alpha)
    deg = [len(x) for x in adj]
    M = 1
    for d in set(deg):
        M = lcm(M, d)
    W = [M // d for d in deg]
    low = None   # (num, den, k, u)
    up = None
    empty_row = None
    viol = []
    best23 = (-1e300, None)
    lcfail = []
    rows = [] if keep_rows else None
    for k in range(lo, hi + 1):
        c1, c2 = k * co(I, k - 1), (k + 1) * co(I, k)
        D = [c1 * co(Jv, k) - c2 * co(Jv, k - 1) for Jv in J]
        sD = sum(D)
        assert sD == k * (k + 1) * (co(I, k - 1) * co(I, k + 1) - co(I, k) ** 2)
        if sD > 0:
            lcfail.append(k)
        WD = [W[v] * D[v] for v in range(n)]
        nrm = 3 * M * k * co(I, k - 1) * co(I, k)
        for u in range(n):
            b = 0
            for v in adj[u]:
                b += WD[v]
            a = M * D[u] - b
            if keep_rows:
                rows.append((k, u, a, b))
            if a > 0:
                cand = (-b, a)
                if up is None or lt(cand, up[:2]):
                    up = (cand[0], cand[1], k, u)
            elif a < 0:
                cand = (b, -a)
                if low is None or lt(low[:2], cand):
                    low = (cand[0], cand[1], k, u)
            elif b > 0 and empty_row is None:
                empty_row = (k, u)
            s23 = 2 * a + 3 * b
            if s23 > 0:
                viol.append((u, k, deg[u]))
            f = s23 / nrm            # float, ranking only
            if f > best23[0]:
                best23 = (f, (u, k, deg[u]))
    out = dict(n=n, alpha=alpha, window=[lo, hi], lcfail=lcfail, viol23=viol[:20], n_viol23=len(viol),
               maxL23=best23[0], maxL23_at=best23[1], empty_row=empty_row)
    if low is not None:
        out['low'] = [str(low[0]), str(low[1])]
        out['low_f'] = low[0] / low[1]
        out['low_at'] = (low[3], low[2], deg[low[3]])
    else:
        out['low_f'] = -1e300
    if up is not None:
        out['up'] = [str(up[0]), str(up[1])]
        out['up_f'] = up[0] / up[1]
        out['up_at'] = (up[3], up[2], deg[up[3]])
    else:
        out['up_f'] = 1e300
    if keep_rows:
        out['rows'] = rows
    return out


def edges_of(adj):
    return [(u, v) for u in range(len(adj)) for v in adj[u] if u < v]


# ---------------- constructors (vertex 0 is the root) ----------------
def add_hubstar(b, parent, m, s):
    """attach a hub to `parent`; hub carries m supports each with s leaves. Returns hub."""
    h = b.add(parent)
    for _ in range(m):
        x = b.add(h)
        for _ in range(s):
            b.add(x)
    return h


def hubstar_rooted(m, s):
    b = B()
    for _ in range(m):
        x = b.add(0)
        for _ in range(s):
            b.add(x)
    return b


def msh(hubs, s=2, centre_leaves=0, centre_path=0, sub=1, svar=None):
    """centre -> (path of `sub` edges) -> hub H(m, s) for m in hubs.
    svar: optional list of s per hub."""
    b = B()
    for i, m in enumerate(hubs):
        p = b.path(0, sub - 1) if sub > 1 else 0
        add_hubstar(b, p, m, svar[i] if svar else s)
    for _ in range(centre_leaves):
        b.add(0)
    if centre_path:
        b.path(0, centre_path)
    return b.adj()
