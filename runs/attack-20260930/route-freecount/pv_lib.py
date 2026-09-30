"""Exact-arithmetic helpers for the per-vertex (PV) free-count route.

All inequality checks use Python integers (cross-multiplied); floats appear
only in clearly labelled diagnostic fields.

Notation (tree T on n vertices, vertex v, level k):
  I      = independence sequence of T            (i_k)
  J_v    = independence sequence of T - N[v]     (j_k)
  E_v    = independence sequence of T - v        (e_k)
  alpha  = len(I) - 1,  q = ceil((2 alpha - 1)/3)
  window = [ceil(n/4), q]   (LC needed there; Fang Prop 8.2 covers below)

Exact identities used (asserted at run time):
  k i_k(T)       = sum_v i_{k-1}(T - N[v])
  (k+1) i_{k+1}  = sum_v i_k(T - N[v])
Hence LC at k,  i_{k-1} i_{k+1} <= i_k^2, is the v-sum of
  PV_v(k):  k * i_{k-1}(T) * j_k  <=  (k+1) * i_k(T) * j_{k-1}.
"""

from fractions import Fraction


def parse_parent_array(tokens):
    """1-indexed parent array (0 = root), as written by gentreeg -p and the
    LC census files. Returns adjacency lists on 0..n-1."""
    par = [int(t) for t in tokens]
    n = len(par)
    adj = [[] for _ in range(n)]
    for i, p in enumerate(par):
        if p != 0:
            a, b = i, p - 1
            adj[a].append(b)
            adj[b].append(a)
    return adj


def polymul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def polyadd(a, b):
    if len(a) < len(b):
        a, b = b, a
    out = list(a)
    for i, y in enumerate(b):
        out[i] += y
    return out


def indep_seq(adj, alive):
    """Independence sequence of the forest induced on `alive` (a set)."""
    seen = set()
    total = [1]
    for root in alive:
        if root in seen:
            continue
        # iterative DFS order
        order = []
        parent = {root: None}
        stack = [root]
        seen.add(root)
        while stack:
            u = stack.pop()
            order.append(u)
            for w in adj[u]:
                if w in alive and w not in seen:
                    seen.add(w)
                    parent[w] = u
                    stack.append(w)
        dp0 = {}
        dp1 = {}
        for u in reversed(order):
            p0 = [1]
            p1 = [0, 1]
            for w in adj[u]:
                if w in alive and parent.get(w) == u:
                    p0 = polymul(p0, polyadd(dp0[w], dp1[w]))
                    p1 = polymul(p1, dp0[w])
            dp0[u] = p0
            dp1[u] = p1
        total = polymul(total, polyadd(dp0[root], dp1[root]))
    while len(total) > 1 and total[-1] == 0:
        total.pop()
    return total


def co(seq, k):
    return seq[k] if 0 <= k < len(seq) else 0


def ceil_div(a, b):
    return -((-a) // b)


def window(n, alpha):
    lo = ceil_div(n, 4)
    q = ceil_div(2 * alpha - 1, 3)
    return lo, q


def tree_data(adj):
    n = len(adj)
    V = set(range(n))
    I = indep_seq(adj, V)
    J = []
    E = []
    for v in range(n):
        J.append(indep_seq(adj, V - {v} - set(adj[v])))
        E.append(indep_seq(adj, V - {v}))
    alpha = len(I) - 1
    # assert the two double counts at every level
    for k in range(1, alpha + 2):
        assert k * co(I, k) == sum(co(Jv, k - 1) for Jv in J), ("dc1", k)
    return I, J, E, alpha


def pv_lhs_rhs(I, Jv, k, slack=True):
    """PV_v(k): lhs <= rhs. slack=False gives the no-slack (Mason-type) form."""
    if slack:
        return k * co(I, k - 1) * co(Jv, k), (k + 1) * co(I, k) * co(Jv, k - 1)
    return co(I, k - 1) * co(Jv, k), co(I, k) * co(Jv, k - 1)


def pv_ratio(I, Jv, k):
    lhs, rhs = pv_lhs_rhs(I, Jv, k)
    if rhs == 0:
        return None if lhs == 0 else float("inf")
    return Fraction(lhs, rhs)


def lc_ok(I, k):
    return co(I, k - 1) * co(I, k + 1) <= co(I, k) ** 2
