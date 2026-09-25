"""Exact free-vertex statistics over uniform k-independent sets of trees.

For an independent set S of a tree T, let
  f(S) = number of free vertices = |V \\ N[S]|,
  e(S) = number of edges of T with both endpoints free,
  c(S) = f(S) - e(S) = number of components of the free forest.

We compute, by a 3-state tree DP, the trivariate generating polynomial
  G(x, y, z) = sum_S x^{|S|} y^{f(S)} z^{e(S)}
truncated to the moments we need (f, f^2, e), and then verify the exact identity

  (k+1)(k+2) i_{k+2} = sum_{|S|=k} [f^2 - f - 2e]

and report, per (T, k), the dispersion diagnostics:
  mu_k   = E_k[f]                     (so i_{k+1}/i_k = mu_k/(k+1))
  Var_k  = Var_k(f)
  Ec_k   = E_k[c]
  valley-at-(k+1) margin:  mu(k+5-mu) - Var - 2Ec  (must be >=0 whenever mu < k+1)
  LC-at-(k+1) margin:      mu^2/(k+1) + 3mu - Var - 2Ec

Everything is exact (Fractions).
"""
from fractions import Fraction
import sys
import networkx as nx

# state of a vertex relative to S: 0 = in S, 1 = neighbour of S (dominated), 2 = free
# We do a rooted DP; for each vertex v we keep, for each state of v, a dict
# k -> (m0, m1, m2, me) where m0 = count, m1 = sum f, m2 = sum f^2, me = sum e
# over independent sets of the subtree rooted at v, restricted to the subtree,
# where the "state" of v is provisional: state 1 (dominated) must be certified
# by some child in S OR by the parent; state 2 (free) requires no child in S and
# parent not in S.  We handle this by letting the subtree DP compute, for v:
#   A[k] : v in S
#   B[k] : v not in S, some child in S   (v dominated from below)
#   C[k] : v not in S, no child in S     (v undetermined: free unless parent in S)
# Free-vertex count contributed by v is decided when we know the parent state,
# so f-moments are accumulated as polynomials in a formal marker; simplest is
# to carry moments of (f, e) *excluding* v's own contribution and add it at the
# parent.  We instead carry full moments but "commit" v as free in C, and when
# the parent is in S we subtract; subtraction of moments needs joint info, so
# we carry the pair (f, e) via moments m1=sum f, m2=sum f^2, me=sum e, and
# also mfe = sum f*e is not needed.  To subtract 1 from f for all sets in C we
# need sum (f-1)^2 = m2 - 2 m1 + m0: fine.  Edges: an edge (v, child) is free
# iff both endpoints are free; child free <=> child in state C and v not in S;
# v free <=> v in C and parent not in S.  So the edge (v,child) is free iff
# v in C, child in C, and parent not in S.  We commit edges (v,child) with
# child in C as free when v is in C, and subtract them if the parent is in S.
# For that subtraction we need, for v in C, the count of free child-edges;
# we carry for state C an extra moment: mq = sum q where q = #children in C,
# and the joint sum m1q = sum f*q is not needed since we only need sums of
# e - q and f - 1 separately (second moment only of f).  Good.

def combine(dicts_list):
    """Convolve a list of per-child dicts k->(m0,m1,m2,me) (moments of f,e) as
    independent products: sizes add, f adds, e adds."""
    acc = {0: (1, 0, 0, 0)}
    for d in dicts_list:
        new = {}
        for k1, (a0, a1, a2, ae) in acc.items():
            for k2, (b0, b1, b2, be) in d.items():
                k = k1 + k2
                # f = f1 + f2 ; sum f = a1*b0 + a0*b1 ; sum f^2 = a2*b0 + 2 a1 b1 + a0 b2
                o = new.get(k, (0, 0, 0, 0))
                new[k] = (o[0] + a0 * b0,
                          o[1] + a1 * b0 + a0 * b1,
                          o[2] + a2 * b0 + 2 * a1 * b1 + a0 * b2,
                          o[3] + ae * b0 + a0 * be)
        acc = new
    return acc

def add_dicts(*ds):
    out = {}
    for d in ds:
        for k, v in d.items():
            o = out.get(k, (0, 0, 0, 0))
            out[k] = tuple(o[i] + v[i] for i in range(4))
    return out

def shift_f(d, delta):
    """add delta to f for every set (delta may be -1, +1)."""
    out = {}
    for k, (m0, m1, m2, me) in d.items():
        out[k] = (m0, m1 + delta * m0, m2 + 2 * delta * m1 + delta * delta * m0, me)
    return out

def shift_e(d, delta_counts):
    # delta_counts: dict k -> sum over sets of delta_e  (we pass sums directly)
    out = {}
    for k, (m0, m1, m2, me) in d.items():
        out[k] = (m0, m1, m2, me + delta_counts.get(k, 0))
    return out

def rooted_dp(T, root):
    order = list(nx.dfs_preorder_nodes(T, root))
    parent = {root: None}
    for u in order:
        for w in T[u]:
            if w != parent[u] and w not in parent:
                parent[w] = u
    children = {u: [w for w in T[u] if w != parent[u]] for u in order}
    # returns per vertex three dicts A (v in S), B (v dominated from below), C (v undetermined)
    # In C, v is committed as free, and edges (v,child in C) committed free; we
    # also store Q[k] = sum over sets of q = #children-in-C so a parent in S can uncommit.
    memo = {}
    for u in reversed(order):
        ch = children[u]
        A_u = combine([add_dicts(memo[w][1], memo[w][2]) for w in ch])  # children not in S: B or C ; u in S
        # children in state C under a parent in S: their f loses 1 and their committed edges (child, grandchild in C) are unaffected (grandchildren still free? grandchild free needs child not in S: yes, child not in S) -- only the child's own freeness is lost, and edges (child, grandchild) need both free: child is no longer free, so those edges lose freeness too!
        # So under parent in S, a child in C loses: f -= 1, e -= q_child.
        # Redo A_u properly: for each child, choose B (unchanged) or C' (C with f-1, e-q).
        parts = []
        for w in ch:
            Cw = memo[w][2]
            Qw = memo[w][3]
            Cprime = shift_f(Cw, -1)
            Cprime = {k: (m0, m1, m2, me - Qw.get(k, 0)) for k, (m0, m1, m2, me) in Cprime.items()}
            parts.append(add_dicts(memo[w][1], Cprime))
        A_u = combine(parts)
        # shift size by 1 (u in S)
        A_u = {k + 1: v for k, v in A_u.items()}
        # u not in S: each child in A, B, or C.
        # B_u: at least one child in A.  C_u: no child in A (all children B or C); u committed free (+1 f), edges (u, child in C) committed free.
        # For C_u we need Q_u[k] = sum over sets of (#children in C).
        allBC = combine([add_dicts(memo[w][1], memo[w][2]) for w in ch])
        allABC = combine([add_dicts(memo[w][0], memo[w][1], memo[w][2]) for w in ch])
        B_u = {}
        for k, v in allABC.items():
            o = allBC.get(k, (0, 0, 0, 0))
            B_u[k] = tuple(v[i] - o[i] for i in range(4))
        B_u = {k: v for k, v in B_u.items() if v[0] != 0}
        # C_u: children in B or C, with q = #children in C.  Need sum q and the e-shift by q.
        # Compute by convolving with a marker: track (m0, m1, m2, me, mq) where mq = sum q.
        acc = {0: (1, 0, 0, 0, 0)}
        for w in ch:
            Bw, Cw = memo[w][1], memo[w][2]
            new = {}
            for k1, (a0, a1, a2, ae, aq) in acc.items():
                for k2, (b0, b1, b2, be) in Bw.items():
                    k = k1 + k2
                    o = new.get(k, (0, 0, 0, 0, 0))
                    new[k] = (o[0] + a0 * b0, o[1] + a1 * b0 + a0 * b1, o[2] + a2 * b0 + 2 * a1 * b1 + a0 * b2,
                              o[3] + ae * b0 + a0 * be, o[4] + aq * b0)
                for k2, (b0, b1, b2, be) in Cw.items():
                    k = k1 + k2
                    o = new.get(k, (0, 0, 0, 0, 0))
                    new[k] = (o[0] + a0 * b0, o[1] + a1 * b0 + a0 * b1, o[2] + a2 * b0 + 2 * a1 * b1 + a0 * b2,
                              o[3] + ae * b0 + a0 * be, o[4] + aq * b0 + a0 * b0)  # q += 1
            acc = new
        C_u = {}
        Q_u = {}
        for k, (m0, m1, m2, me, mq) in acc.items():
            # commit u free: f += 1 ; commit edges: e += q
            C_u[k] = (m0, m1 + m0, m2 + 2 * m1 + m0, me + mq)
            Q_u[k] = mq
        memo[u] = (A_u, B_u, C_u, Q_u)
    A, B, C, Q = memo[root]
    return add_dicts(A, B, C)

def indpoly(T):
    # independent DP for i_k as a check
    root = next(iter(T.nodes))
    order = list(nx.dfs_preorder_nodes(T, root))
    parent = {root: None}
    for u in order:
        for w in T[u]:
            if w != parent[u] and w not in parent:
                parent[w] = u
    P = {}; R = {}
    for u in reversed(order):
        ch = [w for w in T[u] if w != parent[u]]
        p = [0, 1]
        r = [1]
        for w in ch:
            p = polymul(p, R[w])
            r = polymul(r, polyadd(P[w], R[w]))
        P[u] = p; R[u] = r
    return polyadd(P[root], R[root])

def polymul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out

def polyadd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]

def analyze(T):
    G = rooted_dp(T, next(iter(T.nodes)))
    ipoly = indpoly(T)
    alpha = len(ipoly) - 1
    rows = []
    for k in range(alpha + 1):
        m0, m1, m2, me = G.get(k, (0, 0, 0, 0))
        assert m0 == ipoly[k], (k, m0, ipoly[k])
        # identity: (k+1) i_{k+1} = sum f ; (k+1)(k+2) i_{k+2} = sum (f^2 - f - 2e)
        ik1 = ipoly[k + 1] if k + 1 <= alpha else 0
        ik2 = ipoly[k + 2] if k + 2 <= alpha else 0
        assert m1 == (k + 1) * ik1, (k, m1, ik1)
        assert m2 - m1 - 2 * me == (k + 1) * (k + 2) * ik2, (k, m2, m1, me, ik2)
        mu = Fraction(m1, m0)
        var = Fraction(m2, m0) - mu * mu
        Ee = Fraction(me, m0)
        Ec = mu - Ee
        rows.append((k, m0, mu, var, Ec))
    return ipoly, rows

if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    worst_valley = []   # (margin/mu, n, k, ...)
    worst_disp = []
    count = 0
    for n in range(2, nmax + 1):
        for T in nx.nonisomorphic_trees(n):
            count += 1
            ipoly, rows = analyze(T)
            alpha = len(ipoly) - 1
            thr = -(-(2 * alpha - 1) // 3)
            for (k, m0, mu, var, Ec) in rows:
                if k + 2 > alpha:
                    continue
                D = (var + 2 * Ec) / mu if mu else None
                # valley-at-(k+1) region: mu < k+1
                if mu < k + 1:
                    margin = mu * (k + 5 - mu) - var - 2 * Ec
                    worst_valley.append((float(margin / (mu * (k + 5 - mu))), n, k, alpha, thr, float(mu), float(var), float(Ec), nx.to_graph6_bytes(T, header=False).decode().strip()))
                if D is not None:
                    worst_disp.append((float(D), n, k, alpha, thr, float(mu), float(var), float(Ec), float(mu) < k + 1, nx.to_graph6_bytes(T, header=False).decode().strip()))
    print("trees:", count)
    worst_valley.sort()
    print("\nSmallest relative valley margins (mu<k+1 region): margin/(mu(k+5-mu))")
    for r in worst_valley[:12]:
        print(r)
    worst_disp.sort(reverse=True)
    print("\nLargest dispersion index D=(Var+2E[c])/mu overall:")
    for r in worst_disp[:8]:
        print(r)
    print("\nLargest D restricted to mu<k+1:")
    for r in [r for r in worst_disp if r[8]][:8]:
        print(r)
    print("\nLargest D restricted to mu<k+1 and k <= thr-2:")
    for r in [r for r in worst_disp if r[8] and r[2] <= r[4] - 2][:8]:
        print(r)
