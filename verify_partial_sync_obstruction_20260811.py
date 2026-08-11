#!/usr/bin/env python3
"""Locate partial-synchronization failures along the Bautista-Ramos
pattern-graph recurrences (arXiv:2603.14204).

The double-broom log-concavity proof (notes/double_broom_log_concavity.md)
runs a partial-synchronization induction (Hu--Wang--Zhao--Zhao, MIA 20
(2017) 91--103) through the connector recurrence F_t = F_{t-1} + x F_{t-2}.
The known non-log-concave tree families all fit the Bautista-Ramos
pattern-graph framework, which supplies positive two-term decompositions:

  k-direction (their Eq. (3), H = P_2 rooted at a leaf):
      F_k = (H\\w) F_{k-1} + x (H\\N[w]) H^{k-1} G
  m-direction (their Eq. (10)):
      U_m = (G\\v) Z_k(H)^m + x (G\\N[v]) H^{km},  Z_k(H) = H^k + x (H\\w)^k

Both summands are products of log-concave polynomials, hence log-concave.
By HWZZ, if the two summands were partially synchronized the sum would be
log-concave. These families break log-concavity, so summand partial
synchronization MUST fail. This script locates exactly where: first
failing parameter, violating index pairs (m, n), and how the violation
window sits relative to the log-concavity break indices.

Definition (partial synchronization, double-broom note Eq. (3)):
  A ~_p B  iff  a_m b_n + a_n b_m >= a_{m+1} b_{n-1} + a_{n-1} b_{m+1}
  for all m >= n, out-of-support coefficients read as 0.
(The relation is symmetric in A, B: swapping A and B fixes both sides.)

Exact integer arithmetic throughout. Verification per instance:
  1. decomposition polynomial == direct rooted-tree DP on the explicit tree;
  2. reflected head coefficients == closed forms of arXiv:2603.14204
     Corollaries 4.8, 4.12, 4.13, 4.14;
  3. LC break indices == the corollaries' stated indices in range.
Double brooms run as a negative control (zero violations expected).
"""

import json
import time

# ---------- exact polynomial helpers (ascending int lists) ----------

def pmul(a, b):
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                out[i + j] += ai * bj
    return out


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
            for i in range(n)]


def pshift(a, j):
    return [0] * j + list(a)


def ppow(a, e):
    out = [1]
    base = list(a)
    while e:
        if e & 1:
            out = pmul(out, base)
        e >>= 1
        if e:
            base = pmul(base, base)
    return out


# ---------- rooted trees and independence DP ----------

def tree_ipoly(n, adj, root=0):
    parent = [-1] * n
    order = []
    stack = [root]
    seen = [False] * n
    seen[root] = True
    while stack:
        v = stack.pop()
        order.append(v)
        for u in adj[v]:
            if not seen[u]:
                seen[u] = True
                parent[u] = v
                stack.append(u)
    f = [None] * n
    g = [None] * n
    for v in reversed(order):
        fv, gv = [1], [0, 1]
        for u in adj[v]:
            if parent[u] == v:
                fv = pmul(fv, padd(f[u], g[u]))
                gv = pmul(gv, f[u])
        f[v], g[v] = fv, gv
    return padd(f[root], g[root])


class RTree:
    def __init__(self):
        self.adj = []
        self.root = 0

    def new_vertex(self):
        self.adj.append([])
        return len(self.adj) - 1

    def edge(self, u, v):
        self.adj[u].append(v)
        self.adj[v].append(u)

    def graft(self, other):
        base = len(self.adj)
        for row in other.adj:
            self.adj.append([x + base for x in row])
        return base

    @property
    def n(self):
        return len(self.adj)

    def ipoly(self):
        return tree_ipoly(self.n, self.adj, self.root)


def path_tree(t):
    T = RTree()
    for _ in range(t):
        T.new_vertex()
    for i in range(t - 1):
        T.edge(i, i + 1)
    return T


def spider2(nlegs):
    """S_{2,n}: center root, n legs of length 2."""
    T = RTree()
    c = T.new_vertex()
    for _ in range(nlegs):
        a = T.new_vertex()
        b = T.new_vertex()
        T.edge(c, a)
        T.edge(a, b)
    return T


def T1_ell(ell):
    """T_{1,ell}: root v -- hub -- ell pendant cherries (P_2's)."""
    T = RTree()
    v = T.new_vertex()
    hub = T.new_vertex()
    T.edge(v, hub)
    for _ in range(ell):
        w = T.new_vertex()
        u = T.new_vertex()
        T.edge(hub, w)
        T.edge(w, u)
    return T


def attach_unit(T, H_builder, k):
    hub = T.new_vertex()
    T.edge(T.root, hub)
    for _ in range(k):
        H = H_builder()
        base = T.graft(H)
        T.edge(hub, base + H.root)
    return T


def pattern_tree(G_builder, H_builder, ks):
    T = G_builder()
    for k in ks:
        attach_unit(T, H_builder, k)
    return T


# ---------- diagnostics ----------

def coeff(a, i):
    return a[i] if 0 <= i < len(a) else 0


def lc_breaks(a):
    return [j for j in range(1, len(a) - 1)
            if a[j] * a[j] < a[j - 1] * a[j + 1]]


def is_unimodal(a):
    j = 0
    while j + 1 < len(a) and a[j] <= a[j + 1]:
        j += 1
    while j + 1 < len(a) and a[j] >= a[j + 1]:
        j += 1
    return j == len(a) - 1


def psync_violations(A, B):
    """All (m, n) with m >= n >= 0 violating
    a_m b_n + a_n b_m >= a_{m+1} b_{n-1} + a_{n-1} b_{m+1}."""
    D = max(len(A), len(B))
    out = []
    for m in range(D + 1):
        am = coeff(A, m)
        am1 = coeff(A, m + 1)
        bm = coeff(B, m)
        bm1 = coeff(B, m + 1)
        for n_ in range(m + 1):
            lhs = am * coeff(B, n_) + coeff(A, n_) * bm
            rhs = am1 * coeff(B, n_ - 1) + coeff(A, n_ - 1) * bm1
            if lhs < rhs:
                out.append((m, n_))
    return out


def summarize_pairs(pairs, alpha):
    if not pairs:
        return None
    ms = sorted({m for m, _ in pairs})
    ns = sorted({n for _, n in pairs})
    return {
        "count": len(pairs),
        "m_range": [ms[0], ms[-1]],
        "n_range": [ns[0], ns[-1]],
        "m_refl": [alpha - ms[-1], alpha - ms[0]],
        "n_refl": [alpha - ns[-1], alpha - ns[0]],
        "pairs_head": pairs[:12],
    }


# ---------- family instances ----------

P1 = [1, 1]
P2 = [1, 2]


def kl_3kk_instance(k):
    """(T_{1,3}:P_2)_k^(2); Eq.(3) at the second unit."""
    G = lambda: T1_ell(3)
    H = lambda: path_tree(2)
    F = pattern_tree(G, H, [k, k]).ipoly()
    order = pattern_tree(G, H, [k, k]).n
    base = pattern_tree(G, H, [k]).ipoly()
    Fprev = pattern_tree(G, H, [k, k - 1]).ipoly()
    A = pmul(P1, Fprev)                               # (H\w) F_{k-1}
    B = pshift(pmul(ppow(P2, k - 1), base), 1)        # x H^{k-1} G
    assert padd(A, B) == F, f"Eq.(3) mismatch at k={k}"
    return order, F, A, B


def m_family_instance(ell, nlegs, k, m):
    """(T_{1,ell}:S_{2,n})_k^(m); Eq.(10) decomposition."""
    G = lambda: T1_ell(ell)
    H = lambda: spider2(nlegs)
    tree = pattern_tree(G, H, [k] * m)
    F = tree.ipoly()
    Hp = H().ipoly()
    Hw = ppow(P2, nlegs)                              # S_{2,n} minus center
    Z = padd(ppow(Hp, k), pshift(ppow(Hw, k), 1))
    Gv = spider2(ell).ipoly()                         # T_{1,ell}\v = S_{2,ell}
    GNv = ppow(P2, ell)                               # T_{1,ell}\N[v]
    A = pmul(Gv, ppow(Z, m))
    B = pshift(pmul(GNv, ppow(Hp, k * m)), 1)
    assert padd(A, B) == F, f"Eq.(10) mismatch at ({ell},{nlegs},{k},{m})"
    return tree.n, F, A, B


def double_broom_instance(p, q, t):
    def build(tt):
        T = RTree()
        a = T.new_vertex()
        prev = a
        for _ in range(tt):
            v = T.new_vertex()
            T.edge(prev, v)
            prev = v
        b = T.new_vertex()
        T.edge(prev, b)
        for _ in range(p):
            T.edge(a, T.new_vertex())
        for _ in range(q):
            T.edge(b, T.new_vertex())
        return T
    F = build(t).ipoly()
    A = build(t - 1).ipoly()
    B = pshift(build(t - 2).ipoly(), 1)
    assert padd(A, B) == F, f"broom recurrence mismatch at t={t}"
    return build(t).n, F, A, B


# ---------- closed-form verification (arXiv:2603.14204) ----------

def closed_form_check(rec, RF):
    fam = rec["family"]
    if fam == "KL_3kk":
        k = rec["k"]
        return (RF[0] == 1
                and RF[1] == 2 ** (k + 1) + 2 * k + 11
                and RF[2] == (9 * 4 ** k + 11 * 2 ** (k + 1)
                              + 3 * k * 2 ** k + 2 * k * k + 21 * k + 15)
                and rec["alpha"] == 2 * k + 6)
    if fam == "T12_S24_2":
        m = rec["m"]
        return (RF[0] == 5 and RF[1] == 456 * m + 10
                and RF[2] == 47008 * m * m - 41668 * m + 6
                and rec["alpha"] == 10 * m + 3)
    if fam == "T13_S24_2":
        m = rec["m"]
        return (RF[0] == 9 and RF[1] == 616 * m + 23
                and RF[2] == 50208 * m * m - 41164 * m + 21
                and rec["alpha"] == 10 * m + 4)
    if fam == "T17_S25_2":
        m = rec["m"]
        return (RF[0] == 129 and RF[1] == 10570 * m + 583
                and RF[2] == 953266 * m * m - 566943 * m + 1141
                and rec["alpha"] == 12 * m + 8)
    return None


def run_instance(family, params, order, F, A, B):
    alpha = len(F) - 1
    breaks = lc_breaks(F)
    viol = psync_violations(A, B)
    rec = {
        "family": family, **params, "order": order, "alpha": alpha,
        "unimodal": is_unimodal(F),
        "lc_breaks": breaks,
        "lc_breaks_refl": [alpha - j for j in reversed(breaks)],
        "sync": summarize_pairs(viol, alpha),
    }
    RF = list(reversed(F))[:5]
    cf = closed_form_check(rec, RF)
    if cf is not None:
        rec["closed_form_ok"] = cf
    return rec


def main():
    t0 = time.time()
    results = []

    for t in range(2, 26):                       # control
        order, F, A, B = double_broom_instance(5, 5, t)
        results.append(run_instance("double_broom_ctrl",
                                    {"p": 5, "q": 5, "t": t},
                                    order, F, A, B))

    for k in range(2, 25):                       # KL 3,k,k; breaks for k>=4
        order, F, A, B = kl_3kk_instance(k)
        results.append(run_instance("KL_3kk", {"k": k}, order, F, A, B))

    for m in range(2, 61):                       # one break for m>=9
        order, F, A, B = m_family_instance(2, 4, 2, m)
        results.append(run_instance("T12_S24_2", {"m": m}, order, F, A, B))

    for m in range(2, 31):                       # two breaks for m>=16
        order, F, A, B = m_family_instance(3, 4, 2, m)
        results.append(run_instance("T13_S24_2", {"m": m}, order, F, A, B))

    for m in range(2, 25):                       # three breaks for m>=13
        order, F, A, B = m_family_instance(7, 5, 2, m)
        results.append(run_instance("T17_S25_2", {"m": m}, order, F, A, B))

    cfs = [r.get("closed_form_ok") for r in results
           if "closed_form_ok" in r]
    out = {
        "date": "2026-08-11",
        "source_paper": "arXiv:2603.14204",
        "hwzz": "Hu-Wang-Zhao-Zhao, Math. Inequal. Appl. 20 (2017) 91-103",
        "definition": ("A ~_p B iff a_m b_n + a_n b_m >= "
                       "a_{m+1} b_{n-1} + a_{n-1} b_{m+1} for m >= n"),
        "closed_form_checks_passed": sum(bool(x) for x in cfs),
        "closed_form_checks_total": len(cfs),
        "runtime_s": round(time.time() - t0, 2),
        "results": results,
    }
    with open("results/partial_sync_obstruction_20260811.json", "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"closed-form checks: {out['closed_form_checks_passed']}"
          f"/{out['closed_form_checks_total']} passed; "
          f"runtime {out['runtime_s']}s")
    fams = {}
    for r in results:
        fams.setdefault(r["family"], []).append(r)
    for fam, rs in fams.items():
        print(f"\n== {fam} ==")
        for r in rs:
            key = {k: v for k, v in r.items() if k in ("k", "m", "t")}
            sv = r["sync"]
            svs = ("none" if sv is None else
                   f"{sv['count']} viol; m {sv['m_range']} "
                   f"refl {sv['m_refl']}; n {sv['n_range']}")
            print(f"  {key} a={r['alpha']} brk={r['lc_breaks']} "
                  f"reflbrk={r['lc_breaks_refl']} uni={r['unimodal']} "
                  f"cf={r.get('closed_form_ok', '-')} | {svs}")


if __name__ == "__main__":
    main()
