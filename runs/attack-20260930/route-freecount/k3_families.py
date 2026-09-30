"""K3: exact PV kill-test on structured adversarial families at larger n.

For each tree: max over (v, k in window) of the PV ratio
  k i_{k-1}(T) j_k / ((k+1) i_k(T) j_{k-1}),
every PV failure in the window (exact), and every PV failure at any k.
Output: JSON lines to stdout.
"""

import json
import sys
from fractions import Fraction

from pv_lib import co, lc_ok, pv_lhs_rhs, tree_data, window


def edges_to_adj(n, edges):
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    return adj


class B:
    def __init__(self):
        self.n = 1
        self.edges = []

    def add(self, parent):
        v = self.n
        self.n += 1
        self.edges.append((parent, v))
        return v

    def path(self, parent, length):
        cur = parent
        for _ in range(length):
            cur = self.add(cur)
        return cur

    def adj(self):
        return edges_to_adj(self.n, self.edges)


def spider(legs):
    b = B()
    for L in legs:
        b.path(0, L)
    return b.adj()


def hub_stars(m, t):
    """hub -> m children -> each with t leaves"""
    b = B()
    for _ in range(m):
        c = b.add(0)
        for _ in range(t):
            b.add(c)
    return b.adj()


def bouquet_P2(m, t):
    """root with m children, each carrying t pendant P2 legs (mason note T_{m,t})"""
    b = B()
    for _ in range(m):
        c = b.add(0)
        for _ in range(t):
            b.path(c, 2)
    return b.adj()


def broom_Tab(a, bb):
    """SCC-failure broom: root with a children, each with bb children, each with 1 leaf"""
    b = B()
    for _ in range(a):
        w = b.add(0)
        for _ in range(bb):
            x = b.add(w)
            b.add(x)
    return b.adj()


def double_star(a, c):
    b = B()
    u = b.add(0)
    for _ in range(a):
        b.add(0)
    for _ in range(c):
        b.add(u)
    return b.adj()


def caterpillar(spine, legs):
    b = B()
    prev = 0
    for i in range(legs[0]):
        b.add(0)
    for s in range(1, spine):
        cur = b.add(prev)
        for _ in range(legs[s]):
            b.add(cur)
        prev = cur
    return b.adj()


def analyse(adj):
    n = len(adj)
    I, J, E, alpha = tree_data(adj)
    lo, q = window(n, alpha)
    lo = max(lo, 1)
    best = None
    wit = None
    fails_w = []
    fails_any = []
    lc_fail = []
    for k in range(1, alpha + 1):
        if not lc_ok(I, k):
            lc_fail.append(k)
        for v in range(n):
            lhs, rhs = pv_lhs_rhs(I, J[v], k)
            if lhs > rhs:
                fails_any.append((k, v, len(adj[v])))
                if lo <= k <= q:
                    fails_w.append((k, v, len(adj[v])))
            if lo <= k <= q and rhs > 0:
                r = Fraction(lhs, rhs)
                if best is None or r > best:
                    best = r
                    wit = (k, v, len(adj[v]))
    return {
        "n": n, "alpha": alpha, "lo": lo, "q": q,
        "max_pv_ratio_window_float": float(best) if best is not None else None,  # diagnostic
        "one_minus_max_times_n": float((1 - best) * n) if best is not None else None,  # diagnostic
        "witness_k_v_deg": wit,
        "pv_fail_window": fails_w[:10], "n_pv_fail_window": len(fails_w),
        "pv_fail_any_k_minus_alpha": sorted(set(k - alpha for k, _, _ in fails_any)),
        "n_pv_fail_any": len(fails_any),
        "lc_fail_k_minus_alpha": [k - alpha for k in lc_fail],
    }


def main():
    fams = []
    for m in range(2, 41):
        fams.append((f"spider_2^{m}", spider([2] * m)))
    for m in [5, 10, 20, 30]:
        for a1 in [1, 3, 5, 10]:
            fams.append((f"spider_1^{a1}_2^{m}", spider([1] * a1 + [2] * m)))
    for m in range(2, 26):
        for t in [1, 2, 3, 4, 6]:
            if m * (t + 1) + 1 <= 130:
                fams.append((f"hubstars_{m}_{t}", hub_stars(m, t)))
    for m in range(2, 16):
        for t in [1, 2, 3, 4]:
            if 1 + m * (1 + 2 * t) <= 130:
                fams.append((f"bouquetP2_{m}_{t}", bouquet_P2(m, t)))
    for a in range(2, 6):
        for bb in range(2, 7):
            if 1 + a * (1 + 2 * bb) <= 130:
                fams.append((f"broom_{a}_{bb}", broom_Tab(a, bb)))
    for a in [5, 10, 20, 40]:
        for c in [1, 5, 10, 20, 40]:
            fams.append((f"doublestar_{a}_{c}", double_star(a, c)))
    for L in [10, 20, 40]:
        for leg in [1, 2, 3]:
            fams.append((f"caterpillar_{L}x{leg}", caterpillar(L, [leg] * L)))
    for name, adj in fams:
        rec = analyse(adj)
        rec["family"] = name
        print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
