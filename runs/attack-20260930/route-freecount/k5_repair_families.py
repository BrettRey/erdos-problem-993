"""K5: exact test of the repairs R1 (closed-neighbourhood averaging) and R2
(edge averaging) on families built around the pointwise-PV counterexample
(hub-stars), where pointwise PV actually fails in the window.

Output JSON lines: per tree, number of window levels with PV / R1 / R2
failures, and float diagnostics (max normalised loads).
"""

import json
from fractions import Fraction

from k3_families import B
from pv_lib import co, tree_data, window


def hub_mixed(ts, extra_leaves=0, tail=0):
    """hub -> one child per entry of ts, child i carrying ts[i] leaves;
    plus `extra_leaves` pendant leaves at the hub and a pendant path of length `tail`."""
    b = B()
    for t in ts:
        c = b.add(0)
        for _ in range(t):
            b.add(c)
    for _ in range(extra_leaves):
        b.add(0)
    if tail:
        b.path(0, tail)
    return b.adj()


def two_hubs(m1, m2, t, bridge):
    b = B()
    for _ in range(m1):
        c = b.add(0)
        for _ in range(t):
            b.add(c)
    end = b.path(0, bridge)
    for _ in range(m2):
        c = b.add(end)
        for _ in range(t):
            b.add(c)
    return b.adj()


def hub_of_hubs(m, s, t):
    """hub -> m children; each child -> s grandchildren; each grandchild -> t leaves"""
    b = B()
    for _ in range(m):
        c = b.add(0)
        for _ in range(s):
            g = b.add(c)
            for _ in range(t):
                b.add(g)
    return b.adj()


def check(adj):
    n = len(adj)
    I, J, _E, alpha = tree_data(adj)
    lo, q = window(n, alpha)
    lo = max(lo, 1)
    res = {"n": n, "alpha": alpha, "lo": lo, "q": q, "pv_fail_levels": 0,
           "r1_fail_levels": 0, "r2_fail_levels": 0, "lc_fail_levels": 0,
           "max_r1_norm": -1e9, "max_r2_norm": -1e9, "pv_fail_degs": []}
    edges = [(u, v) for u in range(n) for v in adj[u] if u < v]
    for k in range(lo, q + 1):
        D = [k * co(I, k - 1) * co(Jv, k) - (k + 1) * co(I, k) * co(Jv, k - 1) for Jv in J]
        if sum(D) > 0:
            res["lc_fail_levels"] += 1
        scale = sum(abs(d) for d in D) / n
        if any(d > 0 for d in D):
            res["pv_fail_levels"] += 1
            res["pv_fail_degs"] = sorted(set(res["pv_fail_degs"]) | {len(adj[v]) for v in range(n) if D[v] > 0})
        r1bad = False
        for u in range(n):
            L = Fraction(D[u], len(adj[u]) + 1) + sum(Fraction(D[v], len(adj[v]) + 1) for v in adj[u])
            res["max_r1_norm"] = max(res["max_r1_norm"], float(L) / scale)
            r1bad |= L > 0
        r2bad = False
        for a, b in edges:
            L = Fraction(D[a], len(adj[a])) + Fraction(D[b], len(adj[b]))
            res["max_r2_norm"] = max(res["max_r2_norm"], float(L) / scale)
            r2bad |= L > 0
        res["r1_fail_levels"] += r1bad
        res["r2_fail_levels"] += r2bad
    return res


def main():
    fams = []
    for t in [1, 2, 3, 4, 5]:
        for m in range(4, 60):
            if 1 + m * (t + 1) <= 110:
                fams.append((f"hubstars_{m}_{t}", hub_mixed([t] * m)))
    for m in range(8, 30, 3):
        for e in [1, 3, 6]:
            fams.append((f"hubstars_{m}_2+{e}leaves", hub_mixed([2] * m, extra_leaves=e)))
        for tail in [1, 2, 5, 10]:
            fams.append((f"hubstars_{m}_2+tail{tail}", hub_mixed([2] * m, tail=tail)))
        for a in [2, 4, 8]:
            fams.append((f"hubmixed_{m}x2_{a}x1", hub_mixed([2] * m + [1] * a)))
            fams.append((f"hubmixed_{m}x2_{a}x3", hub_mixed([2] * m + [3] * a)))
    for m in [9, 12, 15, 20]:
        for bridge in [1, 2, 3]:
            fams.append((f"twohubs_{m}_{m}_2_br{bridge}", two_hubs(m, m, 2, bridge)))
            fams.append((f"twohubs_{m}_4_2_br{bridge}", two_hubs(m, 4, 2, bridge)))
    for m in [3, 5, 8, 12]:
        for s in [2, 3, 4]:
            for t in [1, 2]:
                if 1 + m * (1 + s * (1 + t)) <= 130:
                    fams.append((f"hubofhubs_{m}_{s}_{t}", hub_of_hubs(m, s, t)))
    for name, adj in fams:
        r = check(adj)
        r["family"] = name
        print(json.dumps(r), flush=True)


if __name__ == "__main__":
    main()
