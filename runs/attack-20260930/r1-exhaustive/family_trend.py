"""Trend of the tightest R1 load along the hub-star family H(m,t) (hub joined
to m centres, each carrying t leaves), which the exhaustive census identifies
as the extremisers for 19 <= n <= 26.  Exact Fractions via pv_lib.
norm1 = L_u / (k i_{k-1} i_k),  norm2 = L_u / (sum_v |D_v| / n).
Usage: python3 family_trend.py NMAX > family_trend.json
"""
import json
import sys
from fractions import Fraction

sys.path.insert(0, "/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/runs/attack-20260930/route-freecount")
from pv_lib import co, tree_data, window  # noqa: E402


def hub_stars(m, t):
    adj = [[]]
    for _ in range(m):
        c = len(adj); adj.append([0]); adj[0].append(c)
        for _ in range(t):
            leaf = len(adj); adj.append([c]); adj[c].append(leaf)
    return adj


def scan(adj):
    n = len(adj)
    I, J, E, alpha = tree_data(adj)
    lo, q = window(n, alpha)
    lo = max(lo, 1)
    best1 = best2 = None
    for k in range(lo, q + 1):
        D = [k * co(I, k - 1) * co(J[v], k) - (k + 1) * co(I, k) * co(J[v], k - 1) for v in range(n)]
        absS = sum(abs(d) for d in D)
        for u in range(n):
            L = Fraction(D[u], len(adj[u]) + 1) + sum(Fraction(D[v], len(adj[v]) + 1) for v in adj[u])
            v1 = L / (k * co(I, k - 1) * co(I, k))
            if best1 is None or v1 > best1[0]:
                best1 = (v1, k, u, len(adj[u]))
            if absS:
                v2 = L * n / absS
                if best2 is None or v2 > best2[0]:
                    best2 = (v2, k, u, len(adj[u]))
    return n, alpha, lo, q, best1, best2


def main():
    nmax = int(sys.argv[1])
    rows = []
    for m in range(2, 16):
        for t in range(1, 16):
            n = 1 + m * (1 + t)
            if n > nmax:
                continue
            n, alpha, lo, q, b1, b2 = scan(hub_stars(m, t))
            row = {"m": m, "t": t, "n": n, "alpha": alpha, "window": [lo, q],
                   "norm1_max": str(b1[0]), "norm1_float_diag": float(b1[0]), "norm1_k": b1[1], "norm1_u": b1[2], "norm1_deg_u": b1[3],
                   "norm2_max": str(b2[0]), "norm2_float_diag": float(b2[0]), "norm2_k": b2[1], "norm2_u": b2[2], "norm2_deg_u": b2[3],
                   "r1_holds": b1[0] <= 0}
            rows.append(row)
            print(f"H({m},{t}) n={n} a={alpha} W=[{lo},{q}] n1={float(b1[0]):.5f}@k{b1[1]},deg{b1[3]} n2={float(b2[0]):.4f}@k{b2[1]},deg{b2[3]}", file=sys.stderr, flush=True)
    json.dump(rows, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
