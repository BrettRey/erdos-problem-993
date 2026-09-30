"""Cross-check r1_census (C) against route-freecount/pv_lib.py (Python,
independent forest DP per vertex) on every tree with n <= NMAX.

For every tree and every window level k, compares exactly:
  - the window [lo, q] itself (C emits one DUMP line per window level),
  - I(T), every J_v = I(T - N[v]),
  - every D_v(k),
  - every LAM * L_u(k) against Fraction L_u from pv_lib-based computation.
Also checks the C STATS tree count and that trees with an empty window emit
nothing.  C vertex v (1-indexed) <-> Python vertex v-1.
Usage: python3 crosscheck_pvlib.py NMAX > crosscheck.json
"""
import json
import subprocess
import sys
from fractions import Fraction

sys.path.insert(0, "/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/runs/attack-20260930/route-freecount")
from pv_lib import co, parse_parent_array, tree_data, window  # noqa: E402

HERE = "/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/runs/attack-20260930/r1-exhaustive"


def main():
    nmax = int(sys.argv[1])
    report = {}
    for n in range(2, nmax + 1):
        trees = subprocess.run(["gentreeg", "-p", "-q", str(n)], capture_output=True, text=True).stdout.split("\n")
        trees = [t.split() for t in trees if t.strip()]
        cout = subprocess.run([f"{HERE}/r1_census", str(n), "dump"], input="\n".join(" ".join(t) for t in trees) + "\n",
                              capture_output=True, text=True).stdout.split("\n")
        dumps = {}
        stats = None
        for line in cout:
            if line.startswith("DUMP"):
                f = dict(x.split("=", 1) for x in line.split()[1:])
                dumps.setdefault(f["par"], []).append(f)
            elif line.startswith("STATS"):
                stats = dict(x.split("=", 1) for x in line.split()[1:])
        assert int(stats["trees"]) == len(trees), (n, stats["trees"], len(trees))
        n_levels = n_D = n_L = 0
        for toks in trees:
            key = ",".join(toks)
            adj = parse_parent_array(toks)
            I, J, E, alpha = tree_data(adj)
            lo, q = window(n, alpha)
            lo = max(lo, 1)
            ks = list(range(lo, q + 1))
            got = dumps.get(key, [])
            assert [int(f["k"]) for f in got] == ks, (n, key, ks, [f["k"] for f in got])
            for f in got:
                k = int(f["k"])
                LAM = int(f["LAM"])
                Ic = [int(x) for x in f["I"].split(",")]
                assert Ic == I, (key, Ic, I)
                Jc = [[int(x) for x in s.split(",")] for s in f["J"].split(";")]
                for v in range(n):
                    a = list(Jc[v]); b = list(J[v])
                    while len(a) > 1 and a[-1] == 0: a.pop()
                    while len(b) > 1 and b[-1] == 0: b.pop()
                    assert a == b, (key, v, a, b)
                D = [k * co(I, k - 1) * co(J[v], k) - (k + 1) * co(I, k) * co(J[v], k - 1) for v in range(n)]
                Dc = [int(x) for x in f["D"].split(",")]
                assert Dc == D, (key, k, Dc, D)
                LLc = [int(x) for x in f["LL"].split(",")]
                for u in range(n):
                    L = Fraction(D[u], len(adj[u]) + 1)
                    for v in adj[u]:
                        L += Fraction(D[v], len(adj[v]) + 1)
                    assert Fraction(LLc[u], LAM) == L, (key, k, u, LLc[u], L)
                    n_L += 1
                n_D += n
                n_levels += 1
        report[n] = {"trees": len(trees), "window_levels": n_levels, "D_values": n_D, "L_values": n_L,
                     "c_stats_r1_viol": int(stats["r1_viol"]), "c_stats_pv_viol": int(stats["pv_viol"])}
        print(f"n={n} trees={len(trees)} levels={n_levels} D={n_D} L={n_L} OK", file=sys.stderr, flush=True)
    json.dump({"status": "all exact agreements", "per_n": report}, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
