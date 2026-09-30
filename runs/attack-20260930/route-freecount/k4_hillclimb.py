"""K4: adversarial hill-climb against PV (pointwise) and its local-averaging
repairs R1 (closed neighbourhood) and R2 (edge), inside the window.

  PV_v(k):  P_k(v free) <= (1+1/k) P_{k-1}(v free)      [D_v <= 0]
  R1_u(k):  L_u = sum_{v in N[u]} D_v/(deg v + 1) <= 0
  R2_uv(k): D_u/deg u + D_v/deg v <= 0   for every edge uv
with D_v = k i_{k-1}(T) j^v_k - (k+1) i_k(T) j^v_{k-1}.

Scores are floats (diagnostic, for steering only); every reported violation
is re-checked exactly with integers/Fractions before it is recorded.

Usage: python3 k4_hillclimb.py MODE NMIN NMAX SECONDS SEED > out.json
MODE in {pv, r1, r2}.
"""

import json
import random
import sys
import time
from fractions import Fraction

from pv_lib import co, indep_seq, window


def adj_from_edges(n, edges):
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    return adj


def seqs(adj):
    n = len(adj)
    V = set(range(n))
    I = indep_seq(adj, V)
    J = [indep_seq(adj, V - {v} - set(adj[v])) for v in range(n)]
    return I, J


def defects(I, J, k):
    return [k * co(I, k - 1) * co(Jv, k) - (k + 1) * co(I, k) * co(Jv, k - 1) for Jv in J]


def evaluate(edges, n, mode):
    adj = adj_from_edges(n, edges)
    I, J = seqs(adj)
    alpha = len(I) - 1
    lo, q = window(n, alpha)
    lo = max(lo, 1)
    best = -1e300
    viol = None
    for k in range(lo, q + 1):
        D = defects(I, J, k)
        scale = sum(abs(d) for d in D) / n
        if scale == 0:
            continue
        if mode == "pv":
            for v in range(n):
                den = (k + 1) * co(I, k) * co(J[v], k - 1)
                if den == 0:
                    continue
                s = (k * co(I, k - 1) * co(J[v], k)) / den  # float, steering only
                if s > best:
                    best = s
                if D[v] > 0:
                    viol = viol or {"k": k, "v": v, "deg": len(adj[v]), "D": D[v]}
        elif mode == "r1":
            for u in range(n):
                L = Fraction(D[u], len(adj[u]) + 1)
                for v in adj[u]:
                    L += Fraction(D[v], len(adj[v]) + 1)
                s = float(L) / scale
                if s > best:
                    best = s
                if L > 0:
                    viol = viol or {"k": k, "u": u, "deg": len(adj[u]), "L": str(L)}
        elif mode == "r2":
            for a, b in edges:
                L = Fraction(D[a], len(adj[a])) + Fraction(D[b], len(adj[b]))
                s = float(L) / scale
                if s > best:
                    best = s
                if L > 0:
                    viol = viol or {"k": k, "edge": [a, b], "L": str(L)}
    return best, viol, alpha, lo, q


def random_tree(n, rng):
    return [(i, rng.randrange(i)) for i in range(1, n)]


def hubstars_edges(m, t):
    edges = []
    nxt = 1
    for _ in range(m):
        c = nxt
        nxt += 1
        edges.append((0, c))
        for _ in range(t):
            edges.append((c, nxt))
            nxt += 1
    return edges, nxt


def mutate(edges, n, rng, nmin, nmax):
    edges = list(edges)
    deg = [0] * n
    for a, b in edges:
        deg[a] += 1
        deg[b] += 1
    leaves = [v for v in range(n) if deg[v] == 1]
    r = rng.random()
    if r < 0.55 or not (nmin < n or n < nmax):
        # move a leaf elsewhere
        x = rng.choice(leaves)
        edges = [e for e in edges if x not in e]
        tgt = rng.randrange(n - 1)
        if tgt >= x:
            tgt += 1
        edges.append((x, tgt))
        return edges, n
    if r < 0.78 and n < nmax:
        # add a leaf somewhere
        edges.append((n, rng.randrange(n)))
        return edges, n + 1
    if n > nmin:
        # remove a leaf (relabel last vertex into its slot)
        x = rng.choice(leaves)
        edges = [e for e in edges if x not in e]
        last = n - 1
        if x != last:
            edges = [tuple(x if w == last else w for w in e) for e in edges]
        return edges, n - 1
    edges.append((n, rng.randrange(n)))
    return edges, n + 1


def main():
    mode, nmin, nmax, secs, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
    rng = random.Random(seed)
    t0 = time.time()
    evals = 0
    found = []
    global_best = (-1e300, None)
    restarts = 0
    while time.time() - t0 < secs:
        restarts += 1
        if restarts % 3 == 1:
            m = rng.randrange(4, 12)
            t = rng.choice([1, 2, 3])
            edges, n = hubstars_edges(m, t)
            while n > nmax:
                m -= 1
                edges, n = hubstars_edges(m, t)
            while n < nmin:
                edges.append((n, rng.randrange(n)))
                n += 1
        else:
            n = rng.randrange(nmin, nmax + 1)
            edges = random_tree(n, rng)
        cur, viol, *_ = evaluate(edges, n, mode)
        evals += 1
        stall = 0
        while stall < 300 and time.time() - t0 < secs:
            e2, n2 = mutate(edges, n, rng, nmin, nmax)
            s2, v2, alpha, lo, q = evaluate(e2, n2, mode)
            evals += 1
            if v2 is not None and len(found) < 50:
                key = (n2, tuple(sorted(tuple(sorted(e)) for e in e2)))
                if all(f["key"] != str(key) for f in found):
                    found.append({"key": str(key), "n": n2, "edges": e2, "alpha": alpha,
                                  "window": [lo, q], "violation": v2, "score": s2})
            if s2 >= cur:
                if s2 > cur:
                    stall = 0
                else:
                    stall += 1
                edges, n, cur = e2, n2, s2
            else:
                stall += 1
            if cur > global_best[0]:
                global_best = (cur, {"n": n, "edges": edges})
    found.sort(key=lambda f: f["n"])
    json.dump({"mode": mode, "nmin": nmin, "nmax": nmax, "seconds": secs, "seed": seed,
               "evaluations": evals, "restarts": restarts,
               "best_score_float": global_best[0], "best_tree": global_best[1],
               "violations_found": len(found), "smallest_violations": found[:10]},
              sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
