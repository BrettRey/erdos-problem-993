"""Simulated annealing against lemma R1 over trees with NMIN <= n <= NMAX.
Objective OBJ in {F, A, E}: maximise max_{u, k in W} norm_OBJ(u,k) (see r1lib).
Annealed on g = -log(-s) (s<0); any s >= 0 is an exact violation (the exact
integer test in r1lib.evaluate decides; floats only steer).
Moves: leaf move, leaf add/remove, cherry add/remove, subdivide/contract,
subtree prune-and-regraft.  Restarts from structured seeds and random trees.
Usage: python sa_r1.py OBJ NMIN NMAX SECONDS SEED OUT.json [SEEDFILE.jsonl]
"""
import json
import math
import os
import random
import sys
import time
from r1lib import adj_from_edges, evaluate, B


def edges_of(adj):
    return [(u, v) for u in range(len(adj)) for v in adj[u] if u < v]


def relabel(adj):
    return adj  # adjacency lists are always on 0..n-1


def remove_vertex(adj, x):
    """remove vertex x (assumed leaf or isolated after edits); relabel n-1 -> x"""
    n = len(adj)
    edges = [(a, b) for (a, b) in edges_of(adj) if a != x and b != x]
    last = n - 1
    if x != last:
        edges = [tuple(x if w == last else w for w in e) for e in edges]
    return adj_from_edges(n - 1, edges)


def mutate(adj, rng, nmin, nmax):
    n = len(adj)
    edges = edges_of(adj)
    leaves = [v for v in range(n) if len(adj[v]) == 1]
    r = rng.random()
    if r < 0.30:  # leaf move
        x = rng.choice(leaves)
        edges = [e for e in edges if x not in e]
        tgt = rng.randrange(n - 1)
        if tgt >= x:
            tgt += 1
        edges.append((x, tgt))
        return adj_from_edges(n, edges)
    if r < 0.40 and n < nmax:  # add leaf
        edges.append((n, rng.randrange(n)))
        return adj_from_edges(n + 1, edges)
    if r < 0.50 and n > nmin:  # remove leaf
        return remove_vertex(adj, rng.choice(leaves))
    if r < 0.58 and n + 3 <= nmax:  # add cherry (vertex with t leaves)
        t = rng.choice([1, 2, 2, 3])
        if n + 1 + t > nmax:
            t = nmax - n - 1
        c = n
        edges.append((rng.randrange(n), c))
        for i in range(t):
            edges.append((c, n + 1 + i))
        return adj_from_edges(n + 1 + t, edges)
    if r < 0.64 and n > nmin + 3:  # remove a support vertex whose other nbrs are leaves
        cands = [v for v in range(n) if len(adj[v]) >= 2 and sum(1 for w in adj[v] if len(adj[w]) != 1) == 1]
        if cands:
            v = rng.choice(cands)
            kill = [w for w in adj[v] if len(adj[w]) == 1]
            if n - 1 - len(kill) >= nmin:
                keep = [w for w in range(n) if w != v and w not in kill]
                idx = {w: i for i, w in enumerate(keep)}
                e2 = [(idx[a], idx[b]) for (a, b) in edges if a in idx and b in idx]
                return adj_from_edges(len(keep), e2)
    if r < 0.72 and n < nmax:  # subdivide an edge
        a, b = rng.choice(edges)
        edges = [e for e in edges if e != (a, b)]
        edges += [(a, n), (n, b)]
        return adj_from_edges(n + 1, edges)
    if r < 0.80 and n > nmin:  # contract a degree-2 vertex
        d2 = [v for v in range(n) if len(adj[v]) == 2]
        if d2:
            x = rng.choice(d2)
            a, b = adj[x]
            edges = [e for e in edges if x not in e] + [(a, b)]
            last = n - 1
            if x != last:
                edges = [tuple(x if w == last else w for w in e) for e in edges]
            return adj_from_edges(n - 1, edges)
    # subtree prune and regraft
    a, b = rng.choice(edges)
    # side of b when edge removed
    side = {b}
    stack = [b]
    while stack:
        u = stack.pop()
        for w in adj[u]:
            if w not in side and not (u == b and w == a):
                side.add(w)
                stack.append(w)
    other = [v for v in range(n) if v not in side]
    root = rng.choice(sorted(side)) if rng.random() < 0.3 else b
    # re-root the moved subtree at `root`: edges inside side unchanged
    tgt = rng.choice(other)
    edges = [e for e in edges if e != (a, b) and e != (b, a)]
    edges.append((root, tgt))
    return adj_from_edges(n, edges)


def hub_seed(m, t):
    b = B()
    for _ in range(m):
        c = b.add(0)
        for _ in range(t):
            b.add(c)
    return b.adj()


def star_hubs_seed(h, m, t):
    b = B()
    for _ in range(h):
        w = b.add(0)
        for _ in range(m):
            c = b.add(w)
            for _ in range(t):
                b.add(c)
    return b.adj()


def random_tree(n, rng):
    return adj_from_edges(n, [(i, rng.randrange(i)) for i in range(1, n)])


def score(r, obj):
    s = r[obj]
    if r["viol"]:
        return 1e9
    if s >= 0:
        return 1e8
    return -math.log(-s)


def main():
    obj, nmin, nmax, secs, seed, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    seeds_file = sys.argv[7] if len(sys.argv) > 7 else None
    rng = random.Random(seed)
    ext = []
    if seeds_file:
        for line in open(seeds_file):
            d = json.loads(line)
            if nmin <= d["n"] <= nmax:
                ext.append(adj_from_edges(d["n"], [tuple(e) for e in d["edges"]]))
    t0 = time.time()
    evals = 0
    restarts = 0
    best_by_n = {}
    violations = []
    global_best = None
    while time.time() - t0 < secs:
        restarts += 1
        kind = restarts % 4
        if os.environ.get("SEEDONLY") and ext:
            kind = 1
        if kind == 1 and ext:
            adj = rng.choice(ext)
        elif kind == 2:
            t = rng.choice([1, 2, 2, 3, 4, 5])
            m = max(3, min((nmax - 1) // (t + 1), rng.randrange(max(3, (nmin - 1) // (t + 1)), (nmax - 1) // (t + 1) + 1)))
            adj = hub_seed(m, t)
        elif kind == 3:
            h = rng.choice([2, 3, 4])
            t = rng.choice([1, 2, 3])
            m = max(2, (rng.randrange(nmin, nmax + 1) // h - 1) // (t + 1))
            adj = star_hubs_seed(h, m, t)
        else:
            adj = random_tree(rng.randrange(nmin, nmax + 1), rng)
        while len(adj) < nmin:
            adj = adj_from_edges(len(adj) + 1, edges_of(adj) + [(len(adj), rng.randrange(len(adj)))])
        if len(adj) > nmax:
            adj = random_tree(rng.randrange(nmin, nmax + 1), rng)
        r = evaluate(adj)
        evals += 1
        cur = score(r, obj)
        steps = 2500
        T0, T1 = float(os.environ.get("T0", 0.08)), 0.001
        for step in range(steps):
            if time.time() - t0 >= secs:
                break
            T = T0 * (T1 / T0) ** (step / steps)
            try:
                adj2 = mutate(adj, rng, nmin, nmax)
            except (IndexError, ValueError):
                continue
            n2 = len(adj2)
            if n2 < nmin or n2 > nmax:
                continue
            r2 = evaluate(adj2)
            evals += 1
            s2 = score(r2, obj)
            if r2["viol"] and len(violations) < 20:
                violations.append({"n": n2, "edges": edges_of(adj2), "viol": r2["viol"][:10], "window": r2["window"], "alpha": r2["alpha"]})
            if s2 >= cur or rng.random() < math.exp((s2 - cur) / T):
                adj, cur, r = adj2, s2, r2
                n = len(adj)
                b = best_by_n.get(n)
                if b is None or cur > b["g"]:
                    best_by_n[n] = {"g": cur, "val": r[obj], "n": n, "edges": edges_of(adj), "window": r["window"], "alpha": r["alpha"],
                                    "at": r[obj + "_at"], "A": r["A"], "B": r["B"], "C": r["C"], "E": r["E"], "F": r["F"],
                                    "degseq": sorted((len(x) for x in adj), reverse=True)[:8]}
                if global_best is None or cur > global_best["g"]:
                    global_best = best_by_n[n]
        with open(out, "w") as f:
            json.dump({"obj": obj, "nmin": nmin, "nmax": nmax, "seed": seed, "elapsed": time.time() - t0,
                       "evaluations": evals, "restarts": restarts, "violations": violations,
                       "global_best": global_best, "best_by_n": [best_by_n[k] for k in sorted(best_by_n)]}, f)


if __name__ == "__main__":
    main()
