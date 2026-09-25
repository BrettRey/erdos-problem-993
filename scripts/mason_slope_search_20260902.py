"""Adversarial hill-climb: maximize the descending-prefix slope

    obj(T) = max { mu_{k+1} - mu_k : k <= thr-2, mu_k < k+1 }

(a valley needs obj > 1; S1-desc says obj <= 1).  Also tracks the plain prefix
slope max { mu_{k+1}-mu_k : k <= thr-2 } (N-type; stars give -2/(s+1)).
Exact integer DP for the polynomial; Fractions for the objective.
Trees as parent arrays over vertices 0..n-1 (root 0).
"""
import random, sys, json, time
from fractions import Fraction

def ipoly_parent(par):
    n = len(par)
    children = [[] for _ in range(n)]
    order = []
    for v in range(1, n):
        children[par[v]].append(v)
    # BFS order from root 0
    stack = [0]
    while stack:
        u = stack.pop(); order.append(u); stack.extend(children[u])
    P = [None] * n; R = [None] * n
    for u in reversed(order):
        p = [0, 1]; r = [1]
        for c in children[u]:
            p = pmul(p, R[c]); r = pmul(r, padd(P[c], R[c]))
        P[u] = p; R[u] = r
    return padd(P[0], R[0])

def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out

def padd(a, b):
    m = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(m)]

def objective(poly):
    alpha = len(poly) - 1
    thr = -(-(2 * alpha - 1) // 3)
    mu = [Fraction((k + 1) * poly[k + 1], poly[k]) for k in range(alpha)]
    best_desc = None; best_pref = None
    for k in range(min(thr - 2, alpha - 2) + 1):
        s = mu[k + 1] - mu[k]
        if best_pref is None or s > best_pref: best_pref = s
        if mu[k] < k + 1 and (best_desc is None or s > best_desc): best_desc = s
    return best_desc, best_pref

def is_unimodal(a):
    i = 0
    while i + 1 < len(a) and a[i] <= a[i + 1]: i += 1
    while i + 1 < len(a) and a[i] >= a[i + 1]: i += 1
    return i + 1 >= len(a)

# ---- tree edits on parent arrays (keep root 0) ----
def canon(par):
    """re-root/relabel not needed for correctness; just keep as is."""
    return par

def subtree(par, v):
    n = len(par); ch = [[] for _ in range(n)]
    for u in range(1, n): ch[par[u]].append(u)
    out = []; st = [v]
    while st:
        u = st.pop(); out.append(u); st.extend(ch[u])
    return out

def mutate(par, rng, ncap):
    n = len(par)
    par = list(par)
    r = rng.random()
    if r < 0.25 and n < ncap:            # add leaf
        par.append(rng.randrange(n)); return par
    if r < 0.45 and n > 4:               # remove a leaf (non-root)
        ch = [0] * n
        for u in range(1, n): ch[par[u]] += 1
        leaves = [v for v in range(1, n) if ch[v] == 0]
        v = rng.choice(leaves)
        # relabel: drop v, shift labels > v
        new = []
        for u in range(n):
            if u == v: continue
            p = par[u]
            new.append(p - 1 if p > v else p)
        return new
    if r < 0.75:                         # regraft: move a non-root vertex's subtree under another vertex
        v = rng.randrange(1, n)
        sub = set(subtree(par, v))
        cand = [u for u in range(n) if u not in sub]
        par[v] = rng.choice(cand); return par
    if r < 0.9 and n < ncap:             # subdivide edge (v, par[v])
        v = rng.randrange(1, n)
        par.append(par[v]); par[v] = n; return par
    # add a pendant P2 at a random vertex
    if n + 2 <= ncap:
        u = rng.randrange(n); par.append(u); par.append(n); return par
    return par

def seeds_from_graph6(g6list):
    import networkx as nx
    out = []
    for g in g6list:
        T = nx.from_graph6_bytes(g.encode())
        # parent array rooted at 0 with BFS relabel
        order = list(nx.bfs_tree(T, 0).nodes)
        idx = {v: i for i, v in enumerate(order)}
        par = [0] * len(order)
        for v in order:
            if v != 0:
                par[idx[v]] = idx[next(iter(nx.bfs_tree(T, 0).predecessors(v)))]
        out.append(par)
    return out

def star(s): return [0] + [0] * s
def spider2(t):
    par = [0]
    for _ in range(t):
        a = len(par); par.append(0); par.append(a)
    return par
def kl(a, m, n):
    par = [0]
    for arms in (a, m, n):
        h = len(par); par.append(0)
        for _ in range(arms):
            x = len(par); par.append(h); par.append(x)
    return par

def run(seed, minutes, ncap, mode, seeds, log):
    rng = random.Random(seed)
    pop = []
    for s in seeds:
        p = ipoly_parent(s); d, pr = objective(p)
        pop.append(((d if mode == 'desc' else pr), s))
    pop = [x for x in pop if x[0] is not None]
    pop.sort(key=lambda x: x[0], reverse=True)
    pop = pop[:40]
    t0 = time.time(); evals = 0; alarms = 0; best_hist = []
    while time.time() - t0 < 60 * minutes:
        base = rng.choice(pop[:20])[1]
        child = base
        for _ in range(rng.randint(1, 3)):
            child = mutate(child, rng, ncap)
        p = ipoly_parent(child); evals += 1
        if not is_unimodal(p):
            alarms += 1
            with open(log + '.ALARM', 'a') as f: f.write(json.dumps(dict(par=child, poly=p)) + '\n')
        d, pr = objective(p)
        val = d if mode == 'desc' else pr
        if val is None: continue
        if val > pop[-1][0] or len(pop) < 60:
            pop.append((val, child)); pop.sort(key=lambda x: x[0], reverse=True); pop = pop[:60]
            if val >= pop[0][0]:
                best_hist.append((round(time.time() - t0), evals, float(val), len(child)))
    best = pop[0]
    p = ipoly_parent(best[1])
    res = dict(seed=seed, mode=mode, ncap=ncap, evals=evals, alarms=alarms, best=float(best[0]), best_n=len(best[1]),
               best_par=best[1], best_alpha=len(p) - 1, best_hist=best_hist[-10:])
    with open(log, 'w') as f: json.dump(res, f)
    return res

if __name__ == '__main__':
    seed = int(sys.argv[1]); minutes = float(sys.argv[2]); ncap = int(sys.argv[3]); mode = sys.argv[4]; log = sys.argv[5]
    seeds = [star(s) for s in (5, 8, 12, 20)] + [spider2(t) for t in (3, 5, 8, 12)] + [kl(3, 4, 4), kl(3, 6, 6), kl(4, 4, 4), kl(2, 3, 5)]
    # random trees
    rng = random.Random(seed)
    for _ in range(20):
        n = rng.randint(10, ncap); seeds.append([0] + [rng.randrange(i) for i in range(1, n)])
    res = run(seed, minutes, ncap, mode, seeds, log)
    print(json.dumps({k: v for k, v in res.items() if k != 'best_par'}))
