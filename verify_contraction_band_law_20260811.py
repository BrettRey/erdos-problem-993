#!/usr/bin/env python3
"""Contraction band law: exhaustive and adversarial evidence.

For a tree T and edge e, let b = coeffs I(T), c = coeffs I(T/e),
N = deg I(T), and define (out-of-range coefficients read 0)

    U_k = c_k b_k - c_{k-1} b_{k+1}
    V_k = c_k b_k - c_{k+1} b_{k-1}
    W_k = 2 c_{k-1} b_k - c_k b_{k-1} - c_{k-2} b_{k+1}

Two universal conjectures about these forms were stated and REFUTED on
2026-08-11 (see notes/contraction_band_law_2026-08-11.md): (1) V, W >= 0
at every index for every tree edge; (2) all failures confined to
reflected depth N - k <= 2. This script regenerates the full evidence,
including the killing sweep, exactly (integer arithmetic):
 1. exhaustive: all non-isomorphic trees n = 3..16 (geng), all edges;
 2. random trees to n = 18, all edges;
 3. adversarial non-LC composites: Kadrawi-Levit trees KL(3,k,k) for
    k = 4, 5, 6 (log-concavity already broken) with a pendant path of
    length 1..3 attached at every vertex, all edges of each composite
    (kills conjecture 1: V fails at depth 1, W at depths 0-1);
 4. multi-break composites: pattern-family trees with 1, 2, 3
    consecutive LC breaks (orders 177-315) plus one pendant leaf at
    eight spread vertices, all edges (kills conjecture 2: U-failures
    to depth 23, V to 5, W to 4).
Records per-form violation depth histograms per sweep. What survives is
the C2-restricted clean laws (path pieces) and the exhaustive small-tree
fact; see the note.

Rooted at the contracted edge, B and C always have the two-hub shape
B = A_uA_v + x(D_uA_v + A_uD_v), C = A_uA_v + xD_uD_v with
A_w = I(T_w - w), D_w = I(T_w - N[w]); the C2 base-facts pair of
notes/c2_bounded_pendant_core_2026-08-08.md is the special case where
all pieces are products of path polynomials.  See
notes/contraction_band_law_2026-08-11.md.
"""
import json
import random
import subprocess
import time
from collections import Counter

GENG = "/opt/homebrew/bin/geng"


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                out[i + j] += ai * bj
    return out


def padd(a, b):
    m = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
            for i in range(m)]


def ipoly(n, adj, root=0):
    parent = [-1] * n
    order = []
    seen = [False] * n
    st = [root]
    seen[root] = True
    while st:
        v = st.pop()
        order.append(v)
        for u in adj[v]:
            if not seen[u]:
                seen[u] = True
                parent[u] = v
                st.append(u)
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


def parse_graph6(line):
    bs = [ord(ch) - 63 for ch in line.strip()]
    n = bs[0]
    bits = []
    for x in bs[1:]:
        for i in range(5, -1, -1):
            bits.append((x >> i) & 1)
    edges = []
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if bits[idx]:
                edges.append((i, j))
            idx += 1
    return n, edges


def contract(n, edges, e):
    u, v = e
    m = {}
    idx = 0
    for w in range(n):
        if w == v:
            continue
        m[w] = idx
        idx += 1
    m[v] = m[u]
    ne = set()
    for a, b in edges:
        x, y = m[a], m[b]
        if x != y:
            ne.add((min(x, y), max(x, y)))
    adj = [[] for _ in range(n - 1)]
    for a, b in ne:
        adj[a].append(b)
        adj[b].append(a)
    return adj


def cf(p, i):
    return p[i] if 0 <= i < len(p) else 0


def tally(B, C, hists):
    N = len(B) - 1
    for k in range(N + 2):
        if cf(C, k) * cf(B, k) - cf(C, k + 1) * cf(B, k - 1) < 0:
            hists["V"][N - k] += 1
        if (2 * cf(C, k - 1) * cf(B, k) - cf(C, k) * cf(B, k - 1)
                - cf(C, k - 2) * cf(B, k + 1)) < 0:
            hists["W"][N - k] += 1
        if cf(C, k) * cf(B, k) - cf(C, k - 1) * cf(B, k + 1) < 0:
            hists["U"][N - k] += 1


def adj_of(n, edges):
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    return adj


def all_edge_pairs(n, edges, hists):
    B = ipoly(n, adj_of(n, edges))
    cnt = 0
    for e in edges:
        C = ipoly(n - 1, contract(n, edges, e))
        tally(B, C, hists)
        cnt += 1
    return cnt


def kl_tree(k):
    """(T_{1,3}:P_2)_k^(2): root v - hub - 3 cherries; two units of hub+k P2s."""
    adj = [[] for _ in range(1)]
    def new():
        adj.append([])
        return len(adj) - 1
    def edge(a, b):
        adj[a].append(b)
        adj[b].append(a)
    hub = new()
    edge(0, hub)
    for _ in range(3):
        w = new(); u2 = new()
        edge(hub, w); edge(w, u2)
    for _ in range(2):
        uh = new()
        edge(0, uh)
        for _ in range(k):
            a = new(); b = new()
            edge(uh, a); edge(a, b)
    return len(adj), [(a, b) for a in range(len(adj)) for b in adj[a] if a < b]


def main():
    t0 = time.time()
    report = {"date": "2026-08-11",
              "status": "both universal conjectures refuted; see note"}

    # 1. exhaustive n <= 16
    hists = {f: Counter() for f in "UVW"}
    ntrees = npairs = 0
    for n in range(3, 17):
        out = subprocess.run([GENG, "-q", "-c", str(n), f"{n-1}:{n-1}"],
                             capture_output=True, text=True).stdout
        for line in out.splitlines():
            nn, edges = parse_graph6(line)
            ntrees += 1
            npairs += all_edge_pairs(nn, edges, hists)
    report["exhaustive_n_le_16"] = {
        "trees": ntrees, "pairs": npairs,
        "hist": {f: dict(hists[f]) for f in "UVW"}}

    # 2. random trees to n = 18
    random.seed(993)
    import heapq
    rh = {f: Counter() for f in "UVW"}
    rpairs = 0
    for _ in range(400):
        n = random.randint(4, 18)
        if n <= 2:
            continue
        pr = [random.randrange(n) for _ in range(n - 2)]
        deg = [1] * n
        for x in pr:
            deg[x] += 1
        leaves = [i for i in range(n) if deg[i] == 1]
        heapq.heapify(leaves)
        edges = []
        for x in pr:
            leaf = heapq.heappop(leaves)
            edges.append((leaf, x))
            deg[x] -= 1
            if deg[x] == 1:
                heapq.heappush(leaves, x)
        edges.append((heapq.heappop(leaves), heapq.heappop(leaves)))
        rpairs += all_edge_pairs(n, edges, rh)
    report["random_n_le_18"] = {
        "pairs": rpairs, "hist": {f: dict(rh[f]) for f in "UVW"}}

    # 3. non-LC composites
    ch = {f: Counter() for f in "UVW"}
    cpairs = 0
    for k in (4, 5, 6):
        n0, edges0 = kl_tree(k)
        for glen in (1, 2, 3):
            for w in range(n0):
                n = n0 + glen
                edges = list(edges0) + [(w, n0)]
                for i in range(glen - 1):
                    edges.append((n0 + i, n0 + i + 1))
                cpairs += all_edge_pairs(n, edges, ch)
    report["nonlc_composites"] = {
        "pairs": cpairs, "hist": {f: dict(ch[f]) for f in "UVW"}}

    # 4. multi-break composites (kills the absolute band law)
    def spider2(nlegs):
        adj = [[]]
        def new():
            adj.append([])
            return len(adj) - 1
        for _ in range(nlegs):
            a = new(); b = new()
            adj[0].append(a); adj[a].append(0)
            adj[a].append(b); adj[b].append(a)
        return adj

    def t1_ell(ell):
        adj = [[], []]
        adj[0].append(1); adj[1].append(0)
        def new():
            adj.append([])
            return len(adj) - 1
        for _ in range(ell):
            w = new(); u2 = new()
            adj[1].append(w); adj[w].append(1)
            adj[w].append(u2); adj[u2].append(w)
        return adj

    def pattern(ell, nlegs, k, mm):
        adj = t1_ell(ell)
        def new():
            adj.append([])
            return len(adj) - 1
        for _ in range(mm):
            hub = new()
            adj[0].append(hub); adj[hub].append(0)
            for _ in range(k):
                sp = spider2(nlegs)
                base = len(adj)
                for row in sp:
                    adj.append([x + base for x in row])
                adj[hub].append(base); adj[base].append(hub)
        return adj

    mh = {f: Counter() for f in "UVW"}
    mpairs = 0
    for ell, nlegs, mm in ((2, 4, 9), (3, 4, 16), (7, 5, 13)):
        adj0 = pattern(ell, nlegs, 2, mm)
        n0 = len(adj0)
        for w in sorted(set([0, 1, 2, n0 // 4, n0 // 2,
                             3 * n0 // 4, n0 - 2, n0 - 1])):
            adj = [list(r) for r in adj0] + [[]]
            n = n0 + 1
            adj[w].append(n0); adj[n0].append(w)
            edges = [(a, b) for a in range(n) for b in adj[a] if a < b]
            mpairs += all_edge_pairs(n, edges, mh)
    report["multibreak_composites"] = {
        "pairs": mpairs, "hist": {f: dict(mh[f]) for f in "UVW"}}

    deepest = max([d for h in (hists, rh, ch, mh) for f in "UVW"
                   for d in h[f]] or [-1])
    report["deepest_violation_depth_any_form"] = deepest
    report["conjecture1_clean_VW_refuted"] = bool(
        ch["V"] or ch["W"] or mh["V"] or mh["W"])
    report["conjecture2_absolute_band_refuted"] = deepest > 2
    report["runtime_s"] = round(time.time() - t0, 2)

    with open("results/contraction_band_law_20260811.json", "w") as fh:
        json.dump(report, fh, indent=1)
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
