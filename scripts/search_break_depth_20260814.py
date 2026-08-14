#!/usr/bin/env python3
"""Break-depth maximization lane (2026-08-14).

Objective, motivated by the Rota flat-count disproof (arXiv:2608.07342, see
notes/matroid_nonunimodality_2026-08-14.md): push tree LC breaks DOWN toward
the Levit-Kadrawi prefix threshold thr = ceil((2*alpha-1)/3). Every known
break sits at dist = k - thr >= 4 (n=26/28 exhaustive data and every known
family). A break with dist < 0 would kill the prefix-LC route to #993.

Metrics per tree (all exact integer/Fraction arithmetic):
  dist      = min over LC breaks k of (k - thr); None if LC
  band_r    = max Turan ratio i_{k-1} i_{k+1} / i_k^2 over k in [thr-2, thr+2]
Any non-unimodal sequence is a counterexample to #993: dumped to an ALARM
file and printed loudly.

Subcommands:
  selftest            validate against stored exhaustive data + TG break sets
  sweep               structured family grid (TG, pattern grammar, KL, brooms)
  evolve --minutes M  multi-island evolutionary search, checkpointed
"""
import argparse
import json
import math
import os
import random
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)

import networkx as nx

from verify_partial_sync_obstruction_20260811 import (
    tree_ipoly, pattern_tree, path_tree, spider2, T1_ell, RTree)
import probe_absorption_margin_20260811 as fam

OUTDIR = os.path.join(ROOT, "results", "break_depth_lane_20260814")


# ---------- core analysis ----------

def thr_of(alpha):
    return (2 * alpha + 1) // 3          # == ceil((2*alpha - 1) / 3)


def analyze(poly):
    alpha = len(poly) - 1
    thr = thr_of(alpha)
    breaks = [k for k in range(1, alpha)
              if poly[k - 1] * poly[k + 1] > poly[k] * poly[k]]
    dist = min((k - thr) for k in breaks) if breaks else None
    return alpha, thr, breaks, dist


def is_unimodal(a):
    rising = True
    for i in range(1, len(a)):
        if rising:
            if a[i] < a[i - 1]:
                rising = False
        elif a[i] > a[i - 1]:
            return False
    return True


def ratio(poly, k):
    return Fraction(poly[k - 1] * poly[k + 1], poly[k] * poly[k])


def band_ratio(poly, lo, hi):
    alpha = len(poly) - 1
    ks = range(max(1, lo), min(alpha - 1, hi) + 1)
    return max((ratio(poly, k) for k in ks), default=Fraction(0))


def g6_of(edges):
    G = nx.Graph()
    G.add_edges_from(edges)
    G = nx.convert_node_labels_to_integers(G)
    return nx.to_graph6_bytes(G, header=False).decode().strip()


def alarm_if_nonunimodal(poly, edges, origin):
    if is_unimodal(poly):
        return False
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, f"ALARM_counterexample_{int(time.time())}.json")
    with open(path, "w") as fh:
        json.dump({"origin": origin, "graph6": g6_of(edges) if edges else None,
                   "edges": edges, "poly": [str(c) for c in poly]}, fh, indent=1)
    print("!" * 70)
    print("NON-UNIMODAL INDEPENDENCE SEQUENCE FOUND — #993 COUNTEREXAMPLE?")
    print("dumped to", path)
    print("!" * 70)
    return True


# ---------- tree representation: edge list on 0..n-1 ----------

def edges_to_adj(edges):
    n = max(max(e) for e in edges) + 1
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return adj


def rtree_edges(T):
    out = []
    for u, row in enumerate(T.adj):
        for v in row:
            if u < v:
                out.append((u, v))
    return out


def relabel(edges):
    verts = sorted({x for e in edges for x in e})
    m = {v: i for i, v in enumerate(verts)}
    return [(m[u], m[v]) for u, v in edges]


def check_tree(edges):
    adj = edges_to_adj(edges)
    n = len(adj)
    if len(edges) != n - 1:
        return False
    seen = [False] * n
    stack = [0]
    seen[0] = True
    c = 1
    while stack:
        v = stack.pop()
        for u in adj[v]:
            if not seen[u]:
                seen[u] = True
                c += 1
                stack.append(u)
    return c == n


def eval_edges(edges):
    adj = edges_to_adj(edges)
    return tree_ipoly(len(adj), adj, 0)


# ---------- seeds ----------

def kl_seed_edges():
    out = []
    d = json.load(open(os.path.join(ROOT, "results", "analysis_n26.json")))
    for f in d["lc_failures"]:
        G = nx.from_graph6_bytes(f["graph6"].encode())
        out.append(("kl26", [tuple(e) for e in G.edges()]))
    d = json.load(open(os.path.join(ROOT, "results",
                                    "analysis_n28_modal_lc_nm.json")))
    for f in d["top_lc_failures"]:
        G = nx.from_graph6_bytes(f["graph6"].encode())
        out.append(("kl28", [tuple(e) for e in G.edges()]))
    return out


def structured_seed_edges(ncap):
    seeds = []
    for (m, t) in [(2, 2), (2, 4), (4, 2), (4, 4), (8, 2), (8, 4), (12, 4)]:
        T = fam.TG_mt(m, t)
        if T.n <= ncap:
            seeds.append((f"TG_{m}_{t}", rtree_edges(T)))
    for (l, legs, k, m) in [(3, 2, 3, 2), (3, 4, 3, 3), (2, 2, 4, 2)]:
        T = pattern_tree(lambda: T1_ell(l), lambda: spider2(legs), [k] * m)
        if T.n <= ncap:
            seeds.append((f"pat_{l}_{legs}_{k}_{m}", rtree_edges(T)))
    # KL pattern grammar (H = P2): the sweep's whole dist frontier lives here
    for l in (2, 3, 4):
        for k in (3, 4, 5):
            for m in (2, 3, 4):
                T = pattern_tree(lambda: T1_ell(l), lambda: path_tree(2),
                                 [k] * m)
                if T.n <= ncap:
                    seeds.append((f"klpat_{l}_{k}_{m}", rtree_edges(T)))
    rng = random.Random(9)
    for i in range(4):
        n = min(ncap, 30 + 10 * i)
        e = [(j, rng.randrange(j)) for j in range(1, n)]  # random recursive tree
        seeds.append((f"rrt{i}", e))
    return seeds


# ---------- mutations ----------

def subtree_below(adj, u, v):
    """Vertices on v's side of edge (u, v)."""
    seen = {u, v}
    stack = [v]
    comp = [v]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                comp.append(y)
                stack.append(y)
    return comp


def mutate(edges, rng, ncap):
    for _ in range(10):
        e = list(edges)
        n = max(max(x) for x in e) + 1
        adj = edges_to_adj(e)
        deg = [len(a) for a in adj]
        op = rng.choices(
            ["leaf", "p2", "delleaf", "regraft", "subdiv", "smooth",
             "arm", "broom", "unit"],
            weights=[2, 3, 2, 4, 2, 2, 2, 2, 3])[0]
        if op == "leaf" and n + 1 <= ncap:
            e.append((rng.randrange(n), n))
        elif op == "p2" and n + 2 <= ncap:
            v = rng.randrange(n)
            e += [(v, n), (n, n + 1)]
        elif op == "delleaf" and n > 12:
            leaves = [v for v in range(n) if deg[v] == 1]
            v = rng.choice(leaves)
            e = [x for x in e if v not in x]
            e = relabel(e)
        elif op == "regraft" and n > 6:
            u, v = rng.choice(e)
            if rng.random() < 0.5:
                u, v = v, u
            comp = set(subtree_below(adj, u, v))
            rest = [w for w in range(n) if w not in comp]
            if len(rest) < 2:
                continue
            w = rng.choice([x for x in rest if x != u] or rest)
            e = [x for x in e if x != (u, v) and x != (v, u)]
            e.append((w, v))
        elif op == "subdiv" and n + 1 <= ncap:
            u, v = rng.choice(e)
            e = [x for x in e if x != (u, v)]
            e += [(u, n), (n, v)]
        elif op == "smooth" and n > 12:
            mids = [v for v in range(n) if deg[v] == 2]
            if not mids:
                continue
            v = rng.choice(mids)
            a, b = adj[v]
            e = [x for x in e if v not in x]
            e.append((a, b))
            e = relabel(e)
        elif op == "arm":
            t = rng.randrange(1, 5)
            if n + 1 + 2 * t > ncap:
                continue
            v = rng.randrange(n)
            hub = n
            e.append((v, hub))
            nx_id = n + 1
            for _ in range(t):
                e += [(hub, nx_id), (nx_id, nx_id + 1)]
                nx_id += 2
        elif op == "broom":
            j = rng.randrange(1, 4)
            if n + j > ncap:
                continue
            v = rng.randrange(n)
            for i in range(j):
                e.append((v, n + i))
        elif op == "unit":
            # graft a T_{1,l} pattern unit (the KL-grammar gadget) at a
            # random vertex, root-to-vertex
            l = rng.randrange(2, 5)
            U = T1_ell(l)
            if n + U.n > ncap:
                continue
            for u, row in enumerate(U.adj):
                for w in row:
                    if u < w:
                        e.append((n + u, n + w))
            e.append((rng.randrange(n), n + U.root))
        else:
            continue
        if check_tree(e):
            return e
    return list(edges)


# ---------- fitness ----------

BIG = 10 ** 9


def fitness(poly, cfg):
    alpha, thr, breaks, dist = analyze(poly)
    mode = max(range(len(poly)), key=lambda i: poly[i])
    if cfg["obj"] == "band":
        primary = band_ratio(poly, thr - 2, thr + 2)
        secondary = Fraction(-dist) if dist is not None else Fraction(-BIG)
        return (primary, secondary), dist
    # obj == 'deepen'
    if breaks:
        kbest = min(breaks)
        lo = mode + 1 if cfg["zone"] == "descending" else 1
        ks = [k for k in range(lo, kbest) if abs(k - mode) > 2]
        proxy = max((ratio(poly, k) for k in ks), default=Fraction(0))
        return (Fraction(-dist), proxy), dist
    head = band_ratio(poly, thr, alpha - 1)
    return (Fraction(-BIG), head), dist


# ---------- evolution ----------

def run_island(args):
    cfg, deadline, seed = args
    rng = random.Random(seed)
    os.makedirs(OUTDIR, exist_ok=True)
    ckpt = os.path.join(OUTDIR, f"island_{cfg['id']}.json")
    pop = []
    seen = set()
    evals = 0
    alarms = 0

    def add(edges, origin):
        nonlocal evals, alarms
        poly = eval_edges(edges)
        evals += 1
        key = hash(tuple(poly))
        if key in seen:
            return
        seen.add(key)
        if alarm_if_nonunimodal(poly, edges, f"island{cfg['id']}:{origin}"):
            alarms += 1
        fit, dist = fitness(poly, cfg)
        pop.append({"fit": fit, "dist": dist, "edges": edges,
                    "n": max(max(e) for e in edges) + 1})

    for name, e in kl_seed_edges() + structured_seed_edges(cfg["ncap"]):
        if max(max(x) for x in e) + 1 <= cfg["ncap"]:
            add(e, f"seed:{name}")
    pop.sort(key=lambda r: r["fit"], reverse=True)
    pop = pop[:cfg["mu"]]

    best_hist = []
    last_ckpt = 0.0
    gen = 0
    while time.monotonic() < deadline:
        gen += 1
        parents = [max(rng.sample(pop, min(3, len(pop))),
                       key=lambda r: r["fit"]) for _ in range(cfg["lam"])]
        for p in parents:
            child = mutate(p["edges"], rng, cfg["ncap"])
            add(child, "mut")
        pop.sort(key=lambda r: r["fit"], reverse=True)
        pop = pop[:cfg["mu"]]
        b = pop[0]
        best_hist.append({"gen": gen, "dist": b["dist"], "n": b["n"],
                          "fit0": float(b["fit"][0]), "fit1": float(b["fit"][1])})
        if time.monotonic() - last_ckpt > 60:
            last_ckpt = time.monotonic()
            snapshot(cfg, ckpt, pop, evals, gen, alarms, best_hist, False)
    snapshot(cfg, ckpt, pop, evals, gen, alarms, best_hist, True)
    return ckpt


def snapshot(cfg, path, pop, evals, gen, alarms, hist, final):
    top = []
    for r in pop[:5]:
        poly = eval_edges(r["edges"])
        alpha, thr, breaks, dist = analyze(poly)
        top.append({
            "n": r["n"], "alpha": alpha, "thr": thr, "dist": dist,
            "breaks": breaks,
            "band_r_at_thr": float(band_ratio(poly, thr - 2, thr + 2)),
            "graph6": g6_of(r["edges"]),
        })
    with open(path, "w") as fh:
        json.dump({"cfg": cfg, "final": final, "evals": evals, "gens": gen,
                   "alarms": alarms, "top": top,
                   "best_hist": hist[-50:]}, fh, indent=1)


# ---------- structured sweep ----------

def sweep():
    rows = []

    def rec(label, n, poly, edges=None):
        alpha, thr, breaks, dist = analyze(poly)
        alarm_if_nonunimodal(poly, edges, f"sweep:{label}")
        rows.append({
            "label": label, "n": n, "alpha": alpha, "thr": thr,
            "breaks": breaks, "dist": dist,
            "max_ratio": float(max((ratio(poly, k) for k in breaks),
                                   default=Fraction(0))),
            "band_r_at_thr": float(band_ratio(poly, thr - 2, thr + 2)),
        })

    t0 = time.time()
    for m in range(1, 21):
        for t in range(1, 9):
            T = fam.TG_mt(m, t)
            if T.n > 1200:
                continue
            rec(f"TG_{m}_{t}", T.n, T.ipoly())
    print(f"TG grid done ({time.time()-t0:.0f}s, {len(rows)} rows)")

    for l in (1, 2, 3, 4):
        for legs in (1, 2, 3, 4):
            for k in (2, 3, 4, 5, 6):
                for m in (1, 2, 3, 4, 6, 8):
                    T = pattern_tree(lambda: T1_ell(l),
                                     lambda: spider2(legs), [k] * m)
                    if T.n > 900:
                        continue
                    rec(f"pat_l{l}_s{legs}_k{k}_m{m}", T.n, T.ipoly())
    print(f"pattern grid done ({time.time()-t0:.0f}s, {len(rows)} rows)")

    for l in (2, 3, 4):
        for k in (2, 3, 4, 5):
            for m in (2, 3, 4, 6):
                T = pattern_tree(lambda: T1_ell(l), lambda: path_tree(2),
                                 [k] * m)
                if T.n > 900:
                    continue
                rec(f"klpat_l{l}_k{k}_m{m}", T.n, T.ipoly())

    import itertools
    for l in (2, 3, 4):
        kss = [list(t) for t in itertools.product((2, 3, 4, 5, 6), repeat=2)]
        kss += [list(t) for t in itertools.product((2, 3, 4, 5), repeat=3)]
        for ks in kss:
            T = pattern_tree(lambda: T1_ell(l), lambda: path_tree(2), ks)
            if T.n > 70:
                continue
            rec(f"klmix_l{l}_{'-'.join(map(str, ks))}", T.n, T.ipoly())
    print(f"mixed-ks grid done ({time.time()-t0:.0f}s, {len(rows)} rows)")

    for p in (3, 6, 9):
        for q in (3, 6, 9):
            for t in (5, 15, 30):
                T = RTree()
                a = T.new_vertex()
                prev = a
                for _ in range(t):
                    v = T.new_vertex()
                    T.edge(prev, v)
                    prev = v
                for _ in range(p):
                    T.edge(a, T.new_vertex())
                for _ in range(q):
                    T.edge(prev, T.new_vertex())
                rec(f"broom_{p}_{q}_{t}", T.n, T.ipoly())

    with_breaks = [r for r in rows if r["dist"] is not None]
    with_breaks.sort(key=lambda r: r["dist"])
    out = os.path.join(ROOT, "results", "break_depth_sweep_20260814.json")
    with open(out, "w") as fh:
        json.dump({"rows": rows, "n_rows": len(rows),
                   "n_with_breaks": len(with_breaks)}, fh, indent=1)
    print(f"\nsweep: {len(rows)} trees, {len(with_breaks)} with breaks")
    print("deepest 10 by dist = k - thr:")
    for r in with_breaks[:10]:
        print(f"  {r['label']:<22} n={r['n']:<5} alpha={r['alpha']:<5} "
              f"thr={r['thr']:<5} dist={r['dist']:<4} "
              f"band_r@thr={r['band_r_at_thr']:.6f}")
    tight = sorted(rows, key=lambda r: -r["band_r_at_thr"])[:5]
    print("tightest Turan ratio at the threshold band:")
    for r in tight:
        print(f"  {r['label']:<22} n={r['n']:<5} band_r@thr="
              f"{r['band_r_at_thr']:.6f} dist={r['dist']}")
    print("wrote", out)


# ---------- deterministic mutation-ball audits ----------

def det_neighbors(edges):
    """Every single application of: pendant leaf, pendant P2, edge
    subdivision, leaf deletion, degree-2 smoothing."""
    n = max(max(e) for e in edges) + 1
    adj = edges_to_adj(edges)
    deg = [len(x) for x in adj]
    out = []
    for v in range(n):
        out.append(edges + [(v, n)])
        out.append(edges + [(v, n), (n, n + 1)])
    for (u, v) in edges:
        out.append([e for e in edges if e != (u, v)] + [(u, n), (n, v)])
    for v in range(n):
        if deg[v] == 1:
            out.append(relabel([e for e in edges if v not in e]))
        if deg[v] == 2:
            a_, b_ = adj[v]
            out.append(relabel([e for e in edges if v not in e] + [(a_, b_)]))
    return out


def balls():
    """1-mutation ball around the frontier trees (KL + klpat grammar) and
    the complete 2-mutation ball around all 21 known LC-failure trees."""
    kl = [e for _, e in kl_seed_edges()]
    klpat = []
    for l in (2, 3, 4):
        for ks in ([4, 4, 4], [3, 4, 4], [4, 5], [4, 6], [5, 4], [4, 4]):
            T = pattern_tree(lambda: T1_ell(l), lambda: path_tree(2), ks)
            if T.n <= 60:
                klpat.append((f"klpat_l{l}_{ks}", rtree_edges(T)))

    def scan(seed_sets, radius):
        seen = set()
        stats = {"distinct_polys": 0, "min_dist": None, "dist_le3_hits": [],
                 "depth2_at_alpha_le18": None}

        def visit(e):
            p = eval_edges(e)
            key = hash(tuple(p))
            if key in seen:
                return
            seen.add(key)
            stats["distinct_polys"] += 1
            alarm_if_nonunimodal(p, e, "balls")
            alpha, thr, breaks, dist = analyze(p)
            if dist is not None:
                if stats["min_dist"] is None or dist < stats["min_dist"]:
                    stats["min_dist"] = dist
                if dist <= 3:
                    stats["dist_le3_hits"].append(
                        {"dist": dist, "alpha": alpha, "graph6": g6_of(e)})
                if alpha <= 18 and max(alpha - k for k in breaks) >= 2:
                    stats["depth2_at_alpha_le18"] = {"alpha": alpha,
                                                     "graph6": g6_of(e)}

        for e0 in seed_sets:
            for e1 in det_neighbors(e0):
                visit(e1)
                if radius >= 2:
                    for e2 in det_neighbors(e1):
                        visit(e2)
        return stats

    r1 = scan([e for _, e in klpat] + kl[:6], 1)
    r2 = scan(kl, 2)
    out = os.path.join(ROOT, "results", "break_depth_balls_20260814.json")
    with open(out, "w") as fh:
        json.dump({"one_ball_frontier_seeds": r1,
                   "two_ball_all_known_lc_failures": r2}, fh, indent=1)
    print("1-ball (frontier seeds):", r1["distinct_polys"], "polys, min dist",
          r1["min_dist"], "| hits:", len(r1["dist_le3_hits"]))
    print("2-ball (all 21 known LC failures):", r2["distinct_polys"],
          "polys, min dist", r2["min_dist"],
          "| hits:", len(r2["dist_le3_hits"]),
          "| depth>=2 at alpha<=18:", r2["depth2_at_alpha_le18"])
    print("wrote", out)


# ---------- selftest ----------

def selftest():
    d = json.load(open(os.path.join(ROOT, "results", "analysis_n26.json")))
    for f in d["lc_failures"]:
        G = nx.from_graph6_bytes(f["graph6"].encode())
        edges = [tuple(e) for e in G.edges()]
        poly = eval_edges(relabel(edges))
        assert poly == f["poly"], "n26 poly mismatch"
        alpha, thr, breaks, dist = analyze(poly)
        assert f["lc_pos"] in breaks and dist == 4
    F = fam.TG_mt(4, 8).ipoly()
    alpha, thr, breaks, dist = analyze(F)
    assert breaks == [102, 104, 106, 108] and thr == 73 and dist == 29, breaks
    assert not is_unimodal([1, 3, 2, 3, 1])
    assert is_unimodal([1, 3, 3, 2, 1])
    rng = random.Random(1)
    e = structured_seed_edges(96)[0][1]
    for _ in range(500):
        e = mutate(e, rng, 96)
        assert check_tree(e)
    cfgA = {"obj": "deepen", "zone": "full"}
    fitA, distA = fitness(F, cfgA)
    assert distA == 29 and fitA[0] == Fraction(-29)
    cfgB = {"obj": "band"}
    fitB, _ = fitness(F, cfgB)
    assert Fraction(0) < fitB[0] < Fraction(1)
    print("selftest: all passed")


# ---------- main ----------

def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    sub.add_parser("sweep")
    sub.add_parser("balls")
    ev = sub.add_parser("evolve")
    ev.add_argument("--minutes", type=float, default=30.0)
    ev.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()

    if args.cmd == "selftest":
        selftest()
    elif args.cmd == "sweep":
        sweep()
    elif args.cmd == "balls":
        balls()
    else:
        os.makedirs(OUTDIR, exist_ok=True)
        deadline = time.monotonic() + 60 * args.minutes
        configs = [
            {"id": 0, "obj": "deepen", "zone": "full", "ncap": 48,
             "mu": 14, "lam": 42},
            {"id": 1, "obj": "deepen", "zone": "descending", "ncap": 96,
             "mu": 14, "lam": 42},
            {"id": 2, "obj": "deepen", "zone": "full", "ncap": 96,
             "mu": 14, "lam": 42},
            {"id": 3, "obj": "deepen", "zone": "descending", "ncap": 200,
             "mu": 14, "lam": 42},
            {"id": 4, "obj": "band", "ncap": 96, "mu": 14, "lam": 42},
            {"id": 5, "obj": "band", "ncap": 200, "mu": 14, "lam": 42},
        ]
        tasks = [(c, deadline, 1000 + c["id"]) for c in configs]
        with Pool(min(args.workers, len(tasks))) as pool:
            paths = pool.map(run_island, tasks)
        merged = []
        for p in paths:
            with open(p) as fh:
                merged.append(json.load(fh))
        best_dist = min((t["dist"] for isl in merged for t in isl["top"]
                         if t["dist"] is not None), default=None)
        out = os.path.join(ROOT, "results", "break_depth_search_20260814.json")
        with open(out, "w") as fh:
            json.dump({"budget_minutes": args.minutes,
                       "total_evals": sum(i["evals"] for i in merged),
                       "total_alarms": sum(i["alarms"] for i in merged),
                       "best_dist_overall": best_dist,
                       "islands": merged}, fh, indent=1)
        print("evolve done:", sum(i["evals"] for i in merged), "evals,",
              "best dist overall:", best_dist, "| alarms:",
              sum(i["alarms"] for i in merged))
        print("wrote", out)


if __name__ == "__main__":
    main()
