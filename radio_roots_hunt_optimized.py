#!/usr/bin/env python3
"""
Optimized Radio Roots Hunt: Search for Erdos #993 counterexample.

Optimization:
  - PRE-FILTER: Only compute roots if Near-Miss Ratio (NM) > 0.9.
  - This avoids O(N^3) root finding for 99% of candidates.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import pickle
import random
import sys
import time
from collections import deque
from typing import List, Tuple, Optional, Dict, Any

# Numpy is REQUIRED
import numpy as np

try:
    from indpoly import (
        is_unimodal, 
        poly_add, 
        poly_mul,
        first_unimodality_break,
        near_miss_ratio
    )
except ImportError:
    sys.path.append(os.getcwd())
    from indpoly import (
        is_unimodal, 
        _polyadd as poly_add, 
        _polymul as poly_mul, 
        is_unimodal as _is_unimodal,
        near_miss_ratio
    )
    def first_unimodality_break(seq):
        down = False
        for i in range(len(seq)-1):
            if seq[i+1] < seq[i]: down = True
            if down and seq[i+1] > seq[i]: return {"down_at": i-1, "up_at": i}
        return None

# -------------------------
# SCORING (OPTIMIZED)
# -------------------------

def radio_score(seq: List[int]) -> Tuple[float, float, float, float]:
    """
    Computes the 'Radio Roots' score.
    Returns tuple for sorting: (Composite_Score, NM_Ratio, Max_Re, -Max_Im).
    """
    if len(seq) < 2:
        return (-999.0, 0.0, -999.0, 0.0)
    
    # 1. Near Miss Ratio (Primary Signal) - CHEAP
    nm, _ = near_miss_ratio(seq)
    
    # OPTIMIZATION: If NM is low, don't bother with expensive roots
    if nm < 0.9:
        # Return score based only on NM, with dummy root values
        # We penalize this so high-NM candidates always win
        return (nm * 20.0 - 1000.0, nm, -999.0, 0.0)

    # 2. Roots Analysis - EXPENSIVE
    try:
        coeffs = np.array(seq[::-1], dtype=float)
        max_coeff = np.max(np.abs(coeffs))
        if max_coeff > 0:
            coeffs /= max_coeff
        roots = np.roots(coeffs)
    except Exception:
        return (-999.0, nm, -999.0, 0.0)

    if len(roots) == 0:
        return (-999.0, nm, -999.0, 0.0)

    max_real = float(np.max(roots.real))
    max_imag = float(np.max(np.abs(roots.imag)))
    
    if max_real > 2.0:
        return (-9999.0, nm, max_real, -max_imag)
    
    score = (nm * 20.0) + (max_real * 5.0) - (max_imag * 0.5)
    return (score, nm, max_real, -max_imag)


# -------------------------
# Graph / Tree Helpers (Reused)
# -------------------------

def check_tree(n: int, edges: List[Tuple[int, int]]) -> bool:
    if n == 0: return True
    if len(edges) != n - 1: return False
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    seen = [False] * n
    q = deque([0])
    seen[0] = True
    while q:
        x = q.popleft()
        for y in adj[x]:
            if not seen[y]:
                seen[y] = True
                q.append(y)
    return all(seen)

def random_tree_edges(n: int, rng: random.Random) -> List[Tuple[int, int]]:
    if n <= 1: return []
    prufer = [rng.randrange(n) for _ in range(n - 2)]
    degree = [1] * n
    for x in prufer: degree[x] += 1
    import heapq
    leaves = [i for i, d in enumerate(degree) if d == 1]
    heapq.heapify(leaves)
    edges = []
    for x in prufer:
        leaf = heapq.heappop(leaves)
        edges.append((leaf, x))
        degree[leaf] -= 1
        degree[x] -= 1
        if degree[x] == 1: heapq.heappush(leaves, x)
    edges.append((heapq.heappop(leaves), heapq.heappop(leaves)))
    return edges

def path_edges(n: int) -> List[Tuple[int, int]]:
    return [(i, i + 1) for i in range(n - 1)]

def star_edges(n: int) -> List[Tuple[int, int]]:
    return [(0, i) for i in range(1, n)]

def broom_edges(handle: int, leaves: int) -> Tuple[int, List[Tuple[int, int]]]:
    n = handle + 1 + leaves
    edges = [(i, i + 1) for i in range(handle)]
    end = handle
    for j in range(leaves):
        edges.append((end, handle + 1 + j))
    return n, edges

def caterpillar_edges(spine_len: int, leaf_counts: List[int]) -> Tuple[int, List[Tuple[int, int]]]:
    n = spine_len + sum(leaf_counts)
    edges = [(i, i + 1) for i in range(spine_len - 1)]
    v = spine_len
    for i, c in enumerate(leaf_counts):
        for _ in range(c):
            edges.append((i, v))
            v += 1
    return n, edges

def spherically_symmetric_tree(branching: List[int]) -> Tuple[int, List[Tuple[int, int]]]:
    edges = []
    current_level = [0]
    next_id = 1
    for b in branching:
        new_level = []
        for v in current_level:
            for _ in range(b):
                child = next_id
                next_id += 1
                edges.append((v, child))
                new_level.append(child)
        current_level = new_level
    return next_id, edges

def tree_T_de(d: int, e: int) -> Tuple[int, List[Tuple[int, int]]]:
    return spherically_symmetric_tree([d, e, 1])

def choose_attach_vertex(n: int, edges: List[Tuple[int, int]]) -> int:
    deg = [0] * n
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1
    maxdeg = max(deg) if n else 0
    for i, d in enumerate(deg):
        if d == maxdeg: return i
    return 0

def independence_fg_tree(n: int, edges: List[Tuple[int, int]], root: int) -> Tuple[List[int], List[int]]:
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    parent = [-1] * n
    parent[root] = root
    order = [root]
    for v in order:
        for w in adj[v]:
            if parent[w] == -1:
                parent[w] = v
                order.append(w)
    post = order[::-1]
    x_poly = [0, 1]
    f = [None] * n
    g = [None] * n
    for v in post:
        children = [w for w in adj[v] if parent[w] == v]
        if not children:
            f[v] = [1]
            g[v] = x_poly
        else:
            prod_h = [1]
            prod_f = [1]
            for w in children:
                h_w = poly_add(f[w], g[w])
                prod_h = poly_mul(prod_h, h_w)
                prod_f = poly_mul(prod_f, f[w])
            f[v] = prod_h
            g[v] = poly_mul(x_poly, prod_f)
    return f[root], g[root]

@dataclasses.dataclass
class TreeEntry:
    name: str
    n: int
    edges: List[Tuple[int, int]]
    attach: int
    f: List[int]
    g: List[int]
    h: List[int]
    radio_score: Tuple[float, float, float, float]

def make_entry(name: str, n: int, edges: List[Tuple[int, int]], attach: Optional[int] = None) -> TreeEntry:
    if attach is None:
        attach = choose_attach_vertex(n, edges)
    f, g = independence_fg_tree(n, edges, attach)
    h = poly_add(f, g)
    score = radio_score(h)
    return TreeEntry(
        name=name, n=n, edges=edges, attach=attach,
        f=f, g=g, h=h, radio_score=score
    )

def build_library(library_size: int, max_n: int, seed: int) -> List[TreeEntry]:
    rng = random.Random(seed)
    lib: List[TreeEntry] = []
    print(f"[library] Building library size={library_size} max_n={max_n}...")

    # Standard families
    for n in range(2, min(60, max_n) + 1):
        lib.append(make_entry(f"path_{n}", n, path_edges(n), attach=0))
        lib.append(make_entry(f"star_{n}", n, star_edges(n), attach=0))

    # Brooms
    for _ in range(int(library_size * 0.2)):
        handle = rng.randrange(5, 60)
        leaves = rng.randrange(5, 120)
        n, edges = broom_edges(handle, leaves)
        if n <= max_n: lib.append(make_entry(f"broom_{handle}_{leaves}", n, edges))

    # Caterpillars
    for _ in range(int(library_size * 0.2)):
        spine = rng.randrange(5, 70)
        leaf_counts = [rng.randrange(0, 10) for _ in range(spine)]
        n, edges = caterpillar_edges(spine, leaf_counts)
        if n <= max_n: lib.append(make_entry(f"cat_{spine}", n, edges))

    # Galvin Trees
    for _ in range(int(library_size * 0.3)):
        d = rng.randrange(2, 90)
        e = rng.randrange(2, 25)
        n, edges = tree_T_de(d, e)
        if n <= max_n: lib.append(make_entry(f"T_{d}_{e}", n, edges, attach=0))

    # Random
    while len(lib) < library_size:
        n = rng.randrange(30, max_n + 1)
        edges = random_tree_edges(n, rng)
        lib.append(make_entry(f"rnd_{n}_{rng.randrange(10**9)}", n, edges))

    lib.sort(key=lambda e: e.radio_score, reverse=True)
    return lib

def _guided_pool(lib: List[TreeEntry], rng: random.Random, pool_size: int) -> List[TreeEntry]:
    pool_size = min(pool_size, len(lib))
    return lib[:pool_size]

def beam_search_forest(
    lib: List[TreeEntry],
    rng: random.Random,
    k_min: int,
    k_max: int,
    beam_width: int,
    pool_size: int,
) -> Optional[Dict[str, Any]]:
    pool = _guided_pool(lib, rng, max(pool_size, beam_width))
    states: List[Tuple[List[int], List[TreeEntry]]] = [([1], [])]
    
    for depth in range(1, k_max + 1):
        new_states = []
        for poly, comps in states:
            for e in pool:
                poly2 = poly_mul(poly, e.h)
                score = radio_score(poly2)
                new_states.append((score, poly2, comps + [e]))
        
        new_states.sort(key=lambda s: s[0], reverse=True)
        states = [(poly, comps) for _, poly, comps in new_states[:beam_width]]
        
        if depth >= k_min:
            for _, poly, comps in new_states[:beam_width]:
                if not is_unimodal(poly):
                    return {
                        "type": "forest",
                        "components": comps,
                        "sequence": poly,
                        "witness": first_unimodality_break(poly)
                    }
    return None

def beam_search_star(
    lib: List[TreeEntry],
    rng: random.Random,
    k_min: int,
    k_max: int,
    beam_width: int,
    pool_size: int,
) -> Optional[Dict[str, Any]]:
    pool = _guided_pool(lib, rng, max(pool_size, beam_width))
    states: List[Tuple[List[int], List[int], List[TreeEntry]]] = [([1], [1], [])]
    
    for depth in range(1, k_max + 1):
        new_states = []
        for prod_h, prod_f, branches in states:
            for br in pool:
                prod_h2 = poly_mul(prod_h, br.h)
                prod_f2 = poly_mul(prod_f, br.f)
                poly = poly_add(prod_h2, [0] + prod_f2)
                score = radio_score(poly)
                new_states.append((score, prod_h2, prod_f2, poly, branches + [br]))
        
        new_states.sort(key=lambda s: s[0], reverse=True)
        states = [(ph, pf, brs) for _, ph, pf, _, brs in new_states[:beam_width]]
        
        if depth >= k_min:
            for _, ph, pf, poly, brs in new_states[:beam_width]:
                if not is_unimodal(poly):
                    return {
                        "type": "tree_star",
                        "branches": brs,
                        "sequence": poly,
                        "witness": first_unimodality_break(poly)
                    }
    return None

def verify_and_write(result: Dict[str, Any], out_dir: str, seed: int) -> str:
    os.makedirs(out_dir, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    path = os.path.join(out_dir, f"certificate_radio_{result['type']}_{ts}.json")
    
    cert = {
        "problem": "Erdos Problem #993",
        "method": "Radio Roots Search (High NM, High Re, Low Im)",
        "found_at": time.strftime("%c"),
        "seed": seed,
    }
    
    if result["type"] == "forest":
        comps = []
        for c in result["components"]:
            comps.append({"name": c.name, "n": c.n, "edges": [[int(u),int(v)] for u,v in c.edges], "attach": int(c.attach)})
        cert["type"] = "forest"
        cert["components"] = comps
        cert["sequence"] = [int(x) for x in result["sequence"]]
        cert["witness"] = result["witness"]
        
    elif result["type"] == "tree_star":
        branches = []
        for c in result["branches"]:
            branches.append({"name": c.name, "n": c.n, "attach": int(c.attach)})
        
        n_total = 0
        offset = 1
        edges_total = []
        for br in result["branches"]:
            edges_total.append((0, offset + br.attach))
            for u, v in br.edges:
                edges_total.append((offset + u, offset + v))
            offset += br.n
        n_total = offset
        
        cert["type"] = "tree"
        cert["construction"] = "star_of_subtrees"
        cert["n"] = int(n_total)
        cert["edges"] = [[int(u),int(v)] for u,v in edges_total]
        cert["branches"] = branches
        cert["sequence"] = [int(x) for x in result["sequence"]]
        cert["witness"] = result["witness"]

    with open(path, "w", encoding="utf-8") as f:
        json.dump(cert, f, indent=2)
    return path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library-size", type=int, default=1000)
    ap.add_argument("--max-n", type=int, default=180)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--beam-width", type=int, default=50)
    ap.add_argument("--pool-size", type=int, default=200)
    ap.add_argument("--time-limit", type=float, default=3600.0)
    ap.add_argument("--out-dir", type=str, default="out_radio_optimized")
    args = ap.parse_args()
    
    print(f"Starting Optimized Radio Roots Hunt (seed={args.seed})...")
    
    lib = build_library(args.library_size, args.max_n, args.seed)
    
    t0 = time.time()
    iteration = 0
    while True:
        elapsed = time.time() - t0
        if elapsed > args.time_limit:
            print("Time limit reached.")
            break
            
        iteration += 1
        best_entry = lib[0]
        b_score, b_nm, b_re, b_im = best_entry.radio_score
        
        print(f"[Iter {iteration}] Best: {best_entry.name} (n={best_entry.n})")
        print(f"    Score={b_score:.2f} (NM={b_nm:.4f}, Re={b_re:.4f}, Im={-b_im:.4f})")
        
        try:
            roots = np.roots(best_entry.h[::-1])
            roots = sorted(roots, key=lambda r: r.real, reverse=True)
            print(f"    Top 5 Roots (by Re): {[f'{r:.2f}' for r in roots[:5]]}")
        except:
            pass

        print(f"[Iter {iteration}] Elapsed: {elapsed:.1f}s. Running Forest Beam Search...")
        
        res = beam_search_forest(lib, random.Random(args.seed + iteration), 2, 5, args.beam_width, args.pool_size)
        if res:
            print("COUNTEREXAMPLE FOUND (Forest)!")
            p = verify_and_write(res, args.out_dir, args.seed)
            print(f"Certificate: {p}")
            break
            
        print(f"[Iter {iteration}] Running Star Beam Search...")
        res = beam_search_star(lib, random.Random(args.seed + iteration*100), 3, 10, args.beam_width, args.pool_size)
        if res:
            print("COUNTEREXAMPLE FOUND (Tree Star)!")
            p = verify_and_write(res, args.out_dir, args.seed)
            print(f"Certificate: {p}")
            break
            
        random.shuffle(lib)
        lib.sort(key=lambda e: e.radio_score, reverse=True) 
        
if __name__ == "__main__":
    main()
