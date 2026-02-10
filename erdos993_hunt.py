#!/usr/bin/env python3
"""
Erdos Problem #993 hunt: find a tree or forest whose independent set sequence is NOT unimodal.

Problem statement (ErdosProblems #993): "The independent set sequence of any tree or forest is unimodal."
One explicit counterexample tree/forest settles it.

This script:
  - Builds (or loads) a library of candidate trees, including structured families and random trees.
  - Searches in parallel in two ways:
      (A) Forest search: disjoint unions of library trees (convolution of coefficient sequences).
      (B) Star-of-subtrees search: a single tree formed by joining a new root to the roots of several library trees.
  - If a non-unimodal sequence is found, it writes a JSON certificate you can post publicly.

Usage (minimal):
  python erdos993_hunt.py

Useful options:
  python erdos993_hunt.py --workers 8 --time-limit 36000
  python erdos993_hunt.py --library-size 3000 --max-n 260
  python erdos993_hunt.py --mode forest
  python erdos993_hunt.py --mode star

Files written:
  - library.pkl : cached library of tree components
  - certificate_*.json : when a counterexample is found
  - progress_worker*.json : periodic progress snapshots

Notes:
  - Requires Python 3.9+. Uses numpy if available (much faster polynomial convolution).
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import math
import os
import pickle
import random
import sys
import time
from collections import deque
from typing import List, Tuple, Optional, Dict, Any

# Optional numpy acceleration for polynomial multiplication.
try:
    import numpy as _np  # type: ignore
    _HAVE_NUMPY = True
except Exception:
    _HAVE_NUMPY = False


# -------------------------
# Polynomial utilities
# -------------------------

def poly_add(a: List[int], b: List[int]) -> List[int]:
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        if i < len(a):
            out[i] += a[i]
        if i < len(b):
            out[i] += b[i]
    # trim
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out

def poly_mul_py(a: List[int], b: List[int]) -> List[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0:
                continue
            out[i + j] += ai * bj
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out

def poly_mul(a: List[int], b: List[int]) -> List[int]:
    # Heuristic: numpy object convolution is often faster for moderate degrees.
    if _HAVE_NUMPY and (len(a) * len(b) >= 20_000):
        return list(_np.convolve(_np.array(a, dtype=object), _np.array(b, dtype=object)))
    return poly_mul_py(a, b)

def poly_pow(p: List[int], e: int) -> List[int]:
    # exponentiation by squaring
    out = [1]
    base = p
    exp = e
    while exp > 0:
        if exp & 1:
            out = poly_mul(out, base)
        exp >>= 1
        if exp:
            base = poly_mul(base, base)
    return out

def is_unimodal(seq: List[int]) -> bool:
    # unimodal allowing plateaus: nondecreasing then nonincreasing
    n = len(seq)
    i = 0
    while i + 1 < n and seq[i] <= seq[i + 1]:
        i += 1
    while i + 1 < n and seq[i] >= seq[i + 1]:
        i += 1
    return i == n - 1

def first_unimodality_break(seq: List[int]) -> Optional[Dict[str, int]]:
    """
    Returns a small witness of non-unimodality:
      an index i where seq[i+1] < seq[i] (a "down"),
      and later an index j where seq[j+1] > seq[j] (an "up after down").
    """
    saw_down_at = None
    for i in range(len(seq) - 1):
        if seq[i + 1] < seq[i] and saw_down_at is None:
            saw_down_at = i
        if saw_down_at is not None and seq[i + 1] > seq[i]:
            return {"down_at": saw_down_at, "up_at": i}
    return None

def log_concavity_break_score(seq: List[int]) -> Tuple[int, int]:
    """
    Returns (count, total_gap) where count is how many k satisfy seq[k]^2 < seq[k-1]*seq[k+1],
    and total_gap is sum(seq[k-1]*seq[k+1] - seq[k]^2) over those k.
    A nonzero score is often a useful "wiggliness" indicator even if unimodality still holds.
    """
    cnt = 0
    tot = 0
    for k in range(1, len(seq) - 1):
        lhs = seq[k] * seq[k]
        rhs = seq[k - 1] * seq[k + 1]
        if lhs < rhs:
            cnt += 1
            tot += (rhs - lhs)
    return cnt, tot


def sequence_score(seq: List[int]) -> Tuple[int, int, float, int]:
    """
    Heuristic score for "wiggliness" to guide search.
    Higher is better; intended for best-first/beam search over polynomials.
    """
    lc_breaks, lc_gap = log_concavity_break_score(seq)
    if len(seq) <= 1:
        return (lc_breaks, lc_gap, 0.0, len(seq))
    tv = 0
    for i in range(len(seq) - 1):
        tv += abs(seq[i + 1] - seq[i])
    total = sum(seq)
    rough = tv / max(1, total)
    return (lc_breaks, lc_gap, rough, len(seq))


# -------------------------
# Graph / tree utilities
# -------------------------

def check_tree(n: int, edges: List[Tuple[int, int]]) -> bool:
    if n == 0:
        return True
    if len(edges) != n - 1:
        return False
    adj = [[] for _ in range(n)]
    for u, v in edges:
        if not (0 <= u < n and 0 <= v < n) or u == v:
            return False
        adj[u].append(v)
        adj[v].append(u)
    # connectivity via BFS
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
    # Prüfer sequence method (uniform labelled tree)
    if n <= 1:
        return []
    prufer = [rng.randrange(n) for _ in range(n - 2)]
    degree = [1] * n
    for x in prufer:
        degree[x] += 1
    import heapq
    leaves = [i for i, d in enumerate(degree) if d == 1]
    heapq.heapify(leaves)
    edges: List[Tuple[int, int]] = []
    for x in prufer:
        leaf = heapq.heappop(leaves)
        edges.append((leaf, x))
        degree[leaf] -= 1
        degree[x] -= 1
        if degree[x] == 1:
            heapq.heappush(leaves, x)
    a = heapq.heappop(leaves)
    b = heapq.heappop(leaves)
    edges.append((a, b))
    return edges

def path_edges(n: int) -> List[Tuple[int, int]]:
    return [(i, i + 1) for i in range(n - 1)]

def star_edges(n: int) -> List[Tuple[int, int]]:
    return [(0, i) for i in range(1, n)]

def broom_edges(handle: int, leaves: int) -> Tuple[int, List[Tuple[int, int]]]:
    # path of length handle (handle+1 vertices), leaves pendant to end
    n = handle + 1 + leaves
    edges = [(i, i + 1) for i in range(handle)]
    end = handle
    for j in range(leaves):
        edges.append((end, handle + 1 + j))
    return n, edges

def caterpillar_edges(spine_len: int, leaf_counts: List[int]) -> Tuple[int, List[Tuple[int, int]]]:
    assert spine_len == len(leaf_counts)
    n = spine_len + sum(leaf_counts)
    edges = [(i, i + 1) for i in range(spine_len - 1)]
    v = spine_len
    for i, c in enumerate(leaf_counts):
        for _ in range(c):
            edges.append((i, v))
            v += 1
    return n, edges

def spherically_symmetric_tree(branching: List[int]) -> Tuple[int, List[Tuple[int, int]]]:
    """
    branching[d] = number of children each vertex at depth d has.
    Depth runs 0..len(branching)-1 (so total height = len(branching)).
    """
    edges: List[Tuple[int, int]] = []
    # start with root 0
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
    # The depth-3 spherically symmetric tree from Galvin's construction: branching [d, e, 1]
    return spherically_symmetric_tree([d, e, 1])

def choose_attach_vertex(n: int, edges: List[Tuple[int, int]]) -> int:
    # pick a high-degree vertex as the attachment root for "star-of-subtrees" constructions
    deg = [0] * n
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1
    maxdeg = max(deg) if n else 0
    for i, d in enumerate(deg):
        if d == maxdeg:
            return i
    return 0


# -------------------------
# Independence polynomial on trees
# -------------------------

def independence_fg_tree(n: int, edges: List[Tuple[int, int]], root: int) -> Tuple[List[int], List[int]]:
    """
    For a rooted tree (n, edges), returns:
      f = coefficient list for independent sets in the tree that EXCLUDE the root
      g = coefficient list for independent sets in the tree that INCLUDE the root
    Then total independence polynomial is h = f + g.
    """
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
    f = [None] * n  # type: ignore
    g = [None] * n  # type: ignore

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

def independence_poly_tree(n: int, edges: List[Tuple[int, int]], root: int = 0) -> List[int]:
    f, g = independence_fg_tree(n, edges, root)
    return poly_add(f, g)


# -------------------------
# Library entries
# -------------------------

@dataclasses.dataclass
class TreeEntry:
    name: str
    n: int
    edges: List[Tuple[int, int]]
    attach: int
    f: List[int]
    g: List[int]
    h: List[int]        # total polynomial at attach
    alpha: int          # independence number (degree of polynomial)
    peak: int           # index of a (first) maximum coefficient
    unimodal: bool
    lc_breaks: int
    lc_gap: int

def _peak_index(seq: List[int]) -> int:
    m = max(seq)
    for i, x in enumerate(seq):
        if x == m:
            return i
    return 0

def make_entry(name: str, n: int, edges: List[Tuple[int, int]], attach: Optional[int] = None) -> TreeEntry:
    if attach is None:
        attach = choose_attach_vertex(n, edges)
    f, g = independence_fg_tree(n, edges, attach)
    h = poly_add(f, g)
    unim = is_unimodal(h)
    alpha = len(h) - 1
    peak = _peak_index(h)
    lc_breaks, lc_gap = log_concavity_break_score(h)
    return TreeEntry(
        name=name,
        n=n,
        edges=edges,
        attach=attach,
        f=f,
        g=g,
        h=h,
        alpha=alpha,
        peak=peak,
        unimodal=unim,
        lc_breaks=lc_breaks,
        lc_gap=lc_gap,
    )

def build_library(library_size: int, max_n: int, seed: int, verbose: bool = True) -> List[TreeEntry]:
    rng = random.Random(seed)
    lib: List[TreeEntry] = []

    def add(entry: TreeEntry):
        lib.append(entry)
        if verbose and len(lib) % 200 == 0:
            print(f"[library] {len(lib)} entries...")

    # 1) Small structured families (cheap, diverse)
    for n in range(2, min(60, max_n) + 1):
        add(make_entry(f"path_{n}", n, path_edges(n), attach=0))
        add(make_entry(f"star_{n}", n, star_edges(n), attach=0))

    # 2) Brooms (handle+leaves)
    for _ in range(300):
        handle = rng.randrange(5, 60)
        leaves = rng.randrange(5, 120)
        n, edges = broom_edges(handle, leaves)
        if n <= max_n:
            add(make_entry(f"broom_{handle}_{leaves}", n, edges))

    # 3) Caterpillars with random leaf counts along spine
    for _ in range(600):
        spine = rng.randrange(5, 70)
        leaf_counts = [rng.randrange(0, 10) for _ in range(spine)]
        n, edges = caterpillar_edges(spine, leaf_counts)
        if n <= max_n:
            add(make_entry(f"cat_{spine}", n, edges))

    # 4) Galvin-style depth-3 symmetric trees T(d,e) = branching [d,e,1]
    #    (from the "non log-concave" constructions; good for "wiggly" coefficients)
    for d in range(2, 90):
        for e in range(2, 25):
            n, edges = tree_T_de(d, e)
            if n <= max_n:
                add(make_entry(f"T_{d}_{e}", n, edges, attach=0))
        if len(lib) >= library_size:
            break

    # 4b) Galvin-style with unary tail: branching [d, e] + [1]*t
    #     This extends the depth-3 construction with a long unary tail.
    for d in range(2, 60):
        for e in range(2, 20):
            for t in range(2, 8):
                # quick size check for [d, e] + [1]*t
                n_est = 1 + d + d * e * (1 + t)
                if n_est > max_n:
                    continue
                branching = [d, e] + [1] * t
                n, edges = spherically_symmetric_tree(branching)
                if n <= max_n:
                    add(make_entry(f"T_{d}_{e}_1^{t}", n, edges, attach=0))
        if len(lib) >= library_size:
            break

    # 5) Deeper symmetric trees: [2]*a + [1]*b
    #    Motivated by examples where log-concavity breaks in multiple places.
    for _ in range(800):
        a = rng.randrange(2, 14)   # number of "binary" levels
        b = rng.randrange(1, 50)   # number of "unary" levels after that
        branching = [2] * a + [1] * b
        n, edges = spherically_symmetric_tree(branching)
        if n <= max_n:
            add(make_entry(f"sym_2^{a}_1^{b}", n, edges, attach=0))

    # 6) Random labelled trees (Prüfer)
    while len(lib) < library_size:
        n = rng.randrange(30, max_n + 1)
        edges = random_tree_edges(n, rng)
        add(make_entry(f"rnd_{n}_{rng.randrange(10**9)}", n, edges))

    # If any entry is already non-unimodal, we've already solved the problem.
    bad = [e for e in lib if not e.unimodal]
    if bad and verbose:
        print(f"[library] WARNING: found {len(bad)} non-unimodal tree(s) during library build!")
        print(f"[library] First example: {bad[0].name} (n={bad[0].n})")

    # Sort so "interesting" (log-concavity-breaking) entries are earlier
    lib.sort(key=lambda e: (e.lc_breaks, e.lc_gap, e.n), reverse=True)

    return lib


# -------------------------
# Deterministic Galvin-style sweep
# -------------------------

def galvin_sweep(
    max_n: int,
    d_max: int,
    e_max: int,
    t_max: int,
    top_k: int,
) -> Dict[str, Any]:
    """
    Deterministic sweep of Galvin-style symmetric trees with a unary tail:
      branching = [d, e] + [1] * t, where t >= 1.
    Returns top_k entries ranked by sequence_score.
    """
    top: List[Tuple[Tuple[int, int, float, int], str, Dict[str, Any]]] = []
    checked = 0
    first_counter: Optional[Dict[str, Any]] = None

    for d in range(2, d_max + 1):
        for e in range(2, e_max + 1):
            for t in range(1, t_max + 1):
                n_est = 1 + d + d * e * (1 + t)
                if n_est > max_n:
                    continue
                branching = [d, e] + [1] * t
                n, edges = spherically_symmetric_tree(branching)
                if n > max_n:
                    continue
                name = f"T_{d}_{e}_1^{t}"
                entry = make_entry(name, n, edges, attach=0)
                score = sequence_score(entry.h)
                data = {
                    "name": name,
                    "branching": branching,
                    "n": entry.n,
                    "alpha": entry.alpha,
                    "peak": entry.peak,
                    "unimodal": entry.unimodal,
                    "lc_breaks": entry.lc_breaks,
                    "lc_gap": entry.lc_gap,
                    "sequence": [int(x) for x in entry.h],
                }
                if not entry.unimodal:
                    data["witness"] = first_unimodality_break(entry.h)
                    if first_counter is None:
                        first_counter = {
                            "name": name,
                            "n": entry.n,
                            "edges": entry.edges,
                            "attach": entry.attach,
                            "sequence": [int(x) for x in entry.h],
                            "witness": data["witness"],
                            "branching": branching,
                        }
                checked += 1
                top.append((score, name, data))
                top.sort(key=lambda x: (x[0], x[1]), reverse=True)
                if len(top) > top_k:
                    top.pop()

    top_out: List[Dict[str, Any]] = []
    for score, _, data in top:
        data["score"] = [int(score[0]), int(score[1]), float(score[2]), int(score[3])]
        top_out.append(data)

    return {
        "checked": checked,
        "max_n": max_n,
        "d_max": d_max,
        "e_max": e_max,
        "t_max": t_max,
        "top_k": top_k,
        "top": top_out,
        "first_counterexample": first_counter,
    }


def load_sweep_report(path: str, max_n: int, limit: int) -> List[TreeEntry]:
    if limit <= 0 or not os.path.exists(path):
        return []
    with open(path, "r", encoding="utf-8") as f:
        report = json.load(f)
    top = report.get("top", [])
    out: List[TreeEntry] = []
    for item in top:
        branching = item.get("branching")
        name = item.get("name")
        if not isinstance(branching, list) or not branching:
            continue
        if not isinstance(name, str) or not name:
            name = "sweep_entry"
        n, edges = spherically_symmetric_tree(branching)
        if n > max_n:
            continue
        out.append(make_entry(name, n, edges, attach=0))
        if len(out) >= limit:
            break
    return out


# -------------------------
# Constructions to search
# -------------------------

def build_star_tree(branches: List[TreeEntry]) -> Tuple[int, List[Tuple[int, int]]]:
    """
    Construct a tree by taking a new root vertex 0 and attaching each branch's attach-vertex to 0.
    The branch trees are vertex-disjoint except for the new root.
    """
    edges: List[Tuple[int, int]] = []
    offset = 1
    # root is 0
    for br in branches:
        # map br vertices to [offset, offset+br.n-1]
        # connect root to mapped attach
        edges.append((0, offset + br.attach))
        for u, v in br.edges:
            edges.append((offset + u, offset + v))
        offset += br.n
    n_total = offset
    return n_total, edges


def build_double_star_tree(
    left: List[TreeEntry], right: List[TreeEntry]
) -> Tuple[int, List[Tuple[int, int]]]:
    """
    Construct a tree with two hubs (0 and 1) connected by an edge.
    Left branches attach to hub 0, right branches attach to hub 1.
    """
    edges: List[Tuple[int, int]] = []
    edges.append((0, 1))
    offset = 2

    for br in left:
        edges.append((0, offset + br.attach))
        for u, v in br.edges:
            edges.append((offset + u, offset + v))
        offset += br.n

    for br in right:
        edges.append((1, offset + br.attach))
        for u, v in br.edges:
            edges.append((offset + u, offset + v))
        offset += br.n

    n_total = offset
    return n_total, edges

def star_poly_from_branches(branches: List[TreeEntry]) -> List[int]:
    """
    Independence polynomial of the star-of-subtrees construction, computed from branch f/g at attach:
      f_root = Π (h_i)
      g_root = x * Π (f_i)
      total = f_root + g_root
    """
    prod_h = [1]
    prod_f = [1]
    for br in branches:
        prod_h = poly_mul(prod_h, br.h)
        prod_f = poly_mul(prod_f, br.f)
    g_root = poly_mul([0, 1], prod_f)  # x * prod_f
    return poly_add(prod_h, g_root)

def forest_poly_from_components(comps: List[TreeEntry]) -> List[int]:
    poly = [1]
    for e in comps:
        poly = poly_mul(poly, e.h)  # component total polynomial
    return poly


def double_star_poly_from_branches(
    left: List[TreeEntry], right: List[TreeEntry]
) -> List[int]:
    """
    Two-hub construction: an edge between hubs u--v.
    Each hub has its own set of branches (star-of-subtrees).

    For each hub:
      f = product(h_i)  (hub excluded)
      g = x * product(f_i) (hub included)

    Total polynomial on edge u--v:
      f_u * f_v + g_u * f_v + f_u * g_v
    """
    # left hub
    f_left = [1]
    g_left_base = [1]
    for br in left:
        f_left = poly_mul(f_left, br.h)
        g_left_base = poly_mul(g_left_base, br.f)
    g_left = poly_mul([0, 1], g_left_base)

    # right hub
    f_right = [1]
    g_right_base = [1]
    for br in right:
        f_right = poly_mul(f_right, br.h)
        g_right_base = poly_mul(g_right_base, br.f)
    g_right = poly_mul([0, 1], g_right_base)

    # combine (u excluded, v excluded) + (u included, v excluded) + (u excluded, v included)
    term1 = poly_mul(f_left, f_right)
    term2 = poly_mul(g_left, f_right)
    term3 = poly_mul(f_left, g_right)
    return poly_add(poly_add(term1, term2), term3)


# -------------------------
# Search loops
# -------------------------

def biased_choice(lib: List[TreeEntry], rng: random.Random) -> TreeEntry:
    """
    Choose an entry with a bias toward "wiggly" ones, but still some randomness.
    We use a simple mixture: with probability ~0.7 choose from top chunk.
    """
    if not lib:
        raise ValueError("empty library")
    r = rng.random()
    if r < 0.70:
        top = min(250, len(lib))
        return lib[rng.randrange(top)]
    if r < 0.95:
        mid = min(1200, len(lib))
        return lib[rng.randrange(mid)]
    return lib[rng.randrange(len(lib))]

def search_forest(
    lib: List[TreeEntry],
    rng: random.Random,
    k_min: int,
    k_max: int,
    max_tries: int,
    end_time: Optional[float] = None,
) -> Optional[Dict[str, Any]]:
    for _ in range(max_tries):
        if end_time is not None and time.time() >= end_time:
            return None
        k = rng.randrange(k_min, k_max + 1)
        comps = [biased_choice(lib, rng) for _ in range(k)]
        poly = forest_poly_from_components(comps)
        if not is_unimodal(poly):
            witness = first_unimodality_break(poly)
            return {
                "type": "forest",
                "components": comps,
                "sequence": poly,
                "witness": witness,
            }
    return None

def search_star(
    lib: List[TreeEntry],
    rng: random.Random,
    k_min: int,
    k_max: int,
    max_tries: int,
    end_time: Optional[float] = None,
) -> Optional[Dict[str, Any]]:
    for _ in range(max_tries):
        if end_time is not None and time.time() >= end_time:
            return None
        k = rng.randrange(k_min, k_max + 1)
        branches = [biased_choice(lib, rng) for _ in range(k)]
        poly = star_poly_from_branches(branches)
        if not is_unimodal(poly):
            witness = first_unimodality_break(poly)
            return {
                "type": "tree_star",
                "branches": branches,
                "sequence": poly,
                "witness": witness,
            }
    return None


def search_double_star(
    lib: List[TreeEntry],
    rng: random.Random,
    left_min: int,
    left_max: int,
    right_min: int,
    right_max: int,
    max_tries: int,
    end_time: Optional[float] = None,
) -> Optional[Dict[str, Any]]:
    for _ in range(max_tries):
        if end_time is not None and time.time() >= end_time:
            return None
        k_left = rng.randrange(left_min, left_max + 1)
        k_right = rng.randrange(right_min, right_max + 1)
        left = [biased_choice(lib, rng) for _ in range(k_left)]
        right = [biased_choice(lib, rng) for _ in range(k_right)]
        poly = double_star_poly_from_branches(left, right)
        if not is_unimodal(poly):
            witness = first_unimodality_break(poly)
            return {
                "type": "tree_double_star",
                "left": left,
                "right": right,
                "sequence": poly,
                "witness": witness,
            }
    return None


def _guided_pool(lib: List[TreeEntry], rng: random.Random, pool_size: int) -> List[TreeEntry]:
    pool_size = min(pool_size, len(lib))
    scored: List[Tuple[Tuple[int, int, float, int], float, TreeEntry]] = []
    for e in lib:
        scored.append((sequence_score(e.h), rng.random(), e))
    scored.sort(key=lambda x: (x[0], x[1]), reverse=True)
    pool = [e for _, _, e in scored[:pool_size]]
    return pool


def beam_search_forest(
    lib: List[TreeEntry],
    rng: random.Random,
    k_min: int,
    k_max: int,
    beam_width: int,
    pool_size: int,
) -> Optional[Dict[str, Any]]:
    if not lib:
        return None
    pool = _guided_pool(lib, rng, max(pool_size, beam_width))

    states: List[Tuple[List[int], List[TreeEntry]]] = [([1], [])]
    for depth in range(1, k_max + 1):
        new_states: List[Tuple[Tuple[int, int, float, int], List[int], List[TreeEntry]]] = []
        for poly, comps in states:
            for e in pool:
                poly2 = poly_mul(poly, e.h)
                score = sequence_score(poly2)
                new_states.append((score, poly2, comps + [e]))
        new_states.sort(key=lambda s: s[0], reverse=True)
        states = [(poly, comps) for _, poly, comps in new_states[:beam_width]]
        if depth >= k_min:
            for poly, comps in states:
                if not is_unimodal(poly):
                    witness = first_unimodality_break(poly)
                    return {
                        "type": "forest",
                        "components": comps,
                        "sequence": poly,
                        "witness": witness,
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
    if not lib:
        return None
    pool = _guided_pool(lib, rng, max(pool_size, beam_width))

    states: List[Tuple[List[int], List[int], List[TreeEntry]]] = [([1], [1], [])]
    for depth in range(1, k_max + 1):
        new_states: List[Tuple[Tuple[int, int, float, int], List[int], List[int], List[int], List[TreeEntry]]] = []
        for prod_h, prod_f, branches in states:
            for br in pool:
                prod_h2 = poly_mul(prod_h, br.h)
                prod_f2 = poly_mul(prod_f, br.f)
                poly = poly_add(prod_h2, [0] + prod_f2)
                score = sequence_score(poly)
                new_states.append((score, prod_h2, prod_f2, poly, branches + [br]))
        new_states.sort(key=lambda s: s[0], reverse=True)
        states = [(ph, pf, brs) for _, ph, pf, _, brs in new_states[:beam_width]]
        if depth >= k_min:
            for _, ph, pf, poly, brs in new_states[:beam_width]:
                if not is_unimodal(poly):
                    witness = first_unimodality_break(poly)
                    return {
                        "type": "tree_star",
                        "branches": brs,
                        "sequence": poly,
                        "witness": witness,
                    }
    return None


def guided_search(lib: List[TreeEntry], rng: random.Random, args: argparse.Namespace) -> Optional[Dict[str, Any]]:
    pool_size = max(10, args.guided_pool)
    beam_width = max(5, min(args.guided_beam, pool_size))
    if args.mode in ("mixed", "forest"):
        res = beam_search_forest(
            lib, rng, args.forest_k_min, args.forest_k_max, beam_width, pool_size
        )
        if res is not None:
            return res
    if args.mode in ("mixed", "star"):
        res = beam_search_star(
            lib, rng, args.star_k_min, args.star_k_max, beam_width, pool_size
        )
        if res is not None:
            return res
    return None


# -------------------------
# Certificate writing + verification
# -------------------------

def verify_and_write_certificate(result: Dict[str, Any], out_dir: str, seed: int, meta: Dict[str, Any]) -> str:
    os.makedirs(out_dir, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out_path = os.path.join(out_dir, f"certificate_{result['type']}_{ts}.json")

    cert: Dict[str, Any] = {
        "problem": "Erdos Problem #993",
        "found_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "seed": seed,
        "meta": meta,
    }

    if result["type"] == "forest":
        comps: List[TreeEntry] = result["components"]
        # recompute each component polynomial independently for safety
        comp_certs = []
        for c in comps:
            poly_c = independence_poly_tree(c.n, c.edges, root=c.attach)
            comp_certs.append({
                "name": c.name,
                "n": c.n,
                "edges": [[int(u), int(v)] for (u, v) in c.edges],
                "attach": int(c.attach),
                "independent_set_sequence": [int(x) for x in poly_c],
                "unimodal": bool(is_unimodal(poly_c)),
            })
        # forest sequence = product of component sequences
        poly = [1]
        for cc in comp_certs:
            poly = poly_mul(poly, cc["independent_set_sequence"])
        cert["type"] = "forest"
        cert["components"] = comp_certs
        cert["forest_sequence"] = [int(x) for x in poly]
        cert["forest_unimodal"] = bool(is_unimodal(poly))
        cert["witness"] = first_unimodality_break(poly)

    elif result["type"] == "tree_star":
        branches: List[TreeEntry] = result["branches"]
        n_total, edges_total = build_star_tree(branches)
        # verify it's a tree and recompute polynomial on the full tree
        if not check_tree(n_total, edges_total):
            raise RuntimeError("Internal error: constructed object is not a tree.")
        poly_full = independence_poly_tree(n_total, edges_total, root=0)
        cert["type"] = "tree"
        cert["construction"] = {
            "kind": "star_of_subtrees",
            "branches": [{"name": b.name, "n": b.n, "attach": int(b.attach)} for b in branches],
        }
        cert["n"] = int(n_total)
        cert["edges"] = [[int(u), int(v)] for (u, v) in edges_total]
        cert["independent_set_sequence"] = [int(x) for x in poly_full]
        cert["unimodal"] = bool(is_unimodal(poly_full))
        cert["witness"] = first_unimodality_break(poly_full)

    elif result["type"] == "tree_double_star":
        left: List[TreeEntry] = result["left"]
        right: List[TreeEntry] = result["right"]
        n_total, edges_total = build_double_star_tree(left, right)
        if not check_tree(n_total, edges_total):
            raise RuntimeError("Internal error: constructed object is not a tree.")
        poly_full = independence_poly_tree(n_total, edges_total, root=0)
        cert["type"] = "tree"
        cert["construction"] = {
            "kind": "double_star",
            "left": [{"name": b.name, "n": b.n, "attach": int(b.attach)} for b in left],
            "right": [{"name": b.name, "n": b.n, "attach": int(b.attach)} for b in right],
        }
        cert["n"] = int(n_total)
        cert["edges"] = [[int(u), int(v)] for (u, v) in edges_total]
        cert["independent_set_sequence"] = [int(x) for x in poly_full]
        cert["unimodal"] = bool(is_unimodal(poly_full))
        cert["witness"] = first_unimodality_break(poly_full)

    else:
        raise ValueError(f"Unknown result type: {result['type']}")

    # sanity: must be non-unimodal
    if cert.get("type") == "forest":
        if cert.get("forest_unimodal", True):
            raise RuntimeError("Internal error: certificate forest is unimodal after verification.")
    else:
        if cert.get("unimodal", True):
            raise RuntimeError("Internal error: certificate tree is unimodal after verification.")

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(cert, f, indent=2)
    return out_path


def write_tree_certificate(
    n: int,
    edges: List[Tuple[int, int]],
    attach: int,
    poly: List[int],
    out_dir: str,
    seed: int,
    meta: Dict[str, Any],
    prefix: str = "certificate_tree",
) -> str:
    os.makedirs(out_dir, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out_path = os.path.join(out_dir, f"{prefix}_{ts}.json")
    cert = {
        "problem": "Erdos Problem #993",
        "found_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "seed": seed,
        "type": "tree",
        "n": int(n),
        "edges": [[int(u), int(v)] for (u, v) in edges],
        "attach": int(attach),
        "independent_set_sequence": [int(x) for x in poly],
        "unimodal": bool(is_unimodal(poly)),
        "witness": first_unimodality_break(poly),
        "meta": meta,
    }
    if cert["unimodal"]:
        raise RuntimeError("Internal error: certificate tree is unimodal after verification.")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(cert, f, indent=2)
    return out_path


# -------------------------
# Worker loop
# -------------------------

def worker_main(worker_id: int, args: argparse.Namespace, lib_path: str, found_flag, queue) -> None:
    # Each worker loads the same library from disk (simpler than cross-process sharing of big ints).
    with open(lib_path, "rb") as f:
        lib: List[TreeEntry] = pickle.load(f)

    rng = random.Random(args.seed + 10_000 * worker_id + int(time.time()) % 10_000)
    t0 = time.time()
    end_time = t0 + args.time_limit if args.time_limit is not None else None
    last_report = time.time()
    last_guided = time.time()

    best = {"score": (0, 0), "kind": None, "name": None}
    tries = 0

    if worker_id == 0 and args.guided_every and args.guided_every > 0:
        res = guided_search(lib, rng, args)
        if res is not None:
            queue.put(res)
            found_flag.set()
            return
        last_guided = time.time()

    while not found_flag.is_set():
        if end_time is not None and time.time() >= end_time:
            return

        if worker_id == 0 and args.guided_every and args.guided_every > 0:
            if (time.time() - last_guided) >= args.guided_every:
                res = guided_search(lib, rng, args)
                if res is not None:
                    queue.put(res)
                    found_flag.set()
                    return
                last_guided = time.time()

        if args.mode in ("mixed", "forest"):
            res = search_forest(
                lib,
                rng,
                args.forest_k_min,
                args.forest_k_max,
                args.batch,
                end_time=end_time,
            )
            if res is not None:
                queue.put(res)
                found_flag.set()
                return

        if args.mode in ("mixed", "star"):
            res = search_star(
                lib,
                rng,
                args.star_k_min,
                args.star_k_max,
                args.batch,
                end_time=end_time,
            )
            if res is not None:
                queue.put(res)
                found_flag.set()
                return

        if args.mode in ("mixed", "double"):
            res = search_double_star(
                lib,
                rng,
                args.double_left_min,
                args.double_left_max,
                args.double_right_min,
                args.double_right_max,
                args.batch,
                end_time=end_time,
            )
            if res is not None:
                queue.put(res)
                found_flag.set()
                return

        tries += args.batch

        # lightweight progress report + checkpoint
        if time.time() - last_report >= args.checkpoint_every:
            # report some "wiggly" library stats
            top = lib[0]
            prog = {
                "worker_id": worker_id,
                "elapsed_sec": round(time.time() - t0, 2),
                "tries": tries,
                "mode": args.mode,
                "top_library_entry": {
                    "name": top.name,
                    "n": top.n,
                    "alpha": top.alpha,
                    "peak": top.peak,
                    "lc_breaks": top.lc_breaks,
                },
            }
            outp = os.path.join(args.out_dir, f"progress_worker{worker_id}.json")
            os.makedirs(args.out_dir, exist_ok=True)
            with open(outp, "w", encoding="utf-8") as f:
                json.dump(prog, f, indent=2)
            last_report = time.time()


# -------------------------
# Main
# -------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["mixed", "forest", "star", "double"], default="mixed",
                    help="Search mode: forests only, star-of-subtrees only, double-star only, or mixed.")
    ap.add_argument("--library-size", type=int, default=2000, help="Number of library trees to build.")
    ap.add_argument("--max-n", type=int, default=260, help="Max vertex count for any library tree.")
    ap.add_argument("--library-file", type=str, default="library.pkl", help="Cache file for library.")
    ap.add_argument("--seed", type=int, default=0, help="Base RNG seed.")
    ap.add_argument("--workers", type=int, default=0, help="Number of worker processes (0 = auto).")
    ap.add_argument("--time-limit", type=float, default=None, help="Time limit in seconds (default: run until found).")
    ap.add_argument("--batch", type=int, default=200, help="How many random candidates per inner batch.")
    ap.add_argument("--checkpoint-every", type=float, default=300.0, help="Seconds between progress checkpoints per worker.")
    ap.add_argument("--guided-every", type=float, default=900.0, help="Seconds between guided beam-search passes (0 to disable).")
    ap.add_argument("--guided-pool", type=int, default=200, help="Library subset size for guided search.")
    ap.add_argument("--guided-beam", type=int, default=40, help="Beam width for guided search.")
    ap.add_argument("--guided-once", action="store_true", help="Run one guided search pass and exit.")
    ap.add_argument("--out-dir", type=str, default="out_erdos993", help="Directory for output files.")
    ap.add_argument("--forest-k-min", type=int, default=2)
    ap.add_argument("--forest-k-max", type=int, default=5)
    ap.add_argument("--star-k-min", type=int, default=3)
    ap.add_argument("--star-k-max", type=int, default=10)
    ap.add_argument("--double-left-min", type=int, default=2)
    ap.add_argument("--double-left-max", type=int, default=5)
    ap.add_argument("--double-right-min", type=int, default=2)
    ap.add_argument("--double-right-max", type=int, default=5)
    ap.add_argument("--galvin-sweep", action="store_true", help="Run deterministic Galvin-style sweep and exit.")
    ap.add_argument("--sweep-max-n", type=int, default=None, help="Max n for the Galvin sweep (default: --max-n).")
    ap.add_argument("--sweep-d-max", type=int, default=80, help="Max d in branching [d, e] + [1]*t.")
    ap.add_argument("--sweep-e-max", type=int, default=30, help="Max e in branching [d, e] + [1]*t.")
    ap.add_argument("--sweep-t-max", type=int, default=10, help="Max unary tail length t in branching [d, e] + [1]*t.")
    ap.add_argument("--sweep-top-k", type=int, default=200, help="Top K sequences to record in sweep report.")
    ap.add_argument("--sweep-out", type=str, default="galvin_sweep.json", help="Output JSON filename for sweep report.")
    ap.add_argument("--sweep-report", type=str, default=None, help="Optional sweep report to seed the library.")
    ap.add_argument("--sweep-include-top", type=int, default=100, help="Max entries to import from sweep report.")
    ap.add_argument("--sweep-only", action="store_true", help="Use only sweep-report entries as the library.")

    args = ap.parse_args()

    if args.sweep_max_n is None:
        args.sweep_max_n = args.max_n

    if args.galvin_sweep:
        os.makedirs(args.out_dir, exist_ok=True)
        report = galvin_sweep(
            max_n=args.sweep_max_n,
            d_max=args.sweep_d_max,
            e_max=args.sweep_e_max,
            t_max=args.sweep_t_max,
            top_k=args.sweep_top_k,
        )
        out_path = os.path.join(args.out_dir, args.sweep_out)
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2)
        print(f"[sweep] checked={report['checked']} wrote {out_path}")
        if report.get("first_counterexample") is not None:
            ce = report["first_counterexample"]
            meta = {
                "platform": sys.platform,
                "python": sys.version,
                "numpy": _HAVE_NUMPY,
                "sweep": {
                    "max_n": args.sweep_max_n,
                    "d_max": args.sweep_d_max,
                    "e_max": args.sweep_e_max,
                    "t_max": args.sweep_t_max,
                },
            }
            cert_path = write_tree_certificate(
                n=ce["n"],
                edges=ce["edges"],
                attach=ce["attach"],
                poly=ce["sequence"],
                out_dir=args.out_dir,
                seed=args.seed,
                meta=meta,
                prefix="certificate_tree_sweep",
            )
            print(f"[sweep] COUNTEREXAMPLE FOUND. Certificate: {cert_path}")
        return

    # build/load library (unless sweep-only)
    lib_path = os.path.join(args.out_dir, args.library_file)
    os.makedirs(args.out_dir, exist_ok=True)

    # optional: load sweep entries (used for seeding or sweep-only mode)
    sweep_path = args.sweep_report
    if sweep_path is None:
        sweep_path = os.path.join(args.out_dir, "galvin_sweep.json")
    sweep_entries = load_sweep_report(sweep_path, args.max_n, args.sweep_include_top)

    lib_dirty = False

    if args.sweep_only:
        if not sweep_entries:
            print(f"[main] ERROR: --sweep-only set but no sweep entries found at {sweep_path}")
            sys.exit(1)
        lib = sweep_entries
        lib.sort(key=lambda e: (e.lc_breaks, e.lc_gap, e.n), reverse=True)
        print(f"[main] Using sweep-only library with {len(lib)} entries from {sweep_path}")
        lib_dirty = True
    else:
        if os.path.exists(lib_path):
            print(f"[main] Loading library from {lib_path}")
            with open(lib_path, "rb") as f:
                lib = pickle.load(f)
        else:
            print(f"[main] Building library (size={args.library_size}, max_n={args.max_n}) ...")
            lib = build_library(args.library_size, args.max_n, args.seed, verbose=True)
            with open(lib_path, "wb") as f:
                pickle.dump(lib, f)
            print(f"[main] Saved library to {lib_path}")

        if sweep_entries:
            existing = {e.name for e in lib}
            added = [e for e in sweep_entries if e.name not in existing]
            if added:
                lib.extend(added)
                lib.sort(key=lambda e: (e.lc_breaks, e.lc_gap, e.n), reverse=True)
                print(f"[main] Added {len(added)} entries from sweep report {sweep_path}")
                lib_dirty = True

    if lib_dirty or not os.path.exists(lib_path):
        with open(lib_path, "wb") as f:
            pickle.dump(lib, f)
        if lib_dirty:
            print(f"[main] Saved updated library to {lib_path}")

    if args.guided_once:
        rng = random.Random(args.seed)
        res = guided_search(lib, rng, args)
        if res is None:
            print("[guided] No counterexample found in guided pass.")
            return
        meta = {
            "platform": sys.platform,
            "python": sys.version,
            "numpy": _HAVE_NUMPY,
            "guided_once": True,
        }
        out_path = verify_and_write_certificate(res, args.out_dir, args.seed, meta)
        print(f"[guided] COUNTEREXAMPLE FOUND. Certificate: {out_path}")
        return

    # immediate win?
    bad = [e for e in lib if not e.unimodal]
    if bad:
        e = bad[0]
        print("[main] Found a non-unimodal TREE already in the library build.")
        print(f"        name={e.name} n={e.n} alpha={e.alpha} peak={e.peak}")
        meta = {
            "platform": sys.platform,
            "python": sys.version,
            "numpy": _HAVE_NUMPY,
        }
        out_path = write_tree_certificate(
            n=e.n,
            edges=e.edges,
            attach=e.attach,
            poly=e.h,
            out_dir=args.out_dir,
            seed=args.seed,
            meta=meta,
            prefix="certificate_tree",
        )
        print(f"[main] Wrote certificate to {out_path}")
        return

    # number of workers
    if args.workers <= 0:
        try:
            import multiprocessing as mp
            args.workers = max(1, mp.cpu_count() - 1)
        except Exception:
            args.workers = 1

    print(f"[main] Starting search: mode={args.mode}, workers={args.workers}, numpy={_HAVE_NUMPY}")
    print(f"[main] Forest components k in [{args.forest_k_min},{args.forest_k_max}]")
    print(f"[main] Star branches    k in [{args.star_k_min},{args.star_k_max}]")
    if args.mode in ("mixed", "double"):
        print(
            f"[main] Double-star left k in [{args.double_left_min},{args.double_left_max}], "
            f"right k in [{args.double_right_min},{args.double_right_max}]"
        )
    if args.guided_every and args.guided_every > 0:
        print(
            f"[main] Guided beam search every {args.guided_every:.0f}s "
            f"(pool={args.guided_pool}, beam={args.guided_beam})"
        )
    else:
        print("[main] Guided beam search disabled")

    import multiprocessing as mp
    found_flag = mp.Event()
    queue = mp.Queue()

    # spawn workers
    procs = []
    for wid in range(args.workers):
        p = mp.Process(target=worker_main, args=(wid, args, lib_path, found_flag, queue), daemon=True)
        p.start()
        procs.append(p)

    try:
        # wait for a result
        while True:
            try:
                res = queue.get(timeout=1.0)
            except Exception:
                # time limit?
                if args.time_limit is not None:
                    # If all workers are dead, stop.
                    if all((not p.is_alive()) for p in procs):
                        print("[main] Time limit reached or workers stopped; no counterexample found.")
                        break
                continue

            # got a candidate
            found_flag.set()
            meta = {
                "platform": sys.platform,
                "python": sys.version,
                "numpy": _HAVE_NUMPY,
                "mode": args.mode,
                "library_size": args.library_size,
                "max_n": args.max_n,
            }
            out_path = verify_and_write_certificate(res, args.out_dir, args.seed, meta)
            print("[main] COUNTEREXAMPLE FOUND!")
            print(f"[main] Certificate written to: {out_path}")
            break

    except KeyboardInterrupt:
        print("\n[main] Interrupted by user.")

    finally:
        found_flag.set()
        for p in procs:
            try:
                p.join(timeout=1.0)
            except Exception:
                pass


if __name__ == "__main__":
    main()
