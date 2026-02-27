#!/usr/bin/env python3
"""Pair-AB push-fast2 exhaustive bounded search (exact rationals).

Base pair:
  A = O??????_A?C?E?@_WG@j?
  B = P????A?OD?E?E?B??o?E?OO?

Search objective:
  - target m == m_target
  - target lambda == lambda_target
  - odd->even adjacent witness: same (m,lambda,rho) with deltaN = +1

Method:
  1) Build full rooted gadget universe (trees up to gadget_max_n) filtered by internal d_leaf<=1,
     dedup by (n,F,G).
  2) Enumerate all gadget multisets with:
       #gadgets <= max_gadgets, total added size <= max_added_size.
  3) Fast ratio gate per side using truncated degree-5 coefficients:
       require 9*i4 == 7*i5 for lambda_target = 7/9.
  4) Full check only on ratio-pass patterns:
       leftmost mode m == m_target and lambda == lambda_target.
  5) Collect A-good and B-good pattern records, then cross-compare to:
       - find first odd->even rho witness (deltaN=+1),
       - produce shared-key set and nearest-miss frontier (general deltaN).

No floating-point arithmetic is used.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import deque
from fractions import Fraction
from heapq import heappush, heappushpop
from typing import Dict, List, Optional, Sequence, Tuple

import networkx as nx
from networkx.generators.nonisomorphic_trees import nonisomorphic_trees


# -----------------------------
# Polynomial utilities
# -----------------------------

def poly_add(a: List[int], b: List[int]) -> List[int]:
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_mul(a: List[int], b: List[int]) -> List[int]:
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


def poly_mul_trunc(a: List[int], b: List[int], kmax: int) -> List[int]:
    out = [0] * (kmax + 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        if i > kmax:
            break
        for j, bj in enumerate(b):
            if bj == 0:
                continue
            ij = i + j
            if ij > kmax:
                break
            out[ij] += ai * bj
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_shift_x(a: List[int], k: int = 1) -> List[int]:
    return [0] * k + a


def poly_eval_frac(p: List[int], x: Fraction) -> Fraction:
    s = Fraction(0, 1)
    powx = Fraction(1, 1)
    for c in p:
        if c:
            s += Fraction(c, 1) * powx
        powx *= x
    return s


def compute_I_from_PQ(P: List[int], Q: List[int]) -> List[int]:
    # I = (1+2x)P + (1+x)Q
    P1 = poly_add(P, poly_shift_x([2 * c for c in P], 1))
    Q1 = poly_add(Q, poly_shift_x(Q, 1))
    return poly_add(P1, Q1)


def compute_I_from_PQ_trunc(P: List[int], Q: List[int], kmax: int) -> List[int]:
    # Truncated variant of I = (1+2x)P + (1+x)Q
    P2 = [2 * c for c in P]
    P1 = poly_add(P, poly_shift_x(P2, 1))
    Q1 = poly_add(Q, poly_shift_x(Q, 1))
    I = poly_add(P1, Q1)
    if len(I) > kmax + 1:
        I = I[: kmax + 1]
    while len(I) > 1 and I[-1] == 0:
        I.pop()
    return I


def leftmost_mode(coeffs: List[int]) -> int:
    mx = max(coeffs)
    for i, v in enumerate(coeffs):
        if v == mx:
            return i
    raise RuntimeError("unreachable")


# -----------------------------
# Graph / canonical DP helpers
# -----------------------------

def d_leaf_from_adj(adj: List[set]) -> int:
    deg = [len(nb) for nb in adj]
    leaves = {i for i, d in enumerate(deg) if d == 1}
    maxc = 0
    for v in range(len(adj)):
        c = sum(1 for w in adj[v] if w in leaves)
        if c > maxc:
            maxc = c
    return maxc


def admissible_triplets(adj: List[set]) -> List[Tuple[int, int, int]]:
    deg = [len(nb) for nb in adj]
    trips: List[Tuple[int, int, int]] = []
    for s in range(len(adj)):
        if deg[s] != 2:
            continue
        nbrs = list(adj[s])
        if len(nbrs) != 2:
            continue
        for l in nbrs:
            if deg[l] == 1 and adj[l] == {s}:
                u = nbrs[0] if nbrs[1] == l else nbrs[1]
                trips.append((l, s, u))
    return trips


def canonical_triplet_min_u(adj: List[set]) -> Optional[Tuple[int, int, int]]:
    trips = admissible_triplets(adj)
    if not trips:
        return None
    return sorted(trips, key=lambda t: (t[2], t[1], t[0]))[0]


def delete_vertices(adj: List[set], to_delete: Sequence[int]) -> Tuple[List[set], Dict[int, int]]:
    delset = set(to_delete)
    mapping: Dict[int, int] = {}
    nxt = 0
    for i in range(len(adj)):
        if i not in delset:
            mapping[i] = nxt
            nxt += 1
    out = [set() for _ in range(nxt)]
    for i in range(len(adj)):
        if i in delset:
            continue
        ni = mapping[i]
        for j in adj[i]:
            if j in delset:
                continue
            out[ni].add(mapping[j])
    return out, mapping


def rooted_parent_children(adj: List[set], root: int) -> Tuple[List[int], List[List[int]], List[int]]:
    n = len(adj)
    parent = [-1] * n
    parent[root] = root
    order = [root]
    q = deque([root])
    while q:
        v = q.popleft()
        for w in adj[v]:
            if parent[w] != -1:
                continue
            parent[w] = v
            order.append(w)
            q.append(w)
    children = [[] for _ in range(n)]
    for v in range(n):
        if v == root:
            continue
        p = parent[v]
        children[p].append(v)
    return parent, children, order


def rooted_dp_polys(adj: List[set], root: int) -> Tuple[List[List[int]], List[List[int]]]:
    _, children, order = rooted_parent_children(adj, root)
    dp0: List[Optional[List[int]]] = [None] * len(adj)
    dp1: List[Optional[List[int]]] = [None] * len(adj)
    for v in reversed(order):
        prod = [1]
        prod0 = [1]
        for c in children[v]:
            sumc = poly_add(dp0[c], dp1[c])  # type: ignore[arg-type]
            prod = poly_mul(prod, sumc)
            prod0 = poly_mul(prod0, dp0[c])  # type: ignore[arg-type]
        dp0[v] = prod
        dp1[v] = poly_shift_x(prod0, 1)
    return dp0, dp1  # type: ignore[return-value]


def canonical_record_from_g6(g6: str) -> dict:
    G = nx.from_graph6_bytes(g6.encode("ascii"))
    n = G.number_of_nodes()
    adj = [set() for _ in range(n)]
    for u, v in G.edges():
        adj[u].add(v)
        adj[v].add(u)

    trip = canonical_triplet_min_u(adj)
    if trip is None:
        raise ValueError("no canonical triplet")
    l, s, u = trip

    B_adj, mapping = delete_vertices(adj, [l, s])
    u_new = mapping[u]

    dp0, dp1 = rooted_dp_polys(B_adj, u_new)
    P = dp0[u_new]
    Q = dp1[u_new]
    I = compute_I_from_PQ(P, Q)
    m = leftmost_mode(I)
    if m == 0:
        raise ValueError("m=0 unsupported")
    lam = Fraction(I[m - 1], I[m])
    rho = poly_eval_frac(Q, lam) / poly_eval_frac(P, lam)

    return {
        "g6": g6,
        "canonical_triplet": [int(l), int(s), int(u)],
        "P": P,
        "Q": Q,
        "I": I,
        "m": int(m),
        "lambda": [int(lam.numerator), int(lam.denominator)],
        "rho": [int(rho.numerator), int(rho.denominator)],
        "N": int(P[1] if len(P) > 1 else 0),
    }


# -----------------------------
# Gadget universe
# -----------------------------

def build_gadget_universe(max_n: int) -> Tuple[List[dict], Dict[int, int], float]:
    t0 = time.perf_counter()
    raw = []
    hist: Dict[int, int] = {}

    for n in range(2, max_n + 1):
        local = set()
        for T in nonisomorphic_trees(n):
            adj = [set() for _ in range(n)]
            for u, v in T.edges():
                adj[u].add(v)
                adj[v].add(u)
            if d_leaf_from_adj(adj) > 1:
                continue
            for root in range(n):
                dp0, dp1 = rooted_dp_polys(adj, root)
                F = tuple(poly_add(dp0[root], dp1[root]))
                G = tuple(dp0[root])
                local.add((n, F, G))
        for x in local:
            raw.append(x)
        hist[n] = len(local)

    raw = sorted(set(raw), key=lambda t: (t[0], t[1], t[2]))
    gadgets = [
        {
            "gid": i,
            "n": n,
            "F": list(F),
            "G": list(G),
        }
        for i, (n, F, G) in enumerate(raw)
    ]
    t1 = time.perf_counter()
    return gadgets, hist, t1 - t0


# -----------------------------
# Main search
# -----------------------------

def parse_fraction(s: str) -> Fraction:
    if "/" in s:
        a, b = s.split("/", 1)
        return Fraction(int(a), int(b))
    return Fraction(int(s), 1)


def coeff_at(p: List[int], i: int) -> int:
    return p[i] if i < len(p) else 0


def cfg_counts(seq: List[int]) -> List[dict]:
    out: Dict[int, int] = {}
    for i in seq:
        out[i] = out.get(i, 0) + 1
    return [{"gid": int(k), "count": int(v)} for k, v in sorted(out.items())]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gadget-max-n", type=int, required=True)
    ap.add_argument("--max-gadgets", type=int, required=True)
    ap.add_argument("--max-added-size", type=int, required=True)
    ap.add_argument("--m-target", type=int, required=True)
    ap.add_argument("--lambda-target", type=str, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--g6A", default="O??????_A?C?E?@_WG@j?")
    ap.add_argument("--g6B", default="P????A?OD?E?E?B??o?E?OO?")
    ap.add_argument("--top-k", type=int, default=50)
    args = ap.parse_args()

    lambda_target = parse_fraction(args.lambda_target)
    if lambda_target <= 0:
        raise ValueError("lambda-target must be positive")

    cmd = " ".join([sys.executable] + sys.argv)

    baseA = canonical_record_from_g6(args.g6A)
    baseB = canonical_record_from_g6(args.g6B)

    # Build gadgets
    gadgets, size_hist, build_s = build_gadget_universe(args.gadget_max_n)

    # Fast arrays
    g_n = [g["n"] for g in gadgets]
    g_F = [g["F"] for g in gadgets]
    g_G = [g["G"] for g in gadgets]

    # Persistent collections
    goodA: List[dict] = []
    goodB: List[dict] = []

    patterns_total = 0
    ratio_pass_A = 0
    ratio_pass_B = 0

    fullcheck_s_A = 0.0
    fullcheck_s_B = 0.0

    enum_t0 = time.perf_counter()

    seq: List[int] = []

    # Initial full and truncated states (degree <=5)
    P_A0 = baseA["P"]
    Q_A0 = baseA["Q"]
    P_B0 = baseB["P"]
    Q_B0 = baseB["Q"]

    P_A0_t = P_A0[:6]
    Q_A0_t = Q_A0[:6]
    P_B0_t = P_B0[:6]
    Q_B0_t = Q_B0[:6]

    def maybe_fullcheck(side: str, P: List[int], Q: List[int], added: int, cfg: List[dict]) -> Optional[dict]:
        nonlocal fullcheck_s_A, fullcheck_s_B
        t0 = time.perf_counter()
        I = compute_I_from_PQ(P, Q)
        m = leftmost_mode(I)
        if m != args.m_target or m == 0:
            dt = time.perf_counter() - t0
            if side == "A":
                fullcheck_s_A += dt
            else:
                fullcheck_s_B += dt
            return None
        lam = Fraction(I[m - 1], I[m])
        if lam != lambda_target:
            dt = time.perf_counter() - t0
            if side == "A":
                fullcheck_s_A += dt
            else:
                fullcheck_s_B += dt
            return None
        rho = poly_eval_frac(Q, lambda_target) / poly_eval_frac(P, lambda_target)
        N = int(P[1] if len(P) > 1 else 0)
        rec = {
            "config": cfg,
            "added_size": int(added),
            "N": N,
            "m": int(m),
            "lambda": [int(lam.numerator), int(lam.denominator)],
            "rho": [int(rho.numerator), int(rho.denominator)],
            "P": P,
            "Q": Q,
            "I": I,
        }
        dt = time.perf_counter() - t0
        if side == "A":
            fullcheck_s_A += dt
        else:
            fullcheck_s_B += dt
        return rec

    def dfs(
        start_idx: int,
        depth: int,
        added_size: int,
        P_A: List[int], Q_A: List[int], P_A_t: List[int], Q_A_t: List[int],
        P_B: List[int], Q_B: List[int], P_B_t: List[int], Q_B_t: List[int],
    ) -> None:
        nonlocal patterns_total, ratio_pass_A, ratio_pass_B

        patterns_total += 1

        # Ratio gate (target 7/9): 9*i4 == 7*i5 on each side using truncated polynomials
        I_A_t = compute_I_from_PQ_trunc(P_A_t, Q_A_t, 5)
        i4A = coeff_at(I_A_t, 4)
        i5A = coeff_at(I_A_t, 5)
        passA = (i5A > 0 and 9 * i4A == 7 * i5A)

        I_B_t = compute_I_from_PQ_trunc(P_B_t, Q_B_t, 5)
        i4B = coeff_at(I_B_t, 4)
        i5B = coeff_at(I_B_t, 5)
        passB = (i5B > 0 and 9 * i4B == 7 * i5B)

        cfg = cfg_counts(seq)

        if passA:
            ratio_pass_A += 1
            recA = maybe_fullcheck("A", P_A, Q_A, added_size, cfg)
            if recA is not None:
                goodA.append(recA)

        if passB:
            ratio_pass_B += 1
            recB = maybe_fullcheck("B", P_B, Q_B, added_size, cfg)
            if recB is not None:
                goodB.append(recB)

        if depth == args.max_gadgets:
            return

        for i in range(start_idx, len(gadgets)):
            ni = g_n[i]
            if added_size + ni > args.max_added_size:
                continue
            seq.append(i)
            dfs(
                i,
                depth + 1,
                added_size + ni,
                poly_mul(P_A, g_F[i]),
                poly_mul(Q_A, g_G[i]),
                poly_mul_trunc(P_A_t, g_F[i], 5),
                poly_mul_trunc(Q_A_t, g_G[i], 5),
                poly_mul(P_B, g_F[i]),
                poly_mul(Q_B, g_G[i]),
                poly_mul_trunc(P_B_t, g_F[i], 5),
                poly_mul_trunc(Q_B_t, g_G[i], 5),
            )
            seq.pop()

    dfs(0, 0, 0, P_A0, Q_A0, P_A0_t, Q_A0_t, P_B0, Q_B0, P_B0_t, Q_B0_t)

    enum_s = time.perf_counter() - enum_t0

    # Build key maps by (m, lambda, rho)
    mapA: Dict[Tuple[int, Tuple[int, int], Tuple[int, int]], List[dict]] = {}
    mapB: Dict[Tuple[int, Tuple[int, int], Tuple[int, int]], List[dict]] = {}

    for r in goodA:
        key = (r["m"], tuple(r["lambda"]), tuple(r["rho"]))
        mapA.setdefault(key, []).append(r)
    for r in goodB:
        key = (r["m"], tuple(r["lambda"]), tuple(r["rho"]))
        mapB.setdefault(key, []).append(r)

    # Witness search + nearest miss frontier
    first_witness = None

    # shared keys in sense used by bounded certificate: (deltaN,m,lambda)
    shared_delta_keys = set()

    # nearest misses across all pairs with same (m,lambda), deltaN != 0
    # heap stores (-gap, entry) so we keep top_k smallest gaps
    heap: List[Tuple[Fraction, dict]] = []

    def push_frontier(gap: Fraction, entry: dict) -> None:
        key = -gap
        if len(heap) < args.top_k:
            heappush(heap, (key, entry))
        else:
            if key > heap[0][0]:
                heappushpop(heap, (key, entry))

    # quick grouping by (m,lambda)
    A_by_ml: Dict[Tuple[int, Tuple[int, int]], List[dict]] = {}
    B_by_ml: Dict[Tuple[int, Tuple[int, int]], List[dict]] = {}
    for r in goodA:
        A_by_ml.setdefault((r["m"], tuple(r["lambda"])), []).append(r)
    for r in goodB:
        B_by_ml.setdefault((r["m"], tuple(r["lambda"])), []).append(r)

    # Witness: same (m,lambda,rho) and deltaN=+1 odd->even
    for key in sorted(set(mapA.keys()) & set(mapB.keys())):
        m, lam, _rho = key
        for a in mapA[key]:
            NA = int(a["N"])
            for b in mapB[key]:
                NB = int(b["N"])
                dN = NB - NA
                shared_delta_keys.add((dN, m, lam))
                if dN == 1 and (NA % 2 == 1) and (NB % 2 == 0):
                    first_witness = {
                        "A": {
                            "g6": baseA["g6"],
                            "canonical_triplet": baseA["canonical_triplet"],
                            "config": a["config"],
                            "P": a["P"],
                            "Q": a["Q"],
                            "I": a["I"],
                            "m": m,
                            "lambda": list(lam),
                            "rho": list(_rho),
                            "N": NA,
                        },
                        "B": {
                            "g6": baseB["g6"],
                            "canonical_triplet": baseB["canonical_triplet"],
                            "config": b["config"],
                            "P": b["P"],
                            "Q": b["Q"],
                            "I": b["I"],
                            "m": m,
                            "lambda": list(lam),
                            "rho": list(_rho),
                            "N": NB,
                        },
                    }
                    break
            if first_witness is not None:
                break
        if first_witness is not None:
            break

    # General frontier on shared (m,lambda)
    for ml in sorted(set(A_by_ml.keys()) & set(B_by_ml.keys())):
        m, lam = ml
        A_list = A_by_ml[ml]
        B_list = B_by_ml[ml]
        for a in A_list:
            NA = int(a["N"])
            rhoA = Fraction(a["rho"][0], a["rho"][1])
            for b in B_list:
                NB = int(b["N"])
                dN = NB - NA
                if dN == 0:
                    continue
                rhoB = Fraction(b["rho"][0], b["rho"][1])
                gap = abs(rhoA - rhoB)
                entry = {
                    "deltaN": int(dN),
                    "m": int(m),
                    "lambda": [int(lam[0]), int(lam[1])],
                    "rhoA": [int(rhoA.numerator), int(rhoA.denominator)],
                    "rhoB": [int(rhoB.numerator), int(rhoB.denominator)],
                    "rho_gap": [int(gap.numerator), int(gap.denominator)],
                    "wA": int(a["added_size"]),
                    "wB": int(b["added_size"]),
                }
                push_frontier(gap, entry)

    top = [(-k, e) for (k, e) in heap]
    top.sort(key=lambda t: t[0])
    top_entries = [e for _, e in top]

    out = {
        "script_path": os.path.abspath(__file__),
        "command": cmd,
        "witness_found": bool(first_witness is not None),
        "first_witness": first_witness,
        "bounded_no_go_certificate": {
            "base_pair": {
                "A": {
                    "g6": baseA["g6"],
                    "canonical_triplet": baseA["canonical_triplet"],
                    "P": baseA["P"],
                    "Q": baseA["Q"],
                    "I": baseA["I"],
                    "m": baseA["m"],
                    "lambda": baseA["lambda"],
                    "N": baseA["N"],
                },
                "B": {
                    "g6": baseB["g6"],
                    "canonical_triplet": baseB["canonical_triplet"],
                    "P": baseB["P"],
                    "Q": baseB["Q"],
                    "I": baseB["I"],
                    "m": baseB["m"],
                    "lambda": baseB["lambda"],
                    "N": baseB["N"],
                },
                "base_deltaN": int(baseB["N"] - baseA["N"]),
            },
            "bounds": {
                "gadget_max_n": int(args.gadget_max_n),
                "max_gadgets": int(args.max_gadgets),
                "max_added_size": int(args.max_added_size),
                "m_target": int(args.m_target),
                "lambda_target": [int(lambda_target.numerator), int(lambda_target.denominator)],
                "adjacent_target": {"deltaN": 1, "odd_even": True},
            },
            "gadget_universe": {
                "generation_rule": {
                    "trees": "all nonisomorphic trees",
                    "n_range": [2, int(args.gadget_max_n)],
                    "filter": "internal d_leaf<=1",
                    "rooting": "all roots, then deduplicate by (F,G,n)",
                },
                "universe_size": int(len(gadgets)),
                "size_hist": {str(k): int(v) for k, v in sorted(size_hist.items())},
            },
            "search_space": {
                "patterns_enumerated_total": int(patterns_total),
                "patterns_passing_ratio_check_A": int(ratio_pass_A),
                "patterns_passing_ratio_check_B": int(ratio_pass_B),
                "good_patterns_A_count": int(len(goodA)),
                "good_patterns_B_count": int(len(goodB)),
                "pruning_rules": [
                    "multisets only (nondecreasing gadget id)",
                    "total added size <= max_added_size",
                    "#gadgets <= max_gadgets",
                    "gadgets restricted to internal d_leaf<=1 rooted tree types (n>=2)",
                    "ratio gate (side-specific): require 9*i4 == 7*i5 (lambda=7/9)",
                    "full check: require leftmost mode m == m_target and lambda == lambda_target",
                ],
                "#gadgets <= max_gadgets": None,
                "total added size <= max_added_size": None,
            },
            "shared_key_set": [
                [int(dN), int(m), [int(lam[0]), int(lam[1])]]
                for (dN, m, lam) in sorted(shared_delta_keys)
            ],
            "timing": {
                "build_gadgets_s": float(build_s),
                "enumeration_s": float(enum_s),
                "fullcheck_s_A": float(fullcheck_s_A),
                "fullcheck_s_B": float(fullcheck_s_B),
                "bottleneck": "enumeration" if enum_s >= max(fullcheck_s_A, fullcheck_s_B, build_s) else "fullcheck",
            },
        },
        "frontier": {
            "all_keys_count": int(len(shared_delta_keys)),
            "shared_key_set": [
                [int(dN), int(m), [int(lam[0]), int(lam[1])]]
                for (dN, m, lam) in sorted(shared_delta_keys)
            ],
            "top_nearest_misses": top_entries,
        },
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
