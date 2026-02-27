#!/usr/bin/env python3
"""Deterministic pair-AB attachment search (multiset subset).

Goal:
  Given fixed base trees A and B (graph6), enumerate bounded multisets of rooted
  d_leaf<=1 gadgets attached at the canonical root u (after deleting canonical
  {leaf,support}).

  For each attachment configuration, compute derived (m,lambda,rho) exactly.

  Early-abort on first true same-(m,lambda,rho) split with different N.
  Always emit a nearest-miss frontier: top-k smallest |rhoA-rhoB| among pairs
  with identical (m,lambda) and different N.

This script is self-contained (uses networkx for graph6 decoding and generating
nonisomorphic trees). It is designed to be drop-in in-repo.

Exact arithmetic:
  - lambda is computed as i_{m-1}/i_m (Fraction)
  - rho is Q(lambda)/P(lambda) (Fraction)
  - All rationals serialized as {"num":..,"den":..}

CLI (required by prompt):
  --gadget-max-n --max-gadgets --max-added-size --m-min --subset-size --seed --out

Optional:
  --g6A --g6B --top-k
"""

from __future__ import annotations

import argparse
import bisect
import heapq
import json
import os
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, List, Optional, Sequence, Tuple

import networkx as nx
from networkx.generators.nonisomorphic_trees import nonisomorphic_trees


def frac_to_numden(x: Fraction) -> dict:
    return {"num": int(x.numerator), "den": int(x.denominator)}


def numden_to_frac(d: dict) -> Fraction:
    return Fraction(int(d["num"]), int(d["den"]))


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


def poly_scale(a: List[int], s: int) -> List[int]:
    return [s * c for c in a]


def poly_shift_x(a: List[int], k: int = 1) -> List[int]:
    return [0] * k + a


def poly_eval(a: List[int], x: Fraction) -> Fraction:
    acc = Fraction(0, 1)
    powx = Fraction(1, 1)
    for c in a:
        if c:
            acc += c * powx
        powx *= x
    return acc


def compute_I_from_PQ(P: List[int], Q: List[int]) -> List[int]:
    # I = (1+2x)P + (1+x)Q
    termP = poly_add(P, poly_shift_x(poly_scale(P, 2), 1))
    termQ = poly_add(Q, poly_shift_x(Q, 1))
    return poly_add(termP, termQ)


def mode_lambda(I: List[int]) -> Tuple[int, Optional[Fraction]]:
    mx = max(I)
    m = I.index(mx)
    if m == 0:
        return m, None
    return m, Fraction(I[m - 1], I[m])


def d_leaf_from_adj(adj: List[set]) -> int:
    deg = [len(nb) for nb in adj]
    leaves = {i for i, d in enumerate(deg) if d == 1}
    maxc = 0
    for v in range(len(adj)):
        c = sum(1 for w in adj[v] if w in leaves)
        if c > maxc:
            maxc = c
    return maxc


def canonical_triplet_min_u(adj: List[set]) -> Optional[Tuple[int, int, int]]:
    deg = [len(nb) for nb in adj]
    leaves = [i for i, d in enumerate(deg) if d == 1]
    if not leaves:
        return None
    best_leaf = None
    best_key = None
    best_parent = None
    for leaf in leaves:
        parent = next(iter(adj[leaf]))
        key = (deg[parent], leaf)
        if best_key is None or key < best_key:
            best_key = key
            best_leaf = leaf
            best_parent = parent
    assert best_leaf is not None and best_parent is not None
    support = best_parent
    if deg[support] != 2:
        return None
    nbs = list(adj[support])
    u = nbs[0] if nbs[1] == best_leaf else nbs[1]
    return best_leaf, support, u


def delete_vertices(adj: List[set], to_delete: Sequence[int]) -> Tuple[List[set], Dict[int, int]]:
    delset = set(to_delete)
    mapping: Dict[int, int] = {}
    new_idx = 0
    for i in range(len(adj)):
        if i not in delset:
            mapping[i] = new_idx
            new_idx += 1
    new_adj = [set() for _ in range(new_idx)]
    for i in range(len(adj)):
        if i in delset:
            continue
        ni = mapping[i]
        for j in adj[i]:
            if j in delset:
                continue
            new_adj[ni].add(mapping[j])
    return new_adj, mapping


def rooted_parent_children(adj: List[set], root: int) -> Tuple[List[int], List[List[int]], List[int]]:
    n = len(adj)
    parent = [-1] * n
    order = [root]
    parent[root] = root
    for v in order:
        for w in adj[v]:
            if parent[w] == -1:
                parent[w] = v
                order.append(w)
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


def canonical_PQ_I_from_g6(g6: str) -> dict:
    G = nx.from_graph6_bytes(g6.encode("ascii"))
    n = G.number_of_nodes()
    adj = [set() for _ in range(n)]
    for u, v in G.edges():
        adj[u].add(v)
        adj[v].add(u)

    trip = canonical_triplet_min_u(adj)
    if trip is None:
        raise ValueError("canonical triplet does not exist")
    leaf, support, u = trip

    B_adj, mapping = delete_vertices(adj, [leaf, support])
    u_new = mapping[u]
    dp0, dp1 = rooted_dp_polys(B_adj, u_new)
    P = dp0[u_new]
    Q = dp1[u_new]
    I = compute_I_from_PQ(P, Q)

    m, lam = mode_lambda(I)
    if lam is None:
        raise ValueError("mode m=0 encountered")
    rho = poly_eval(Q, lam) / poly_eval(P, lam)

    return {
        "g6": g6,
        "triplet": {"leaf": int(leaf), "support": int(support), "u": int(u)},
        "P": P,
        "Q": Q,
        "I": I,
        "m": int(m),
        "lambda": frac_to_numden(lam),
        "rho": frac_to_numden(rho),
        "N": int(P[1] if len(P) > 1 else 0),
    }


@dataclass(frozen=True)
class Gadget:
    # index in the *full* sorted gadget list (so subset can be reproduced)
    gid: int
    n: int
    F: Tuple[int, ...]
    G: Tuple[int, ...]


def build_gadget_library(max_n: int) -> List[Gadget]:
    """All rooted gadget types (dedup by (n,F,G)) from nonisomorphic trees with d_leaf<=1."""
    raw = []

    # n=1
    if max_n >= 1:
        adj = [set()]
        dp0, dp1 = rooted_dp_polys(adj, 0)
        F = poly_add(dp0[0], dp1[0])
        G = dp0[0]
        raw.append((1, tuple(F), tuple(G)))

    for n in range(2, max_n + 1):
        for T in nonisomorphic_trees(n):
            # adjacency
            adj = [set() for _ in range(n)]
            for u, v in T.edges():
                adj[u].add(v)
                adj[v].add(u)
            if d_leaf_from_adj(adj) > 1:
                continue
            for root in range(n):
                dp0, dp1 = rooted_dp_polys(adj, root)
                F = poly_add(dp0[root], dp1[root])
                G = dp0[root]
                raw.append((n, tuple(F), tuple(G)))

    # dedup and sort
    raw = sorted(set(raw), key=lambda t: (t[0], t[1], t[2]))
    out: List[Gadget] = []
    for gid, (n, F, G) in enumerate(raw):
        out.append(Gadget(gid=gid, n=n, F=F, G=G))
    return out


def select_subset(gadgets: List[Gadget], subset_size: int, seed: int) -> List[Gadget]:
    import random

    rng = random.Random(seed)
    idxs = list(range(len(gadgets)))
    rng.shuffle(idxs)
    chosen = set(idxs[: min(subset_size, len(idxs))])
    return [gadgets[i] for i in sorted(chosen)]


def counts_from_seq(seq: List[int]) -> List[dict]:
    # seq entries are subset indices (not gids)
    counts: Dict[int, int] = {}
    for i in seq:
        counts[i] = counts.get(i, 0) + 1
    return [{"subset_index": int(i), "count": int(c)} for i, c in sorted(counts.items())]


def apply_config(P0: List[int], Q0: List[int], subset: List[Gadget], config: List[dict]) -> Tuple[List[int], List[int], int]:
    P = P0
    Q = Q0
    added = 0
    for item in config:
        i = int(item["subset_index"])
        c = int(item["count"])
        g = subset[i]
        for _ in range(c):
            P = poly_mul(P, list(g.F))
            Q = poly_mul(Q, list(g.G))
            added += g.n
    return P, Q, added


def scan_side_records(
    base: dict,
    subset: List[Gadget],
    max_gadgets: int,
    max_added_size: int,
    m_min: int,
) -> Tuple[Dict[Tuple[int, Fraction], List[Tuple[Fraction, dict]]], Dict[Tuple[int, Fraction, Fraction], Dict[int, dict]], dict]:
    """Enumerate all attachment configs for one base tree.

    Returns:
      by_mlam: (m,lam)-> sorted list of (rho, rec)
      by_key3: (m,lam,rho)-> {N: rec}
      totals: checked, kept

    rec contains: config (counts), added_size, N, m, lambda(Fraction), rho(Fraction)
    """

    P0: List[int] = base["P"]
    Q0: List[int] = base["Q"]
    baseN = int(base["N"])

    checked = 0
    kept = 0

    by_mlam: Dict[Tuple[int, Fraction], List[Tuple[Fraction, dict]]] = {}
    by_key3: Dict[Tuple[int, Fraction, Fraction], Dict[int, dict]] = {}

    # arrays for fast access
    g_ns = [g.n for g in subset]
    g_F = [list(g.F) for g in subset]
    g_G = [list(g.G) for g in subset]

    seq: List[int] = []

    def dfs(start_pos: int, depth: int, cur_P: List[int], cur_Q: List[int], added: int) -> None:
        nonlocal checked, kept

        checked += 1
        I = compute_I_from_PQ(cur_P, cur_Q)
        m, lam = mode_lambda(I)
        if lam is not None and m >= m_min:
            rho = poly_eval(cur_Q, lam) / poly_eval(cur_P, lam)
            N = baseN + added
            if len(cur_P) > 1 and cur_P[1] != N:
                N = cur_P[1]
            rec = {
                "config": counts_from_seq(seq),
                "added_size": int(added),
                "N": int(N),
                "m": int(m),
                "lambda": lam,
                "rho": rho,
            }
            kept += 1
            key2 = (m, lam)
            by_mlam.setdefault(key2, []).append((rho, rec))
            key3 = (m, lam, rho)
            by_key3.setdefault(key3, {})[N] = rec

        if depth == max_gadgets:
            return

        for i in range(start_pos, len(subset)):
            n_i = g_ns[i]
            if added + n_i > max_added_size:
                continue
            seq.append(i)
            dfs(i, depth + 1, poly_mul(cur_P, g_F[i]), poly_mul(cur_Q, g_G[i]), added + n_i)
            seq.pop()

    dfs(0, 0, P0, Q0, 0)

    for key2, lst in by_mlam.items():
        lst.sort(key=lambda t: t[0])

    return by_mlam, by_key3, {"checked": int(checked), "kept": int(kept)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gadget-max-n", type=int, required=True)
    ap.add_argument("--max-gadgets", type=int, required=True)
    ap.add_argument("--max-added-size", type=int, required=True)
    ap.add_argument("--m-min", type=int, required=True)
    ap.add_argument("--subset-size", type=int, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", required=True)

    ap.add_argument("--g6A", default="O??????_A?C?E?@_WG@j?")
    ap.add_argument("--g6B", default="P????A?OD?E?E?B??o?E?OO?")
    ap.add_argument("--top-k", type=int, default=50)

    args = ap.parse_args()

    cmd = " ".join([sys.executable] + sys.argv)
    t0 = time.perf_counter()

    baseA = canonical_PQ_I_from_g6(args.g6A)
    baseB = canonical_PQ_I_from_g6(args.g6B)

    gadgets = build_gadget_library(args.gadget_max_n)
    subset = select_subset(gadgets, args.subset_size, args.seed)

    # scan A fully
    A_mlam, A_key3, A_tot = scan_side_records(
        baseA, subset, args.max_gadgets, args.max_added_size, args.m_min
    )

    # Prepare A lookup per (m,lam): sorted rhos list and recs list
    A_lookup: Dict[Tuple[int, Fraction], Tuple[List[Fraction], List[dict]]] = {}
    for key2, lst in A_mlam.items():
        rhos = [x[0] for x in lst]
        recs = [x[1] for x in lst]
        A_lookup[key2] = (rhos, recs)

    # scan B with early-abort witness and keep nearest misses
    P0B: List[int] = baseB["P"]
    Q0B: List[int] = baseB["Q"]
    baseNB = int(baseB["N"])

    g_ns = [g.n for g in subset]
    g_F = [list(g.F) for g in subset]
    g_G = [list(g.G) for g in subset]

    checkedB = 0
    keptB = 0
    B_key2_seen = set()

    witness = None

    # heap keeps the *largest* diff at top via key=-diff
    heap: List[Tuple[Fraction, dict]] = []

    def push_miss(diff: Fraction, entry: dict) -> None:
        key = -diff
        if len(heap) < args.top_k:
            heapq.heappush(heap, (key, entry))
        else:
            # heap[0] is largest diff (most negative key)
            if key > heap[0][0]:
                heapq.heapreplace(heap, (key, entry))

    seq: List[int] = []

    def dfsB(start_pos: int, depth: int, cur_P: List[int], cur_Q: List[int], added: int) -> None:
        nonlocal checkedB, keptB, witness
        if witness is not None:
            return

        checkedB += 1
        I = compute_I_from_PQ(cur_P, cur_Q)
        m, lam = mode_lambda(I)
        if lam is not None and m >= args.m_min:
            rho = poly_eval(cur_Q, lam) / poly_eval(cur_P, lam)
            N = baseNB + added
            if len(cur_P) > 1 and cur_P[1] != N:
                N = cur_P[1]
            recB = {
                "config": counts_from_seq(seq),
                "added_size": int(added),
                "N": int(N),
                "m": int(m),
                "lambda": lam,
                "rho": rho,
            }
            keptB += 1

            key2 = (m, lam)
            B_key2_seen.add(key2)

            key3 = (m, lam, rho)
            if key3 in A_key3:
                for NA, recA in A_key3[key3].items():
                    if int(NA) != int(N):
                        # rebuild full witness records with P,Q,I
                        PA, QA, _ = apply_config(baseA["P"], baseA["Q"], subset, recA["config"])
                        IA = compute_I_from_PQ(PA, QA)
                        PB, QB, _ = apply_config(baseB["P"], baseB["Q"], subset, recB["config"])
                        IB = compute_I_from_PQ(PB, QB)
                        witness = {
                            "A": {
                                "g6": baseA["g6"],
                                "triplet": baseA["triplet"],
                                "config": recA["config"],
                                "P": PA,
                                "Q": QA,
                                "I": IA,
                                "m": int(m),
                                "lambda": frac_to_numden(lam),
                                "rho": frac_to_numden(rho),
                                "N": int(NA),
                            },
                            "B": {
                                "g6": baseB["g6"],
                                "triplet": baseB["triplet"],
                                "config": recB["config"],
                                "P": PB,
                                "Q": QB,
                                "I": IB,
                                "m": int(m),
                                "lambda": frac_to_numden(lam),
                                "rho": frac_to_numden(rho),
                                "N": int(N),
                            },
                        }
                        return

            # nearest miss among same (m,lam)
            if key2 in A_lookup:
                rhosA, recsA = A_lookup[key2]
                pos = bisect.bisect_left(rhosA, rho)
                cand_idxs = []
                if pos < len(rhosA):
                    cand_idxs.append(pos)
                if pos > 0:
                    cand_idxs.append(pos - 1)
                for j in cand_idxs:
                    rhoA = rhosA[j]
                    recA = recsA[j]
                    NA = int(recA["N"])
                    if NA == int(N):
                        continue
                    diff = abs(rhoA - rho)
                    push_miss(
                        diff,
                        {
                            "m": int(m),
                            "lambda": frac_to_numden(lam),
                            "rhoA": frac_to_numden(rhoA),
                            "rhoB": frac_to_numden(rho),
                            "abs_diff_rho": frac_to_numden(diff),
                            "NA": int(NA),
                            "NB": int(N),
                            "deltaN": int(N) - int(NA),
                            "configA": recA["config"],
                            "configB": recB["config"],
                        },
                    )

        if depth == args.max_gadgets:
            return
        for i in range(start_pos, len(subset)):
            n_i = g_ns[i]
            if added + n_i > args.max_added_size:
                continue
            seq.append(i)
            dfsB(i, depth + 1, poly_mul(cur_P, g_F[i]), poly_mul(cur_Q, g_G[i]), added + n_i)
            seq.pop()
            if witness is not None:
                return

    dfsB(0, 0, P0B, Q0B, 0)

    shared_keys = len(set(A_lookup.keys()) & B_key2_seen)

    misses = [(-k, e) for (k, e) in heap]  # convert key=-diff back to diff
    misses.sort(key=lambda t: t[0])
    frontier = [e for _, e in misses]

    out = {
        "command": cmd,
        "params": {
            "gadget_max_n": int(args.gadget_max_n),
            "max_gadgets": int(args.max_gadgets),
            "max_added_size": int(args.max_added_size),
            "m_min": int(args.m_min),
            "subset_size": int(args.subset_size),
            "seed": int(args.seed),
            "top_k": int(args.top_k),
            "exact_rationals": True,
            "early_abort_on_witness": True,
        },
        "base": {
            "A": baseA,
            "B": baseB,
        },
        "gadget_subset": [
            {"subset_index": i, "gid": g.gid, "n": g.n, "F": list(g.F), "G": list(g.G)}
            for i, g in enumerate(subset)
        ],
        "totals": {
            "checked_A": int(A_tot["checked"]),
            "kept_A": int(A_tot["kept"]),
            "checked_B": int(checkedB),
            "kept_B": int(keptB),
            "checked_total": int(A_tot["checked"] + checkedB),
            "shared_keys": int(shared_keys),
            "witness_found": bool(witness is not None),
        },
        "witness": witness,
        "nearest_miss_frontier": frontier,
        "runtime_s": float(time.perf_counter() - t0),
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
