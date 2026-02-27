#!/usr/bin/env python3

"""Scan all unlabeled trees for n in [min_n,max_n] for canonical lock (m,lambda).

Canonicalization:
  - require d_leaf<=1
  - admissible triplets (ℓ,s,u): ℓ leaf, deg(s)=2, N(s)={ℓ,u}
  - choose min (u,s,ℓ) in current labeling (min-u)
  - delete {ℓ,s} and run rooted DP at u to get P=dp0[u], Q=dp1[u]
  - I=(1+2x)P+(1+x)Q, m=leftmost mode, lambda=i_{m-1}/i_m, rho=Q(lambda)/P(lambda)

Outputs JSON with per-n counts and lists of matches.
Also checks whether any rho value appears at two different n in the scanned range.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from collections import deque, defaultdict

import networkx as nx
from networkx.generators.nonisomorphic_trees import nonisomorphic_trees


def poly_add(a, b):
    n = max(len(a), len(b))
    r = [0] * n
    for i in range(n):
        if i < len(a):
            r[i] += a[i]
        if i < len(b):
            r[i] += b[i]
    while len(r) > 1 and r[-1] == 0:
        r.pop()
    return r


def poly_mul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0:
                continue
            r[i + j] += ai * bj
    while len(r) > 1 and r[-1] == 0:
        r.pop()
    return r


def poly_shift_x(a):
    return [0] + a


def poly_eval_frac(p, x: Fraction) -> Fraction:
    s = Fraction(0, 1)
    powx = Fraction(1, 1)
    for c in p:
        s += Fraction(c, 1) * powx
        powx *= x
    return s


def is_dleaf_le_1(G: nx.Graph) -> bool:
    deg = dict(G.degree())
    leaves = {v for v, d in deg.items() if d == 1}
    for v in G.nodes():
        if sum(1 for w in G.neighbors(v) if w in leaves) > 1:
            return False
    return True


def admissible_triplets(G: nx.Graph):
    deg = dict(G.degree())
    trips = []
    for s in G.nodes():
        if deg[s] != 2:
            continue
        nbrs = list(G.neighbors(s))
        for l in nbrs:
            if deg[l] == 1 and set(G.neighbors(l)) == {s}:
                u = nbrs[0] if nbrs[1] == l else nbrs[1]
                trips.append((l, s, u))
    return trips


def canonical_triplet_min_u(G: nx.Graph):
    trips = admissible_triplets(G)
    if not trips:
        return None
    return sorted(trips, key=lambda t: (t[2], t[1], t[0]))[0]


def dp_root_polys(G: nx.Graph, root: int):
    parent = {root: None}
    order = [root]
    dq = deque([root])
    while dq:
        v = dq.popleft()
        for w in G.neighbors(v):
            if w == parent.get(v):
                continue
            if w in parent:
                continue
            parent[w] = v
            order.append(w)
            dq.append(w)

    children = {v: [] for v in order}
    for v in order[1:]:
        children[parent[v]].append(v)

    dp0 = {}
    dp1 = {}
    for v in reversed(order):
        prod0 = [1]
        prodS = [1]
        for c in children[v]:
            prodS = poly_mul(prodS, poly_add(dp0[c], dp1[c]))
            prod0 = poly_mul(prod0, dp0[c])
        dp0[v] = prodS
        dp1[v] = poly_shift_x(prod0)
    return dp0[root], dp1[root]


def canonical_PQ(G: nx.Graph):
    trip = canonical_triplet_min_u(G)
    if trip is None:
        return None
    l, s, u = trip
    B = G.copy()
    B.remove_nodes_from([l, s])
    P, Q = dp_root_polys(B, u)
    return trip, P, Q


def compute_I(P, Q):
    # (1+2x)P + (1+x)Q
    P1 = poly_add(P, [0] + [2 * c for c in P])
    Q1 = poly_add(Q, [0] + Q)
    return poly_add(P1, Q1)


def leftmost_mode(coeffs):
    mx = max(coeffs)
    for k, v in enumerate(coeffs):
        if v == mx:
            return k
    raise RuntimeError("no mode")


def frac_pair(f: Fraction):
    return [f.numerator, f.denominator]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, required=True)
    ap.add_argument("--max-n", type=int, required=True)
    ap.add_argument("--m-target", type=int, required=True)
    ap.add_argument("--lambda-target", type=str, required=True)
    ap.add_argument("--out", type=str, required=True)
    args = ap.parse_args()

    ln, ld = args.lambda_target.split("/")
    lam_target = Fraction(int(ln), int(ld))

    G_A = nx.from_graph6_bytes(b"O??????_A?C?E?@_WG@j?")
    G_B = nx.from_graph6_bytes(b"P????A?OD?E?E?B??o?E?OO?")

    out = {
        "params": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "m_target": args.m_target,
            "lambda_target": frac_pair(lam_target),
        },
        "per_n": {},
        "rho_cross_matches": [],
        "verdict": None,
    }

    rho_by_n = defaultdict(lambda: defaultdict(list))

    for n in range(args.min_n, args.max_n + 1):
        total = 0
        gated = 0
        lock_ok = 0
        matches = []

        for T in nonisomorphic_trees(n):
            total += 1
            if not is_dleaf_le_1(T):
                continue
            cpq = canonical_PQ(T)
            if cpq is None:
                continue
            trip, P, Q = cpq
            gated += 1

            I = compute_I(P, Q)
            m = leftmost_mode(I)
            if m != args.m_target:
                continue
            lam = Fraction(I[m - 1], I[m])
            if lam != lam_target:
                continue

            lock_ok += 1
            rho = poly_eval_frac(Q, lam) / poly_eval_frac(P, lam)
            g6 = nx.to_graph6_bytes(T, header=False).decode().strip()
            rec = {
                "n": n,
                "N": (P[1] if len(P) > 1 else 0),
                "g6": g6,
                "canonical_triplet": list(trip),
                "m": m,
                "lambda": frac_pair(lam),
                "rho": frac_pair(rho),
                "isomorphic_to_A": nx.is_isomorphic(T, G_A),
                "isomorphic_to_B": nx.is_isomorphic(T, G_B),
            }
            matches.append(rec)
            rho_by_n[n][rho].append(rec)

        out["per_n"][str(n)] = {
            "total_trees": total,
            "gated_ok": gated,
            "lock_ok": lock_ok,
            "matches": matches,
        }

    # cross-match rho across different n
    ns = list(range(args.min_n, args.max_n + 1))
    for i in range(len(ns)):
        for j in range(i + 1, len(ns)):
            n1, n2 = ns[i], ns[j]
            common = set(rho_by_n[n1].keys()) & set(rho_by_n[n2].keys())
            for rho in sorted(common):
                r1 = rho_by_n[n1][rho][0]
                r2 = rho_by_n[n2][rho][0]
                out["rho_cross_matches"].append({
                    "n1": n1,
                    "n2": n2,
                    "rho": frac_pair(rho),
                    "tree1": {k: r1[k] for k in ["g6", "N", "isomorphic_to_A", "isomorphic_to_B"]},
                    "tree2": {k: r2[k] for k in ["g6", "N", "isomorphic_to_A", "isomorphic_to_B"]},
                })

    out["verdict"] = "outside-model witness found" if out["rho_cross_matches"] else "outside-model witness not found in tested bounds"

    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)

    print(json.dumps({
        "verdict": out["verdict"],
        "out": args.out,
        "per_n_lock_ok": {n: out["per_n"][str(n)]["lock_ok"] for n in ns},
        "rho_cross_matches": len(out["rho_cross_matches"]),
    }, indent=2))


if __name__ == "__main__":
    main()
