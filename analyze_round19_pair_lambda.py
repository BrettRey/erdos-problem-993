#!/usr/bin/env python3
"""
Compute per-pair lambda* bounds for shifted reserve in Round 19 regime.

Condition:
    sum_err <= D + lambda * (C10 + C01 + C11)

Scope:
- exhaustive trees via geng
- support-root processing
- child order by rooted subtree size (then id)
- step t >= 2, prefix k < mode(I_new)
- collect only X_k < 0 cases
- report pairwise stats for step-2 pairs (a,b)
"""

from __future__ import annotations

import argparse
import json
import subprocess
from collections import defaultdict
from dataclasses import dataclass

GENG = "/opt/homebrew/bin/geng"


def poly_mul(a: list[int], b: list[int]) -> list[int]:
    c = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            c[i + j] += ai * bj
    return c


def coeff(p: list[int], i: int) -> int:
    if 0 <= i < len(p):
        return p[i]
    return 0


def mode_smallest(p: list[int]) -> int:
    m = max(p)
    for i, v in enumerate(p):
        if v == m:
            return i
    return 0


def graph6_to_adj(s: str) -> list[list[int]]:
    s = s.strip()
    n = ord(s[0]) - 63
    bits = []
    for ch in s[1:]:
        v = ord(ch) - 63
        for b in range(5, -1, -1):
            bits.append((v >> b) & 1)
    adj = [[] for _ in range(n)]
    k = 0
    for j in range(1, n):
        for i in range(j):
            if k < len(bits) and bits[k]:
                adj[i].append(j)
                adj[j].append(i)
            k += 1
    return adj


def tree_dp_with_sizes(adj: list[list[int]], root: int):
    n = len(adj)
    children = [[] for _ in range(n)]
    order = []
    vis = [False] * n
    st = [root]
    vis[root] = True
    while st:
        v = st.pop()
        order.append(v)
        for u in adj[v]:
            if not vis[u]:
                vis[u] = True
                children[v].append(u)
                st.append(u)

    I = [None] * n
    E = [None] * n
    J = [None] * n
    subsize = [0] * n
    for v in reversed(order):
        if not children[v]:
            E[v] = [1]
            J[v] = [1]
            I[v] = [1, 1]
            subsize[v] = 1
            continue
        ev = [1]
        jv = [1]
        s = 1
        for c in children[v]:
            ev = poly_mul(ev, I[c])
            jv = poly_mul(jv, E[c])
            s += subsize[c]
        E[v] = ev
        J[v] = jv
        deg = max(len(ev), len(jv) + 1)
        iv = [0] * deg
        for i, val in enumerate(ev):
            iv[i] += val
        for i, val in enumerate(jv):
            iv[i + 1] += val
        I[v] = iv
        subsize[v] = s
    return I, E, J, children, subsize


@dataclass
class PairAgg:
    count: int = 0
    lambda_star: float = 0.0
    witness: dict | None = None


def update_lambda(agg: PairAgg, need: int, reserve: int, wit: dict):
    if need <= 0:
        return
    if reserve <= 0:
        cand = float("inf")
    else:
        cand = need / reserve
    if cand > agg.lambda_star:
        agg.lambda_star = cand
        agg.witness = wit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=18)
    ap.add_argument(
        "--out-json",
        type=str,
        default="results/round19_pair_lambda_profile_n18.json",
    )
    args = ap.parse_args()

    total_xneg = 0
    total_xneg_step2 = 0
    global_lambda = 0.0
    global_witness = None
    pair_aggs: dict[tuple[int, int], PairAgg] = defaultdict(PairAgg)

    for n in range(3, args.max_n + 1):
        g6_lines = subprocess.run(
            [GENG, "-cq", str(n), f"{n-1}:{n-1}"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip().splitlines()

        n_xneg = 0
        for g6 in g6_lines:
            adj = graph6_to_adj(g6)
            for root in range(n):
                if not any(len(adj[nb]) == 1 for nb in adj[root]):
                    continue
                I, E, J, children, subsize = tree_dp_with_sizes(adj, root)

                leaf_kids = []
                nonleaf_kids = []
                for c in children[root]:
                    if not children[c]:
                        leaf_kids.append(c)
                    else:
                        nonleaf_kids.append(c)
                if len(nonleaf_kids) < 2:
                    continue

                nonleaf_kids.sort(key=lambda c: (subsize[c], c))

                e_acc = [1]
                for _ in range(len(leaf_kids)):
                    e_acc = poly_mul(e_acc, [1, 1])
                j_acc = [1]

                for t_idx, c in enumerate(nonleaf_kids):
                    P = I[c]
                    Q = E[c]
                    e_new = poly_mul(e_acc, P)
                    j_new = poly_mul(j_acc, Q)
                    step = t_idx + 1
                    if step < 2:
                        e_acc = e_new
                        j_acc = j_new
                        continue

                    A = e_acc
                    B = j_acc
                    max_new = max(len(e_new), len(j_new))
                    max_ab = max(len(A), len(B))

                    lam_old = [
                        coeff(A, k) * coeff(B, k) - coeff(A, k - 1) * coeff(B, k + 1)
                        for k in range(max_ab + 1)
                    ]
                    D = [
                        sum(lam_old[i] * coeff(P, k - i) * coeff(Q, k - i) for i in range(len(lam_old)))
                        for k in range(max_new + 1)
                    ]
                    lam_new = [
                        coeff(e_new, k) * coeff(j_new, k) - coeff(e_new, k - 1) * coeff(j_new, k + 1)
                        for k in range(max_new + 1)
                    ]

                    I_new = [0] * max(len(e_new), len(j_new) + 1)
                    for i, v in enumerate(e_new):
                        I_new[i] += v
                    for i, v in enumerate(j_new):
                        I_new[i + 1] += v
                    mode = mode_smallest(I_new)

                    lo = -1
                    hi = max(len(A), len(B))
                    for k in range(min(mode, max_new + 1)):
                        X = lam_new[k] - D[k]
                        if X >= 0:
                            continue

                        total_xneg += 1
                        n_xneg += 1
                        if step == 2:
                            total_xneg_step2 += 1

                        sum_err = 0
                        for s in range(2 * lo, 2 * hi + 1):
                            vals_u = []
                            vals_D = []
                            for i in range(lo, hi + 1):
                                u = coeff(A, i) * coeff(B, s - i)
                                wi = coeff(P, k - i) * coeff(Q, k - s + i)
                                wip1 = coeff(P, k - (i + 1)) * coeff(Q, k - s + (i + 1))
                                vals_u.append(u)
                                vals_D.append(wi - wip1)
                            m_idx = None
                            for off, d_i in enumerate(vals_D):
                                if d_i >= 0:
                                    m_idx = off
                                    break
                            if m_idx is None:
                                continue
                            u_m = vals_u[m_idx]
                            for off in range(m_idx, len(vals_D)):
                                d_i = vals_D[off]
                                if d_i <= 0:
                                    continue
                                u_i = vals_u[off]
                                if u_i < u_m:
                                    sum_err += (u_m - u_i) * d_i

                        C10 = sum(
                            lam_old[i] * coeff(P, k - i - 1) * coeff(Q, k - i)
                            for i in range(len(lam_old))
                        )
                        C01 = sum(
                            lam_old[i] * coeff(P, k - i) * coeff(Q, k - i - 1)
                            for i in range(len(lam_old))
                        )
                        C11 = sum(
                            lam_old[i] * coeff(P, k - i - 1) * coeff(Q, k - i - 1)
                            for i in range(len(lam_old))
                        )
                        reserve = C10 + C01 + C11
                        need = max(0, sum_err - D[k])

                        wit = {
                            "n": n,
                            "g6": g6,
                            "root": root,
                            "step": step,
                            "k": k,
                            "mode": mode,
                            "X": X,
                            "D": D[k],
                            "sum_err": sum_err,
                            "need": need,
                            "C10": C10,
                            "C01": C01,
                            "C11": C11,
                            "reserve": reserve,
                        }

                        if reserve <= 0 and need > 0:
                            global_lambda = float("inf")
                            global_witness = wit
                        elif reserve > 0:
                            cand = need / reserve
                            if cand > global_lambda:
                                global_lambda = cand
                                global_witness = wit

                        if step == 2:
                            a = subsize[nonleaf_kids[0]]
                            b = subsize[nonleaf_kids[1]]
                            key = (a, b)
                            agg = pair_aggs[key]
                            agg.count += 1
                            update_lambda(agg, need, reserve, wit)

                    e_acc = e_new
                    j_acc = j_new

        print(f"n={n:2d}: X<0={n_xneg}", flush=True)

    pair_rows = []
    for (a, b), agg in sorted(pair_aggs.items(), key=lambda kv: (-kv[1].count, kv[0])):
        pair_rows.append(
            {
                "a": a,
                "b": b,
                "count": agg.count,
                "lambda_star": agg.lambda_star,
                "witness": agg.witness,
            }
        )

    out = {
        "scope": "support-step prefix X<0, n<=max_n",
        "max_n": args.max_n,
        "total_xneg": total_xneg,
        "total_xneg_step2": total_xneg_step2,
        "global_lambda_star": global_lambda,
        "global_witness": global_witness,
        "pair_rows": pair_rows,
    }

    with open(args.out_json, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, sort_keys=True)

    print(json.dumps(out, indent=2, sort_keys=True))
    print(f"\nWrote {args.out_json}")


if __name__ == "__main__":
    main()
