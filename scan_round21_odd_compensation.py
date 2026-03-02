#!/usr/bin/env python3
"""
Round 21 diagnostic: odd/even compensation profile on X<0 cases.

Conventions match Round 19 scripts:
- support-root processing
- non-leaf children sorted by rooted subtree size, then id
- boundary-correct indexing
- prefix k < mode(I_new), smallest-index mode

For each X<0 case we compute:
- even_err, odd_err, sum_err from diagonal err_s definition
- D, R_shift = C10 + C01 + C11
- compensation quantities relative to lambda0
"""

from __future__ import annotations

import argparse
import json
import subprocess
from collections import defaultdict
from dataclasses import dataclass

GENG = "/opt/homebrew/bin/geng"
LAMBDA0 = 0.05201381704686925


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
    sub = [0] * n
    for v in reversed(order):
        if not children[v]:
            E[v] = [1]
            J[v] = [1]
            I[v] = [1, 1]
            sub[v] = 1
            continue
        ev = [1]
        jv = [1]
        sz = 1
        for c in children[v]:
            ev = poly_mul(ev, I[c])
            jv = poly_mul(jv, E[c])
            sz += sub[c]
        E[v] = ev
        J[v] = jv
        deg = max(len(ev), len(jv) + 1)
        iv = [0] * deg
        for i, val in enumerate(ev):
            iv[i] += val
        for i, val in enumerate(jv):
            iv[i + 1] += val
        I[v] = iv
        sub[v] = sz
    return I, E, J, children, sub


@dataclass
class Witness:
    value: float
    n: int
    g6: str
    root: int
    step: int
    k: int
    a: int
    b: int
    extra: dict


def update_witness(w: Witness | None, value: float, n: int, g6: str, root: int, step: int, k: int, a: int, b: int, extra: dict):
    if w is None or value > w.value:
        return Witness(value=value, n=n, g6=g6, root=root, step=step, k=k, a=a, b=b, extra=extra)
    return w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=18)
    ap.add_argument("--out-json", type=str, default="results/round21_odd_compensation_n18.json")
    args = ap.parse_args()

    total_xneg = 0
    xneg_step_not2 = 0

    even_bound_fail = 0
    global_bound_fail = 0

    odd_over_r_fail = 0
    odd_over_r_fail_with_comp = 0

    w_max_odd_over_r = None
    w_max_needed_global_lambda = None
    w_max_odd_surplus_over_even_slack = None

    pair_stats = defaultdict(lambda: {
        "count": 0,
        "odd_over_r_fail": 0,
        "odd_over_r_fail_with_comp": 0,
        "max_needed_global_lambda": 0.0,
        "max_needed_odd_lambda": 0.0,
        "max_odd_surplus_over_even_slack": 0.0,
        "witness_global_lambda": None,
        "witness_odd_lambda": None,
        "witness_comp_ratio": None,
    })

    for n in range(3, args.max_n + 1):
        lines = subprocess.run(
            [GENG, "-cq", str(n), f"{n-1}:{n-1}"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip().splitlines()

        n_xneg = 0
        for g6 in lines:
            adj = graph6_to_adj(g6)
            for root in range(n):
                if not any(len(adj[nb]) == 1 for nb in adj[root]):
                    continue
                I, E, J, children, sub = tree_dp_with_sizes(adj, root)

                leaf = []
                nonleaf = []
                for c in children[root]:
                    if not children[c]:
                        leaf.append(c)
                    else:
                        nonleaf.append(c)
                if len(nonleaf) < 2:
                    continue

                nonleaf.sort(key=lambda c: (sub[c], c))

                e_acc = [1]
                for _ in range(len(leaf)):
                    e_acc = poly_mul(e_acc, [1, 1])
                j_acc = [1]

                for t_idx, c in enumerate(nonleaf):
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
                    D_vals = [
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
                    a = sub[nonleaf[0]]
                    b = sub[nonleaf[1]]

                    for k in range(min(mode, max_new + 1)):
                        X = lam_new[k] - D_vals[k]
                        if X >= 0:
                            continue

                        total_xneg += 1
                        n_xneg += 1
                        if step != 2:
                            xneg_step_not2 += 1

                        odd_err = 0
                        even_err = 0
                        for s in range(2 * lo, 2 * hi + 1):
                            vals_u = []
                            vals_D = []
                            for i in range(lo, hi + 1):
                                u_i = coeff(A, i) * coeff(B, s - i)
                                w_i = coeff(P, k - i) * coeff(Q, k - s + i)
                                w_ip1 = coeff(P, k - (i + 1)) * coeff(Q, k - s + (i + 1))
                                vals_u.append(u_i)
                                vals_D.append(w_i - w_ip1)

                            m_idx = None
                            for off, d_i in enumerate(vals_D):
                                if d_i >= 0:
                                    m_idx = off
                                    break

                            err_s = 0
                            if m_idx is not None:
                                u_m = vals_u[m_idx]
                                for off in range(m_idx, len(vals_D)):
                                    d_i = vals_D[off]
                                    if d_i <= 0:
                                        continue
                                    u_i = vals_u[off]
                                    if u_i < u_m:
                                        err_s += (u_m - u_i) * d_i

                            if s % 2 == 0:
                                even_err += err_s
                            else:
                                odd_err += err_s

                        sum_err = odd_err + even_err

                        C10 = sum(lam_old[i] * coeff(P, k - i - 1) * coeff(Q, k - i) for i in range(len(lam_old)))
                        C01 = sum(lam_old[i] * coeff(P, k - i) * coeff(Q, k - i - 1) for i in range(len(lam_old)))
                        C11 = sum(lam_old[i] * coeff(P, k - i - 1) * coeff(Q, k - i - 1) for i in range(len(lam_old)))
                        R_shift = C10 + C01 + C11

                        Dk = D_vals[k]
                        rhs_global = Dk + LAMBDA0 * R_shift

                        if even_err > Dk:
                            even_bound_fail += 1
                        if sum_err > rhs_global + 1e-12:
                            global_bound_fail += 1

                        odd_rhs = LAMBDA0 * R_shift
                        odd_over_r = odd_err - odd_rhs
                        odd_surplus = max(0.0, odd_over_r)
                        even_slack = max(0.0, Dk - even_err)

                        if odd_over_r > 1e-12:
                            odd_over_r_fail += 1
                            if odd_surplus <= even_slack + 1e-12:
                                odd_over_r_fail_with_comp += 1

                        needed_odd_lambda = (odd_err / R_shift) if R_shift > 0 else (float("inf") if odd_err > 0 else 0.0)
                        needed_global_lambda = (max(0.0, sum_err - Dk) / R_shift) if R_shift > 0 else (float("inf") if sum_err > Dk else 0.0)
                        comp_ratio = (odd_surplus / even_slack) if even_slack > 0 else (float("inf") if odd_surplus > 0 else 0.0)

                        info = {
                            "X": X,
                            "D": Dk,
                            "odd_err": odd_err,
                            "even_err": even_err,
                            "sum_err": sum_err,
                            "R_shift": R_shift,
                            "needed_odd_lambda": needed_odd_lambda,
                            "needed_global_lambda": needed_global_lambda,
                            "odd_surplus": odd_surplus,
                            "even_slack": even_slack,
                            "comp_ratio": comp_ratio,
                            "C10": C10,
                            "C01": C01,
                            "C11": C11,
                        }

                        w_max_odd_over_r = update_witness(
                            w_max_odd_over_r,
                            needed_odd_lambda,
                            n,
                            g6,
                            root,
                            step,
                            k,
                            a,
                            b,
                            info,
                        )
                        w_max_needed_global_lambda = update_witness(
                            w_max_needed_global_lambda,
                            needed_global_lambda,
                            n,
                            g6,
                            root,
                            step,
                            k,
                            a,
                            b,
                            info,
                        )
                        w_max_odd_surplus_over_even_slack = update_witness(
                            w_max_odd_surplus_over_even_slack,
                            comp_ratio,
                            n,
                            g6,
                            root,
                            step,
                            k,
                            a,
                            b,
                            info,
                        )

                        ps = pair_stats[(a, b)]
                        ps["count"] += 1
                        if odd_over_r > 1e-12:
                            ps["odd_over_r_fail"] += 1
                            if odd_surplus <= even_slack + 1e-12:
                                ps["odd_over_r_fail_with_comp"] += 1

                        if needed_global_lambda > ps["max_needed_global_lambda"]:
                            ps["max_needed_global_lambda"] = needed_global_lambda
                            ps["witness_global_lambda"] = {
                                "n": n,
                                "g6": g6,
                                "root": root,
                                "step": step,
                                "k": k,
                                **info,
                            }
                        if needed_odd_lambda > ps["max_needed_odd_lambda"]:
                            ps["max_needed_odd_lambda"] = needed_odd_lambda
                            ps["witness_odd_lambda"] = {
                                "n": n,
                                "g6": g6,
                                "root": root,
                                "step": step,
                                "k": k,
                                **info,
                            }
                        if comp_ratio > ps["max_odd_surplus_over_even_slack"]:
                            ps["max_odd_surplus_over_even_slack"] = comp_ratio
                            ps["witness_comp_ratio"] = {
                                "n": n,
                                "g6": g6,
                                "root": root,
                                "step": step,
                                "k": k,
                                **info,
                            }

                    e_acc = e_new
                    j_acc = j_new

        print(f"n={n:2d}: X<0={n_xneg}", flush=True)

    pair_rows = []
    for (a, b), ps in sorted(pair_stats.items(), key=lambda kv: (-kv[1]["count"], kv[0])):
        pair_rows.append({"a": a, "b": b, **ps})

    summary = {
        "scope": "support-root step-prefix X<0 cases, boundary-correct indexing",
        "max_n": args.max_n,
        "lambda0": LAMBDA0,
        "total_xneg": total_xneg,
        "xneg_step_not2": xneg_step_not2,
        "even_bound_fail": even_bound_fail,
        "global_bound_fail": global_bound_fail,
        "odd_over_r_fail": odd_over_r_fail,
        "odd_over_r_fail_with_comp": odd_over_r_fail_with_comp,
        "w_max_odd_over_r": None if w_max_odd_over_r is None else w_max_odd_over_r.__dict__,
        "w_max_needed_global_lambda": None if w_max_needed_global_lambda is None else w_max_needed_global_lambda.__dict__,
        "w_max_odd_surplus_over_even_slack": None if w_max_odd_surplus_over_even_slack is None else w_max_odd_surplus_over_even_slack.__dict__,
        "pair_rows": pair_rows,
    }

    with open(args.out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote {args.out_json}")


if __name__ == "__main__":
    main()
