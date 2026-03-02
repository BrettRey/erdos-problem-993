#!/usr/bin/env python3
"""
Round 21 neighbor-compensation scan.

For X<0 step-prefix cases (support-root, boundary-correct indexing), define:

- odd local deficit at odd diagonal s=2t+1:
    d_t := max(0, err_{2t+1} - lambda0 * Lambda_old(t) * S_t)
  where S_t = 10 + 01 + 11 shifted kernel sum.

- even local slack at even diagonal s=2t:
    e_t := max(0, Lambda_old(t)*P(k-t)Q(k-t) - err_{2t})

We test:
  (A) d_t <= e_t
  (B) d_t <= e_t + e_{t+1}
  (C) sum_t d_t <= sum_t e_t   (per case)
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
    k: int
    step: int
    t: int
    a: int
    b: int
    extra: dict


def update_witness(w: Witness | None, value: float, n: int, g6: str, root: int, k: int, step: int, t: int, a: int, b: int, extra: dict):
    if w is None or value > w.value:
        return Witness(value=value, n=n, g6=g6, root=root, k=k, step=step, t=t, a=a, b=b, extra=extra)
    return w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=18)
    ap.add_argument("--out-json", type=str, default="results/round21_neighbor_compensation_n18.json")
    args = ap.parse_args()

    total_xneg = 0
    xneg_step_not2 = 0

    odd_slots = 0
    fail_d_le_e = 0
    fail_d_le_eep = 0
    fail_sumd_lesume = 0

    w_max_ratio_d_over_e = None
    w_max_ratio_d_over_eep = None
    w_max_ratio_sumd_over_sume = None

    pair_stats = defaultdict(lambda: {
        "odd_slots": 0,
        "fail_d_le_e": 0,
        "fail_d_le_eep": 0,
        "max_ratio_d_over_e": 0.0,
        "max_ratio_d_over_eep": 0.0,
        "witness_d_over_e": None,
        "witness_d_over_eep": None,
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

                    if step != 2:
                        e_acc = e_new
                        j_acc = j_new
                        continue

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

                        odd_err_by_t = defaultdict(int)
                        even_err_by_t = defaultdict(int)
                        even_rhs_by_t = defaultdict(int)
                        odd_local_rhs_by_t = defaultdict(float)

                        t_min = 10**9
                        t_max = -10**9

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
                                t = s // 2
                                even_err_by_t[t] += err_s
                                even_rhs_by_t[t] += coeff(lam_old, t) * coeff(P, k - t) * coeff(Q, k - t)
                            else:
                                t = s // 2
                                odd_err_by_t[t] += err_s
                                odd_local_rhs_by_t[t] += LAMBDA0 * coeff(lam_old, t) * (
                                    coeff(P, k - t - 1) * coeff(Q, k - t)
                                    + coeff(P, k - t) * coeff(Q, k - t - 1)
                                    + coeff(P, k - t - 1) * coeff(Q, k - t - 1)
                                )

                            if t < t_min:
                                t_min = t
                            if t > t_max:
                                t_max = t

                        sum_d = 0.0
                        sum_e = 0.0
                        for t in range(t_min, t_max + 1):
                            d_t = max(0.0, odd_err_by_t[t] - odd_local_rhs_by_t[t])
                            e_t = max(0.0, even_rhs_by_t[t] - even_err_by_t[t])
                            e_tp1 = max(0.0, even_rhs_by_t[t + 1] - even_err_by_t[t + 1])
                            e_t_plus = e_t + e_tp1

                            sum_d += d_t
                            sum_e += e_t

                            odd_slots += 1
                            ps = pair_stats[(a, b)]
                            ps["odd_slots"] += 1

                            if d_t > e_t + 1e-12:
                                fail_d_le_e += 1
                                ps["fail_d_le_e"] += 1
                                ratio = d_t / e_t if e_t > 0 else float("inf")
                                if ratio > ps["max_ratio_d_over_e"]:
                                    ps["max_ratio_d_over_e"] = ratio
                                    ps["witness_d_over_e"] = {
                                        "n": n,
                                        "g6": g6,
                                        "root": root,
                                        "step": step,
                                        "k": k,
                                        "t": t,
                                        "d_t": d_t,
                                        "e_t": e_t,
                                        "e_t_plus": e_t_plus,
                                        "odd_err_t": odd_err_by_t[t],
                                        "odd_local_rhs_t": odd_local_rhs_by_t[t],
                                        "even_err_t": even_err_by_t[t],
                                        "even_rhs_t": even_rhs_by_t[t],
                                        "even_err_tp1": even_err_by_t[t + 1],
                                        "even_rhs_tp1": even_rhs_by_t[t + 1],
                                    }
                                w_max_ratio_d_over_e = update_witness(
                                    w_max_ratio_d_over_e,
                                    ratio,
                                    n,
                                    g6,
                                    root,
                                    k,
                                    step,
                                    t,
                                    a,
                                    b,
                                    {
                                        "d_t": d_t,
                                        "e_t": e_t,
                                        "e_t_plus": e_t_plus,
                                    },
                                )

                            if d_t > e_t_plus + 1e-12:
                                fail_d_le_eep += 1
                                ps["fail_d_le_eep"] += 1
                                ratio = d_t / e_t_plus if e_t_plus > 0 else float("inf")
                                if ratio > ps["max_ratio_d_over_eep"]:
                                    ps["max_ratio_d_over_eep"] = ratio
                                    ps["witness_d_over_eep"] = {
                                        "n": n,
                                        "g6": g6,
                                        "root": root,
                                        "step": step,
                                        "k": k,
                                        "t": t,
                                        "d_t": d_t,
                                        "e_t": e_t,
                                        "e_t_plus": e_t_plus,
                                        "odd_err_t": odd_err_by_t[t],
                                        "odd_local_rhs_t": odd_local_rhs_by_t[t],
                                        "even_err_t": even_err_by_t[t],
                                        "even_rhs_t": even_rhs_by_t[t],
                                        "even_err_tp1": even_err_by_t[t + 1],
                                        "even_rhs_tp1": even_rhs_by_t[t + 1],
                                    }
                                w_max_ratio_d_over_eep = update_witness(
                                    w_max_ratio_d_over_eep,
                                    ratio,
                                    n,
                                    g6,
                                    root,
                                    k,
                                    step,
                                    t,
                                    a,
                                    b,
                                    {
                                        "d_t": d_t,
                                        "e_t": e_t,
                                        "e_t_plus": e_t_plus,
                                    },
                                )

                        if sum_d > sum_e + 1e-12:
                            fail_sumd_lesume += 1
                            ratio = sum_d / sum_e if sum_e > 0 else float("inf")
                            w_max_ratio_sumd_over_sume = update_witness(
                                w_max_ratio_sumd_over_sume,
                                ratio,
                                n,
                                g6,
                                root,
                                k,
                                step,
                                -1,
                                a,
                                b,
                                {"sum_d": sum_d, "sum_e": sum_e},
                            )

                    e_acc = e_new
                    j_acc = j_new

        print(f"n={n:2d}: X<0={n_xneg}", flush=True)

    pair_rows = []
    for (a, b), ps in sorted(pair_stats.items(), key=lambda kv: (-kv[1]["odd_slots"], kv[0])):
        pair_rows.append({"a": a, "b": b, **ps})

    out = {
        "scope": "support-root step-prefix X<0, boundary-correct indexing",
        "max_n": args.max_n,
        "lambda0": LAMBDA0,
        "total_xneg": total_xneg,
        "xneg_step_not2": xneg_step_not2,
        "odd_slots": odd_slots,
        "fail_d_le_e": fail_d_le_e,
        "fail_d_le_eep": fail_d_le_eep,
        "fail_sumd_lesume": fail_sumd_lesume,
        "w_max_ratio_d_over_e": None if w_max_ratio_d_over_e is None else w_max_ratio_d_over_e.__dict__,
        "w_max_ratio_d_over_eep": None if w_max_ratio_d_over_eep is None else w_max_ratio_d_over_eep.__dict__,
        "w_max_ratio_sumd_over_sume": None if w_max_ratio_sumd_over_sume is None else w_max_ratio_sumd_over_sume.__dict__,
        "pair_rows": pair_rows,
    }

    with open(args.out_json, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, sort_keys=True)

    print(json.dumps(out, indent=2, sort_keys=True))
    print(f"\nWrote {args.out_json}")


if __name__ == "__main__":
    main()
