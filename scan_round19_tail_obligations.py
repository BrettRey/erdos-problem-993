#!/usr/bin/env python3
"""
Round 19 diagnostic scan for diagonal-tail obligations.

Boundary-correct support-step prefix regime:
- index range includes -1
- prefix means k < mode(I_new), smallest-index mode

Focus on X_k < 0 instances:
1) Sign-change structure of D_i^(s) := W_i^(s)-W_{i+1}^(s)
2) Tail error term Err_s from Abel-shift baseline at turning point m
3) Candidate obligations:
   - odd s: Err_s == 0?
   - even s=2t: Err_s <= Lambda_old[t] * P(k-t) * Q(k-t)?
   - global: sum_s Err_s <= D_k?
"""

from __future__ import annotations

import argparse
import json
import subprocess
from collections import Counter
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
    idx = 0
    n = ord(s[idx]) - 63
    idx += 1
    bits = []
    for ch in s[idx:]:
        val = ord(ch) - 63
        for b in range(5, -1, -1):
            bits.append((val >> b) & 1)
    adj = [[] for _ in range(n)]
    k = 0
    for j in range(1, n):
        for i in range(j):
            if k < len(bits) and bits[k]:
                adj[i].append(j)
                adj[j].append(i)
            k += 1
    return adj


def tree_dp(adj: list[list[int]], root: int):
    n = len(adj)
    parent = [-1] * n
    children = [[] for _ in range(n)]
    order = []
    vis = [False] * n
    stack = [root]
    vis[root] = True
    while stack:
        v = stack.pop()
        order.append(v)
        for u in adj[v]:
            if not vis[u]:
                vis[u] = True
                parent[u] = v
                children[v].append(u)
                stack.append(u)

    I = [None] * n
    E = [None] * n
    J = [None] * n
    for v in reversed(order):
        if not children[v]:
            E[v] = [1]
            J[v] = [1]
            I[v] = [1, 1]
            continue
        ev = [1]
        jv = [1]
        for c in children[v]:
            ev = poly_mul(ev, I[c])
            jv = poly_mul(jv, E[c])
        E[v] = ev
        J[v] = jv
        deg = max(len(ev), len(jv) + 1)
        iv = [0] * deg
        for i, val in enumerate(ev):
            iv[i] += val
        for i, val in enumerate(jv):
            iv[i + 1] += val
        I[v] = iv
    return I, E, J, children


def sign_changes_nonzero(vals: list[int]) -> int:
    s = []
    for v in vals:
        if v > 0:
            s.append(1)
        elif v < 0:
            s.append(-1)
    if not s:
        return 0
    c = 0
    last = s[0]
    for x in s[1:]:
        if x != last:
            c += 1
            last = x
    return c


@dataclass
class Witness:
    value: float
    n: int
    g6: str
    root: int
    step: int
    k: int
    extra: dict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=18)
    ap.add_argument("--out-json", type=str, default="results/round19_tail_obligations_n18.json")
    args = ap.parse_args()

    trees = 0
    support_rootings = 0
    steps = 0
    prefix_checks = 0

    xneg_cases = 0
    xneg_by_step = Counter()

    diag_total = 0
    diag_changes_gt1 = 0
    max_changes = 0
    xneg_with_any_diag_gt1 = 0

    odd_diag_total = 0
    odd_diag_err_positive = 0
    odd_diag_max_err = Witness(
        value=-1.0, n=0, g6="", root=0, step=0, k=0, extra={}
    )

    even_diag_total = 0
    even_diag_bound_fail = 0
    even_diag_max_ratio = Witness(
        value=-1.0, n=0, g6="", root=0, step=0, k=0, extra={}
    )

    global_sumerr_fail = 0
    global_sumerr_max_ratio = Witness(
        value=-1.0, n=0, g6="", root=0, step=0, k=0, extra={}
    )

    first_global_fail = None
    first_odd_err_fail = None
    first_even_bound_fail = None
    first_diag_changes_fail = None

    for n in range(3, args.max_n + 1):
        lines = subprocess.run(
            [GENG, "-cq", str(n), f"{n-1}:{n-1}"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip().splitlines()

        n_prefix = 0
        n_xneg = 0
        for g6 in lines:
            trees += 1
            adj = graph6_to_adj(g6)
            for root in range(n):
                if not any(len(adj[nb]) == 1 for nb in adj[root]):
                    continue
                support_rootings += 1
                I, E, J, children = tree_dp(adj, root)

                leaf_kids = []
                nonleaf_kids = []
                for c in children[root]:
                    if not children[c]:
                        leaf_kids.append(c)
                    else:
                        nonleaf_kids.append(c)
                if len(nonleaf_kids) < 2:
                    continue
                nonleaf_kids.sort(key=lambda c: len(I[c]) - 1)

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

                    steps += 1
                    A = e_acc
                    B = j_acc
                    max_new = max(len(e_new), len(j_new))
                    max_ab = max(len(A), len(B))
                    lambda_old = []
                    for k in range(max_ab + 1):
                        lambda_old.append(
                            coeff(A, k) * coeff(B, k) - coeff(A, k - 1) * coeff(B, k + 1)
                        )
                    D = [0] * (max_new + 1)
                    for k in range(max_new + 1):
                        s = 0
                        for i in range(len(lambda_old)):
                            s += lambda_old[i] * coeff(P, k - i) * coeff(Q, k - i)
                        D[k] = s

                    lambda_new = []
                    for k in range(max_new + 1):
                        lambda_new.append(
                            coeff(e_new, k) * coeff(j_new, k)
                            - coeff(e_new, k - 1) * coeff(j_new, k + 1)
                        )

                    I_new = [0] * max(len(e_new), len(j_new) + 1)
                    for i, v in enumerate(e_new):
                        I_new[i] += v
                    for i, v in enumerate(j_new):
                        I_new[i + 1] += v
                    mode = mode_smallest(I_new)

                    min_idx = -1
                    max_idx = max(len(A), len(B))

                    for k in range(min(mode, max_new + 1)):
                        prefix_checks += 1
                        n_prefix += 1
                        xk = lambda_new[k] - D[k]
                        if xk >= 0:
                            continue
                        xneg_cases += 1
                        n_xneg += 1
                        xneg_by_step[step] += 1

                        sum_err = 0
                        any_diag_gt1 = False
                        for s in range(2 * min_idx, 2 * max_idx + 1):
                            # Sequence along i for this diagonal.
                            vals_D = []
                            vals_u = []
                            for i in range(min_idx, max_idx + 1):
                                u_i = coeff(A, i) * coeff(B, s - i)
                                w_i = coeff(P, k - i) * coeff(Q, k - s + i)
                                w_ip1 = coeff(P, k - (i + 1)) * coeff(Q, k - s + (i + 1))
                                d_i = w_i - w_ip1
                                vals_u.append(u_i)
                                vals_D.append(d_i)

                            diag_total += 1
                            changes = sign_changes_nonzero(vals_D)
                            if changes > max_changes:
                                max_changes = changes
                            if changes > 1:
                                diag_changes_gt1 += 1
                                any_diag_gt1 = True
                                if first_diag_changes_fail is None:
                                    first_diag_changes_fail = (
                                        n,
                                        g6,
                                        root,
                                        step,
                                        k,
                                        s,
                                        changes,
                                    )

                            # m := first i with D_i >= 0
                            m_idx = None
                            for off, d_i in enumerate(vals_D):
                                if d_i >= 0:
                                    m_idx = off
                                    break
                            if m_idx is None:
                                # all D<0: no positive tail -> Err=0
                                err_s = 0
                                m_i = max_idx + 1
                                u_m = 0
                            else:
                                m_i = min_idx + m_idx
                                u_m = vals_u[m_idx]
                                err_s = 0
                                # i >= m, and effectively only D_i>0 matter
                                for off in range(m_idx, len(vals_D)):
                                    d_i = vals_D[off]
                                    if d_i <= 0:
                                        continue
                                    u_i = vals_u[off]
                                    if u_i < u_m:
                                        err_s += (u_m - u_i) * d_i

                            sum_err += err_s

                            if s % 2 != 0:
                                odd_diag_total += 1
                                if err_s > 0:
                                    odd_diag_err_positive += 1
                                    if first_odd_err_fail is None:
                                        first_odd_err_fail = (
                                            n,
                                            g6,
                                            root,
                                            step,
                                            k,
                                            s,
                                            err_s,
                                            m_i,
                                        )
                                if err_s > odd_diag_max_err.value:
                                    odd_diag_max_err = Witness(
                                        value=float(err_s),
                                        n=n,
                                        g6=g6,
                                        root=root,
                                        step=step,
                                        k=k,
                                        extra={"s": s, "err_s": err_s, "m_i": m_i},
                                    )
                            else:
                                even_diag_total += 1
                                t = s // 2
                                rhs = coeff(lambda_old, t) * coeff(P, k - t) * coeff(Q, k - t)
                                if err_s > rhs:
                                    even_diag_bound_fail += 1
                                    if first_even_bound_fail is None:
                                        first_even_bound_fail = (
                                            n,
                                            g6,
                                            root,
                                            step,
                                            k,
                                            s,
                                            err_s,
                                            rhs,
                                            m_i,
                                        )
                                if rhs > 0:
                                    ratio = err_s / rhs
                                    if ratio > even_diag_max_ratio.value:
                                        even_diag_max_ratio = Witness(
                                            value=ratio,
                                            n=n,
                                            g6=g6,
                                            root=root,
                                            step=step,
                                            k=k,
                                            extra={
                                                "s": s,
                                                "t": t,
                                                "err_s": err_s,
                                                "rhs": rhs,
                                                "m_i": m_i,
                                            },
                                        )

                        if any_diag_gt1:
                            xneg_with_any_diag_gt1 += 1

                        if sum_err > D[k]:
                            global_sumerr_fail += 1
                            if first_global_fail is None:
                                first_global_fail = (
                                    n,
                                    g6,
                                    root,
                                    step,
                                    k,
                                    sum_err,
                                    D[k],
                                    xk,
                                )

                        if D[k] > 0:
                            ratio = sum_err / D[k]
                            if ratio > global_sumerr_max_ratio.value:
                                global_sumerr_max_ratio = Witness(
                                    value=ratio,
                                    n=n,
                                    g6=g6,
                                    root=root,
                                    step=step,
                                    k=k,
                                    extra={
                                        "sum_err": sum_err,
                                        "D": D[k],
                                        "X": xk,
                                    },
                                )

                    e_acc = e_new
                    j_acc = j_new

        print(
            f"n={n:2d}: trees={len(lines):8d} prefix={n_prefix:9d} X<0={n_xneg:8d}",
            flush=True,
        )

    summary = {
        "max_n": args.max_n,
        "totals": {
            "trees": trees,
            "support_rootings": support_rootings,
            "steps_t_ge_2": steps,
            "prefix_checks": prefix_checks,
            "xneg_cases": xneg_cases,
        },
        "diag_sign_changes_on_xneg": {
            "diag_total": diag_total,
            "diag_changes_gt1": diag_changes_gt1,
            "max_changes": max_changes,
            "xneg_with_any_diag_gt1": xneg_with_any_diag_gt1,
            "first_fail": first_diag_changes_fail,
        },
        "odd_err": {
            "odd_diag_total": odd_diag_total,
            "odd_diag_err_positive": odd_diag_err_positive,
            "rate": (odd_diag_err_positive / odd_diag_total if odd_diag_total else None),
            "max_err_witness": odd_diag_max_err.__dict__,
            "first_fail": first_odd_err_fail,
        },
        "even_channel_bound": {
            "even_diag_total": even_diag_total,
            "even_diag_bound_fail": even_diag_bound_fail,
            "rate": (even_diag_bound_fail / even_diag_total if even_diag_total else None),
            "max_ratio_witness": even_diag_max_ratio.__dict__,
            "first_fail": first_even_bound_fail,
        },
        "global_sumerr_vs_D": {
            "fail_count": global_sumerr_fail,
            "max_ratio_witness": global_sumerr_max_ratio.__dict__,
            "first_fail": first_global_fail,
        },
        "xneg_by_step": dict(xneg_by_step),
    }

    with open(args.out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)

    print("\n=== ROUND19 TAIL OBLIGATION SUMMARY ===")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote {args.out_json}")


if __name__ == "__main__":
    main()

