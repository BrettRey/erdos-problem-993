#!/usr/bin/env python3
"""
Fit shifted-channel reserve coefficients for Round 19.

Model:
  denom = D + alpha*C10 + beta*C01 + gamma*C11
  score = max over X<0 cases of (sum_err / denom)

All cases use boundary-correct support-step prefix conventions.
"""

from __future__ import annotations

import argparse
import json
import subprocess
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


@dataclass
class Row:
    D: int
    sum_err: int
    C10: int
    C01: int
    C11: int
    n: int
    g6: str
    root: int
    step: int
    k: int
    X: int


def collect_rows(max_n: int) -> list[Row]:
    rows: list[Row] = []
    for n in range(3, max_n + 1):
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

                    A = e_acc
                    B = j_acc
                    max_new = max(len(e_new), len(j_new))
                    max_ab = max(len(A), len(B))
                    lam_old = []
                    for k in range(max_ab + 1):
                        lam_old.append(
                            coeff(A, k) * coeff(B, k) - coeff(A, k - 1) * coeff(B, k + 1)
                        )
                    D = [0] * (max_new + 1)
                    for k in range(max_new + 1):
                        s = 0
                        for i in range(len(lam_old)):
                            s += lam_old[i] * coeff(P, k - i) * coeff(Q, k - i)
                        D[k] = s
                    lam_new = []
                    for k in range(max_new + 1):
                        lam_new.append(
                            coeff(e_new, k) * coeff(j_new, k)
                            - coeff(e_new, k - 1) * coeff(j_new, k + 1)
                        )
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
                        n_xneg += 1

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
                            err_s = 0
                            for off in range(m_idx, len(vals_D)):
                                d_i = vals_D[off]
                                if d_i <= 0:
                                    continue
                                u_i = vals_u[off]
                                if u_i < u_m:
                                    err_s += (u_m - u_i) * d_i
                            sum_err += err_s

                        C10 = 0
                        C01 = 0
                        C11 = 0
                        for i, l_i in enumerate(lam_old):
                            C10 += l_i * coeff(P, k - i - 1) * coeff(Q, k - i)
                            C01 += l_i * coeff(P, k - i) * coeff(Q, k - i - 1)
                            C11 += l_i * coeff(P, k - i - 1) * coeff(Q, k - i - 1)

                        rows.append(
                            Row(
                                D=D[k],
                                sum_err=sum_err,
                                C10=C10,
                                C01=C01,
                                C11=C11,
                                n=n,
                                g6=g6,
                                root=root,
                                step=step,
                                k=k,
                                X=X,
                            )
                        )

                    e_acc = e_new
                    j_acc = j_new
        print(f"n={n:2d}: X<0 cases={n_xneg}", flush=True)
    return rows


def eval_params(rows: list[Row], a: float, b: float, g: float):
    worst = -1.0
    worst_row = None
    for r in rows:
        denom = r.D + a * r.C10 + b * r.C01 + g * r.C11
        if denom <= 0:
            return float("inf"), r
        ratio = r.sum_err / denom
        if ratio > worst:
            worst = ratio
            worst_row = r
    return worst, worst_row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=18)
    ap.add_argument("--out-json", type=str, default="results/round19_shifted_reserve_fit_n18.json")
    args = ap.parse_args()

    rows = collect_rows(args.max_n)
    baseline = max(r.sum_err / r.D for r in rows)

    coarse_vals = [i / 10 for i in range(0, 11)]  # 0.0 .. 1.0
    best = (float("inf"), None, None, None, None)
    for a in coarse_vals:
        for b in coarse_vals:
            for g in coarse_vals:
                w, wr = eval_params(rows, a, b, g)
                if w < best[0]:
                    best = (w, a, b, g, wr)

    _, a0, b0, g0, _ = best
    # Refine around coarse best on ±0.2 window with step 0.05
    def local_grid(center: float):
        vals = []
        start = max(0.0, center - 0.2)
        end = center + 0.2
        x = start
        while x <= end + 1e-12:
            vals.append(round(x, 4))
            x += 0.05
        return vals

    fine_a = local_grid(a0)
    fine_b = local_grid(b0)
    fine_g = local_grid(g0)
    best_f = (float("inf"), None, None, None, None)
    for a in fine_a:
        for b in fine_b:
            for g in fine_g:
                w, wr = eval_params(rows, a, b, g)
                if w < best_f[0]:
                    best_f = (w, a, b, g, wr)

    def row_to_dict(r: Row):
        if r is None:
            return None
        return {
            "n": r.n,
            "g6": r.g6,
            "root": r.root,
            "step": r.step,
            "k": r.k,
            "X": r.X,
            "D": r.D,
            "sum_err": r.sum_err,
            "C10": r.C10,
            "C01": r.C01,
            "C11": r.C11,
        }

    best_c_w, best_c_a, best_c_b, best_c_g, best_c_row = best
    best_f_w, best_f_a, best_f_b, best_f_g, best_f_row = best_f

    summary = {
        "max_n": args.max_n,
        "rows_xneg": len(rows),
        "baseline_max_sumerr_over_D": baseline,
        "coarse_grid": {
            "values": coarse_vals,
            "best": {
                "alpha": best_c_a,
                "beta": best_c_b,
                "gamma": best_c_g,
                "worst_ratio": best_c_w,
                "witness": row_to_dict(best_c_row),
            },
        },
        "fine_grid": {
            "alpha_values": fine_a,
            "beta_values": fine_b,
            "gamma_values": fine_g,
            "best": {
                "alpha": best_f_a,
                "beta": best_f_b,
                "gamma": best_f_g,
                "worst_ratio": best_f_w,
                "witness": row_to_dict(best_f_row),
            },
        },
    }

    with open(args.out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote {args.out_json}")


if __name__ == "__main__":
    main()

