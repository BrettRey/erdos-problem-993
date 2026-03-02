#!/usr/bin/env python3
"""
Round 18 separator scan (boundary-correct CB indexing, prefix regime).

Focus:
1) Negative gap-mass budget vs diagonal reserve D_k.
2) Diagonal-Abel difference sign-change structure.
3) Extremal witnesses for X_k < 0 and NegGap/D.

Conventions:
- Trees from geng.
- Support-root processing: absorb leaf children, then non-leaf children by
  increasing rooted-subtree size.
- Prefix at step: k < mode(I_new), smallest-index mode.
- CB cross indices include boundary -1.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction


GENG = "/opt/homebrew/bin/geng"


def poly_mul(a: list[int], b: list[int]) -> list[int]:
    if not a or not b:
        return [0]
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
    visited = [False] * n
    stack = [root]
    visited[root] = True
    while stack:
        v = stack.pop()
        order.append(v)
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
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
        else:
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


def sign_changes_nonzero(seq: list[int]) -> int:
    signs = []
    for x in seq:
        if x > 0:
            signs.append(1)
        elif x < 0:
            signs.append(-1)
    if not signs:
        return 0
    changes = 0
    last = signs[0]
    for s in signs[1:]:
        if s != last:
            changes += 1
            last = s
    return changes


@dataclass
class Witness:
    value: int
    n: int
    g6: str
    root: int
    step: int
    k: int
    D: int
    X: int


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=15)
    ap.add_argument("--out-json", type=str, default="results/round18_separator_scan.json")
    args = ap.parse_args()

    total_trees = 0
    support_rootings = 0
    steps = 0
    prefix_checks = 0
    x_neg = 0

    neg_gap_le_D = 0
    neg_gap_total = 0
    neg_gap_positive = 0
    neg_gap_over_D = []
    neg_gap_over_D_max = None

    x_neg_ratio_D_over_absX = []
    min_x_witness = None

    # Diagonal-Abel sign-change stats for X<0 cases.
    diag_total = 0
    diag_gt1_changes = 0
    max_changes_seen = 0
    xneg_cases_with_any_diag_gt1 = 0

    by_step = Counter()
    by_mode_dist = Counter()

    for n in range(3, args.max_n + 1):
        cmd = [GENG, "-cq", str(n), f"{n-1}:{n-1}"]
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
        lines = [l for l in out.strip().split("\n") if l]

        n_prefix = 0
        n_xneg = 0
        for g6 in lines:
            total_trees += 1
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
                    P_c = I[c]
                    Q_c = E[c]
                    e_new = poly_mul(e_acc, P_c)
                    j_new = poly_mul(j_acc, Q_c)
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
                    lambda_new = []
                    for k in range(max_new + 1):
                        lambda_new.append(
                            coeff(e_new, k) * coeff(j_new, k)
                            - coeff(e_new, k - 1) * coeff(j_new, k + 1)
                        )
                    D = [0] * (max_new + 1)
                    for k in range(max_new + 1):
                        s = 0
                        for i in range(len(lambda_old)):
                            s += lambda_old[i] * coeff(P_c, k - i) * coeff(Q_c, k - i)
                        D[k] = s

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
                        if xk < 0:
                            x_neg += 1
                            n_xneg += 1
                            if min_x_witness is None or xk < min_x_witness.value:
                                min_x_witness = Witness(
                                    value=xk,
                                    n=n,
                                    g6=g6,
                                    root=root,
                                    step=step,
                                    k=k,
                                    D=D[k],
                                    X=xk,
                                )
                            if xk != 0:
                                x_neg_ratio_D_over_absX.append(Fraction(D[k], -xk))

                        # Gap sums, grouped by g = j-i (boundary-inclusive indices).
                        gap_sums = defaultdict(int)
                        for i in range(min_idx, max_idx + 1):
                            for j in range(i + 1, max_idx + 1):
                                g = j - i
                                delta_ij = coeff(A, i) * coeff(B, j) - coeff(A, i - 1) * coeff(B, j + 1)
                                delta_ji = coeff(A, j) * coeff(B, i) - coeff(A, j - 1) * coeff(B, i + 1)
                                term = (
                                    delta_ij * coeff(P_c, k - i) * coeff(Q_c, k - j)
                                    + delta_ji * coeff(P_c, k - j) * coeff(Q_c, k - i)
                                )
                                gap_sums[g] += term

                        neg_gap = sum(-v for v in gap_sums.values() if v < 0)
                        neg_gap_total += 1
                        if neg_gap > 0:
                            neg_gap_positive += 1
                        if neg_gap <= D[k]:
                            neg_gap_le_D += 1

                        if D[k] > 0:
                            ratio = Fraction(neg_gap, D[k])
                            neg_gap_over_D.append(ratio)
                            if neg_gap_over_D_max is None or ratio > neg_gap_over_D_max[0]:
                                neg_gap_over_D_max = (
                                    ratio,
                                    n,
                                    g6,
                                    root,
                                    step,
                                    k,
                                    neg_gap,
                                    D[k],
                                    xk,
                                )

                        if xk < 0:
                            any_diag_gt1 = False
                            # Diagonal-Abel D_i^{(s)} sign-change profile.
                            s_min = 2 * min_idx
                            s_max = 2 * max_idx
                            for s in range(s_min, s_max + 1):
                                vals = []
                                i_lo = min_idx
                                i_hi = max_idx
                                for i in range(i_lo, i_hi + 1):
                                    w_i = coeff(P_c, k - i) * coeff(Q_c, k - s + i)
                                    w_ip1 = coeff(P_c, k - (i + 1)) * coeff(Q_c, k - s + (i + 1))
                                    vals.append(w_i - w_ip1)
                                chg = sign_changes_nonzero(vals)
                                diag_total += 1
                                max_changes_seen = max(max_changes_seen, chg)
                                if chg > 1:
                                    diag_gt1_changes += 1
                                    any_diag_gt1 = True
                            if any_diag_gt1:
                                xneg_cases_with_any_diag_gt1 += 1

                            by_step[step] += 1
                            by_mode_dist[mode - k] += 1

                    e_acc = e_new
                    j_acc = j_new

        print(
            f"n={n:2d}: trees={len(lines):8d} prefix_checks={n_prefix:10d} X<0={n_xneg:8d}",
            flush=True,
        )

    summary = {
        "max_n": args.max_n,
        "totals": {
            "trees": total_trees,
            "support_rootings": support_rootings,
            "steps_t_ge_2": steps,
            "prefix_checks": prefix_checks,
            "x_neg_count": x_neg,
            "x_neg_rate": (x_neg / prefix_checks if prefix_checks else None),
        },
        "neg_gap_budget": {
            "checks": neg_gap_total,
            "neg_gap_positive": neg_gap_positive,
            "neg_gap_le_D_count": neg_gap_le_D,
            "neg_gap_le_D_rate": (neg_gap_le_D / neg_gap_total if neg_gap_total else None),
            "neg_gap_over_D_max": (
                {
                    "ratio": str(neg_gap_over_D_max[0]),
                    "ratio_float": float(neg_gap_over_D_max[0]),
                    "n": neg_gap_over_D_max[1],
                    "g6": neg_gap_over_D_max[2],
                    "root": neg_gap_over_D_max[3],
                    "step": neg_gap_over_D_max[4],
                    "k": neg_gap_over_D_max[5],
                    "neg_gap": neg_gap_over_D_max[6],
                    "D": neg_gap_over_D_max[7],
                    "X": neg_gap_over_D_max[8],
                }
                if neg_gap_over_D_max
                else None
            ),
        },
        "xneg_D_over_absX": {
            "count": len(x_neg_ratio_D_over_absX),
            "min": (str(min(x_neg_ratio_D_over_absX)) if x_neg_ratio_D_over_absX else None),
            "min_float": (float(min(x_neg_ratio_D_over_absX)) if x_neg_ratio_D_over_absX else None),
            "median": (
                str(sorted(x_neg_ratio_D_over_absX)[len(x_neg_ratio_D_over_absX) // 2])
                if x_neg_ratio_D_over_absX
                else None
            ),
        },
        "diag_abel_sign_changes_on_xneg": {
            "diag_total": diag_total,
            "diag_gt1_changes": diag_gt1_changes,
            "diag_gt1_rate": (diag_gt1_changes / diag_total if diag_total else None),
            "max_changes_seen": max_changes_seen,
            "xneg_cases_with_any_diag_gt1": xneg_cases_with_any_diag_gt1,
        },
        "witnesses": {
            "min_X": (
                {
                    "X": min_x_witness.value,
                    "n": min_x_witness.n,
                    "g6": min_x_witness.g6,
                    "root": min_x_witness.root,
                    "step": min_x_witness.step,
                    "k": min_x_witness.k,
                    "D": min_x_witness.D,
                }
                if min_x_witness
                else None
            ),
        },
        "xneg_by_step": dict(by_step),
        "xneg_by_mode_distance": dict(by_mode_dist),
    }

    with open(args.out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)

    print("\n=== ROUND18 SEPARATOR SUMMARY ===")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote {args.out_json}")


if __name__ == "__main__":
    main()

