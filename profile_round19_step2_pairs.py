#!/usr/bin/env python3
"""
Profile Round 19 step-2 negative cases by child-subtree-size pair.

Conventions:
- Exhaustive trees via geng.
- Support-root processing.
- Non-leaf children sorted by rooted subtree size (vertex count), then id.
- Prefix regime: k < mode(I_new), smallest-index mode.
- Boundary-correct CB indexing (includes -1).
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
class PairStats:
    count: int = 0
    worst_x: int = 0
    worst_x_witness: dict | None = None
    min_d_over_absx: float = float("inf")
    min_d_over_absx_witness: dict | None = None
    max_sumerr_over_d: float = -1.0
    max_sumerr_over_d_witness: dict | None = None
    max_odd_over_d: float = -1.0
    max_odd_over_d_witness: dict | None = None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=18)
    ap.add_argument("--out-json", type=str, default="results/round19_step2_pair_profile_n18.json")
    args = ap.parse_args()

    pair_stats: dict[tuple[int, int], PairStats] = defaultdict(PairStats)
    total_xneg = 0
    xneg_step_not2 = 0

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
                    lam_old = []
                    for k in range(max_ab + 1):
                        lam_old.append(
                            coeff(A, k) * coeff(B, k) - coeff(A, k - 1) * coeff(B, k + 1)
                        )
                    D = [0] * (max_new + 1)
                    for k in range(max_new + 1):
                        D[k] = sum(
                            lam_old[i] * coeff(P, k - i) * coeff(Q, k - i)
                            for i in range(len(lam_old))
                        )
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
                        total_xneg += 1
                        n_xneg += 1
                        if step != 2:
                            xneg_step_not2 += 1

                        # Only profile step-2 by child-size pair.
                        if step != 2:
                            continue
                        a = subsize[nonleaf_kids[0]]
                        b = subsize[nonleaf_kids[1]]
                        key = (a, b)
                        ps = pair_stats[key]
                        ps.count += 1

                        # Error splits.
                        sum_err = 0
                        odd_err = 0
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
                            if s % 2 != 0:
                                odd_err += err_s

                        wit = {
                            "n": n,
                            "g6": g6,
                            "root": root,
                            "step": step,
                            "k": k,
                            "mode": mode,
                            "a": a,
                            "b": b,
                            "X": X,
                            "D": D[k],
                            "sum_err": sum_err,
                            "odd_err": odd_err,
                        }

                        if X < ps.worst_x:
                            ps.worst_x = X
                            ps.worst_x_witness = wit

                        d_over_absx = D[k] / (-X)
                        if d_over_absx < ps.min_d_over_absx:
                            ps.min_d_over_absx = d_over_absx
                            ps.min_d_over_absx_witness = wit

                        if D[k] > 0:
                            r_sum = sum_err / D[k]
                            if r_sum > ps.max_sumerr_over_d:
                                ps.max_sumerr_over_d = r_sum
                                ps.max_sumerr_over_d_witness = wit
                            r_odd = odd_err / D[k]
                            if r_odd > ps.max_odd_over_d:
                                ps.max_odd_over_d = r_odd
                                ps.max_odd_over_d_witness = wit

                    e_acc = e_new
                    j_acc = j_new

        print(f"n={n:2d}: X<0 cases={n_xneg}", flush=True)

    pair_table = []
    for (a, b), ps in sorted(pair_stats.items(), key=lambda kv: (-kv[1].count, kv[0])):
        pair_table.append(
            {
                "a": a,
                "b": b,
                "count": ps.count,
                "worst_x": ps.worst_x,
                "min_d_over_absx": ps.min_d_over_absx,
                "max_sumerr_over_d": ps.max_sumerr_over_d,
                "max_odd_over_d": ps.max_odd_over_d,
                "worst_x_witness": ps.worst_x_witness,
                "min_d_over_absx_witness": ps.min_d_over_absx_witness,
                "max_sumerr_over_d_witness": ps.max_sumerr_over_d_witness,
                "max_odd_over_d_witness": ps.max_odd_over_d_witness,
            }
        )

    summary = {
        "max_n": args.max_n,
        "total_xneg": total_xneg,
        "xneg_step_not2": xneg_step_not2,
        "unique_pairs": len(pair_table),
        "pair_table": pair_table,
    }

    with open(args.out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote {args.out_json}")


if __name__ == "__main__":
    main()

