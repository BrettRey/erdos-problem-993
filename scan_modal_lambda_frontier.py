#!/usr/bin/env python3
"""Modal scan for drift-aware lambda frontier on X<0 corpus.

Target quantity (per X<0 prefix case):
    lambda_needed = max(0, sum_all - D) / R_shift

where
    R_shift = C10 + C01 + C11
and conventions handle R_shift = 0.
"""

from __future__ import annotations

import json
import math
import subprocess
import time
from typing import Any

import modal

app = modal.App("erdos-993-lambda-frontier")

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("build-essential", "curl")
    .run_commands(
        "curl -sL https://pallini.di.uniroma1.it/nauty2_8_9.tar.gz | tar xz",
        "cd nauty2_8_9 && ./configure --quiet && make -j$(nproc) geng",
        "cp nauty2_8_9/geng /usr/local/bin/geng",
        "rm -rf nauty2_8_9",
    )
    .add_local_file("graph6.py", "/root/src/graph6.py")
)


def _dict_name(min_n: int, max_n: int) -> str:
    return f"erdos-993-lambda-frontier-n{min_n}-n{max_n}"


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
    bits: list[int] = []
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


def _pair_key(a: int, b: int) -> str:
    return f"{a},{b}"


@app.function(image=image, timeout=7200, cpu=1)
def scan_partition(n: int, res: int, mod: int, dict_name: str) -> dict[str, Any]:
    """Scan one geng partition for lambda frontier stats on X<0 cases."""
    results_dict = modal.Dict.from_name(dict_name, create_if_missing=True)

    lines = subprocess.run(
        ["geng", "-cq", str(n), f"{n-1}:{n-1}", f"{res}/{mod}"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip().splitlines()

    xneg_total = 0
    xneg_step2 = 0

    lambda_max = 0.0
    lambda_wit: dict[str, Any] | None = None
    impossible_lambda = False

    pair_stats: dict[str, dict[str, Any]] = {}

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
            a = sub[nonleaf[0]]
            b = sub[nonleaf[1]]
            pkey = _pair_key(a, b)
            if pkey not in pair_stats:
                pair_stats[pkey] = {
                    "count": 0,
                    "lambda_max": 0.0,
                    "impossible": False,
                    "wit": None,
                }

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

                for k in range(min(mode, max_new + 1)):
                    X = lam_new[k] - D_vals[k]
                    if X >= 0:
                        continue

                    xneg_total += 1
                    if step == 2:
                        xneg_step2 += 1
                        pair_stats[pkey]["count"] += 1

                    sum_all = 0
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
                        if m_idx is None:
                            continue

                        err_s = 0
                        u_m = vals_u[m_idx]
                        for off in range(m_idx, len(vals_D)):
                            d_i = vals_D[off]
                            if d_i <= 0:
                                continue
                            u_i = vals_u[off]
                            if u_i < u_m:
                                err_s += (u_m - u_i) * d_i
                        sum_all += err_s

                    R_shift = sum(
                        lam_old[i] * coeff(P, k - i - 1) * coeff(Q, k - i)
                        + lam_old[i] * coeff(P, k - i) * coeff(Q, k - i - 1)
                        + lam_old[i] * coeff(P, k - i - 1) * coeff(Q, k - i - 1)
                        for i in range(len(lam_old))
                    )

                    need = max(0, sum_all - D_vals[k])
                    wit = {
                        "n": n,
                        "g6": g6,
                        "root": root,
                        "step": step,
                        "k": k,
                        "a": a,
                        "b": b,
                        "X": X,
                        "R": R_shift,
                        "sum_all": sum_all,
                        "D": D_vals[k],
                        "need": need,
                    }

                    if R_shift > 0:
                        lam_need = need / R_shift
                    elif need == 0:
                        lam_need = 0.0
                    else:
                        lam_need = math.inf

                    if math.isinf(lam_need):
                        impossible_lambda = True
                        lambda_max = math.inf
                        if lambda_wit is None:
                            lambda_wit = {**wit, "lambda_needed": "inf"}
                    elif lam_need > lambda_max:
                        lambda_max = lam_need
                        lambda_wit = {**wit, "lambda_needed": lam_need}

                    if step == 2:
                        row = pair_stats[pkey]
                        if math.isinf(lam_need):
                            row["impossible"] = True
                            row["lambda_max"] = math.inf
                            if row["wit"] is None:
                                row["wit"] = {**wit, "lambda_needed": "inf"}
                        elif (not math.isinf(row["lambda_max"])) and lam_need > row["lambda_max"]:
                            row["lambda_max"] = lam_need
                            row["wit"] = {**wit, "lambda_needed": lam_need}

                e_acc = e_new
                j_acc = j_new

    result = {
        "n": n,
        "res": res,
        "mod": mod,
        "xneg_total": xneg_total,
        "xneg_step2": xneg_step2,
        "lambda_max": lambda_max,
        "impossible_lambda": impossible_lambda,
        "lambda_wit": lambda_wit,
        "pair_stats": pair_stats,
    }
    results_dict[f"{n}/{res}/{mod}"] = result
    return result


@app.function(image=image, timeout=300)
def collect_results(min_n: int, max_n: int, workers: int, dict_name: str) -> dict[str, Any]:
    results_dict = modal.Dict.from_name(dict_name, create_if_missing=True)
    total_expected = (max_n - min_n + 1) * workers
    completed = 0
    xneg = 0
    for n in range(min_n, max_n + 1):
        for res in range(workers):
            key = f"{n}/{res}/{workers}"
            try:
                row = results_dict[key]
                completed += 1
                xneg += row["xneg_total"]
            except KeyError:
                pass
    return {
        "min_n": min_n,
        "max_n": max_n,
        "workers": workers,
        "completed": completed,
        "expected": total_expected,
        "xneg_total_so_far": xneg,
    }


@app.function(image=image, timeout=900)
def launch_partitions(min_n: int, max_n: int, workers: int, dict_name: str) -> dict[str, Any]:
    """Server-side launcher: spawn scan shards from inside Modal."""
    total = (max_n - min_n + 1) * workers
    for n in range(min_n, max_n + 1):
        for res in range(workers):
            scan_partition.spawn(n, res, workers, dict_name)
    return {
        "min_n": min_n,
        "max_n": max_n,
        "workers": workers,
        "total_tasks": total,
        "dict_name": dict_name,
    }


@app.local_entrypoint()
def main(
    min_n: int = 19,
    max_n: int = 22,
    workers: int = 512,
    out_json: str = "",
):
    if out_json == "":
        out_json = f"results/lambda_frontier_modal_n{min_n}_n{max_n}.json"

    dict_name = _dict_name(min_n, max_n)
    print(
        f"Modal lambda-frontier scan: n={min_n}..{max_n}, workers={workers}, "
        f"dict={dict_name}"
    )

    t0 = time.time()
    total_tasks = (max_n - min_n + 1) * workers
    completed = 0

    global_xneg = 0
    global_xneg_step2 = 0
    lambda_max = 0.0
    lambda_wit = None
    impossible_lambda = False

    by_n: dict[int, dict[str, Any]] = {}
    pair_aggs: dict[str, dict[str, Any]] = {}

    for n in range(min_n, max_n + 1):
        by_n[n] = {
            "xneg_total": 0,
            "xneg_step2": 0,
            "lambda_max": 0.0,
            "impossible_lambda": False,
            "lambda_wit": None,
        }

        for result in scan_partition.starmap(
            [(n, res, workers, dict_name) for res in range(workers)]
        ):
            completed += 1
            global_xneg += result["xneg_total"]
            global_xneg_step2 += result["xneg_step2"]

            row_n = by_n[n]
            row_n["xneg_total"] += result["xneg_total"]
            row_n["xneg_step2"] += result["xneg_step2"]

            if result["impossible_lambda"]:
                row_n["impossible_lambda"] = True
                impossible_lambda = True

            if (not math.isinf(row_n["lambda_max"])) and result["lambda_max"] > row_n["lambda_max"]:
                row_n["lambda_max"] = result["lambda_max"]
                row_n["lambda_wit"] = result["lambda_wit"]

            if (not math.isinf(lambda_max)) and result["lambda_max"] > lambda_max:
                lambda_max = result["lambda_max"]
                lambda_wit = result["lambda_wit"]

            for pkey, pstat in result["pair_stats"].items():
                if pkey not in pair_aggs:
                    pair_aggs[pkey] = {
                        "count": 0,
                        "lambda_max": 0.0,
                        "impossible": False,
                        "wit": None,
                    }
                pa = pair_aggs[pkey]
                pa["count"] += pstat["count"]
                if pstat["impossible"]:
                    pa["impossible"] = True
                    pa["lambda_max"] = math.inf
                    if pa["wit"] is None:
                        pa["wit"] = pstat["wit"]
                elif (not math.isinf(pa["lambda_max"])) and pstat["lambda_max"] > pa["lambda_max"]:
                    pa["lambda_max"] = pstat["lambda_max"]
                    pa["wit"] = pstat["wit"]

            if completed % 50 == 0 or completed == total_tasks:
                elapsed = time.time() - t0
                print(
                    f"  {completed}/{total_tasks} tasks done, "
                    f"xneg={global_xneg:,}, elapsed={elapsed:.0f}s"
                )

    elapsed = time.time() - t0

    def norm(x: float):
        if x == math.inf:
            return "inf"
        if x == -math.inf:
            return "-inf"
        return x

    by_n_json: dict[str, Any] = {}
    for n in range(min_n, max_n + 1):
        row = by_n[n]
        by_n_json[str(n)] = {
            "xneg_total": row["xneg_total"],
            "xneg_step2": row["xneg_step2"],
            "lambda_max": norm(row["lambda_max"]),
            "impossible_lambda": row["impossible_lambda"],
            "lambda_wit": row["lambda_wit"],
        }

    pair_json: dict[str, Any] = {}
    for pkey, row in sorted(pair_aggs.items()):
        pair_json[pkey] = {
            "count": row["count"],
            "lambda_max": norm(row["lambda_max"]),
            "impossible": row["impossible"],
            "wit": row["wit"],
        }

    summary = {
        "min_n": min_n,
        "max_n": max_n,
        "workers": workers,
        "xneg_total": global_xneg,
        "xneg_step2_total": global_xneg_step2,
        "lambda_max": norm(lambda_max),
        "impossible_lambda": impossible_lambda,
        "lambda_wit": lambda_wit,
        "by_n": by_n_json,
        "by_pair": pair_json,
        "wall_time_s": round(elapsed, 1),
        "platform": "Modal",
        "dict_name": dict_name,
        "date": time.strftime("%Y-%m-%d"),
    }

    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\nDone")
    print(
        json.dumps(
            {
                "min_n": min_n,
                "max_n": max_n,
                "xneg_total": global_xneg,
                "xneg_step2_total": global_xneg_step2,
                "lambda_max": summary["lambda_max"],
                "wall_time_s": round(elapsed, 1),
            },
            indent=2,
        )
    )
    print(f"Saved: {out_json}")


@app.local_entrypoint()
def dispatch(min_n: int = 19, max_n: int = 22, workers: int = 512):
    """Fire-and-forget launch: spawn one lambda scan task per (n,res)."""
    dict_name = _dict_name(min_n, max_n)
    print(f"Dispatching lambda-frontier tasks via server launcher, dict={dict_name}")
    info = launch_partitions.remote(min_n, max_n, workers, dict_name)
    total = info["total_tasks"]
    print(
        f"Dispatching lambda-frontier tasks for n={min_n}..{max_n}, "
        f"workers={workers}, total_tasks={total}, dict={dict_name}"
    )
    print("Dispatch submitted.")
