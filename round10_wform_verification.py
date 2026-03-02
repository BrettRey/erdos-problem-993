#!/usr/bin/env python3
"""Round 10 prompt 1: W-form verification and profiling.

Tasks covered:
1) Exhaustive W-form verification on support rootings for n=21..22.
2) s>=2 extremal profiling (default n<=20).
3) Factor-level invariant checks on T_{3,4} and T_{3,3} brooms.
4) E>=J at all rootings of T_{3,4}.

All polynomial arithmetic uses indpoly._polymul/_polyadd.
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import subprocess
import time
from collections import defaultdict
from dataclasses import dataclass

from indpoly import _polyadd, _polymul


GENG = "/opt/homebrew/bin/geng"


def xshift(p: list[int]) -> list[int]:
    return [0] + list(p)


def coeff(p: list[int], k: int) -> int:
    return p[k] if 0 <= k < len(p) else 0


def delta_k(f: list[int], g: list[int], k: int) -> int:
    return coeff(f, k + 1) * coeff(g, k) - coeff(f, k) * coeff(g, k + 1)


def d_km1(u: list[int], v: list[int], k: int) -> int:
    return coeff(u, k) * coeff(v, k - 1) - coeff(u, k - 1) * coeff(v, k)


def c_k(p: list[int], k: int) -> int:
    return coeff(p, k) * coeff(p, k) - coeff(p, k - 1) * coeff(p, k + 1)


def parse_g6(g6: str) -> tuple[int, list[list[int]]]:
    s = g6.strip()
    n = ord(s[0]) - 63
    adj = [[] for _ in range(n)]

    bits: list[int] = []
    for ch in s[1:]:
        v = ord(ch) - 63
        for sh in range(5, -1, -1):
            bits.append((v >> sh) & 1)

    k = 0
    for j in range(n):
        for i in range(j):
            if k < len(bits) and bits[k]:
                adj[i].append(j)
                adj[j].append(i)
            k += 1
    return n, adj


def dp_rooted(n: int, adj: list[list[int]], root: int) -> tuple[list[int], list[list[int]], list[list[int]]]:
    """Return (children_of_root, dp0, dp1s).

    dp0[v] = E_v, dp1s[v] = J_v where I_v = E_v + x*J_v.
    """
    children = [[] for _ in range(n)]
    visited = [False] * n
    visited[root] = True
    queue = [root]
    head = 0
    while head < len(queue):
        v = queue[head]
        head += 1
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                children[v].append(u)
                queue.append(u)

    order: list[int] = []
    stack: list[tuple[int, bool]] = [(root, False)]
    while stack:
        v, done = stack.pop()
        if done:
            order.append(v)
            continue
        stack.append((v, True))
        for c in children[v]:
            stack.append((c, False))

    dp0: list[list[int]] = [None] * n  # type: ignore[assignment]
    dp1s: list[list[int]] = [None] * n  # type: ignore[assignment]
    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1s[v] = [1]
        else:
            p_s = [1]
            p_e = [1]
            for c in children[v]:
                s_c = _polyadd(dp0[c], xshift(dp1s[c]))
                p_s = _polymul(p_s, s_c)
                p_e = _polymul(p_e, dp0[c])
            dp0[v] = p_s
            dp1s[v] = p_e
    return children[root], dp0, dp1s


def degree_sequence(adj: list[list[int]]) -> list[int]:
    return sorted((len(nei) for nei in adj), reverse=True)


@dataclass
class RootProfile:
    ratio: float
    abs_margin: int
    step: int
    k: int
    term1: int
    term2: int


def make_extremal_record(
    n: int,
    root: int,
    s_total: int,
    g6: str,
    deg_seq: list[int],
    prof: RootProfile,
) -> dict:
    return {
        "n": n,
        "root": root,
        "s": s_total,
        "k": prof.k,
        "step": prof.step,
        "ratio": prof.ratio,
        "abs_margin": prof.abs_margin,
        "term1": prof.term1,
        "term2": prof.term2,
        "degree_sequence": deg_seq,
        "g6": g6,
    }


def scan_partition(args: tuple[int, int, int, int, int, int]) -> dict:
    """Scan a geng residue class for one n.

    Args tuple:
      (n, res, mod, task1_min_n, task1_max_n, profile_max_n)
    """
    n, res, mod, task1_min_n, task1_max_n, profile_max_n = args
    task1_active = task1_min_n <= n <= task1_max_n
    profile_active = n <= profile_max_n

    cmd = [GENG, "-q", str(n), f"{n-1}:{n-1}", "-c", f"{res}/{mod}"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

    trees = 0
    support_vertices = 0
    incremental_steps = 0
    w_checks = 0
    w_fail_count = 0
    min_w_margin = None
    failures: list[dict] = []

    # Task 2 aggregations
    root_count_by_s: dict[int, int] = defaultdict(int)
    roots_with_neg_term2_by_s: dict[int, int] = defaultdict(int)
    best_ratio_by_s: dict[int, dict] = {}
    best_abs_by_s: dict[int, dict] = {}
    best_ratio_by_sn: dict[tuple[int, int], dict] = {}
    best_abs_by_sn: dict[tuple[int, int], dict] = {}

    assert proc.stdout is not None
    for line in proc.stdout:
        g6 = line.strip()
        if not g6:
            continue
        trees += 1

        nn, adj = parse_g6(g6)
        assert nn == n
        deg = [len(adj[v]) for v in range(n)]

        for root in range(n):
            if not any(deg[u] == 1 for u in adj[root]):
                continue

            support_vertices += 1
            children_root, dp0, dp1s = dp_rooted(n, adj, root)

            leaf_children: list[int] = []
            nonleaf_children: list[int] = []
            for c in children_root:
                if len([u for u in adj[c] if u != root]) == 0:
                    leaf_children.append(c)
                else:
                    nonleaf_children.append(c)

            s_total = len(nonleaf_children)
            do_profile_root = profile_active and s_total >= 2

            # Stage 0: process leaves only.
            e_old = [1]
            for _ in leaf_children:
                e_old = _polymul(e_old, [1, 1])
            j_old = [1]

            root_min_ratio = math.inf
            root_ratio_prof: RootProfile | None = None
            root_min_abs = math.inf
            root_abs_prof: RootProfile | None = None
            saw_neg_term2 = False

            step = 0
            for c in nonleaf_children:
                step += 1
                e_c = dp0[c]
                j_c = dp1s[c]

                a = _polymul(e_old, e_c)
                b = _polymul(e_old, j_c)
                c_poly = _polymul(j_old, e_c)
                j_poly = c_poly

                max_k = max(len(a), len(b), len(c_poly)) + 1
                incremental_steps += 1
                w_checks += max_k

                for k in range(max_k):
                    dk = delta_k(a, c_poly, k)
                    dkm1 = d_km1(b, c_poly, k)
                    term1 = coeff(j_poly, k) * dk
                    term2 = coeff(j_poly, k + 1) * dkm1
                    w_val = term1 + term2

                    if min_w_margin is None or w_val < min_w_margin:
                        min_w_margin = w_val

                    if task1_active and w_val < 0:
                        w_fail_count += 1
                        if len(failures) < 50:
                            failures.append(
                                {
                                    "n": n,
                                    "root": root,
                                    "step": step,
                                    "k": k,
                                    "W": w_val,
                                    "term1": term1,
                                    "term2": term2,
                                    "g6": g6,
                                }
                            )

                    if do_profile_root:
                        if term2 < 0:
                            saw_neg_term2 = True
                            ratio = term1 / abs(term2)
                            if ratio < root_min_ratio:
                                root_min_ratio = ratio
                                root_ratio_prof = RootProfile(
                                    ratio=ratio,
                                    abs_margin=w_val,
                                    step=step,
                                    k=k,
                                    term1=term1,
                                    term2=term2,
                                )

                        if w_val < root_min_abs:
                            root_min_abs = w_val
                            root_abs_prof = RootProfile(
                                ratio=(term1 / abs(term2)) if term2 != 0 else math.inf,
                                abs_margin=w_val,
                                step=step,
                                k=k,
                                term1=term1,
                                term2=term2,
                            )

                e_old = _polyadd(a, xshift(b))
                j_old = c_poly

            if do_profile_root:
                root_count_by_s[s_total] += 1
                key_sn = (s_total, n)
                if saw_neg_term2 and root_ratio_prof is not None:
                    roots_with_neg_term2_by_s[s_total] += 1
                    deg_seq = degree_sequence(adj)
                    rec_ratio = make_extremal_record(n, root, s_total, g6, deg_seq, root_ratio_prof)

                    cur_s = best_ratio_by_s.get(s_total)
                    if cur_s is None or rec_ratio["ratio"] < cur_s["ratio"]:
                        best_ratio_by_s[s_total] = rec_ratio

                    cur_sn = best_ratio_by_sn.get(key_sn)
                    if cur_sn is None or rec_ratio["ratio"] < cur_sn["ratio"]:
                        best_ratio_by_sn[key_sn] = rec_ratio

                if root_abs_prof is not None:
                    deg_seq = degree_sequence(adj)
                    rec_abs = make_extremal_record(n, root, s_total, g6, deg_seq, root_abs_prof)

                    cur_abs_s = best_abs_by_s.get(s_total)
                    if cur_abs_s is None or rec_abs["abs_margin"] < cur_abs_s["abs_margin"]:
                        best_abs_by_s[s_total] = rec_abs

                    cur_abs_sn = best_abs_by_sn.get(key_sn)
                    if cur_abs_sn is None or rec_abs["abs_margin"] < cur_abs_sn["abs_margin"]:
                        best_abs_by_sn[key_sn] = rec_abs

    proc.wait()

    return {
        "n": n,
        "trees": trees,
        "support_vertices": support_vertices,
        "incremental_steps": incremental_steps,
        "w_checks": w_checks,
        "w_fail_count": w_fail_count,
        "min_w_margin": min_w_margin,
        "failures": failures,
        "root_count_by_s": dict(root_count_by_s),
        "roots_with_neg_term2_by_s": dict(roots_with_neg_term2_by_s),
        "best_ratio_by_s": best_ratio_by_s,
        "best_abs_by_s": best_abs_by_s,
        "best_ratio_by_sn": {f"{s}|{nn}": rec for (s, nn), rec in best_ratio_by_sn.items()},
        "best_abs_by_sn": {f"{s}|{nn}": rec for (s, nn), rec in best_abs_by_sn.items()},
    }


def merge_partition_results(parts: list[dict]) -> dict:
    out = {
        "trees": 0,
        "support_vertices": 0,
        "incremental_steps": 0,
        "w_checks": 0,
        "w_fail_count": 0,
        "min_w_margin": None,
        "failures": [],
        "root_count_by_s": defaultdict(int),
        "roots_with_neg_term2_by_s": defaultdict(int),
        "best_ratio_by_s": {},
        "best_abs_by_s": {},
        "best_ratio_by_sn": {},
        "best_abs_by_sn": {},
    }

    for p in parts:
        out["trees"] += p["trees"]
        out["support_vertices"] += p["support_vertices"]
        out["incremental_steps"] += p["incremental_steps"]
        out["w_checks"] += p["w_checks"]
        out["w_fail_count"] += p["w_fail_count"]

        if p["min_w_margin"] is not None:
            if out["min_w_margin"] is None or p["min_w_margin"] < out["min_w_margin"]:
                out["min_w_margin"] = p["min_w_margin"]

        if len(out["failures"]) < 50:
            room = 50 - len(out["failures"])
            out["failures"].extend(p["failures"][:room])

        for s, c in p["root_count_by_s"].items():
            out["root_count_by_s"][int(s)] += c
        for s, c in p["roots_with_neg_term2_by_s"].items():
            out["roots_with_neg_term2_by_s"][int(s)] += c

        for s, rec in p["best_ratio_by_s"].items():
            ss = int(s)
            cur = out["best_ratio_by_s"].get(ss)
            if cur is None or rec["ratio"] < cur["ratio"]:
                out["best_ratio_by_s"][ss] = rec

        for s, rec in p["best_abs_by_s"].items():
            ss = int(s)
            cur = out["best_abs_by_s"].get(ss)
            if cur is None or rec["abs_margin"] < cur["abs_margin"]:
                out["best_abs_by_s"][ss] = rec

        for key, rec in p["best_ratio_by_sn"].items():
            cur = out["best_ratio_by_sn"].get(key)
            if cur is None or rec["ratio"] < cur["ratio"]:
                out["best_ratio_by_sn"][key] = rec

        for key, rec in p["best_abs_by_sn"].items():
            cur = out["best_abs_by_sn"].get(key)
            if cur is None or rec["abs_margin"] < cur["abs_margin"]:
                out["best_abs_by_sn"][key] = rec

    out["root_count_by_s"] = dict(out["root_count_by_s"])
    out["roots_with_neg_term2_by_s"] = dict(out["roots_with_neg_term2_by_s"])
    return out


def build_broom(m: int, t: int) -> tuple[int, list[list[int]], dict[int, str], list[int], list[int], list[int]]:
    n = 1 + m + m * t + m * t
    adj = [[] for _ in range(n)]
    labels: dict[int, str] = {}

    v = 0
    labels[v] = "v"

    nxt = 1
    w_ids: list[int] = []
    for i in range(m):
        w = nxt
        nxt += 1
        w_ids.append(w)
        labels[w] = f"w_{i}"
        adj[v].append(w)
        adj[w].append(v)

    x_ids: list[int] = []
    for i in range(m):
        for j in range(t):
            x = nxt
            nxt += 1
            x_ids.append(x)
            labels[x] = f"x_{i}_{j}"
            adj[w_ids[i]].append(x)
            adj[x].append(w_ids[i])

    y_ids: list[int] = []
    for i in range(m):
        for j in range(t):
            y = nxt
            nxt += 1
            y_ids.append(y)
            labels[y] = f"y_{i}_{j}"
            xi = i * t + j
            adj[x_ids[xi]].append(y)
            adj[y].append(x_ids[xi])

    assert nxt == n
    return n, adj, labels, x_ids, w_ids, y_ids


def check_poly_dominance(f: list[int], g: list[int]) -> tuple[bool, int]:
    min_margin = None
    max_k = max(len(f), len(g)) + 1
    for k in range(max_k):
        d = delta_k(f, g, k)
        if min_margin is None or d < min_margin:
            min_margin = d
        if d < 0:
            return False, min_margin
    assert min_margin is not None
    return True, min_margin


def check_poly_lc(p: list[int]) -> tuple[bool, int]:
    min_gap = None
    for k in range(len(p)):
        g = c_k(p, k)
        if min_gap is None or g < min_gap:
            min_gap = g
        if g < 0:
            return False, min_gap
    assert min_gap is not None
    return True, min_gap


def factor_checks_for_broom(m: int, t: int) -> dict:
    n, adj, labels, _, _, _ = build_broom(m, t)
    deg = [len(adj[v]) for v in range(n)]

    total_factors = 0
    fail_counts = {
        "E_lc": 0,
        "SCC": 0,
        "E_ge_J": 0,
        "leaf_aug": 0,
    }
    roots: list[dict] = []

    for root in range(n):
        if not any(deg[u] == 1 for u in adj[root]):
            continue

        children_root, dp0, dp1s = dp_rooted(n, adj, root)
        nonleaf_children = [c for c in children_root if len([u for u in adj[c] if u != root]) > 0]

        root_rec = {
            "root": root,
            "label": labels[root],
            "degree": deg[root],
            "nonleaf_children": len(nonleaf_children),
            "factors": [],
        }

        for c in nonleaf_children:
            total_factors += 1
            e_t = dp0[c]
            j_t = dp1s[c]
            i_t = _polyadd(e_t, xshift(j_t))
            scc_poly = _polyadd(i_t, xshift(i_t))
            leaf_aug_poly = _polyadd(i_t, xshift(e_t))

            ok_lc, lc_margin = check_poly_lc(e_t)
            ok_scc, scc_margin = check_poly_dominance(scc_poly, e_t)
            ok_ej, ej_margin = check_poly_dominance(e_t, j_t)
            ok_la, la_margin = check_poly_dominance(leaf_aug_poly, e_t)

            if not ok_lc:
                fail_counts["E_lc"] += 1
            if not ok_scc:
                fail_counts["SCC"] += 1
            if not ok_ej:
                fail_counts["E_ge_J"] += 1
            if not ok_la:
                fail_counts["leaf_aug"] += 1

            root_rec["factors"].append(
                {
                    "child": c,
                    "label": labels[c],
                    "E_lc": ok_lc,
                    "E_lc_min_gap": lc_margin,
                    "SCC": ok_scc,
                    "SCC_min_margin": scc_margin,
                    "E_ge_J": ok_ej,
                    "E_ge_J_min_margin": ej_margin,
                    "leaf_aug": ok_la,
                    "leaf_aug_min_margin": la_margin,
                }
            )

        roots.append(root_rec)

    return {
        "tree": f"T_{{{m},{t}}}",
        "n": n,
        "support_roots": len(roots),
        "total_factors": total_factors,
        "fail_counts": fail_counts,
        "roots": roots,
    }


def all_rootings_egej_broom(m: int, t: int) -> dict:
    n, adj, labels, _, _, _ = build_broom(m, t)
    deg = [len(adj[v]) for v in range(n)]
    records: list[dict] = []
    fail_roots: list[dict] = []

    for root in range(n):
        _, dp0, dp1s = dp_rooted(n, adj, root)
        e = dp0[root]
        j = dp1s[root]

        ok, min_margin = check_poly_dominance(e, j)
        is_support = any(deg[u] == 1 for u in adj[root])
        label = labels[root]
        kind = label.split("_", 1)[0]
        rec = {
            "root": root,
            "label": label,
            "kind": kind,
            "degree": deg[root],
            "is_support": is_support,
            "E_ge_J": ok,
            "min_margin": min_margin,
        }
        records.append(rec)
        if not ok:
            fail_roots.append(rec)

    by_kind = defaultdict(lambda: {"total": 0, "fail": 0, "min_margin": None})
    for r in records:
        k = r["kind"]
        by_kind[k]["total"] += 1
        if not r["E_ge_J"]:
            by_kind[k]["fail"] += 1
        mm = by_kind[k]["min_margin"]
        if mm is None or r["min_margin"] < mm:
            by_kind[k]["min_margin"] = r["min_margin"]

    return {
        "tree": f"T_{{{m},{t}}}",
        "n": n,
        "all_rootings": len(records),
        "fail_count": len(fail_roots),
        "fail_roots": fail_roots,
        "by_kind": dict(by_kind),
        "records": records,
    }


def run_scan(
    task1_min_n: int,
    task1_max_n: int,
    profile_max_n: int,
    workers: int,
) -> dict:
    n_max = max(task1_max_n, profile_max_n)

    task1_by_n: dict[int, dict] = {}

    # Task 2 accumulators (global across n<=profile_max_n)
    global_root_count_by_s: dict[int, int] = defaultdict(int)
    global_neg_by_s: dict[int, int] = defaultdict(int)
    global_best_ratio_by_s: dict[int, dict] = {}
    global_best_abs_by_s: dict[int, dict] = {}
    global_best_ratio_by_sn: dict[str, dict] = {}
    global_best_abs_by_sn: dict[str, dict] = {}

    t0 = time.time()
    with mp.Pool(processes=workers) as pool:
        for n in range(3, n_max + 1):
            if not ((task1_min_n <= n <= task1_max_n) or (n <= profile_max_n)):
                continue

            args_list = [(n, r, workers, task1_min_n, task1_max_n, profile_max_n) for r in range(workers)]
            parts = pool.map(scan_partition, args_list)
            merged = merge_partition_results(parts)

            if task1_min_n <= n <= task1_max_n:
                task1_by_n[n] = {
                    "trees": merged["trees"],
                    "support_vertices": merged["support_vertices"],
                    "incremental_steps": merged["incremental_steps"],
                    "w_checks": merged["w_checks"],
                    "w_fail_count": merged["w_fail_count"],
                    "min_w_margin": merged["min_w_margin"],
                    "failures": merged["failures"],
                }

            if n <= profile_max_n:
                for s, c in merged["root_count_by_s"].items():
                    global_root_count_by_s[int(s)] += c
                for s, c in merged["roots_with_neg_term2_by_s"].items():
                    global_neg_by_s[int(s)] += c

                for s, rec in merged["best_ratio_by_s"].items():
                    ss = int(s)
                    cur = global_best_ratio_by_s.get(ss)
                    if cur is None or rec["ratio"] < cur["ratio"]:
                        global_best_ratio_by_s[ss] = rec

                for s, rec in merged["best_abs_by_s"].items():
                    ss = int(s)
                    cur = global_best_abs_by_s.get(ss)
                    if cur is None or rec["abs_margin"] < cur["abs_margin"]:
                        global_best_abs_by_s[ss] = rec

                for key, rec in merged["best_ratio_by_sn"].items():
                    cur = global_best_ratio_by_sn.get(key)
                    if cur is None or rec["ratio"] < cur["ratio"]:
                        global_best_ratio_by_sn[key] = rec

                for key, rec in merged["best_abs_by_sn"].items():
                    cur = global_best_abs_by_sn.get(key)
                    if cur is None or rec["abs_margin"] < cur["abs_margin"]:
                        global_best_abs_by_sn[key] = rec

            elapsed = time.time() - t0
            print(
                f"n={n}: trees={merged['trees']:,} support={merged['support_vertices']:,} "
                f"steps={merged['incremental_steps']:,} Wchecks={merged['w_checks']:,} "
                f"Wfail={merged['w_fail_count']:,} minW={merged['min_w_margin']} "
                f"elapsed={elapsed:.1f}s",
                flush=True,
            )

    # Task 1 total
    task1_total = {
        "trees": sum(v["trees"] for v in task1_by_n.values()),
        "support_vertices": sum(v["support_vertices"] for v in task1_by_n.values()),
        "incremental_steps": sum(v["incremental_steps"] for v in task1_by_n.values()),
        "w_checks": sum(v["w_checks"] for v in task1_by_n.values()),
        "w_fail_count": sum(v["w_fail_count"] for v in task1_by_n.values()),
        "min_w_margin": min((v["min_w_margin"] for v in task1_by_n.values() if v["min_w_margin"] is not None), default=None),
        "failures": [],
    }
    for n in sorted(task1_by_n):
        task1_total["failures"].extend(task1_by_n[n]["failures"])
        if len(task1_total["failures"]) >= 50:
            task1_total["failures"] = task1_total["failures"][:50]
            break

    # Task 2 by-s table
    task2_by_s: dict[int, dict] = {}
    all_s = sorted(set(global_root_count_by_s.keys()) | set(global_best_ratio_by_s.keys()) | set(global_best_abs_by_s.keys()))
    for s in all_s:
        task2_by_s[s] = {
            "roots_total": global_root_count_by_s.get(s, 0),
            "roots_with_neg_term2": global_neg_by_s.get(s, 0),
            "best_ratio_record": global_best_ratio_by_s.get(s),
            "best_abs_margin_record": global_best_abs_by_s.get(s),
        }

    # Task 2 by-(s,n)
    task2_by_sn: dict[str, dict] = {}
    for key in sorted(set(global_best_ratio_by_sn.keys()) | set(global_best_abs_by_sn.keys())):
        task2_by_sn[key] = {
            "best_ratio_record": global_best_ratio_by_sn.get(key),
            "best_abs_margin_record": global_best_abs_by_sn.get(key),
        }

    return {
        "task1": {
            "range": [task1_min_n, task1_max_n],
            "by_n": task1_by_n,
            "total": task1_total,
        },
        "task2": {
            "profile_max_n": profile_max_n,
            "by_s": task2_by_s,
            "by_s_n": task2_by_sn,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task1-min-n", type=int, default=21)
    parser.add_argument("--task1-max-n", type=int, default=22)
    parser.add_argument("--profile-max-n", type=int, default=20)
    parser.add_argument("--workers", type=int, default=max(1, min(8, (os.cpu_count() or 8))))
    parser.add_argument("--out", type=str, default="results/round10_wform_verification.json")
    args = parser.parse_args()

    t0 = time.time()

    scan = run_scan(
        task1_min_n=args.task1_min_n,
        task1_max_n=args.task1_max_n,
        profile_max_n=args.profile_max_n,
        workers=args.workers,
    )

    factor_34 = factor_checks_for_broom(3, 4)
    factor_33 = factor_checks_for_broom(3, 3)
    rootings_34 = all_rootings_egej_broom(3, 4)

    out = {
        "meta": {
            "script": "round10_wform_verification.py",
            "geng": GENG,
            "workers": args.workers,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "elapsed_seconds": time.time() - t0,
        },
        "scan": scan,
        "task3": {
            "T_3_4": factor_34,
            "T_3_3": factor_33,
        },
        "task4": {
            "T_3_4_all_rootings_E_ge_J": rootings_34,
        },
    }

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)

    print("\nSaved:", args.out)
    print("Elapsed: %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
