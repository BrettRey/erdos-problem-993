#!/usr/bin/env python3
"""Scan Lamport slacks at k=d and k=d+1 by child arity.

Also scans the leaf-child quantitative bound:
  Delta(PV)_{d+1} <= (1/8) * (-Delta(IU)_{d+1}).
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import shutil
import sys
import time
from typing import Any

# Allow imports from repository root when running as `python notes/...`.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from trees import trees, trees_geng_raw


def _trim(poly: list[int]) -> list[int]:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def poly_add(a: list[int], b: list[int]) -> list[int]:
    out = [0] * max(len(a), len(b))
    for i, v in enumerate(a):
        out[i] += v
    for i, v in enumerate(b):
        out[i] += v
    return _trim(out)


def poly_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0:
                continue
            out[i + j] += ai * bj
    return _trim(out)


def deltas(seq: list[int]) -> list[int]:
    return [seq[i + 1] - seq[i] for i in range(len(seq) - 1)]


def first_descent(seq: list[int]) -> int | None:
    for k in range(len(seq) - 1):
        if seq[k + 1] < seq[k]:
            return k
    return None


def rooted_dp_all(adj: list[list[int]], root: int) -> tuple[list[list[int]], list[list[int]], list[list[int]], list[int]]:
    n = len(adj)
    parent = [-1] * n
    parent[root] = root
    order = [root]
    for v in order:
        for w in adj[v]:
            if parent[w] == -1:
                parent[w] = v
                order.append(w)

    children: list[list[int]] = [[] for _ in range(n)]
    for v in order:
        if v != root:
            children[parent[v]].append(v)

    f: list[list[int]] = [[] for _ in range(n)]
    g: list[list[int]] = [[] for _ in range(n)]
    sz = [0] * n

    for v in reversed(order):
        if not children[v]:
            f[v] = [1]
            g[v] = [0, 1]
            sz[v] = 1
            continue

        prod_h = [1]
        prod_f = [1]
        size_v = 1
        for c in children[v]:
            prod_h = poly_mul(prod_h, poly_add(f[c], g[c]))
            prod_f = poly_mul(prod_f, f[c])
            size_v += sz[c]

        f[v] = prod_h
        g[v] = [0] + prod_f
        sz[v] = size_v

    return children, f, g, sz


def _arity_record_init() -> dict[str, Any]:
    return {
        "count": 0,
        "min_s0": None,
        "min_s0_witness": None,
        "min_s1": None,
        "min_s1_witness": None,
    }


def _arity_update_min(
    rec: dict[str, Any],
    field: str,
    witness_field: str,
    value: int,
    witness: dict[str, Any],
) -> None:
    cur = rec[field]
    if cur is None or value < cur:
        rec[field] = value
        rec[witness_field] = witness


def _empty_scan_record(label: str) -> dict[str, Any]:
    return {
        "label": label,
        "trees_scanned": 0,
        "steps_with_descent": 0,
        "by_child_arity": {},
        "leaf_q18": {
            "checks": 0,
            "violations": 0,
            "equalities": 0,
            "max_ratio": -1.0,
            "max_ratio_witness": None,
            "min_margin": None,
            "min_margin_witness": None,
        },
    }


def _update_one_step(
    out: dict[str, Any],
    *,
    n: int,
    tree_index: int,
    graph6: str | None,
    root: int,
    parent_v: int,
    child_u: int,
    child_pos: int,
    child_arity: int,
    child_size: int,
    P: list[int],
    Q: list[int],
    U: list[int],
    V: list[int],
) -> None:
    I = poly_add(P, Q)
    IU = poly_mul(I, U)
    PV = poly_mul(P, V)

    d = first_descent(IU)
    if d is None:
        return

    out["steps_with_descent"] += 1

    dIU = deltas(IU)
    dPV = deltas(PV)
    L = max(len(dIU), len(dPV))
    dIU += [0] * (L - len(dIU))
    dPV += [0] * (L - len(dPV))

    k = str(child_arity)
    by_ar = out["by_child_arity"]
    if k not in by_ar:
        by_ar[k] = _arity_record_init()
    rec = by_ar[k]
    rec["count"] += 1

    s0 = (-dIU[d]) - dPV[d]
    w0 = {
        "n": n,
        "tree_index": tree_index,
        "graph6": graph6,
        "root": root,
        "parent": parent_v,
        "child": child_u,
        "child_pos": child_pos,
        "child_arity": child_arity,
        "child_size": child_size,
        "d": d,
        "s0": s0,
        "delta_IU_d": dIU[d],
        "delta_PV_d": dPV[d],
    }
    _arity_update_min(rec, "min_s0", "min_s0_witness", s0, w0)

    if d + 1 < L:
        s1 = (-dIU[d + 1]) - dPV[d + 1]
        w1 = {
            "n": n,
            "tree_index": tree_index,
            "graph6": graph6,
            "root": root,
            "parent": parent_v,
            "child": child_u,
            "child_pos": child_pos,
            "child_arity": child_arity,
            "child_size": child_size,
            "d": d,
            "s1": s1,
            "delta_IU_dplus1": dIU[d + 1],
            "delta_PV_dplus1": dPV[d + 1],
        }
        _arity_update_min(rec, "min_s1", "min_s1_witness", s1, w1)

        # Leaf-child quantitative scan (Q1/8).
        if child_arity == 0:
            leaf = out["leaf_q18"]
            num = dPV[d + 1]
            den = -dIU[d + 1]
            if den > 0:
                leaf["checks"] += 1
                ratio = num / den
                if ratio > leaf["max_ratio"]:
                    leaf["max_ratio"] = ratio
                    leaf["max_ratio_witness"] = {
                        "n": n,
                        "tree_index": tree_index,
                        "graph6": graph6,
                        "root": root,
                        "parent": parent_v,
                        "child": child_u,
                        "child_pos": child_pos,
                        "d": d,
                        "num": num,
                        "den": den,
                        "ratio": ratio,
                    }
                margin = den - 8 * num
                if leaf["min_margin"] is None or margin < leaf["min_margin"]:
                    leaf["min_margin"] = margin
                    leaf["min_margin_witness"] = {
                        "n": n,
                        "tree_index": tree_index,
                        "graph6": graph6,
                        "root": root,
                        "parent": parent_v,
                        "child": child_u,
                        "child_pos": child_pos,
                        "d": d,
                        "num": num,
                        "den": den,
                        "margin": margin,
                    }
                if margin < 0:
                    leaf["violations"] += 1
                if margin == 0:
                    leaf["equalities"] += 1


def _scan_tree(out: dict[str, Any], n: int, adj: list[list[int]], tree_index: int, graph6: str | None) -> None:
    out["trees_scanned"] += 1

    for root in range(n):
        children, f, g, sz = rooted_dp_all(adj, root)

        for v in range(n):
            P = [1]
            Q = [0, 1]
            for child_pos, c in enumerate(children[v]):
                U = f[c]
                V = g[c]
                _update_one_step(
                    out,
                    n=n,
                    tree_index=tree_index,
                    graph6=graph6,
                    root=root,
                    parent_v=v,
                    child_u=c,
                    child_pos=child_pos,
                    child_arity=len(children[c]),
                    child_size=sz[c],
                    P=P,
                    Q=Q,
                    U=U,
                    V=V,
                )
                P = poly_mul(P, poly_add(U, V))
                Q = poly_mul(Q, U)


def _iter_trees_exhaustive(n: int, backend: str):
    if n == 1:
        yield 1, [[]], None
        return

    resolved = backend
    if resolved == "auto":
        resolved = "geng" if shutil.which("geng") else "networkx"

    if resolved == "geng":
        for n_out, adj, raw in trees_geng_raw(n):
            yield n_out, adj, raw.decode("ascii")
    else:
        for n_out, adj in trees(n, backend=resolved):
            yield n_out, adj, None


def _iter_trees_sampled(n: int, mod: int, per_partition: int):
    for res in range(mod):
        count = 0
        for n_out, adj, raw in trees_geng_raw(n, res=res, mod=mod):
            yield n_out, adj, raw.decode("ascii")
            count += 1
            if count >= per_partition:
                break


def scan_exhaustive(min_n: int, max_n: int, backend: str) -> dict[str, Any]:
    out = _empty_scan_record("exhaustive")
    out["range"] = {"min_n": min_n, "max_n": max_n}
    out["backend"] = backend

    t0 = time.time()
    by_n: list[dict[str, Any]] = []
    for n in range(min_n, max_n + 1):
        nrec = _empty_scan_record(f"n={n}")
        n_t0 = time.time()
        for tree_idx, (n_out, adj, g6) in enumerate(_iter_trees_exhaustive(n, backend), start=1):
            _scan_tree(nrec, n_out, adj, tree_idx, g6)

        by_n.append({
            "n": n,
            "elapsed_s": round(time.time() - n_t0, 3),
            "trees_scanned": nrec["trees_scanned"],
            "steps_with_descent": nrec["steps_with_descent"],
            "leaf_q18": nrec["leaf_q18"],
            "by_child_arity": nrec["by_child_arity"],
        })

        # merge nrec into out
        out["trees_scanned"] += nrec["trees_scanned"]
        out["steps_with_descent"] += nrec["steps_with_descent"]

        # merge leaf stats
        ldst = out["leaf_q18"]
        lsrc = nrec["leaf_q18"]
        ldst["checks"] += lsrc["checks"]
        ldst["violations"] += lsrc["violations"]
        ldst["equalities"] += lsrc["equalities"]
        if lsrc["max_ratio"] > ldst["max_ratio"]:
            ldst["max_ratio"] = lsrc["max_ratio"]
            ldst["max_ratio_witness"] = lsrc["max_ratio_witness"]
        if ldst["min_margin"] is None or (lsrc["min_margin"] is not None and lsrc["min_margin"] < ldst["min_margin"]):
            ldst["min_margin"] = lsrc["min_margin"]
            ldst["min_margin_witness"] = lsrc["min_margin_witness"]

        # merge arity dict
        for k, rec in nrec["by_child_arity"].items():
            if k not in out["by_child_arity"]:
                out["by_child_arity"][k] = _arity_record_init()
            dst = out["by_child_arity"][k]
            dst["count"] += rec["count"]
            if dst["min_s0"] is None or (rec["min_s0"] is not None and rec["min_s0"] < dst["min_s0"]):
                dst["min_s0"] = rec["min_s0"]
                dst["min_s0_witness"] = rec["min_s0_witness"]
            if dst["min_s1"] is None or (rec["min_s1"] is not None and rec["min_s1"] < dst["min_s1"]):
                dst["min_s1"] = rec["min_s1"]
                dst["min_s1_witness"] = rec["min_s1_witness"]

        print(
            f"[exhaustive] n={n} trees={nrec['trees_scanned']} steps={nrec['steps_with_descent']} "
            f"leaf_max_ratio={nrec['leaf_q18']['max_ratio']:.6f} elapsed={time.time()-n_t0:.2f}s",
            flush=True,
        )

    out["elapsed_s"] = round(time.time() - t0, 3)
    out["by_n"] = by_n
    return out


def scan_sample(ns: list[int], mod: int, per_partition: int) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []

    for n in ns:
        rec = _empty_scan_record(f"sample_n={n}")
        t0 = time.time()
        for tree_idx, (n_out, adj, g6) in enumerate(_iter_trees_sampled(n, mod, per_partition), start=1):
            _scan_tree(rec, n_out, adj, tree_idx, g6)

        rec["n"] = n
        rec["sampling"] = {
            "mod": mod,
            "per_partition": per_partition,
        }
        rec["elapsed_s"] = round(time.time() - t0, 3)
        out.append(rec)

        print(
            f"[sample] n={n} trees={rec['trees_scanned']} steps={rec['steps_with_descent']} "
            f"leaf_max_ratio={rec['leaf_q18']['max_ratio']:.6f} elapsed={time.time()-t0:.2f}s",
            flush=True,
        )

    return out


def parse_int_list(s: str) -> list[int]:
    s = s.strip()
    if not s:
        return []
    return [int(x) for x in s.split(",") if x.strip()]


def main() -> None:
    ap = argparse.ArgumentParser(description="Shannon slack scan")
    ap.add_argument("--min-n", type=int, default=1)
    ap.add_argument("--max-n", type=int, default=17)
    ap.add_argument("--backend", choices=["auto", "geng", "networkx"], default="auto")
    ap.add_argument("--sample-n", type=str, default="")
    ap.add_argument("--sample-mod", type=int, default=16)
    ap.add_argument("--sample-per-partition", type=int, default=200)
    ap.add_argument("--out", type=str, required=True)
    args = ap.parse_args()

    report: dict[str, Any] = {
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "description": "Lamport slack scan: boundary slack s0, next-step slack s1, and leaf Q1/8 ratio",
        "definitions": {
            "d": "first descent index of IU",
            "s0": "(-Delta(IU)_d) - Delta(PV)_d",
            "s1": "(-Delta(IU)_{d+1}) - Delta(PV)_{d+1}",
            "leaf_q18_margin": "(-Delta(IU)_{d+1}) - 8*Delta(PV)_{d+1} for child arity 0",
        },
    }

    if args.min_n <= args.max_n:
        report["exhaustive"] = scan_exhaustive(args.min_n, args.max_n, args.backend)

    sample_ns = parse_int_list(args.sample_n)
    if sample_ns:
        report["stratified_samples"] = scan_sample(sample_ns, args.sample_mod, args.sample_per_partition)

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
