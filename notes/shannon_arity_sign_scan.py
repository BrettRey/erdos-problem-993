#!/usr/bin/env python3
"""Scan Lamport transition signs at k=d and k=d+1 by child arity threshold.

For each rooted composition step with
  I' = IU + PV,  d = first descent index of IU,
this script records sign counts of Delta(PV)_d and Delta(PV)_{d+1}
restricted to steps where child-root arity >= threshold.
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


def rooted_dp_all(
    adj: list[list[int]], root: int
) -> tuple[list[list[int]], list[list[int]], list[list[int]]]:
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

    for v in reversed(order):
        if not children[v]:
            f[v] = [1]
            g[v] = [0, 1]
            continue

        prod_h = [1]
        prod_f = [1]
        for c in children[v]:
            prod_h = poly_mul(prod_h, poly_add(f[c], g[c]))
            prod_f = poly_mul(prod_f, f[c])

        f[v] = prod_h
        g[v] = [0] + prod_f

    return children, f, g


def _empty_stats(threshold: int) -> dict[str, Any]:
    return {
        "arity_threshold": threshold,
        "trees_scanned": 0,
        "steps_with_descent": 0,
        "steps_meeting_threshold": 0,
        "dPV_d": {"positive": 0, "zero": 0, "negative": 0},
        "dPV_dplus1": {"positive": 0, "zero": 0, "negative": 0, "missing": 0},
        "s0_fail_count": 0,
        "s1_fail_count": 0,
        "first_positive_d_witness": None,
        "first_positive_dplus1_witness": None,
        "first_s0_fail_witness": None,
        "first_s1_fail_witness": None,
    }


def _scan_step(
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
    if child_arity < out["arity_threshold"]:
        return
    out["steps_meeting_threshold"] += 1

    dIU = deltas(IU)
    dPV = deltas(PV)
    L = max(len(dIU), len(dPV))
    dIU += [0] * (L - len(dIU))
    dPV += [0] * (L - len(dPV))

    witness_base = {
        "n": n,
        "tree_index": tree_index,
        "graph6": graph6,
        "root": root,
        "parent": parent_v,
        "child": child_u,
        "child_pos": child_pos,
        "child_arity": child_arity,
        "d": d,
    }

    x0 = dPV[d]
    if x0 > 0:
        out["dPV_d"]["positive"] += 1
        if out["first_positive_d_witness"] is None:
            w = dict(witness_base)
            w["delta_PV_d"] = x0
            w["delta_IU_d"] = dIU[d]
            out["first_positive_d_witness"] = w
    elif x0 == 0:
        out["dPV_d"]["zero"] += 1
    else:
        out["dPV_d"]["negative"] += 1

    if x0 > -dIU[d]:
        out["s0_fail_count"] += 1
        if out["first_s0_fail_witness"] is None:
            w = dict(witness_base)
            w["delta_PV_d"] = x0
            w["delta_IU_d"] = dIU[d]
            out["first_s0_fail_witness"] = w

    if d + 1 >= L:
        out["dPV_dplus1"]["missing"] += 1
        return

    x1 = dPV[d + 1]
    if x1 > 0:
        out["dPV_dplus1"]["positive"] += 1
        if out["first_positive_dplus1_witness"] is None:
            w = dict(witness_base)
            w["delta_PV_dplus1"] = x1
            w["delta_IU_dplus1"] = dIU[d + 1]
            out["first_positive_dplus1_witness"] = w
    elif x1 == 0:
        out["dPV_dplus1"]["zero"] += 1
    else:
        out["dPV_dplus1"]["negative"] += 1

    if x1 > -dIU[d + 1]:
        out["s1_fail_count"] += 1
        if out["first_s1_fail_witness"] is None:
            w = dict(witness_base)
            w["delta_PV_dplus1"] = x1
            w["delta_IU_dplus1"] = dIU[d + 1]
            out["first_s1_fail_witness"] = w


def _scan_tree(out: dict[str, Any], n: int, adj: list[list[int]], tree_index: int, graph6: str | None) -> None:
    out["trees_scanned"] += 1
    for root in range(n):
        children, f, g = rooted_dp_all(adj, root)
        for v in range(n):
            P = [1]
            Q = [0, 1]
            for child_pos, c in enumerate(children[v]):
                U = f[c]
                V = g[c]
                _scan_step(
                    out,
                    n=n,
                    tree_index=tree_index,
                    graph6=graph6,
                    root=root,
                    parent_v=v,
                    child_u=c,
                    child_pos=child_pos,
                    child_arity=len(children[c]),
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


def parse_int_list(s: str) -> list[int]:
    s = s.strip()
    if not s:
        return []
    return [int(x) for x in s.split(",") if x.strip()]


def scan_exhaustive(min_n: int, max_n: int, backend: str, threshold: int) -> dict[str, Any]:
    out = _empty_stats(threshold)
    out["range"] = {"min_n": min_n, "max_n": max_n}
    out["backend"] = backend

    by_n: list[dict[str, Any]] = []
    t0 = time.time()

    for n in range(min_n, max_n + 1):
        nrec = _empty_stats(threshold)
        n_t0 = time.time()
        for tree_idx, (n_out, adj, g6) in enumerate(_iter_trees_exhaustive(n, backend), start=1):
            _scan_tree(nrec, n_out, adj, tree_idx, g6)

        by_n.append(
            {
                "n": n,
                "elapsed_s": round(time.time() - n_t0, 3),
                "trees_scanned": nrec["trees_scanned"],
                "steps_meeting_threshold": nrec["steps_meeting_threshold"],
                "dPV_d": nrec["dPV_d"],
                "dPV_dplus1": nrec["dPV_dplus1"],
                "s0_fail_count": nrec["s0_fail_count"],
                "s1_fail_count": nrec["s1_fail_count"],
            }
        )

        out["trees_scanned"] += nrec["trees_scanned"]
        out["steps_with_descent"] += nrec["steps_with_descent"]
        out["steps_meeting_threshold"] += nrec["steps_meeting_threshold"]
        for key in ("positive", "zero", "negative"):
            out["dPV_d"][key] += nrec["dPV_d"][key]
        for key in ("positive", "zero", "negative", "missing"):
            out["dPV_dplus1"][key] += nrec["dPV_dplus1"][key]
        out["s0_fail_count"] += nrec["s0_fail_count"]
        out["s1_fail_count"] += nrec["s1_fail_count"]

        if out["first_positive_d_witness"] is None and nrec["first_positive_d_witness"] is not None:
            out["first_positive_d_witness"] = nrec["first_positive_d_witness"]
        if out["first_positive_dplus1_witness"] is None and nrec["first_positive_dplus1_witness"] is not None:
            out["first_positive_dplus1_witness"] = nrec["first_positive_dplus1_witness"]
        if out["first_s0_fail_witness"] is None and nrec["first_s0_fail_witness"] is not None:
            out["first_s0_fail_witness"] = nrec["first_s0_fail_witness"]
        if out["first_s1_fail_witness"] is None and nrec["first_s1_fail_witness"] is not None:
            out["first_s1_fail_witness"] = nrec["first_s1_fail_witness"]

        print(
            f"[exhaustive] n={n} trees={nrec['trees_scanned']} "
            f"steps_thr={nrec['steps_meeting_threshold']} "
            f"s0_fail={nrec['s0_fail_count']} s1_fail={nrec['s1_fail_count']} "
            f"elapsed={time.time()-n_t0:.2f}s",
            flush=True,
        )

    out["elapsed_s"] = round(time.time() - t0, 3)
    out["by_n"] = by_n
    return out


def scan_sample(ns: list[int], mod: int, per_partition: int, threshold: int) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []

    for n in ns:
        rec = _empty_stats(threshold)
        t0 = time.time()
        for tree_idx, (n_out, adj, g6) in enumerate(_iter_trees_sampled(n, mod, per_partition), start=1):
            _scan_tree(rec, n_out, adj, tree_idx, g6)

        rec["n"] = n
        rec["sampling"] = {"mod": mod, "per_partition": per_partition}
        rec["elapsed_s"] = round(time.time() - t0, 3)
        out.append(rec)

        print(
            f"[sample] n={n} trees={rec['trees_scanned']} "
            f"steps_thr={rec['steps_meeting_threshold']} "
            f"s0_fail={rec['s0_fail_count']} s1_fail={rec['s1_fail_count']} "
            f"elapsed={time.time()-t0:.2f}s",
            flush=True,
        )

    return out


def _combine_summary(exhaustive: dict[str, Any] | None, sampled: list[dict[str, Any]]) -> dict[str, Any]:
    out = _empty_stats(exhaustive["arity_threshold"] if exhaustive is not None else sampled[0]["arity_threshold"])
    out.pop("arity_threshold")

    def merge(src: dict[str, Any]) -> None:
        out["trees_scanned"] += src["trees_scanned"]
        out["steps_with_descent"] += src["steps_with_descent"]
        out["steps_meeting_threshold"] += src["steps_meeting_threshold"]
        for key in ("positive", "zero", "negative"):
            out["dPV_d"][key] += src["dPV_d"][key]
        for key in ("positive", "zero", "negative", "missing"):
            out["dPV_dplus1"][key] += src["dPV_dplus1"][key]
        out["s0_fail_count"] += src["s0_fail_count"]
        out["s1_fail_count"] += src["s1_fail_count"]
        for w in (
            "first_positive_d_witness",
            "first_positive_dplus1_witness",
            "first_s0_fail_witness",
            "first_s1_fail_witness",
        ):
            if out[w] is None and src[w] is not None:
                out[w] = src[w]

    if exhaustive is not None:
        merge(exhaustive)
    for rec in sampled:
        merge(rec)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Arity-threshold sign scan for Lamport transitions")
    ap.add_argument("--arity-threshold", type=int, default=6)
    ap.add_argument("--min-n", type=int, default=1)
    ap.add_argument("--max-n", type=int, default=17)
    ap.add_argument("--backend", choices=["auto", "geng", "networkx"], default="auto")
    ap.add_argument("--sample-n", type=str, default="")
    ap.add_argument("--sample-mod", type=int, default=16)
    ap.add_argument("--sample-per-partition", type=int, default=200)
    ap.add_argument("--out", type=str, required=True)
    args = ap.parse_args()

    sample_ns = parse_int_list(args.sample_n)

    report: dict[str, Any] = {
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "description": (
            "Lamport transition sign scan at k=d and k=d+1 "
            "for steps with child-root arity >= threshold"
        ),
        "definitions": {
            "d": "first descent index of IU",
            "s0_fail": "Delta(PV)_d > -Delta(IU)_d",
            "s1_fail": "Delta(PV)_{d+1} > -Delta(IU)_{d+1}",
        },
        "parameters": {
            "arity_threshold": args.arity_threshold,
            "min_n": args.min_n,
            "max_n": args.max_n,
            "backend": args.backend,
            "sample_n": sample_ns,
            "sample_mod": args.sample_mod,
            "sample_per_partition": args.sample_per_partition,
        },
    }

    exhaustive = None
    if args.min_n <= args.max_n:
        exhaustive = scan_exhaustive(args.min_n, args.max_n, args.backend, args.arity_threshold)
        report["exhaustive"] = exhaustive

    sampled: list[dict[str, Any]] = []
    if sample_ns:
        if not shutil.which("geng"):
            raise RuntimeError("sample scan requires geng backend availability")
        sampled = scan_sample(sample_ns, args.sample_mod, args.sample_per_partition, args.arity_threshold)
        report["stratified_samples"] = sampled

    if exhaustive is not None or sampled:
        report["combined_summary"] = _combine_summary(exhaustive, sampled)

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
