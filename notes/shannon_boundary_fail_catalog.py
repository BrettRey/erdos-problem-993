#!/usr/bin/env python3
"""Catalog Lamport boundary (k=d) failures by child arity.

A boundary failure is:
  Delta(PV)_d > -Delta(IU)_d,
where d is the first descent index of IU.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import shutil
import sys
import time
from collections import Counter, defaultdict
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
) -> tuple[list[list[int]], list[list[int]], list[list[int]], list[int]]:
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
    subtree_size = [0] * n

    for v in reversed(order):
        if not children[v]:
            f[v] = [1]
            g[v] = [0, 1]
            subtree_size[v] = 1
            continue

        prod_h = [1]
        prod_f = [1]
        size_v = 1
        for c in children[v]:
            prod_h = poly_mul(prod_h, poly_add(f[c], g[c]))
            prod_f = poly_mul(prod_f, f[c])
            size_v += subtree_size[c]

        f[v] = prod_h
        g[v] = [0] + prod_f
        subtree_size[v] = size_v

    return children, f, g, subtree_size


def _iter_trees(n: int, backend: str):
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


def _default_arity_record() -> dict[str, Any]:
    return {
        "failure_count": 0,
        "by_n": {},
        "unique_trees": {},
        "first_witness": None,
        "max_gap": None,
        "max_gap_witness": None,
    }


def _bump_int_key(d: dict[str, int], key: int, inc: int = 1) -> None:
    k = str(key)
    d[k] = d.get(k, 0) + inc


def scan(max_n: int, backend: str, focus_arities: set[int], top_k: int) -> dict[str, Any]:
    totals = {
        "trees_scanned": 0,
        "steps_with_descent": 0,
        "boundary_failures_total": 0,
        "boundary_failures_by_arity": {},
    }

    focus: dict[int, dict[str, Any]] = {a: _default_arity_record() for a in focus_arities}

    t0 = time.time()
    for n in range(1, max_n + 1):
        n_failures = 0
        n_t0 = time.time()
        for tree_index, (n_out, adj, graph6) in enumerate(_iter_trees(n, backend), start=1):
            totals["trees_scanned"] += 1
            for root in range(n_out):
                children, f, g, subtree_size = rooted_dp_all(adj, root)
                for v in range(n_out):
                    P = [1]
                    Q = [0, 1]
                    for child_pos, u in enumerate(children[v]):
                        U = f[u]
                        V = g[u]
                        I = poly_add(P, Q)
                        IU = poly_mul(I, U)
                        PV = poly_mul(P, V)

                        d = first_descent(IU)
                        if d is not None:
                            totals["steps_with_descent"] += 1
                            dIU = deltas(IU)
                            dPV = deltas(PV)
                            L = max(len(dIU), len(dPV))
                            dIU += [0] * (L - len(dIU))
                            dPV += [0] * (L - len(dPV))

                            if dPV[d] > -dIU[d]:
                                n_failures += 1
                                totals["boundary_failures_total"] += 1
                                ar = len(children[u])
                                _bump_int_key(totals["boundary_failures_by_arity"], ar)

                                if ar in focus:
                                    rec = focus[ar]
                                    rec["failure_count"] += 1
                                    _bump_int_key(rec["by_n"], n_out)
                                    tree_key = f"{n_out}:{tree_index}:{graph6}"
                                    rec["unique_trees"][tree_key] = rec["unique_trees"].get(tree_key, 0) + 1

                                    witness = {
                                        "n": n_out,
                                        "tree_index": tree_index,
                                        "graph6": graph6,
                                        "root": root,
                                        "parent": v,
                                        "child": u,
                                        "child_pos": child_pos,
                                        "child_arity": ar,
                                        "child_size": subtree_size[u],
                                        "d": d,
                                        "delta_PV_d": dPV[d],
                                        "delta_IU_d": dIU[d],
                                    }

                                    if rec["first_witness"] is None:
                                        rec["first_witness"] = witness

                                    gap = dPV[d] + dIU[d]
                                    if rec["max_gap"] is None or gap > rec["max_gap"]:
                                        rec["max_gap"] = gap
                                        w = dict(witness)
                                        w["gap"] = gap
                                        rec["max_gap_witness"] = w

                        P = poly_mul(P, poly_add(U, V))
                        Q = poly_mul(Q, U)

        print(
            f"[n={n}] failures={n_failures} elapsed={time.time()-n_t0:.2f}s",
            flush=True,
        )

    out_focus: dict[str, Any] = {}
    for ar in sorted(focus):
        rec = focus[ar]
        uniq_counter = Counter(rec["unique_trees"])
        top_trees = []
        for key, count in uniq_counter.most_common(top_k):
            n_s, ti_s, g6 = key.split(":", 2)
            top_trees.append(
                {
                    "n": int(n_s),
                    "tree_index": int(ti_s),
                    "graph6": g6,
                    "failure_count": count,
                }
            )
        out_focus[str(ar)] = {
            "failure_count": rec["failure_count"],
            "by_n": rec["by_n"],
            "unique_tree_count": len(rec["unique_trees"]),
            "top_trees": top_trees,
            "first_witness": rec["first_witness"],
            "max_gap": rec["max_gap"],
            "max_gap_witness": rec["max_gap_witness"],
        }

    return {
        "elapsed_s": round(time.time() - t0, 3),
        "totals": totals,
        "focus_arities": out_focus,
    }


def parse_int_list(s: str) -> list[int]:
    s = s.strip()
    if not s:
        return []
    return [int(x) for x in s.split(",") if x.strip()]


def main() -> None:
    ap = argparse.ArgumentParser(description="Catalog Lamport boundary failures by child arity")
    ap.add_argument("--max-n", type=int, default=17)
    ap.add_argument("--backend", choices=["auto", "geng", "networkx"], default="auto")
    ap.add_argument("--focus-arities", type=str, default="4,5")
    ap.add_argument("--top-k", type=int, default=12)
    ap.add_argument("--out", type=str, required=True)
    args = ap.parse_args()

    focus_arities = set(parse_int_list(args.focus_arities))
    report = {
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "description": "Boundary failure catalog (k=d): Delta(PV)_d > -Delta(IU)_d",
        "definitions": {
            "d": "first descent index of IU",
            "gap": "Delta(PV)_d + Delta(IU)_d (positive iff boundary fails)",
        },
        "parameters": {
            "max_n": args.max_n,
            "backend": args.backend,
            "focus_arities": sorted(focus_arities),
            "top_k": args.top_k,
        },
    }

    report["scan"] = scan(args.max_n, args.backend, focus_arities, args.top_k)

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
