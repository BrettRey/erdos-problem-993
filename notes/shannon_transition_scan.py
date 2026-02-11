#!/usr/bin/env python3
"""Shannon-style scan of Lamport transition first-difference behavior.

For rooted composition steps
  P' = P (U + V),
  Q' = Q U,
  I' = IU + PV,
this script measures where tail failures can occur relative to
  d = first descent index of IU.

Primary empirical target:
  - Are there any k > d with Delta(I')_k > 0 ?
  - Are Lamport difference-dominance failures confined to k = d ?
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import shutil
import sys
import time
from dataclasses import dataclass
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
    """Return smallest k with seq[k+1] < seq[k], else None."""
    for k in range(len(seq) - 1):
        if seq[k + 1] < seq[k]:
            return k
    return None


def _padded(*arrs: list[int]) -> list[list[int]]:
    m = max((len(a) for a in arrs), default=0)
    out: list[list[int]] = []
    for a in arrs:
        out.append(a + [0] * (m - len(a)))
    return out


def rooted_dp_all(
    adj: list[list[int]], root: int
) -> tuple[list[list[int]], list[list[int]], list[list[int]], list[int]]:
    """Return (children, f, g, subtree_size) for rooted tree at root.

    f[v]: independent-set polynomial in rooted subtree of v with v excluded.
    g[v]: same with v included.
    """
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


@dataclass
class Stats:
    trees_scanned: int = 0
    lamport_steps: int = 0
    steps_no_descent: int = 0
    steps_with_descent: int = 0
    dd_all_hold: int = 0
    dd_failures: int = 0
    dd_failures_boundary_only: int = 0
    dd_failures_interior: int = 0
    positive_tail_after_boundary: int = 0
    pv_positive_after_d: int = 0
    pv_positive_after_d_plus_1: int = 0
    dd_boundary_fail: int = 0
    dd_boundary_hold: int = 0
    dd_boundary_fail_and_pv_positive_after_d: int = 0
    dd_boundary_hold_and_pv_positive_after_d: int = 0
    pv_positive_offset_hist: dict[str, int] | None = None
    profile_by_child_size: dict[str, dict[str, int]] | None = None
    profile_by_child_arity: dict[str, dict[str, int]] | None = None
    profile_by_deg_u: dict[str, dict[str, int]] | None = None
    dplus1_positive_num: int = 0
    dplus1_zero_den_positive_num: int = 0
    dplus1_ratio_max: float = -1.0
    dplus1_ratio_max_witness: dict[str, Any] | None = None
    first_boundary_witness: dict[str, Any] | None = None
    first_interior_witness: dict[str, Any] | None = None
    first_positive_tail_witness: dict[str, Any] | None = None
    first_pv_positive_after_d_witness: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        if self.pv_positive_offset_hist is None:
            self.pv_positive_offset_hist = {}
        if self.profile_by_child_size is None:
            self.profile_by_child_size = {}
        if self.profile_by_child_arity is None:
            self.profile_by_child_arity = {}
        if self.profile_by_deg_u is None:
            self.profile_by_deg_u = {}
        return {
            "trees_scanned": self.trees_scanned,
            "lamport_steps": self.lamport_steps,
            "steps_no_descent": self.steps_no_descent,
            "steps_with_descent": self.steps_with_descent,
            "dd_all_hold": self.dd_all_hold,
            "dd_failures": self.dd_failures,
            "dd_failures_boundary_only": self.dd_failures_boundary_only,
            "dd_failures_interior": self.dd_failures_interior,
            "positive_tail_after_boundary": self.positive_tail_after_boundary,
            "pv_positive_after_d": self.pv_positive_after_d,
            "pv_positive_after_d_plus_1": self.pv_positive_after_d_plus_1,
            "dd_boundary_fail": self.dd_boundary_fail,
            "dd_boundary_hold": self.dd_boundary_hold,
            "dd_boundary_fail_and_pv_positive_after_d": self.dd_boundary_fail_and_pv_positive_after_d,
            "dd_boundary_hold_and_pv_positive_after_d": self.dd_boundary_hold_and_pv_positive_after_d,
            "pv_positive_offset_hist": self.pv_positive_offset_hist,
            "profile_by_child_size": self.profile_by_child_size,
            "profile_by_child_arity": self.profile_by_child_arity,
            "profile_by_deg_u": self.profile_by_deg_u,
            "dplus1_positive_num": self.dplus1_positive_num,
            "dplus1_zero_den_positive_num": self.dplus1_zero_den_positive_num,
            "dplus1_ratio_max": self.dplus1_ratio_max,
            "dplus1_ratio_max_witness": self.dplus1_ratio_max_witness,
            "first_boundary_witness": self.first_boundary_witness,
            "first_interior_witness": self.first_interior_witness,
            "first_positive_tail_witness": self.first_positive_tail_witness,
            "first_pv_positive_after_d_witness": self.first_pv_positive_after_d_witness,
        }


def _profile_bump(
    profile: dict[str, dict[str, int]] | None,
    key: int,
    field: str,
) -> dict[str, dict[str, int]]:
    if profile is None:
        profile = {}
    k = str(key)
    if k not in profile:
        profile[k] = {"steps": 0, "boundary_fail": 0, "pv_positive_after_d": 0}
    profile[k][field] = profile[k].get(field, 0) + 1
    return profile


def _scan_step(
    stats: Stats,
    *,
    n: int,
    tree_index: int,
    graph6: str | None,
    root: int,
    parent_v: int,
    child_u: int,
    child_pos: int,
    child_size: int,
    child_arity: int,
    P: list[int],
    Q: list[int],
    U: list[int],
    V: list[int],
) -> None:
    I = poly_add(P, Q)
    IU = poly_mul(I, U)
    PV = poly_mul(P, V)
    Iprime = poly_add(IU, PV)

    stats.lamport_steps += 1

    d = first_descent(IU)
    if d is None:
        stats.steps_no_descent += 1
        return

    stats.steps_with_descent += 1

    dIU, dPV, dIp = _padded(deltas(IU), deltas(PV), deltas(Iprime))
    L = len(dIU)

    deg_u = len(U) - 1
    stats.profile_by_child_size = _profile_bump(stats.profile_by_child_size, child_size, "steps")
    stats.profile_by_child_arity = _profile_bump(stats.profile_by_child_arity, child_arity, "steps")
    stats.profile_by_deg_u = _profile_bump(stats.profile_by_deg_u, deg_u, "steps")

    if d + 1 < L:
        num = dPV[d + 1]
        den = -dIU[d + 1]
        if num > 0:
            stats.dplus1_positive_num += 1
            if den == 0:
                stats.dplus1_zero_den_positive_num += 1
        if den > 0:
            ratio = num / den
            if ratio > stats.dplus1_ratio_max:
                stats.dplus1_ratio_max = ratio
                stats.dplus1_ratio_max_witness = {
                    "n": n,
                    "tree_index": tree_index,
                    "graph6": graph6,
                    "root": root,
                    "parent": parent_v,
                    "child": child_u,
                    "child_pos": child_pos,
                    "child_size": child_size,
                    "child_arity": child_arity,
                    "deg_u": deg_u,
                    "d": d,
                    "num": num,
                    "den": den,
                    "ratio": ratio,
                    "delta_IU_dplus1": dIU[d + 1],
                    "delta_PV_dplus1": dPV[d + 1],
                }

    boundary_fails = dPV[d] > -dIU[d]
    if boundary_fails:
        stats.dd_boundary_fail += 1
        stats.profile_by_child_size = _profile_bump(stats.profile_by_child_size, child_size, "boundary_fail")
        stats.profile_by_child_arity = _profile_bump(stats.profile_by_child_arity, child_arity, "boundary_fail")
        stats.profile_by_deg_u = _profile_bump(stats.profile_by_deg_u, deg_u, "boundary_fail")
    else:
        stats.dd_boundary_hold += 1

    pv_pos = [k for k in range(d + 1, L) if dPV[k] > 0]
    if pv_pos:
        stats.pv_positive_after_d += 1
        stats.profile_by_child_size = _profile_bump(stats.profile_by_child_size, child_size, "pv_positive_after_d")
        stats.profile_by_child_arity = _profile_bump(stats.profile_by_child_arity, child_arity, "pv_positive_after_d")
        stats.profile_by_deg_u = _profile_bump(stats.profile_by_deg_u, deg_u, "pv_positive_after_d")
        if any(k > d + 1 for k in pv_pos):
            stats.pv_positive_after_d_plus_1 += 1
        if stats.pv_positive_offset_hist is None:
            stats.pv_positive_offset_hist = {}
        for k in pv_pos:
            off = str(k - d)
            stats.pv_positive_offset_hist[off] = stats.pv_positive_offset_hist.get(off, 0) + 1
        if boundary_fails:
            stats.dd_boundary_fail_and_pv_positive_after_d += 1
        else:
            stats.dd_boundary_hold_and_pv_positive_after_d += 1
        if stats.first_pv_positive_after_d_witness is None:
            stats.first_pv_positive_after_d_witness = {
                "n": n,
                "tree_index": tree_index,
                "graph6": graph6,
                "root": root,
                "parent": parent_v,
                "child": child_u,
                "child_pos": child_pos,
                "d": d,
                "positive_k": pv_pos,
                "delta_IU": dIU,
                "delta_PV": dPV,
                "delta_Iprime": dIp,
                "boundary_fails": boundary_fails,
            }

    violations = [k for k in range(d, L) if dPV[k] > -dIU[k]]

    if not violations:
        stats.dd_all_hold += 1
    else:
        stats.dd_failures += 1
        witness = {
            "n": n,
            "tree_index": tree_index,
            "graph6": graph6,
            "root": root,
            "parent": parent_v,
            "child": child_u,
            "child_pos": child_pos,
            "d": d,
            "violations": violations,
            "delta_IU": dIU,
            "delta_PV": dPV,
            "delta_Iprime": dIp,
        }
        if all(k == d for k in violations):
            stats.dd_failures_boundary_only += 1
            if stats.first_boundary_witness is None:
                stats.first_boundary_witness = witness
        else:
            stats.dd_failures_interior += 1
            if stats.first_interior_witness is None:
                stats.first_interior_witness = witness

    pos_after = [k for k in range(d + 1, L) if dIp[k] > 0]
    if pos_after:
        stats.positive_tail_after_boundary += 1
        if stats.first_positive_tail_witness is None:
            stats.first_positive_tail_witness = {
                "n": n,
                "tree_index": tree_index,
                "graph6": graph6,
                "root": root,
                "parent": parent_v,
                "child": child_u,
                "child_pos": child_pos,
                "d": d,
                "positive_k": pos_after,
                "delta_Iprime": dIp,
            }


def _scan_tree(stats: Stats, n: int, adj: list[list[int]], tree_index: int, graph6: str | None) -> None:
    stats.trees_scanned += 1

    for root in range(n):
        children, f, g, subtree_size = rooted_dp_all(adj, root)

        for v in range(n):
            P = [1]
            Q = [0, 1]

            for child_pos, c in enumerate(children[v]):
                U = f[c]
                V = g[c]

                _scan_step(
                    stats,
                    n=n,
                    tree_index=tree_index,
                    graph6=graph6,
                    root=root,
                    parent_v=v,
                    child_u=c,
                    child_pos=child_pos,
                    child_size=subtree_size[c],
                    child_arity=len(children[c]),
                    P=P,
                    Q=Q,
                    U=U,
                    V=V,
                )

                P = poly_mul(P, poly_add(U, V))
                Q = poly_mul(Q, U)

            # Consistency check for the rooted DP state.
            if P != f[v] or Q != g[v]:
                raise RuntimeError("Lamport composition replay mismatch")


def _tree_iter_exhaustive(n: int, backend: str):
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


def scan_exhaustive(min_n: int, max_n: int, backend: str) -> dict[str, Any]:
    by_n: list[dict[str, Any]] = []
    total = Stats()
    t_start = time.time()

    for n in range(min_n, max_n + 1):
        n_stats = Stats()
        n_start = time.time()

        for tree_index, (n_out, adj, g6) in enumerate(_tree_iter_exhaustive(n, backend), start=1):
            _scan_tree(n_stats, n_out, adj, tree_index, g6)

        elapsed = time.time() - n_start

        by_n.append(
            {
                "n": n,
                "elapsed_s": round(elapsed, 3),
                **n_stats.to_dict(),
            }
        )

        total.trees_scanned += n_stats.trees_scanned
        total.lamport_steps += n_stats.lamport_steps
        total.steps_no_descent += n_stats.steps_no_descent
        total.steps_with_descent += n_stats.steps_with_descent
        total.dd_all_hold += n_stats.dd_all_hold
        total.dd_failures += n_stats.dd_failures
        total.dd_failures_boundary_only += n_stats.dd_failures_boundary_only
        total.dd_failures_interior += n_stats.dd_failures_interior
        total.positive_tail_after_boundary += n_stats.positive_tail_after_boundary
        total.pv_positive_after_d += n_stats.pv_positive_after_d
        total.pv_positive_after_d_plus_1 += n_stats.pv_positive_after_d_plus_1
        total.dd_boundary_fail += n_stats.dd_boundary_fail
        total.dd_boundary_hold += n_stats.dd_boundary_hold
        total.dd_boundary_fail_and_pv_positive_after_d += n_stats.dd_boundary_fail_and_pv_positive_after_d
        total.dd_boundary_hold_and_pv_positive_after_d += n_stats.dd_boundary_hold_and_pv_positive_after_d
        total.dplus1_positive_num += n_stats.dplus1_positive_num
        total.dplus1_zero_den_positive_num += n_stats.dplus1_zero_den_positive_num
        if n_stats.dplus1_ratio_max > total.dplus1_ratio_max:
            total.dplus1_ratio_max = n_stats.dplus1_ratio_max
            total.dplus1_ratio_max_witness = n_stats.dplus1_ratio_max_witness
        if total.pv_positive_offset_hist is None:
            total.pv_positive_offset_hist = {}
        if n_stats.pv_positive_offset_hist:
            for k, v in n_stats.pv_positive_offset_hist.items():
                total.pv_positive_offset_hist[k] = total.pv_positive_offset_hist.get(k, 0) + v
        if total.profile_by_child_size is None:
            total.profile_by_child_size = {}
        if total.profile_by_child_arity is None:
            total.profile_by_child_arity = {}
        if total.profile_by_deg_u is None:
            total.profile_by_deg_u = {}
        for attr in ("profile_by_child_size", "profile_by_child_arity", "profile_by_deg_u"):
            src = getattr(n_stats, attr)
            dst = getattr(total, attr)
            if not src:
                continue
            for k, rec in src.items():
                if k not in dst:
                    dst[k] = {"steps": 0, "boundary_fail": 0, "pv_positive_after_d": 0}
                dst[k]["steps"] += rec.get("steps", 0)
                dst[k]["boundary_fail"] += rec.get("boundary_fail", 0)
                dst[k]["pv_positive_after_d"] += rec.get("pv_positive_after_d", 0)
        if total.first_boundary_witness is None and n_stats.first_boundary_witness is not None:
            total.first_boundary_witness = n_stats.first_boundary_witness
        if total.first_interior_witness is None and n_stats.first_interior_witness is not None:
            total.first_interior_witness = n_stats.first_interior_witness
        if total.first_positive_tail_witness is None and n_stats.first_positive_tail_witness is not None:
            total.first_positive_tail_witness = n_stats.first_positive_tail_witness
        if total.first_pv_positive_after_d_witness is None and n_stats.first_pv_positive_after_d_witness is not None:
            total.first_pv_positive_after_d_witness = n_stats.first_pv_positive_after_d_witness

        print(
            f"[exhaustive] n={n} trees={n_stats.trees_scanned} "
            f"steps={n_stats.lamport_steps} dd_interior={n_stats.dd_failures_interior} "
            f"pos_after_d={n_stats.positive_tail_after_boundary} elapsed={elapsed:.2f}s",
            flush=True,
        )

    return {
        "range": {"min_n": min_n, "max_n": max_n},
        "backend": backend,
        "elapsed_total_s": round(time.time() - t_start, 3),
        "by_n": by_n,
        "totals": total.to_dict(),
    }


def scan_stratified_samples(ns: list[int], mod: int, per_partition: int) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []

    for n in ns:
        n_stats = Stats()
        t_start = time.time()
        partition_counts: list[int] = []

        for res in range(mod):
            count = 0
            for tree_index, (n_out, adj, raw) in enumerate(trees_geng_raw(n, res=res, mod=mod), start=1):
                _scan_tree(n_stats, n_out, adj, tree_index, raw.decode("ascii"))
                count += 1
                if count >= per_partition:
                    break
            partition_counts.append(count)

        elapsed = time.time() - t_start
        rec = {
            "n": n,
            "sampling": {
                "type": "geng_partition_prefix",
                "mod": mod,
                "per_partition": per_partition,
                "partition_counts": partition_counts,
            },
            "elapsed_s": round(elapsed, 3),
            **n_stats.to_dict(),
        }
        out.append(rec)

        print(
            f"[sample] n={n} trees={n_stats.trees_scanned} "
            f"steps={n_stats.lamport_steps} dd_interior={n_stats.dd_failures_interior} "
            f"pos_after_d={n_stats.positive_tail_after_boundary} elapsed={elapsed:.2f}s",
            flush=True,
        )

    return out


def parse_int_list(s: str) -> list[int]:
    s = s.strip()
    if not s:
        return []
    return [int(x) for x in s.split(",") if x.strip()]


def main() -> None:
    ap = argparse.ArgumentParser(description="Shannon-style Lamport transition scan")
    ap.add_argument("--min-n", type=int, default=1)
    ap.add_argument("--max-n", type=int, default=17)
    ap.add_argument("--backend", choices=["auto", "geng", "networkx"], default="auto")
    ap.add_argument("--sample-n", type=str, default="")
    ap.add_argument("--sample-mod", type=int, default=16)
    ap.add_argument("--sample-per-partition", type=int, default=0)
    ap.add_argument("--out", type=str, required=True)
    args = ap.parse_args()

    started = time.time()

    report: dict[str, Any] = {
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "description": (
            "Lamport transition scan for rooted tree composition. "
            "Tracks where first-difference tail violations occur relative to d(IU)."
        ),
        "definitions": {
            "transition": "P'=P(U+V), Q'=QU, I'=IU+PV with I=P+Q",
            "d": "first descent index of IU (smallest k with Delta(IU)_k<0)",
            "difference_dominance": "Delta(PV)_k <= -Delta(IU)_k for k>=d",
            "positive_tail_after_boundary": "exists k>d with Delta(I')_k>0",
        },
    }

    if args.min_n <= args.max_n:
        report["exhaustive"] = scan_exhaustive(args.min_n, args.max_n, args.backend)

    sample_ns = parse_int_list(args.sample_n)
    if sample_ns:
        if args.sample_per_partition <= 0:
            raise ValueError("--sample-per-partition must be >0 when --sample-n is set")
        report["stratified_samples"] = scan_stratified_samples(
            sample_ns, mod=args.sample_mod, per_partition=args.sample_per_partition
        )

    report["elapsed_total_s"] = round(time.time() - started, 3)

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
