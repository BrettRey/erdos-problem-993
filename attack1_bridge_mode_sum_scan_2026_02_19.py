#!/usr/bin/env python3
"""Attack 1: bridge decomposition scan for child-mode sums.

For d_leaf <= 1 trees with canonical degree-2 bridge decomposition:
  T -> B rooted at u, children c_1,...,c_d, P = prod_i I(T_{c_i}), m = mode(I(T)).

This script checks:
  1) mode(P) >= m - 1
  2) sum_i mode(I(T_{c_i})) >= m - 1
  3) mode(P) >= sum_i mode(I(T_{c_i}))   (superadditivity over many factors)
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from typing import Any

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from graph6 import parse_graph6
from indpoly import _polyadd, _polymul, independence_poly


def mode_left(poly: list[int]) -> int:
    return max(range(len(poly)), key=lambda i: poly[i])


def choose_min_support_leaf(adj: list[list[int]]) -> tuple[int, int]:
    deg = [len(nb) for nb in adj]
    leaves = [v for v, d in enumerate(deg) if d == 1]
    parent = {l: adj[l][0] for l in leaves}
    # Minimum support degree, tie by largest leaf id.
    leaf = max(leaves, key=lambda l: (-deg[parent[l]], -l))
    return leaf, parent[leaf]


def remove_vertices(adj: list[list[int]], remove_set: set[int]) -> tuple[list[list[int]], dict[int, int]]:
    keep = [v for v in range(len(adj)) if v not in remove_set]
    idx = {v: i for i, v in enumerate(keep)}
    out = [[] for _ in keep]
    for v in keep:
        vv = idx[v]
        for u in adj[v]:
            if u in idx:
                out[vv].append(idx[u])
    return out, idx


def rooted_dp(
    adj: list[list[int]],
    root: int,
) -> tuple[list[list[int]], list[list[int]], list[list[int]], list[int]]:
    """Return (dp0, dp1, children, postorder)."""
    n = len(adj)
    parent = [-1] * n
    children = [[] for _ in range(n)]
    parent[root] = root
    queue = [root]
    for v in queue:
        for w in adj[v]:
            if parent[w] == -1:
                parent[w] = v
                children[v].append(w)
                queue.append(w)

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

    dp0: list[list[int]] = [[] for _ in range(n)]
    dp1: list[list[int]] = [[] for _ in range(n)]
    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
            continue

        p0 = [1]
        for c in children[v]:
            p0 = _polymul(p0, _polyadd(dp0[c], dp1[c]))
        dp0[v] = p0

        p1 = [1]
        for c in children[v]:
            p1 = _polymul(p1, dp0[c])
        dp1[v] = [0] + p1

    return dp0, dp1, children, order


def fresh_stats() -> dict[str, Any]:
    return {
        "seen": 0,
        "considered": 0,
        "checked": 0,
        "no_deg2_support": 0,
        "identity_fail": 0,
        "modeP_ge_m_minus_1_fail": 0,
        "sum_modes_ge_m_minus_1_fail": 0,
        "modeP_ge_sum_modes_fail": 0,
        "min_modeP_minus_m_minus_1": None,
        "min_modeP_minus_m_minus_1_witness": None,
        "min_sum_modes_minus_m_minus_1": None,
        "min_sum_modes_minus_m_minus_1_witness": None,
        "min_modeP_minus_sum_modes": None,
        "min_modeP_minus_sum_modes_witness": None,
        "by_k_children": {},
        "wall_s": 0.0,
    }


def maybe_update_min(
    stats: dict[str, Any],
    key: str,
    witness_key: str,
    value: int,
    witness: dict[str, Any],
) -> None:
    cur = stats[key]
    if cur is None or value < cur:
        stats[key] = value
        stats[witness_key] = witness


def merge_stats(dst: dict[str, Any], src: dict[str, Any]) -> None:
    for k in [
        "seen",
        "considered",
        "checked",
        "no_deg2_support",
        "identity_fail",
        "modeP_ge_m_minus_1_fail",
        "sum_modes_ge_m_minus_1_fail",
        "modeP_ge_sum_modes_fail",
    ]:
        dst[k] += int(src.get(k, 0))

    for key, wkey in [
        ("min_modeP_minus_m_minus_1", "min_modeP_minus_m_minus_1_witness"),
        ("min_sum_modes_minus_m_minus_1", "min_sum_modes_minus_m_minus_1_witness"),
        ("min_modeP_minus_sum_modes", "min_modeP_minus_sum_modes_witness"),
    ]:
        if src.get(key) is not None:
            maybe_update_min(dst, key, wkey, int(src[key]), src[wkey])

    for k, v in src.get("by_k_children", {}).items():
        if k not in dst["by_k_children"]:
            dst["by_k_children"][k] = {
                "checked": 0,
                "modeP_ge_m_minus_1_fail": 0,
                "sum_modes_ge_m_minus_1_fail": 0,
                "modeP_ge_sum_modes_fail": 0,
            }
        for key in [
            "checked",
            "modeP_ge_m_minus_1_fail",
            "sum_modes_ge_m_minus_1_fail",
            "modeP_ge_sum_modes_fail",
        ]:
            dst["by_k_children"][k][key] += int(v.get(key, 0))


def aggregate(per_n: dict[str, dict[str, Any]]) -> dict[str, Any]:
    out = fresh_stats()
    for n_key in sorted(per_n, key=lambda s: int(s)):
        merge_stats(out, per_n[n_key])
    out.pop("wall_s", None)
    return out


def write_payload(out_path: str, params: dict[str, Any], per_n: dict[str, dict[str, Any]]) -> None:
    payload = {
        "params": params,
        "summary": aggregate(per_n),
        "per_n": {k: per_n[k] for k in sorted(per_n, key=lambda s: int(s))},
    }
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def main() -> None:
    ap = argparse.ArgumentParser(description="Attack 1 bridge scan for child-mode sums.")
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--res", type=int, default=None, help="Optional split residue r in r/m.")
    ap.add_argument("--mod", type=int, default=None, help="Optional split modulus m in r/m.")
    ap.add_argument(
        "--out",
        default="results/attack1_bridge_mode_sum_scan_n23_2026_02_19.json",
    )
    ap.add_argument("--no-resume", action="store_true")
    args = ap.parse_args()

    if (args.res is None) != (args.mod is None):
        raise ValueError("Specify both --res and --mod together, or neither.")
    if args.mod is not None and not (0 <= args.res < args.mod):
        raise ValueError("Require 0 <= res < mod.")

    params = {
        "min_n": args.min_n,
        "max_n": args.max_n,
        "res": args.res,
        "mod": args.mod,
    }

    per_n: dict[str, dict[str, Any]] = {}
    if (not args.no_resume) and args.out and os.path.exists(args.out):
        with open(args.out, "r", encoding="utf-8") as f:
            old = json.load(f)
        per_n = old.get("per_n", {})
        done = ", ".join(sorted(per_n, key=lambda s: int(s)))
        print(f"Resuming from {args.out}; completed n: {done}", flush=True)

    split_desc = f", split={args.res}/{args.mod}" if args.mod is not None else ""
    print(
        f"attack1_bridge_mode_sum_scan on d_leaf<=1, n={args.min_n}..{args.max_n}{split_desc}",
        flush=True,
    )
    print("Canonical leaf: minimum support degree (tie: largest leaf id)", flush=True)
    print(f"Output: {args.out}", flush=True)
    print("-" * 96, flush=True)

    t_all = time.time()
    for n in range(args.min_n, args.max_n + 1):
        n_key = str(n)
        if n_key in per_n:
            print(f"n={n:2d}: already complete, skipping", flush=True)
            continue

        t0 = time.time()
        stats = fresh_stats()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        if args.mod is not None:
            cmd.append(f"{args.res}/{args.mod}")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        assert proc.stdout is not None

        for raw in proc.stdout:
            stats["seen"] += 1
            nn, adj = parse_graph6(raw)
            if not is_dleaf_le_1(nn, adj):
                continue
            stats["considered"] += 1

            poly_t = independence_poly(nn, adj)
            m = mode_left(poly_t)
            if m == 0:
                continue

            leaf, support = choose_min_support_leaf(adj)
            if len(adj[support]) != 2:
                stats["no_deg2_support"] += 1
                continue

            u = adj[support][0] if adj[support][1] == leaf else adj[support][1]
            b_adj, idx = remove_vertices(adj, {leaf, support})
            if not b_adj:
                continue
            u_in_b = idx[u]

            dp0, dp1, children, _ = rooted_dp(b_adj, u_in_b)
            p_poly = dp0[u_in_b]

            child_polys: list[list[int]] = []
            for c in children[u_in_b]:
                child_polys.append(_polyadd(dp0[c], dp1[c]))

            # Identity check: P must equal product of child subtree polynomials.
            prod_children = [1]
            for f_c in child_polys:
                prod_children = _polymul(prod_children, f_c)
            if prod_children != p_poly:
                stats["identity_fail"] += 1

            mode_p = mode_left(p_poly)
            sum_modes = sum(mode_left(f_c) for f_c in child_polys)

            margin_mode_p = mode_p - (m - 1)
            margin_sum_modes = sum_modes - (m - 1)
            margin_super = mode_p - sum_modes

            stats["checked"] += 1
            k_children = len(child_polys)
            k_key = str(k_children)
            if k_key not in stats["by_k_children"]:
                stats["by_k_children"][k_key] = {
                    "checked": 0,
                    "modeP_ge_m_minus_1_fail": 0,
                    "sum_modes_ge_m_minus_1_fail": 0,
                    "modeP_ge_sum_modes_fail": 0,
                }
            stats["by_k_children"][k_key]["checked"] += 1

            g6 = raw.decode("ascii").strip()
            witness = {
                "n": nn,
                "g6": g6,
                "m": m,
                "leaf": leaf,
                "support": support,
                "u_in_B": u_in_b,
                "k_children": k_children,
                "mode_P": mode_p,
                "sum_modes_children": sum_modes,
                "margin_modeP_minus_m_minus_1": margin_mode_p,
                "margin_sum_modes_minus_m_minus_1": margin_sum_modes,
                "margin_modeP_minus_sum_modes": margin_super,
                "poly_T": poly_t,
                "poly_P": p_poly,
                "child_polys": child_polys,
            }

            maybe_update_min(
                stats,
                "min_modeP_minus_m_minus_1",
                "min_modeP_minus_m_minus_1_witness",
                margin_mode_p,
                witness,
            )
            maybe_update_min(
                stats,
                "min_sum_modes_minus_m_minus_1",
                "min_sum_modes_minus_m_minus_1_witness",
                margin_sum_modes,
                witness,
            )
            maybe_update_min(
                stats,
                "min_modeP_minus_sum_modes",
                "min_modeP_minus_sum_modes_witness",
                margin_super,
                witness,
            )

            if margin_mode_p < 0:
                stats["modeP_ge_m_minus_1_fail"] += 1
                stats["by_k_children"][k_key]["modeP_ge_m_minus_1_fail"] += 1
            if margin_sum_modes < 0:
                stats["sum_modes_ge_m_minus_1_fail"] += 1
                stats["by_k_children"][k_key]["sum_modes_ge_m_minus_1_fail"] += 1
            if margin_super < 0:
                stats["modeP_ge_sum_modes_fail"] += 1
                stats["by_k_children"][k_key]["modeP_ge_sum_modes_fail"] += 1

        proc.wait()
        stats["wall_s"] = time.time() - t0
        per_n[n_key] = stats
        write_payload(args.out, params, per_n)

        print(
            f"n={n:2d}: seen={stats['seen']:9d} considered={stats['considered']:8d} "
            f"checked={stats['checked']:8d} "
            f"modeP_fail={stats['modeP_ge_m_minus_1_fail']:6d} "
            f"sum_fail={stats['sum_modes_ge_m_minus_1_fail']:6d} "
            f"super_fail={stats['modeP_ge_sum_modes_fail']:6d} "
            f"min(modeP-(m-1))={stats['min_modeP_minus_m_minus_1']} "
            f"min(sum-(m-1))={stats['min_sum_modes_minus_m_minus_1']} "
            f"min(modeP-sum)={stats['min_modeP_minus_sum_modes']} "
            f"({stats['wall_s']:.1f}s)",
            flush=True,
        )

    summary = aggregate(per_n)
    wall = time.time() - t_all
    print("-" * 96, flush=True)
    print(
        f"TOTAL seen={summary['seen']:,} considered={summary['considered']:,} "
        f"checked={summary['checked']:,} identity_fail={summary['identity_fail']:,} "
        f"modeP_fail={summary['modeP_ge_m_minus_1_fail']:,} "
        f"sum_fail={summary['sum_modes_ge_m_minus_1_fail']:,} "
        f"super_fail={summary['modeP_ge_sum_modes_fail']:,} wall={wall:.1f}s",
        flush=True,
    )
    print(
        "Extrema: "
        f"min(modeP-(m-1))={summary['min_modeP_minus_m_minus_1']}, "
        f"min(sum-(m-1))={summary['min_sum_modes_minus_m_minus_1']}, "
        f"min(modeP-sum)={summary['min_modeP_minus_sum_modes']}",
        flush=True,
    )
    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
