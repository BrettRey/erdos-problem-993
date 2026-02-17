#!/usr/bin/env python3
"""Scan marginal-removal structure for WHNC submodularity attacks.

For heavy set H = {h : P(h) > 1/3}, define:
  demand(h) = P(h) - 1/3
  supply(u) = 1/3 - P(u), u in N(H)
  F(S)      = supply(N(S)) - demand(S), for non-empty S subseteq H

Removal marginal:
  M(h, S) = F(S) - F(S \\ {h})
          = supply(N(h) \\ N(S \\ {h})) - demand(h).

This script checks whether there exist subsets S with all removal marginals negative
(M(h,S) < 0 for all h in S), and measures their gap above singleton slack:
  gap(S) = F(S) - min_{h in S} F({h}).
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time

from conjecture_a_hall_subset_scan import hard_core_probs, is_dleaf_le_1, parse_graph6


def decode_subset(mask: int, h_nodes: list[int]) -> list[int]:
    out: list[int] = []
    rem = mask
    while rem:
        lsb = rem & -rem
        i = lsb.bit_length() - 1
        rem ^= lsb
        out.append(h_nodes[i])
    return out


def scan_tree(
    adj: list[list[int]],
    probs: list[float],
    tol: float,
) -> dict:
    one_third = 1.0 / 3.0
    h_nodes = [v for v, pv in enumerate(probs) if pv > one_third + tol]
    m = len(h_nodes)
    if m <= 1:
        return {
            "allneg_count": 0,
            "allneg_with_leaf_count": 0,
            "allneg_without_leaf_count": 0,
            "best_allneg_gap": None,
            "best_allneg_subset_mask": 0,
            "best_allneg_subset_has_leaf": None,
        }

    h_set = set(h_nodes)
    u_nodes = sorted({u for h in h_nodes for u in adj[h] if u not in h_set})
    u_idx = {u: i for i, u in enumerate(u_nodes)}

    dem = [probs[h] - one_third for h in h_nodes]
    sup = [one_third - probs[u] for u in u_nodes]
    heavy_is_leaf = [len(adj[h]) == 1 for h in h_nodes]

    h_masks: list[int] = []
    for h in h_nodes:
        msk = 0
        for u in adj[h]:
            j = u_idx.get(u)
            if j is not None:
                msk |= 1 << j
        h_masks.append(msk)

    u_pow = 1 << len(u_nodes)
    sup_sum = [0.0] * u_pow
    for s in range(1, u_pow):
        lsb = s & -s
        b = lsb.bit_length() - 1
        sup_sum[s] = sup_sum[s ^ lsb] + sup[b]

    h_pow = 1 << m
    union_mask = [0] * h_pow
    dem_sum = [0.0] * h_pow
    min_single = [float("inf")] * h_pow
    singleton = [sup_sum[h_masks[i]] - dem[i] for i in range(m)]

    for s in range(1, h_pow):
        lsb = s & -s
        b = lsb.bit_length() - 1
        prev = s ^ lsb
        union_mask[s] = union_mask[prev] | h_masks[b]
        dem_sum[s] = dem_sum[prev] + dem[b]
        prev_min = min_single[prev]
        si = singleton[b]
        min_single[s] = si if si < prev_min else prev_min

    allneg_count = 0
    allneg_with_leaf_count = 0
    allneg_without_leaf_count = 0
    best_allneg_gap = None
    best_allneg_subset_mask = 0
    best_allneg_subset_has_leaf = None

    for s in range(1, h_pow):
        if s & (s - 1) == 0:
            continue

        all_negative = True
        has_leaf = False

        rem = s
        while rem:
            lsb = rem & -rem
            i = lsb.bit_length() - 1
            rem ^= lsb

            if heavy_is_leaf[i]:
                has_leaf = True

            priv = h_masks[i] & ~union_mask[s ^ lsb]
            marg = sup_sum[priv] - dem[i]
            if marg >= -tol:
                all_negative = False
                break

        if not all_negative:
            continue

        allneg_count += 1
        if has_leaf:
            allneg_with_leaf_count += 1
        else:
            allneg_without_leaf_count += 1

        f_s = sup_sum[union_mask[s]] - dem_sum[s]
        gap = f_s - min_single[s]

        if (best_allneg_gap is None) or (gap < best_allneg_gap):
            best_allneg_gap = gap
            best_allneg_subset_mask = s
            best_allneg_subset_has_leaf = has_leaf

    return {
        "allneg_count": allneg_count,
        "allneg_with_leaf_count": allneg_with_leaf_count,
        "allneg_without_leaf_count": allneg_without_leaf_count,
        "best_allneg_gap": best_allneg_gap,
        "best_allneg_subset_mask": best_allneg_subset_mask,
        "best_allneg_subset_has_leaf": best_allneg_subset_has_leaf,
        "h_nodes": h_nodes,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Marginal-removal scan for WHNC subset structure."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=23)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    total_seen = 0
    total_considered = 0

    trees_with_allneg = 0
    trees_with_allneg_and_leaf = 0
    trees_with_allneg_without_leaf = 0
    total_allneg_subsets = 0
    total_allneg_with_leaf = 0
    total_allneg_without_leaf = 0

    global_best = None
    per_n: dict[str, dict] = {}

    t_all = time.time()
    print("Marginal-removal scan")
    print("Checking subsets with all M(h,S) < 0")
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_trees_with_allneg = 0
        n_trees_with_allneg_leaf = 0
        n_trees_with_allneg_nonleaf = 0
        n_allneg_subsets = 0
        n_allneg_with_leaf = 0
        n_allneg_without_leaf = 0
        n_best_gap = None

        assert proc.stdout is not None
        for line in proc.stdout:
            n0, adj = parse_graph6(line)
            n_seen += 1
            total_seen += 1

            if not is_dleaf_le_1(n0, adj):
                continue

            n_considered += 1
            total_considered += 1
            probs = hard_core_probs(n0, adj)
            res = scan_tree(adj, probs, args.tol)

            allneg_count = res["allneg_count"]
            if allneg_count == 0:
                continue

            n_trees_with_allneg += 1
            trees_with_allneg += 1
            n_allneg_subsets += allneg_count
            total_allneg_subsets += allneg_count

            allneg_with_leaf_count = res["allneg_with_leaf_count"]
            allneg_without_leaf_count = res["allneg_without_leaf_count"]
            n_allneg_with_leaf += allneg_with_leaf_count
            n_allneg_without_leaf += allneg_without_leaf_count
            total_allneg_with_leaf += allneg_with_leaf_count
            total_allneg_without_leaf += allneg_without_leaf_count

            if allneg_with_leaf_count > 0:
                n_trees_with_allneg_leaf += 1
                trees_with_allneg_and_leaf += 1
            if allneg_without_leaf_count > 0:
                n_trees_with_allneg_nonleaf += 1
                trees_with_allneg_without_leaf += 1

            best_gap = res["best_allneg_gap"]
            if (n_best_gap is None) or (
                best_gap is not None and best_gap < n_best_gap
            ):
                n_best_gap = best_gap

            if best_gap is not None and (
                global_best is None or best_gap < global_best["gap"]
            ):
                h_nodes = res["h_nodes"]
                s_mask = res["best_allneg_subset_mask"]
                global_best = {
                    "n": n0,
                    "g6": line.decode("ascii").strip(),
                    "H": h_nodes,
                    "S": decode_subset(s_mask, h_nodes),
                    "gap": best_gap,
                    "subset_has_leaf": res["best_allneg_subset_has_leaf"],
                }

        proc.wait()
        per_n[str(n)] = {
            "seen": n_seen,
            "considered": n_considered,
            "trees_with_allneg": n_trees_with_allneg,
            "trees_with_allneg_leaf": n_trees_with_allneg_leaf,
            "trees_with_allneg_nonleaf": n_trees_with_allneg_nonleaf,
            "allneg_subsets": n_allneg_subsets,
            "allneg_subsets_with_leaf": n_allneg_with_leaf,
            "allneg_subsets_without_leaf": n_allneg_without_leaf,
            "best_gap": n_best_gap,
        }
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
            f"trees_with_allneg={n_trees_with_allneg:7d} "
            f"allneg_subsets={n_allneg_subsets:10d} "
            f"best_gap={('N/A' if n_best_gap is None else f'{n_best_gap:.12f}'):>16} "
            f"({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} "
        f"trees_with_allneg={trees_with_allneg:,} "
        f"allneg_subsets={total_allneg_subsets:,} "
        f"allneg_with_leaf={total_allneg_with_leaf:,} "
        f"allneg_without_leaf={total_allneg_without_leaf:,} "
        f"wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    print(
        f"Tree-level: with_leaf={trees_with_allneg_and_leaf:,} "
        f"without_leaf={trees_with_allneg_without_leaf:,}",
        flush=True,
    )
    print(f"Global best all-negative subset gap witness: {global_best}", flush=True)

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        payload = {
            "params": {
                "min_n": args.min_n,
                "max_n": args.max_n,
                "tol": args.tol,
            },
            "summary": {
                "seen": total_seen,
                "considered": total_considered,
                "trees_with_allneg": trees_with_allneg,
                "trees_with_allneg_leaf": trees_with_allneg_and_leaf,
                "trees_with_allneg_nonleaf": trees_with_allneg_without_leaf,
                "allneg_subsets": total_allneg_subsets,
                "allneg_subsets_with_leaf": total_allneg_with_leaf,
                "allneg_subsets_without_leaf": total_allneg_without_leaf,
                "global_best": global_best,
                "wall_s": time.time() - t_all,
            },
            "per_n": per_n,
        }
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
