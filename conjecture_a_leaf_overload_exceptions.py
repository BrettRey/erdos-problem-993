#!/usr/bin/env python3
"""Scan overload structure inside all-negative leaf-heavy subsets.

For a subset S subseteq H (heavy vertices), define:
  demand(h) = P(h) - 1/3,
  supply(u) = 1/3 - P(u), u in N(H),
  F(S) = supply(N(S)) - demand(S).

Within N(S), define local load at u:
  load_S(u) = sum_{h in N(u) cap S} demand(h).

We call u overloaded if:
  load_S(u) > supply(u).

This script scans all S with:
  - |S| >= 2
  - S contains at least one heavy leaf
  - all removal marginals negative (all-negative subset)

and tracks overload keys (deg_S(u), leaf_S(u)), where leaf_S(u) is the number
of heavy leaves in S adjacent to u.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from collections import Counter

from conjecture_a_hall_subset_scan import hard_core_probs, is_dleaf_le_1, parse_graph6


def decode_subset(mask: int, nodes: list[int]) -> list[int]:
    out: list[int] = []
    rem = mask
    while rem:
        lsb = rem & -rem
        i = lsb.bit_length() - 1
        rem ^= lsb
        out.append(nodes[i])
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Overload-key scan for all-negative leaf-heavy subsets."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=21)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument(
        "--max-exceptions",
        type=int,
        default=40,
        help="Max number of non-leaf-overload exception witnesses to store.",
    )
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    one_third = 1.0 / 3.0
    total_seen = 0
    total_considered = 0
    total_subsets = 0

    subset_with_over = 0
    subset_with_leaf1_over = 0
    subset_with_nonleaf_over = 0
    subset_only_leaf1_over = 0

    by_key: Counter[tuple[int, int]] = Counter()
    over_key: Counter[tuple[int, int]] = Counter()
    exceptions: list[dict] = []

    best_any = None
    best_with_leaf1_over = None
    best_with_nonleaf_over = None
    per_n: dict[str, dict] = {}

    t_all = time.time()
    print("Leaf overload exception scan", flush=True)
    print("Scanning all-negative leaf-heavy subsets and overload-key patterns", flush=True)
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_subsets = 0
        n_over = 0
        n_nonleaf_over = 0

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
            h_nodes = [v for v, pv in enumerate(probs) if pv > one_third + args.tol]
            m = len(h_nodes)
            if m <= 1:
                continue

            h_set = set(h_nodes)
            u_nodes = sorted({u for h in h_nodes for u in adj[h] if u not in h_set})
            if not u_nodes:
                continue
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
                si = singleton[b]
                prev_min = min_single[prev]
                min_single[s] = si if si < prev_min else prev_min

            g6 = line.decode("ascii").strip()
            for s in range(1, h_pow):
                if s & (s - 1) == 0:
                    continue

                rem = s
                subset_bits: list[int] = []
                has_leaf = False
                all_negative = True
                while rem:
                    lsb = rem & -rem
                    i = lsb.bit_length() - 1
                    rem ^= lsb
                    subset_bits.append(i)
                    if heavy_is_leaf[i]:
                        has_leaf = True
                    priv = h_masks[i] & ~union_mask[s ^ lsb]
                    marg = sup_sum[priv] - dem[i]
                    if marg >= -args.tol:
                        all_negative = False
                        break

                if (not all_negative) or (not has_leaf):
                    continue

                n_subsets += 1
                total_subsets += 1

                u_mask = union_mask[s]
                gap = (sup_sum[u_mask] - dem_sum[s]) - min_single[s]
                subset_vertices = [h_nodes[i] for i in subset_bits]

                if (best_any is None) or (gap < best_any["gap"] - 1e-15):
                    best_any = {
                        "n": n0,
                        "g6": g6,
                        "H": h_nodes,
                        "S": subset_vertices,
                        "gap": gap,
                    }

                has_over = False
                has_leaf1_over = False
                has_nonleaf_over = False
                nonleaf_over_details: list[dict] = []

                uu = u_mask
                while uu:
                    ul = uu & -uu
                    j = ul.bit_length() - 1
                    uu ^= ul
                    bit = 1 << j

                    deg_s = 0
                    leaf_s = 0
                    load = 0.0
                    adjacent_heavy: list[int] = []
                    for i in subset_bits:
                        if h_masks[i] & bit:
                            deg_s += 1
                            load += dem[i]
                            adjacent_heavy.append(h_nodes[i])
                            if heavy_is_leaf[i]:
                                leaf_s += 1

                    key = (deg_s, leaf_s)
                    by_key[key] += 1

                    overload = load - sup[j]
                    if overload > args.tol:
                        over_key[key] += 1
                        has_over = True
                        n_over += 1
                        if leaf_s == 1:
                            has_leaf1_over = True
                        else:
                            has_nonleaf_over = True
                            nonleaf_over_details.append(
                                {
                                    "u": u_nodes[j],
                                    "deg_S": deg_s,
                                    "leaf_S": leaf_s,
                                    "overload": overload,
                                    "adjacent_heavy_in_S": adjacent_heavy,
                                }
                            )

                if has_over:
                    subset_with_over += 1
                if has_leaf1_over:
                    subset_with_leaf1_over += 1
                    if (best_with_leaf1_over is None) or (
                        gap < best_with_leaf1_over["gap"] - 1e-15
                    ):
                        best_with_leaf1_over = {
                            "n": n0,
                            "g6": g6,
                            "H": h_nodes,
                            "S": subset_vertices,
                            "gap": gap,
                        }
                if has_nonleaf_over:
                    subset_with_nonleaf_over += 1
                    n_nonleaf_over += 1
                    if (best_with_nonleaf_over is None) or (
                        gap < best_with_nonleaf_over["gap"] - 1e-15
                    ):
                        best_with_nonleaf_over = {
                            "n": n0,
                            "g6": g6,
                            "H": h_nodes,
                            "S": subset_vertices,
                            "gap": gap,
                        }
                    if len(exceptions) < args.max_exceptions:
                        exceptions.append(
                            {
                                "n": n0,
                                "g6": g6,
                                "H": h_nodes,
                                "S": subset_vertices,
                                "gap": gap,
                                "nonleaf_overloads": nonleaf_over_details,
                            }
                        )
                if has_over and has_leaf1_over and (not has_nonleaf_over):
                    subset_only_leaf1_over += 1

        proc.wait()
        per_n[str(n)] = {
            "seen": n_seen,
            "considered": n_considered,
            "leaf_allneg_subsets": n_subsets,
            "overloaded_u_count": n_over,
            "subset_with_nonleaf_over": n_nonleaf_over,
            "wall_s": time.time() - t0,
        }
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
            f"leaf_allneg_subsets={n_subsets:8d} overloaded_u={n_over:8d} "
            f"nonleaf_over_subsets={n_nonleaf_over:6d} ({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} "
        f"leaf_allneg_subsets={total_subsets:,} "
        f"subset_with_over={subset_with_over:,} "
        f"subset_with_leaf1_over={subset_with_leaf1_over:,} "
        f"subset_with_nonleaf_over={subset_with_nonleaf_over:,} "
        f"subset_only_leaf1_over={subset_only_leaf1_over:,} "
        f"wall={time.time() - t_all:.1f}s",
        flush=True,
    )

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

        def encode_counter(src: Counter[tuple[int, int]]) -> list[dict]:
            out = []
            for (deg_s, leaf_s), c in sorted(src.items()):
                out.append(
                    {
                        "deg_S": deg_s,
                        "leaf_S": leaf_s,
                        "count": c,
                    }
                )
            return out

        payload = {
            "params": {
                "min_n": args.min_n,
                "max_n": args.max_n,
                "tol": args.tol,
                "max_exceptions": args.max_exceptions,
            },
            "summary": {
                "seen": total_seen,
                "considered": total_considered,
                "leaf_allneg_subsets": total_subsets,
                "subset_with_over": subset_with_over,
                "subset_with_leaf1_over": subset_with_leaf1_over,
                "subset_with_nonleaf_over": subset_with_nonleaf_over,
                "subset_only_leaf1_over": subset_only_leaf1_over,
                "best_any": best_any,
                "best_with_leaf1_over": best_with_leaf1_over,
                "best_with_nonleaf_over": best_with_nonleaf_over,
                "exceptions_stored": len(exceptions),
                "wall_s": time.time() - t_all,
            },
            "all_u_key_counts": encode_counter(by_key),
            "overloaded_u_key_counts": encode_counter(over_key),
            "exceptions": exceptions,
            "per_n": per_n,
        }
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
