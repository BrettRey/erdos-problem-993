#!/usr/bin/env python3
"""Scan all-negative leaf-heavy subsets by structural signature.

Subset setup:
  H = {h : P(h) > 1/3}
  demand(h) = P(h) - 1/3
  supply(u) = 1/3 - P(u), u in N(H)
  F(S) = supply(N(S)) - demand(S)

For non-empty S subseteq H, define removal marginal:
  M(h,S) = F(S) - F(S\\{h}) = supply(N_priv(h,S)) - demand(h).

This script scans subsets with:
  - |S| >= 2
  - at least one heavy leaf in S
  - M(h,S) < 0 for all h in S  (all-negative marginals)

For each such S, it records a signature:
  (k, l, r, deg_profile)
where
  k = |S|,
  l = # heavy leaves in S,
  r = k - l,
  deg_profile = sorted tuple of deg_S(u) over u in N(S).

It also tracks:
  gap(S) = F(S) - min_{h in S} F({h}),
global minimum gap, counts per signature, and minimum-gap witness per signature.
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
        description="Signature census of all-negative leaf-heavy subsets."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=21)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument(
        "--top-k",
        type=int,
        default=40,
        help="How many signatures to keep in top-by-count and top-by-gap lists.",
    )
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    one_third = 1.0 / 3.0
    total_seen = 0
    total_considered = 0
    total_subsets = 0
    per_n: dict[str, dict] = {}

    sig_count: Counter[tuple[int, int, int, tuple[int, ...]]] = Counter()
    sig_best: dict[tuple[int, int, int, tuple[int, ...]], dict] = {}
    global_best = None

    t_all = time.time()
    print("Leaf signature scan", flush=True)
    print("Collecting all-negative leaf-heavy subsets and their signatures", flush=True)
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_subsets = 0

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

                k = len(subset_bits)
                l = sum(1 for i in subset_bits if heavy_is_leaf[i])
                r = k - l

                u_mask = union_mask[s]
                deg_profile: list[int] = []
                uu = u_mask
                while uu:
                    ul = uu & -uu
                    j = ul.bit_length() - 1
                    uu ^= ul
                    bit = 1 << j
                    deg = 0
                    for i in subset_bits:
                        if h_masks[i] & bit:
                            deg += 1
                    deg_profile.append(deg)
                deg_profile.sort()

                sig = (k, l, r, tuple(deg_profile))
                sig_count[sig] += 1

                f_s = sup_sum[u_mask] - dem_sum[s]
                gap = f_s - min_single[s]

                subset_vertices = [h_nodes[i] for i in subset_bits]
                witness = {
                    "n": n0,
                    "g6": g6,
                    "H": h_nodes,
                    "S": subset_vertices,
                    "gap": gap,
                    "signature": {
                        "k": k,
                        "l": l,
                        "r": r,
                        "deg_profile": deg_profile,
                    },
                }

                prev = sig_best.get(sig)
                if (prev is None) or (gap < prev["gap"] - 1e-15):
                    sig_best[sig] = witness

                if (global_best is None) or (gap < global_best["gap"] - 1e-15):
                    global_best = witness

        proc.wait()
        per_n[str(n)] = {
            "seen": n_seen,
            "considered": n_considered,
            "leaf_allneg_subsets": n_subsets,
            "wall_s": time.time() - t0,
        }
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
            f"leaf_allneg_subsets={n_subsets:9d} ({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} "
        f"leaf_allneg_subsets={total_subsets:,} wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    print(f"Global best: {global_best}", flush=True)

    top_count = sig_count.most_common(args.top_k)
    top_gap = sorted(sig_best.items(), key=lambda kv: kv[1]["gap"])[: args.top_k]

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

        sig_count_out = []
        for sig, c in top_count:
            sig_count_out.append(
                {
                    "signature": {
                        "k": sig[0],
                        "l": sig[1],
                        "r": sig[2],
                        "deg_profile": list(sig[3]),
                    },
                    "count": c,
                    "best_gap": sig_best[sig]["gap"],
                    "best_witness": sig_best[sig],
                }
            )

        sig_gap_out = []
        for sig, witness in top_gap:
            sig_gap_out.append(
                {
                    "signature": {
                        "k": sig[0],
                        "l": sig[1],
                        "r": sig[2],
                        "deg_profile": list(sig[3]),
                    },
                    "count": sig_count[sig],
                    "best_gap": witness["gap"],
                    "best_witness": witness,
                }
            )

        payload = {
            "params": {
                "min_n": args.min_n,
                "max_n": args.max_n,
                "tol": args.tol,
                "top_k": args.top_k,
            },
            "summary": {
                "seen": total_seen,
                "considered": total_considered,
                "leaf_allneg_subsets": total_subsets,
                "signature_count": len(sig_count),
                "global_best": global_best,
                "wall_s": time.time() - t_all,
            },
            "top_by_count": sig_count_out,
            "top_by_best_gap": sig_gap_out,
            "per_n": per_n,
        }
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
