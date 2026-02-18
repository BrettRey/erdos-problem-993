#!/usr/bin/env python3
"""Marginal-removal scan on decimated weighted core model.

For d_leaf<=1 tree T, decimate leaves to core C with activities:
  lambda_v = 1      for v in C\\A,
  lambda_v = 1/2    for v in A (supports adjacent to removed leaves).

Heavy set on core non-support vertices:
  H = {h in C\\A : P(h) > 1/3}.

Weighted supplies on neighbors of H:
  s(u) = 1/3 - P(u)                 for u in C\\A
       = (1/2)*(2/3 - P(u))         for u in A

Demands on H:
  d(h) = P(h) - 1/3.

For S subseteq H, removal marginal is
  M(h,S) = s(N_priv(h,S)) - d(h)
where N_priv(h,S)=N(h)\\N(S\\{h}).

This script checks whether all-negative subsets exist:
  M(h,S) < 0 for all h in S.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time

from conjecture_a_decimation_core_model import decimate_tree, weighted_hard_core_probs
from conjecture_a_hall_subset_scan import is_dleaf_le_1, parse_graph6


def decode_subset(mask: int, nodes: list[int]) -> list[int]:
    out: list[int] = []
    rem = mask
    while rem:
        lsb = rem & -rem
        i = lsb.bit_length() - 1
        rem ^= lsb
        out.append(nodes[i])
    return out


def scan_tree(
    core_adj: list[list[int]],
    core_to_orig: list[int],
    supports_orig: set[int],
    probs: list[float],
    tol: float,
) -> dict:
    one_third = 1.0 / 3.0
    a_idx = {i for i, v in enumerate(core_to_orig) if v in supports_orig}
    h_nodes = [i for i in range(len(core_adj)) if i not in a_idx and probs[i] > one_third + tol]
    m = len(h_nodes)
    if m == 0:
        return {
            "m": 0,
            "allneg_subsets": 0,
            "allneg_with_leaf": 0,
            "allneg_without_leaf": 0,
            "hardest_subset": None,
        }

    h_set = set(h_nodes)
    u_nodes = sorted({u for h in h_nodes for u in core_adj[h] if u not in h_set})
    u_idx = {u: i for i, u in enumerate(u_nodes)}

    dem = [probs[h] - one_third for h in h_nodes]
    h_deg_core = [len(core_adj[h]) for h in h_nodes]

    sup = []
    for u in u_nodes:
        if u in a_idx:
            sup.append(0.5 * (2.0 * one_third - probs[u]))
        else:
            sup.append(one_third - probs[u])

    h_masks: list[int] = []
    for h in h_nodes:
        msk = 0
        for u in core_adj[h]:
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
    dem_sum = [0.0] * h_pow
    nbr_mask = [0] * h_pow

    allneg_subsets = 0
    allneg_with_leaf = 0
    allneg_without_leaf = 0
    hardest_subset = None

    for s in range(1, h_pow):
        lsb = s & -s
        b = lsb.bit_length() - 1
        prev = s ^ lsb
        dem_sum[s] = dem_sum[prev] + dem[b]
        nbr_mask[s] = nbr_mask[prev] | h_masks[b]

        rem = s
        has_nonneg = False
        max_marg = float("-inf")
        min_marg = float("inf")
        while rem:
            bit = rem & -rem
            i = bit.bit_length() - 1
            rem ^= bit
            other = s ^ bit
            priv_mask = h_masks[i] & ~nbr_mask[other]
            marg = sup_sum[priv_mask] - dem[i]
            if marg >= -tol:
                has_nonneg = True
                break
            if marg > max_marg:
                max_marg = marg
            if marg < min_marg:
                min_marg = marg

        if has_nonneg:
            continue

        allneg_subsets += 1
        subset_ids = decode_subset(s, h_nodes)
        has_leaf = any(h_deg_core[i] == 1 for i in range(m) if (s >> i) & 1)
        if has_leaf:
            allneg_with_leaf += 1
        else:
            allneg_without_leaf += 1

        if (hardest_subset is None) or (max_marg > hardest_subset["max_marginal"]):
            hardest_subset = {
                "mask": s,
                "subset": subset_ids,
                "size": len(subset_ids),
                "max_marginal": max_marg,
                "min_marginal": min_marg,
                "has_leaf": has_leaf,
                "slack": sup_sum[nbr_mask[s]] - dem_sum[s],
            }

    return {
        "m": m,
        "allneg_subsets": allneg_subsets,
        "allneg_with_leaf": allneg_with_leaf,
        "allneg_without_leaf": allneg_without_leaf,
        "hardest_subset": hardest_subset,
        "h_nodes": h_nodes,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Scan all-negative removal marginals on decimated weighted model."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=23)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--max-heavy", type=int, default=24)
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    total_seen = 0
    total_considered = 0
    trees_with_h = 0
    trees_with_allneg = 0
    allneg_subsets = 0
    allneg_with_leaf = 0
    allneg_without_leaf = 0
    skipped_heavy = 0

    worst = None
    per_n: dict[str, dict] = {}

    t_all = time.time()
    print("Decimated marginal-removal scan", flush=True)
    print("Check: existence of subsets with M(h,S)<0 for all h in S", flush=True)
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_trees_h = 0
        n_trees_allneg = 0
        n_allneg_subsets = 0
        n_allneg_with_leaf = 0
        n_allneg_without_leaf = 0
        n_skipped_heavy = 0

        assert proc.stdout is not None
        for line in proc.stdout:
            n0, adj = parse_graph6(line)
            n_seen += 1
            total_seen += 1
            if not is_dleaf_le_1(n0, adj):
                continue

            n_considered += 1
            total_considered += 1

            core_adj, core_to_orig, _leaves, supports, lam = decimate_tree(adj)
            probs = weighted_hard_core_probs(core_adj, lam)
            g6 = line.decode("ascii").strip()

            a_idx = {i for i, v in enumerate(core_to_orig) if v in supports}
            m = sum(1 for i in range(len(core_adj)) if i not in a_idx and probs[i] > 1.0 / 3.0 + args.tol)
            if m > args.max_heavy:
                skipped_heavy += 1
                n_skipped_heavy += 1
                continue

            res = scan_tree(core_adj, core_to_orig, set(supports), probs, args.tol)
            if res["m"] > 0:
                trees_with_h += 1
                n_trees_h += 1
            if res["allneg_subsets"] > 0:
                trees_with_allneg += 1
                n_trees_allneg += 1

            allneg_subsets += res["allneg_subsets"]
            allneg_with_leaf += res["allneg_with_leaf"]
            allneg_without_leaf += res["allneg_without_leaf"]

            n_allneg_subsets += res["allneg_subsets"]
            n_allneg_with_leaf += res["allneg_with_leaf"]
            n_allneg_without_leaf += res["allneg_without_leaf"]

            hs = res["hardest_subset"]
            if hs is not None:
                if (worst is None) or (hs["max_marginal"] > worst["max_marginal"]):
                    worst = {
                        "n": n0,
                        "g6": g6,
                        "H_core_indices": res["h_nodes"],
                        "subset": hs["subset"],
                        "subset_size": hs["size"],
                        "max_marginal": hs["max_marginal"],
                        "min_marginal": hs["min_marginal"],
                        "has_leaf": hs["has_leaf"],
                        "subset_slack": hs["slack"],
                    }

        proc.wait()

        per_n[str(n)] = {
            "seen": n_seen,
            "considered": n_considered,
            "trees_with_h": n_trees_h,
            "trees_with_allneg": n_trees_allneg,
            "allneg_subsets": n_allneg_subsets,
            "allneg_with_leaf": n_allneg_with_leaf,
            "allneg_without_leaf": n_allneg_without_leaf,
            "skipped_heavy": n_skipped_heavy,
            "wall_s": time.time() - t0,
        }
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
            f"trees_h={n_trees_h:7d} trees_allneg={n_trees_allneg:7d} "
            f"allneg_subsets={n_allneg_subsets:10d} "
            f"skip_heavy={n_skipped_heavy:5d} ({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} "
        f"trees_h={trees_with_h:,} trees_allneg={trees_with_allneg:,} "
        f"allneg_subsets={allneg_subsets:,} "
        f"allneg_with_leaf={allneg_with_leaf:,} allneg_without_leaf={allneg_without_leaf:,} "
        f"skipped_heavy={skipped_heavy:,} wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    print(f"Hardest all-negative witness: {worst}", flush=True)

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        payload = {
            "params": {
                "min_n": args.min_n,
                "max_n": args.max_n,
                "tol": args.tol,
                "max_heavy": args.max_heavy,
            },
            "summary": {
                "seen": total_seen,
                "considered": total_considered,
                "trees_with_h": trees_with_h,
                "trees_with_allneg": trees_with_allneg,
                "allneg_subsets": allneg_subsets,
                "allneg_with_leaf": allneg_with_leaf,
                "allneg_without_leaf": allneg_without_leaf,
                "skipped_heavy": skipped_heavy,
                "hardest_allneg": worst,
                "wall_s": time.time() - t_all,
            },
            "per_n": per_n,
        }
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
