#!/usr/bin/env python3
"""Exact Steiner-gap scan on decimated core model.

For d_leaf<=1 trees, this checks the gap-formula lane:

  F_gap(S) = sum_{u in N_C(S) cap (C\\A)} (1/3 - P(u))
           + 1/2 * sum_{u in N_C(S) cap A} (1/3 - P(u))
           - sum_{h in S} (P(h) - 1/3),

where H_core = {h in C\\A : P(h) > 1/3}.

Checks:
  (A) F_gap(S) >= 0 for all non-empty S subseteq H_core.
  (B) For each non-singleton S, there exists a Steiner leaf h in S with
      M_gap(h,S) = F_gap(S)-F_gap(S\\{h}) >= 0.
  (C) Track strictness margins and witnesses.

All probabilities are computed exactly via tree DP on the decimated core model
using `weighted_hard_core_probs` (core marginals equal original marginals on C).
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from collections import Counter
from typing import Any

from conjecture_a_decimation_core_model import decimate_tree, weighted_hard_core_probs
from conjecture_a_hall_subset_scan import is_dleaf_le_1, parse_graph6


def decode_subset(mask: int, nodes: list[int]) -> list[int]:
    out: list[int] = []
    rem = mask
    while rem:
        bit = rem & -rem
        i = bit.bit_length() - 1
        rem ^= bit
        out.append(nodes[i])
    return out


def steiner_direction_masks(core_adj: list[list[int]], h_nodes: list[int]) -> list[list[int]]:
    """For each heavy vertex h_i, return masks of heavy vertices by first-hop direction.

    For fixed heavy h_i and each neighbor u of h_i in the tree C, the mask records
    which heavy vertices lie in the component of C-h_i containing u.
    A heavy h_i is a Steiner leaf of subset S iff among its direction masks,
    at most one intersects S\\{h_i}.
    """
    h_idx = {h: i for i, h in enumerate(h_nodes)}
    out: list[list[int]] = []
    for i, h in enumerate(h_nodes):
        dir_masks: list[int] = []
        for start in core_adj[h]:
            mask = 0
            stack: list[tuple[int, int]] = [(start, h)]
            while stack:
                v, parent = stack.pop()
                j = h_idx.get(v)
                if j is not None and j != i:
                    mask |= 1 << j
                for w in core_adj[v]:
                    if w != parent:
                        stack.append((w, v))
            dir_masks.append(mask)
        out.append(dir_masks)
    return out


def scan_tree(
    core_adj: list[list[int]],
    core_to_orig: list[int],
    supports_orig: set[int],
    probs: list[float],
    tol: float,
) -> dict[str, Any]:
    one_third = 1.0 / 3.0
    a_idx = {i for i, v in enumerate(core_to_orig) if v in supports_orig}

    # H_core = heavy non-support core vertices
    h_nodes = [i for i in range(len(core_adj)) if i not in a_idx and probs[i] > one_third + tol]
    m = len(h_nodes)
    if m == 0:
        return {
            "m": 0,
            "subset_count": 0,
            "steiner_subset_count": 0,
            "fgap_ok": True,
            "steiner_ok": True,
            "min_fgap": 0.0,
            "argmin_fgap_mask": 0,
            "argmin_fgap_size": 0,
            "min_steiner_best_m": 0.0,
            "argmin_steiner_mask": 0,
            "fgap_zero_count": 0,
            "steiner_zero_count": 0,
            "h_nodes": [],
        }

    h_set = set(h_nodes)
    u_nodes = sorted({u for h in h_nodes for u in core_adj[h] if u not in h_set})
    u_idx = {u: i for i, u in enumerate(u_nodes)}

    # Gap supplies and heavy demands.
    dem = [probs[h] - one_third for h in h_nodes]
    sup = [((0.5 if u in a_idx else 1.0) * (one_third - probs[u])) for u in u_nodes]

    h_masks: list[int] = []
    for h in h_nodes:
        mask = 0
        for u in core_adj[h]:
            j = u_idx.get(u)
            if j is not None:
                mask |= 1 << j
        h_masks.append(mask)

    # Direction masks for Steiner-leaf test.
    dir_masks = steiner_direction_masks(core_adj, h_nodes)

    u_pow = 1 << len(u_nodes)
    sup_sum = [0.0] * u_pow
    for s in range(1, u_pow):
        bit = s & -s
        b = bit.bit_length() - 1
        sup_sum[s] = sup_sum[s ^ bit] + sup[b]

    h_pow = 1 << m
    dem_sum = [0.0] * h_pow
    nbr_mask = [0] * h_pow

    subset_count = 0
    steiner_subset_count = 0
    fgap_zero_count = 0
    steiner_zero_count = 0
    fgap_ok = True
    steiner_ok = True

    min_fgap = float("inf")
    argmin_fgap_mask = 0
    argmin_fgap_size = 0
    min_steiner_best_m = float("inf")
    argmin_steiner_mask = 0

    for s in range(1, h_pow):
        bit = s & -s
        b = bit.bit_length() - 1
        prev = s ^ bit
        dem_sum[s] = dem_sum[prev] + dem[b]
        nbr_mask[s] = nbr_mask[prev] | h_masks[b]

        subset_count += 1
        fg = sup_sum[nbr_mask[s]] - dem_sum[s]
        if fg < min_fgap:
            min_fgap = fg
            argmin_fgap_mask = s
            argmin_fgap_size = s.bit_count()
        if fg < -tol:
            fgap_ok = False
        if abs(fg) <= tol:
            fgap_zero_count += 1

        if s.bit_count() <= 1:
            continue

        steiner_subset_count += 1
        rem = s
        best_m = float("-inf")
        while rem:
            hbit = rem & -rem
            i = hbit.bit_length() - 1
            rem ^= hbit
            others = s ^ hbit

            # h_i is Steiner leaf iff all other S-vertices lie in at most one
            # first-hop component at h_i.
            nonempty_dirs = 0
            for comp in dir_masks[i]:
                if comp & others:
                    nonempty_dirs += 1
                    if nonempty_dirs > 1:
                        break
            if nonempty_dirs > 1:
                continue

            priv = h_masks[i] & ~nbr_mask[others]
            mval = sup_sum[priv] - dem[i]
            if mval > best_m:
                best_m = mval

        if best_m == float("-inf"):
            # Should not occur for non-singleton subsets; treat as failure.
            steiner_ok = False
            best_m = float("-inf")
        else:
            if best_m < -tol:
                steiner_ok = False
            if abs(best_m) <= tol:
                steiner_zero_count += 1
            if best_m < min_steiner_best_m:
                min_steiner_best_m = best_m
                argmin_steiner_mask = s

    if min_steiner_best_m == float("inf"):
        min_steiner_best_m = 0.0

    return {
        "m": m,
        "subset_count": subset_count,
        "steiner_subset_count": steiner_subset_count,
        "fgap_ok": fgap_ok,
        "steiner_ok": steiner_ok,
        "min_fgap": min_fgap,
        "argmin_fgap_mask": argmin_fgap_mask,
        "argmin_fgap_size": argmin_fgap_size,
        "min_steiner_best_m": min_steiner_best_m,
        "argmin_steiner_mask": argmin_steiner_mask,
        "fgap_zero_count": fgap_zero_count,
        "steiner_zero_count": steiner_zero_count,
        "h_nodes": h_nodes,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Exact full Steiner-gap scan on d_leaf<=1 trees."
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
    total_with_h = 0
    skipped_heavy = 0

    subset_count = 0
    steiner_subset_count = 0
    fgap_fail = 0
    steiner_fail = 0
    fgap_zero_count = 0
    steiner_zero_count = 0

    min_fgap = float("inf")
    min_steiner_best_m = float("inf")
    worst_fgap: dict[str, Any] | None = None
    worst_steiner: dict[str, Any] | None = None

    argmin_fgap_size_dist: Counter[int] = Counter()
    per_n: dict[str, dict[str, Any]] = {}

    t_all = time.time()
    print("Steiner-gap exact scan", flush=True)
    print("Checks: F_gap(S)>=0 and existence of Steiner leaf with M_gap>=0", flush=True)
    print("-" * 88, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_with_h = 0
        n_skipped_heavy = 0
        n_subset_count = 0
        n_steiner_subset_count = 0
        n_fgap_fail = 0
        n_steiner_fail = 0

        assert proc.stdout is not None
        for line in proc.stdout:
            n0, adj = parse_graph6(line)
            n_seen += 1
            total_seen += 1
            if not is_dleaf_le_1(n0, adj):
                continue

            n_considered += 1
            total_considered += 1
            g6 = line.decode("ascii").strip()

            core_adj, core_to_orig, _leaves, supports_orig, lam = decimate_tree(adj)
            probs = weighted_hard_core_probs(core_adj, lam)

            a_idx = {i for i, v in enumerate(core_to_orig) if v in supports_orig}
            m = sum(1 for i in range(len(core_adj)) if i not in a_idx and probs[i] > 1.0 / 3.0 + args.tol)
            if m > args.max_heavy:
                skipped_heavy += 1
                n_skipped_heavy += 1
                continue

            res = scan_tree(core_adj, core_to_orig, set(supports_orig), probs, args.tol)
            if res["m"] > 0:
                total_with_h += 1
                n_with_h += 1
            subset_count += res["subset_count"]
            steiner_subset_count += res["steiner_subset_count"]
            fgap_zero_count += res["fgap_zero_count"]
            steiner_zero_count += res["steiner_zero_count"]
            n_subset_count += res["subset_count"]
            n_steiner_subset_count += res["steiner_subset_count"]

            argmin_fgap_size_dist[res["argmin_fgap_size"]] += 1

            if not res["fgap_ok"]:
                fgap_fail += 1
                n_fgap_fail += 1
            if not res["steiner_ok"]:
                steiner_fail += 1
                n_steiner_fail += 1

            if res["subset_count"] > 0 and res["min_fgap"] < min_fgap:
                min_fgap = res["min_fgap"]
                worst_fgap = {
                    "n": n0,
                    "g6": g6,
                    "H_core_indices": res["h_nodes"],
                    "argmin_size": res["argmin_fgap_size"],
                    "argmin_subset": decode_subset(res["argmin_fgap_mask"], res["h_nodes"]),
                    "min_fgap": res["min_fgap"],
                }

            if res["steiner_subset_count"] > 0 and res["min_steiner_best_m"] < min_steiner_best_m:
                min_steiner_best_m = res["min_steiner_best_m"]
                worst_steiner = {
                    "n": n0,
                    "g6": g6,
                    "H_core_indices": res["h_nodes"],
                    "argmin_subset": decode_subset(res["argmin_steiner_mask"], res["h_nodes"]),
                    "min_best_steiner_marginal": res["min_steiner_best_m"],
                }

        proc.wait()

        per_n[str(n)] = {
            "seen": n_seen,
            "considered": n_considered,
            "with_H": n_with_h,
            "skipped_heavy": n_skipped_heavy,
            "subset_count": n_subset_count,
            "steiner_subset_count": n_steiner_subset_count,
            "fgap_fail": n_fgap_fail,
            "steiner_fail": n_steiner_fail,
            "wall_s": time.time() - t0,
        }
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} with_H={n_with_h:7d} "
            f"subsets={n_subset_count:10d} steiner_subsets={n_steiner_subset_count:10d} "
            f"fgap_fail={n_fgap_fail:4d} steiner_fail={n_steiner_fail:4d} "
            f"skip={n_skipped_heavy:4d} ({time.time()-t0:.1f}s)",
            flush=True,
        )

    print("-" * 88, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} with_H={total_with_h:,} "
        f"subsets={subset_count:,} steiner_subsets={steiner_subset_count:,} "
        f"fgap_fail={fgap_fail:,} steiner_fail={steiner_fail:,} "
        f"fgap_zero={fgap_zero_count:,} steiner_zero={steiner_zero_count:,} "
        f"skipped_heavy={skipped_heavy:,} wall={time.time()-t_all:.1f}s",
        flush=True,
    )
    print(f"Minimum F_gap over non-empty subsets: {min_fgap}", flush=True)
    print(f"Minimum best-Steiner marginal: {min_steiner_best_m}", flush=True)
    print(f"Worst F_gap witness: {worst_fgap}", flush=True)
    print(f"Worst Steiner witness: {worst_steiner}", flush=True)

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
            "with_H": total_with_h,
            "subsets": subset_count,
            "steiner_subsets": steiner_subset_count,
            "fgap_fail": fgap_fail,
            "steiner_fail": steiner_fail,
            "fgap_zero": fgap_zero_count,
            "steiner_zero": steiner_zero_count,
            "skipped_heavy": skipped_heavy,
            "minimum_fgap": min_fgap,
            "minimum_best_steiner_marginal": min_steiner_best_m,
            "worst_fgap": worst_fgap,
            "worst_steiner": worst_steiner,
            "argmin_fgap_size_dist": {str(k): v for k, v in sorted(argmin_fgap_size_dist.items())},
            "wall_s": time.time() - t_all,
        },
        "per_n": per_n,
    }

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
