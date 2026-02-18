#!/usr/bin/env python3
"""Hall-subset scan on the decimated weighted core model.

For d_leaf<=1 tree T, decimate leaves to core C with heterogeneous activities:
  lambda_v = 1      for v in C\\A,
  lambda_v = 1/2    for v in A (supports that had leaves).

Let P(v) be marginals on C under this model and define heavy set:
  H = {h in C\\A : P(h) > 1/3}.

Demands:
  d(h) = P(h) - 1/3,   h in H.

Weighted supplies:
  s(u) = 1/3 - P(u),                 u in C\\A
       = (1/2)*(2/3 - P(u)),         u in A
for u in N(H).

For non-empty S subseteq H:
  slack(S) = s(N(S)) - d(S).

This scan checks:
  (A) Hall feasibility: min_{S!=empty} slack(S) >= 0
  (B) singleton argmin: min_{S!=empty} slack(S) = min_{h in H} slack({h})
  (C) local singleton dominance: slack(S) >= min_{h in S} slack({h}) for all non-empty S.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from collections import Counter

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
            "hall_ok": True,
            "singleton_ok": True,
            "local_ok": True,
            "min_slack": 0.0,
            "min_singleton": 0.0,
            "min_global_gap": 0.0,
            "min_local_gap": 0.0,
            "argmin_mask": 0,
            "argmin_size": 0,
            "local_bad_mask": 0,
            "h_nodes": [],
        }

    h_set = set(h_nodes)
    u_nodes = sorted({u for h in h_nodes for u in core_adj[h] if u not in h_set})
    u_idx = {u: i for i, u in enumerate(u_nodes)}

    dem = [probs[h] - one_third for h in h_nodes]
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
    min_single_in = [float("inf")] * h_pow

    singleton = [sup_sum[h_masks[i]] - dem[i] for i in range(m)]
    min_singleton = min(singleton)

    min_slack = float("inf")
    argmin_mask = 0
    argmin_size = 0
    min_global_gap = float("inf")
    min_local_gap = float("inf")
    local_bad_mask = 0
    hall_ok = True
    local_ok = True

    for s in range(1, h_pow):
        lsb = s & -s
        b = lsb.bit_length() - 1
        prev = s ^ lsb
        dem_sum[s] = dem_sum[prev] + dem[b]
        nbr_mask[s] = nbr_mask[prev] | h_masks[b]
        min_single_in[s] = min(min_single_in[prev], singleton[b])

        sl = sup_sum[nbr_mask[s]] - dem_sum[s]
        if sl < min_slack:
            min_slack = sl
            argmin_mask = s
            argmin_size = s.bit_count()

        if sl < -tol:
            hall_ok = False

        ggap = sl - min_singleton
        if ggap < min_global_gap:
            min_global_gap = ggap

        lgap = sl - min_single_in[s]
        if lgap < min_local_gap:
            min_local_gap = lgap
            if lgap < -tol:
                local_bad_mask = s
        if lgap < -tol:
            local_ok = False

    singleton_ok = min_global_gap >= -tol
    return {
        "m": m,
        "hall_ok": hall_ok,
        "singleton_ok": singleton_ok,
        "local_ok": local_ok,
        "min_slack": min_slack,
        "min_singleton": min_singleton,
        "min_global_gap": min_global_gap,
        "min_local_gap": min_local_gap,
        "argmin_mask": argmin_mask,
        "argmin_size": argmin_size,
        "local_bad_mask": local_bad_mask,
        "h_nodes": h_nodes,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Hall-subset scan on decimated weighted core model."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=21)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--max-heavy", type=int, default=24)
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    total_seen = 0
    total_considered = 0
    hall_fail = 0
    singleton_fail = 0
    local_fail = 0
    skipped_heavy = 0
    nonempty_h_trees = 0

    argmin_size_dist: Counter[int] = Counter()
    worst_hall = None
    worst_singleton = None
    worst_local = None
    per_n: dict[str, dict] = {}

    t_all = time.time()
    print("Decimated Hall subset scan", flush=True)
    print("Checks: Hall feasibility, singleton argmin, local singleton dominance", flush=True)
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_hall_fail = 0
        n_singleton_fail = 0
        n_local_fail = 0
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
            argmin_size_dist[res["argmin_size"]] += 1
            if res["m"] > 0:
                nonempty_h_trees += 1

            if not res["hall_ok"]:
                hall_fail += 1
                n_hall_fail += 1
            if not res["singleton_ok"]:
                singleton_fail += 1
                n_singleton_fail += 1
            if not res["local_ok"]:
                local_fail += 1
                n_local_fail += 1

            if res["m"] > 0 and (
                (worst_hall is None) or (res["min_slack"] < worst_hall["min_slack"])
            ):
                worst_hall = {
                    "n": n0,
                    "g6": g6,
                    "H_core_indices": res["h_nodes"],
                    "argmin_size": res["argmin_size"],
                    "argmin_subset": decode_subset(res["argmin_mask"], res["h_nodes"]),
                    "min_slack": res["min_slack"],
                    "min_singleton": res["min_singleton"],
                }
            if res["m"] > 0 and (
                (worst_singleton is None)
                or (res["min_global_gap"] < worst_singleton["min_global_gap"])
            ):
                worst_singleton = {
                    "n": n0,
                    "g6": g6,
                    "H_core_indices": res["h_nodes"],
                    "argmin_size": res["argmin_size"],
                    "argmin_subset": decode_subset(res["argmin_mask"], res["h_nodes"]),
                    "min_slack": res["min_slack"],
                    "min_singleton": res["min_singleton"],
                    "min_global_gap": res["min_global_gap"],
                }
            if res["m"] > 0 and (
                (worst_local is None) or (res["min_local_gap"] < worst_local["min_local_gap"])
            ):
                worst_local = {
                    "n": n0,
                    "g6": g6,
                    "H_core_indices": res["h_nodes"],
                    "argmin_size": res["argmin_size"],
                    "argmin_subset": decode_subset(res["argmin_mask"], res["h_nodes"]),
                    "local_bad_subset": decode_subset(res["local_bad_mask"], res["h_nodes"]),
                    "min_local_gap": res["min_local_gap"],
                }

        proc.wait()
        per_n[str(n)] = {
            "seen": n_seen,
            "considered": n_considered,
            "hall_fail": n_hall_fail,
            "singleton_fail": n_singleton_fail,
            "local_fail": n_local_fail,
            "skipped_heavy": n_skipped_heavy,
            "wall_s": time.time() - t0,
        }
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
            f"hall_fail={n_hall_fail:5d} singleton_fail={n_singleton_fail:5d} "
            f"local_fail={n_local_fail:5d} skipped_heavy={n_skipped_heavy:4d} "
            f"({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} "
        f"hall_fail={hall_fail:,} singleton_fail={singleton_fail:,} "
        f"local_fail={local_fail:,} skipped_heavy={skipped_heavy:,} "
        f"wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    print(f"Worst Hall slack witness: {worst_hall}", flush=True)
    print(f"Worst singleton gap witness: {worst_singleton}", flush=True)
    print(f"Worst local gap witness: {worst_local}", flush=True)
    print(f"argmin_size_dist={dict(argmin_size_dist)}", flush=True)

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
                "hall_fail": hall_fail,
                "singleton_fail": singleton_fail,
                "local_fail": local_fail,
                "skipped_heavy": skipped_heavy,
                "nonempty_h_trees": nonempty_h_trees,
                "worst_hall": worst_hall,
                "worst_singleton": worst_singleton,
                "worst_local": worst_local,
                "argmin_size_dist": dict(argmin_size_dist),
                "wall_s": time.time() - t_all,
            },
            "per_n": per_n,
        }
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
