#!/usr/bin/env python3
"""Weighted WHNC scan on decimated core model.

Decimation (for d_leaf<=1 tree T):
  - Remove leaves to get core C.
  - Supports A (core vertices adjacent to deleted leaves) get lambda_v = 1/2.
  - Other core vertices get lambda_v = 1.

On C with these lambda_v, define:
  H = {v in C\\A : P(v) > 1/3}.

Weighted neighbor-supply:
  supply_w(u) = (1/3 - P(u))                 if u in C\\A
              = (1/2)*(2/3 - P(u))           if u in A

Checks:
  (A) Global weighted neighbor inequality (decimated WHNC):
      sum_{h in H} (P(h)-1/3) <= sum_{u in N(H)} supply_w(u)

  (B) Local weighted overlap condition:
      for every u in N(H):
      sum_{h in N(u) cap H} (P(h)-1/3) <= supply_w(u)
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time

from conjecture_a_decimation_core_model import decimate_tree, weighted_hard_core_probs
from conjecture_a_hall_subset_scan import is_dleaf_le_1, parse_graph6


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Scan decimated weighted WHNC (global + local overlap)."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=21)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    one_third = 1.0 / 3.0
    total_seen = 0
    total_considered = 0
    global_fail = 0
    local_fail = 0
    min_global_margin = None
    min_local_margin = None
    worst_global = None
    worst_local = None
    per_n: dict[str, dict] = {}

    t_all = time.time()
    print("Decimated weighted WHNC scan", flush=True)
    print("Checks: global weighted neighbor inequality + local weighted overlap", flush=True)
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_global_fail = 0
        n_local_fail = 0

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
            core_adj, core_to_orig, _leaves, supports, lam = decimate_tree(adj)
            probs = weighted_hard_core_probs(core_adj, lam)

            a_idx = {i for i, v in enumerate(core_to_orig) if v in supports}
            h_set = {i for i in range(len(core_adj)) if i not in a_idx and probs[i] > one_third + args.tol}
            if not h_set:
                continue

            excess = sum(probs[h] - one_third for h in h_set)
            n_h = {u for h in h_set for u in core_adj[h] if u not in h_set}
            supply = 0.0
            for u in n_h:
                if u in a_idx:
                    supply += 0.5 * (2.0 * one_third - probs[u])
                else:
                    supply += one_third - probs[u]
            g_margin = supply - excess
            if (min_global_margin is None) or (g_margin < min_global_margin):
                min_global_margin = g_margin
                worst_global = {
                    "n": n0,
                    "g6": g6,
                    "core_size": len(core_adj),
                    "A_size": len(a_idx),
                    "H_size": len(h_set),
                    "excess": excess,
                    "supply": supply,
                    "margin": g_margin,
                }
            if g_margin < -args.tol:
                global_fail += 1
                n_global_fail += 1

            for u in n_h:
                hs = [h for h in core_adj[u] if h in h_set]
                local_excess = sum(probs[h] - one_third for h in hs)
                if u in a_idx:
                    local_supply = 0.5 * (2.0 * one_third - probs[u])
                else:
                    local_supply = one_third - probs[u]
                l_margin = local_supply - local_excess
                if (min_local_margin is None) or (l_margin < min_local_margin):
                    min_local_margin = l_margin
                    worst_local = {
                        "n": n0,
                        "g6": g6,
                        "u_core_index": u,
                        "u_is_A": u in a_idx,
                        "degH_u": len(hs),
                        "heavy_neighbors": hs,
                        "local_excess": local_excess,
                        "local_supply": local_supply,
                        "margin": l_margin,
                    }
                if l_margin < -args.tol:
                    local_fail += 1
                    n_local_fail += 1

        proc.wait()
        per_n[str(n)] = {
            "seen": n_seen,
            "considered": n_considered,
            "global_fail": n_global_fail,
            "local_fail": n_local_fail,
            "wall_s": time.time() - t0,
        }
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
            f"global_fail={n_global_fail:5d} local_fail={n_local_fail:5d} "
            f"({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} "
        f"global_fail={global_fail:,} local_fail={local_fail:,} "
        f"wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    print(f"Minimum global margin witness: {worst_global}", flush=True)
    print(f"Minimum local margin witness: {worst_local}", flush=True)

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
                "global_fail": global_fail,
                "local_fail": local_fail,
                "minimum_global_margin": min_global_margin,
                "minimum_local_margin": min_local_margin,
                "worst_global": worst_global,
                "worst_local": worst_local,
                "wall_s": time.time() - t_all,
            },
            "per_n": per_n,
        }
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
