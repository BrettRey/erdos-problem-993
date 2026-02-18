#!/usr/bin/env python3
"""Verify analytic slack identity for signature (3,2,1,(2,2)).

Signature definition for subset S ⊆ H:
  - |S| = 3
  - exactly two heavy leaves and one heavy non-leaf
  - N(S) has two vertices, each with deg_S = 2

For such S = {h,a,b} (h non-leaf, a,b leaves) and N(S)={u,v}, checks:
  F(S) = 2/3 - P(h) - (P(u)+P(v))/2
       = 1/2 * ((2/3 - P(h)-P(u)) + (2/3 - P(h)-P(v))).
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from itertools import combinations

from conjecture_a_hall_subset_scan import hard_core_probs, is_dleaf_le_1, parse_graph6


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Signature (3,2,1,(2,2)) slack identity/positivity verifier."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=23)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    one_third = 1.0 / 3.0
    total_seen = 0
    total_considered = 0
    total_sig = 0
    id_fail = 0
    pos_fail = 0
    struct_fail = 0
    per_n: dict[str, dict] = {}

    worst_diff = None
    min_slack = None
    min_avg_edge_surplus = None

    t_all = time.time()
    print("Signature (3,2,1,(2,2)) verifier", flush=True)
    print("Checking analytic identity and strict positivity", flush=True)
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_sig = 0
        n_id_fail = 0
        n_pos_fail = 0
        n_struct_fail = 0

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
            if len(h_nodes) < 3:
                continue
            h_set = set(h_nodes)
            g6 = line.decode("ascii").strip()

            for triple in combinations(h_nodes, 3):
                leaves = [x for x in triple if len(adj[x]) == 1]
                nonleaves = [x for x in triple if len(adj[x]) > 1]
                if len(leaves) != 2 or len(nonleaves) != 1:
                    continue
                h = nonleaves[0]

                u_nodes = sorted({u for x in triple for u in adj[x] if u not in h_set})
                if len(u_nodes) != 2:
                    continue

                deg_s = []
                bad_struct = False
                for u in u_nodes:
                    ds = sum(1 for x in triple if u in adj[x])
                    deg_s.append(ds)
                    if h not in adj[u]:
                        bad_struct = True
                deg_s.sort()
                if deg_s != [2, 2]:
                    continue

                if bad_struct:
                    # For this signature, each u should be adjacent to the unique non-leaf h.
                    n_struct_fail += 1
                    struct_fail += 1
                    continue

                n_sig += 1
                total_sig += 1

                # WHNC slack F(S)
                supply = sum(one_third - probs[u] for u in u_nodes)
                demand = sum(probs[x] - one_third for x in triple)
                f_s = supply - demand

                rhs = one_third * 2.0 - probs[h] - 0.5 * (probs[u_nodes[0]] + probs[u_nodes[1]])
                diff = abs(f_s - rhs)
                if diff > args.tol:
                    n_id_fail += 1
                    id_fail += 1

                edge_surplus_1 = one_third * 2.0 - probs[h] - probs[u_nodes[0]]
                edge_surplus_2 = one_third * 2.0 - probs[h] - probs[u_nodes[1]]
                avg_edge_surplus = 0.5 * (edge_surplus_1 + edge_surplus_2)

                if f_s <= args.tol:
                    n_pos_fail += 1
                    pos_fail += 1

                if (worst_diff is None) or (diff > worst_diff["diff"] + 1e-18):
                    worst_diff = {
                        "n": n0,
                        "g6": g6,
                        "S": list(triple),
                        "h_nonleaf": h,
                        "u_nodes": u_nodes,
                        "F_S": f_s,
                        "rhs": rhs,
                        "diff": diff,
                    }
                if (min_slack is None) or (f_s < min_slack["F_S"] - 1e-18):
                    min_slack = {
                        "n": n0,
                        "g6": g6,
                        "S": list(triple),
                        "h_nonleaf": h,
                        "u_nodes": u_nodes,
                        "F_S": f_s,
                        "edge_surplus_hu1": edge_surplus_1,
                        "edge_surplus_hu2": edge_surplus_2,
                        "avg_edge_surplus": avg_edge_surplus,
                    }
                if (min_avg_edge_surplus is None) or (
                    avg_edge_surplus < min_avg_edge_surplus["avg_edge_surplus"] - 1e-18
                ):
                    min_avg_edge_surplus = {
                        "n": n0,
                        "g6": g6,
                        "S": list(triple),
                        "h_nonleaf": h,
                        "u_nodes": u_nodes,
                        "F_S": f_s,
                        "edge_surplus_hu1": edge_surplus_1,
                        "edge_surplus_hu2": edge_surplus_2,
                        "avg_edge_surplus": avg_edge_surplus,
                    }

        proc.wait()
        per_n[str(n)] = {
            "seen": n_seen,
            "considered": n_considered,
            "signature_count": n_sig,
            "identity_fail": n_id_fail,
            "positivity_fail": n_pos_fail,
            "structure_fail": n_struct_fail,
            "wall_s": time.time() - t0,
        }
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
            f"signature={n_sig:8d} id_fail={n_id_fail:4d} pos_fail={n_pos_fail:4d} "
            f"struct_fail={n_struct_fail:4d} ({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} "
        f"signature={total_sig:,} id_fail={id_fail:,} pos_fail={pos_fail:,} "
        f"struct_fail={struct_fail:,} wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    print(f"Worst identity diff witness: {worst_diff}", flush=True)
    print(f"Minimum F(S) witness: {min_slack}", flush=True)
    print(f"Minimum avg-edge-surplus witness: {min_avg_edge_surplus}", flush=True)

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
                "signature_count": total_sig,
                "identity_fail": id_fail,
                "positivity_fail": pos_fail,
                "structure_fail": struct_fail,
                "worst_identity_diff": worst_diff,
                "minimum_slack": min_slack,
                "minimum_avg_edge_surplus": min_avg_edge_surplus,
                "wall_s": time.time() - t_all,
            },
            "per_n": per_n,
        }
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
