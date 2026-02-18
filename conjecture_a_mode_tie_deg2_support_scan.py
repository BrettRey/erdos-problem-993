#!/usr/bin/env python3
"""Scan the degree-2-support leaf bridge condition for focused mode-tie.

For each d_leaf<=1 tree T, let m be leftmost mode of I_T and
lambda_m = i_{m-1}/i_m. For each leaf l whose support s has degree 2,
check:

    Phi_m(T-l; lambda_m) >= 0
    Phi_{m-1}(T-{l,s}; lambda_m) >= 0

where Phi_q(F;lambda) = lambda I'_F(lambda) - (q-1) I_F(lambda).
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
from indpoly import independence_poly


def mode_index_leftmost(poly: list[int]) -> int:
    return max(range(len(poly)), key=lambda i: poly[i])


def phi_q(poly: list[int], q: int, lam: float) -> float:
    z = 0.0
    mu_num = 0.0
    p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p
        z += w
        mu_num += k * w
        p *= lam
    mu = mu_num / z
    return z * (mu - (q - 1))


def remove_vertices(adj: list[list[int]], remove_set: set[int]) -> list[list[int]]:
    keep = [v for v in range(len(adj)) if v not in remove_set]
    idx = {v: i for i, v in enumerate(keep)}
    out = [[] for _ in keep]
    for v in keep:
        vv = idx[v]
        for u in adj[v]:
            if u in idx:
                out[vv].append(idx[u])
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Degree-2-support mode-tie bridge scan")
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--stop-on-first", action="store_true")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    total_seen = 0
    total_considered = 0
    total_checked_trees = 0
    total_checked_deg2_leaves = 0
    fail = 0
    first_fail: dict[str, Any] | None = None

    per_n: dict[str, dict[str, Any]] = {}

    t_all = time.time()
    print("Degree-2-support bridge scan on d_leaf<=1", flush=True)
    print("Check: Phi_m(T-l)>=0 and Phi_{m-1}(T-{l,s})>=0 at lambda_m(T)", flush=True)
    print("-" * 92, flush=True)

    should_stop = False

    for n in range(args.min_n, args.max_n + 1):
        if should_stop:
            break

        t0 = time.time()
        n_seen = 0
        n_considered = 0
        n_checked_trees = 0
        n_checked_deg2_leaves = 0
        n_fail = 0

        proc = subprocess.Popen(
            [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"],
            stdout=subprocess.PIPE,
        )
        assert proc.stdout is not None

        for line in proc.stdout:
            nn, adj = parse_graph6(line)
            n_seen += 1
            total_seen += 1

            if not is_dleaf_le_1(nn, adj):
                continue

            n_considered += 1
            total_considered += 1

            poly = independence_poly(nn, adj)
            m = mode_index_leftmost(poly)
            if m == 0 or poly[m - 1] <= 0 or poly[m] <= 0:
                continue

            lam = poly[m - 1] / poly[m]
            n_checked_trees += 1
            total_checked_trees += 1

            deg = [len(nb) for nb in adj]
            for l, dl in enumerate(deg):
                if dl != 1:
                    continue
                s = adj[l][0]
                if deg[s] != 2:
                    continue

                n_checked_deg2_leaves += 1
                total_checked_deg2_leaves += 1

                a_adj = remove_vertices(adj, {l})
                b_adj = remove_vertices(adj, {l, s})
                a_poly = independence_poly(len(a_adj), a_adj)
                b_poly = independence_poly(len(b_adj), b_adj)

                phi_a = phi_q(a_poly, m, lam)
                phi_b = phi_q(b_poly, m - 1, lam)

                if phi_a < -args.tol or phi_b < -args.tol:
                    fail += 1
                    n_fail += 1
                    if first_fail is None:
                        first_fail = {
                            "n": nn,
                            "g6": line.decode("ascii").strip(),
                            "mode": m,
                            "lambda_mode": lam,
                            "leaf": l,
                            "support": s,
                            "phi_a": phi_a,
                            "phi_b": phi_b,
                        }
                    if args.stop_on_first:
                        should_stop = True
                        proc.kill()
                        break

        proc.wait()

        per_n[str(n)] = {
            "seen": n_seen,
            "considered": n_considered,
            "checked_trees": n_checked_trees,
            "checked_deg2_leaves": n_checked_deg2_leaves,
            "fail": n_fail,
            "wall_s": time.time() - t0,
        }

        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
            f"trees={n_checked_trees:8d} deg2_leaves={n_checked_deg2_leaves:9d} "
            f"fail={n_fail:4d} ({time.time()-t0:.1f}s)",
            flush=True,
        )

    print("-" * 92, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} trees={total_checked_trees:,} "
        f"deg2_leaves={total_checked_deg2_leaves:,} fail={fail:,} wall={time.time()-t_all:.1f}s",
        flush=True,
    )
    if first_fail:
        print(f"First failure: {first_fail}", flush=True)

    payload = {
        "params": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "tol": args.tol,
            "stop_on_first": args.stop_on_first,
        },
        "summary": {
            "seen": total_seen,
            "considered": total_considered,
            "checked_trees": total_checked_trees,
            "checked_deg2_leaves": total_checked_deg2_leaves,
            "fail": fail,
            "first_fail": first_fail,
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
