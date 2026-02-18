#!/usr/bin/env python3
"""Scan leaf-decomposition bridge candidates for focused mode-tie inequality.

For each tree T, with m = leftmost mode(I_T) at lambda=1 and
lambda_m = i_{m-1}/i_m, we inspect each leaf l with support s and define:

- A = T - l
- B = T - {l,s}
- Phi_q(F; lambda) := lambda I'_F(lambda) - (q-1) I_F(lambda)

Identity:
  Phi_m(T;lambda) = Phi_m(A;lambda) + lambda * Phi_{m-1}(B;lambda).

This scanner checks whether there exists a leaf with both
Phi_m(A;lambda_m(T)) >= 0 and Phi_{m-1}(B;lambda_m(T)) >= 0.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from collections import Counter
from typing import Any

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from graph6 import parse_graph6
from indpoly import independence_poly


def mode_index_leftmost(poly: list[int]) -> int:
    return max(range(len(poly)), key=lambda i: poly[i])


def remove_vertices(adj: list[list[int]], remove_set: set[int]) -> list[list[int]]:
    n = len(adj)
    keep = [v for v in range(n) if v not in remove_set]
    idx = {v: i for i, v in enumerate(keep)}
    out = [[] for _ in keep]
    for v in keep:
        vv = idx[v]
        for u in adj[v]:
            if u in idx:
                out[vv].append(idx[u])
    return out


def phi_q(poly: list[int], q: int, lam: float) -> float:
    # Phi_q = lam I' - (q-1) I = I * (mu - (q-1)).
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


def choose_leaf_by_parent_degree(
    leaves: list[int],
    parent: dict[int, int],
    deg: list[int],
    pick_max: bool,
) -> int:
    key = (lambda l: (deg[parent[l]], -l)) if pick_max else (lambda l: (-deg[parent[l]], -l))
    # Because we negated leaf index in key, max() keeps smallest leaf on tie.
    return max(leaves, key=key)


def choose_leaf_parent_deg2_else_max(
    leaves: list[int],
    parent: dict[int, int],
    deg: list[int],
) -> int:
    cand = [l for l in leaves if deg[parent[l]] == 2]
    if cand:
        return min(cand)
    best_deg = max(deg[parent[l]] for l in leaves)
    return min(l for l in leaves if deg[parent[l]] == best_deg)


def evaluate_leaf(adj: list[list[int]], leaf: int, m: int, lam: float) -> tuple[float, float, int, int]:
    s = adj[leaf][0]
    a_adj = remove_vertices(adj, {leaf})
    b_adj = remove_vertices(adj, {leaf, s})

    a_poly = independence_poly(len(a_adj), a_adj)
    b_poly = independence_poly(len(b_adj), b_adj)

    phi_a = phi_q(a_poly, m, lam)
    phi_b = phi_q(b_poly, m - 1, lam)
    mode_a = mode_index_leftmost(a_poly)
    mode_b = mode_index_leftmost(b_poly) if b_poly else 0
    return phi_a, phi_b, mode_a, mode_b


def main() -> None:
    ap = argparse.ArgumentParser(description="Leaf-bridge scan for focused mode-tie proof.")
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--all-trees", action="store_true")
    ap.add_argument("--stop-on-first", action="store_true")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    total_seen = 0
    total_considered = 0
    total_checked = 0

    exists_both_nonneg_fail = 0
    all_both_nonneg_fail = 0

    heuristic_names = [
        "first_leaf",
        "max_parent_deg",
        "min_parent_deg",
        "parent_deg2_else_max",
    ]
    heuristic_fail = {k: 0 for k in heuristic_names}

    first_exists_fail: dict[str, Any] | None = None
    first_heuristic_fail: dict[str, dict[str, Any]] = {}

    worst_best_margin = float("inf")
    worst_best_witness: dict[str, Any] | None = None

    per_n: dict[str, dict[str, Any]] = {}

    t_all = time.time()
    scope = "all trees" if args.all_trees else "d_leaf<=1"
    print(f"Leaf-bridge scan on {scope}", flush=True)
    print("Property: exists leaf with Phi_m(A)>=0 and Phi_{m-1}(B)>=0 at lambda_m(T)", flush=True)
    print("-" * 92, flush=True)

    should_stop = False

    for n in range(args.min_n, args.max_n + 1):
        if should_stop:
            break

        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        assert proc.stdout is not None

        n_seen = 0
        n_considered = 0
        n_checked = 0
        n_exists_fail = 0
        n_all_fail = 0
        n_heur_fail = {k: 0 for k in heuristic_names}

        for line in proc.stdout:
            n0, adj = parse_graph6(line)
            n_seen += 1
            total_seen += 1

            if (not args.all_trees) and (not is_dleaf_le_1(n0, adj)):
                continue

            n_considered += 1
            total_considered += 1

            poly = independence_poly(n0, adj)
            m = mode_index_leftmost(poly)
            if m == 0 or poly[m - 1] <= 0 or poly[m] <= 0:
                continue

            lam = poly[m - 1] / poly[m]
            n_checked += 1
            total_checked += 1

            deg = [len(nb) for nb in adj]
            leaves = [v for v in range(n0) if deg[v] == 1]
            parent = {l: adj[l][0] for l in leaves}

            leaf_data: dict[int, tuple[float, float, int, int]] = {}
            any_both = False
            all_both = True
            best_margin = float("-inf")
            best_leaf = -1

            for l in leaves:
                phi_a, phi_b, mode_a, mode_b = evaluate_leaf(adj, l, m, lam)
                leaf_data[l] = (phi_a, phi_b, mode_a, mode_b)
                if phi_a >= -args.tol and phi_b >= -args.tol:
                    any_both = True
                else:
                    all_both = False
                mm = min(phi_a, phi_b)
                if mm > best_margin:
                    best_margin = mm
                    best_leaf = l

            if best_margin < worst_best_margin:
                worst_best_margin = best_margin
                phi_a, phi_b, mode_a, mode_b = leaf_data[best_leaf]
                worst_best_witness = {
                    "n": n0,
                    "g6": line.decode("ascii").strip(),
                    "mode": m,
                    "lambda_mode": lam,
                    "best_leaf": best_leaf,
                    "best_support": parent[best_leaf],
                    "best_margin": best_margin,
                    "phi_a": phi_a,
                    "phi_b": phi_b,
                    "mode_a": mode_a,
                    "mode_b": mode_b,
                    "degree_signature": dict(sorted(Counter(deg).items())),
                }

            if not any_both:
                n_exists_fail += 1
                exists_both_nonneg_fail += 1
                if first_exists_fail is None:
                    first_exists_fail = {
                        "n": n0,
                        "g6": line.decode("ascii").strip(),
                        "mode": m,
                        "lambda_mode": lam,
                        "degree_signature": dict(sorted(Counter(deg).items())),
                    }
                    if args.stop_on_first:
                        should_stop = True
                        proc.kill()
                        break

            if not all_both:
                n_all_fail += 1
                all_both_nonneg_fail += 1

            # Heuristics
            choices = {
                "first_leaf": min(leaves),
                "max_parent_deg": choose_leaf_by_parent_degree(leaves, parent, deg, True),
                "min_parent_deg": choose_leaf_by_parent_degree(leaves, parent, deg, False),
                "parent_deg2_else_max": choose_leaf_parent_deg2_else_max(leaves, parent, deg),
            }
            for name, l in choices.items():
                phi_a, phi_b, mode_a, mode_b = leaf_data[l]
                ok = (phi_a >= -args.tol and phi_b >= -args.tol)
                if not ok:
                    n_heur_fail[name] += 1
                    heuristic_fail[name] += 1
                    if name not in first_heuristic_fail:
                        first_heuristic_fail[name] = {
                            "n": n0,
                            "g6": line.decode("ascii").strip(),
                            "mode": m,
                            "lambda_mode": lam,
                            "leaf": l,
                            "support": parent[l],
                            "phi_a": phi_a,
                            "phi_b": phi_b,
                            "mode_a": mode_a,
                            "mode_b": mode_b,
                            "degree_signature": dict(sorted(Counter(deg).items())),
                        }

        proc.wait()

        per_n[str(n)] = {
            "seen": n_seen,
            "considered": n_considered,
            "checked": n_checked,
            "exists_both_nonneg_fail": n_exists_fail,
            "all_both_nonneg_fail": n_all_fail,
            "heuristic_fail": n_heur_fail,
            "wall_s": time.time() - t0,
        }

        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} checked={n_checked:8d} "
            f"exists_fail={n_exists_fail:4d} all_fail={n_all_fail:4d} "
            f"h_max_fail={n_heur_fail['max_parent_deg']:4d} h_deg2_fail={n_heur_fail['parent_deg2_else_max']:4d} "
            f"({time.time()-t0:.1f}s)",
            flush=True,
        )

    print("-" * 92, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} checked={total_checked:,} "
        f"exists_fail={exists_both_nonneg_fail:,} all_fail={all_both_nonneg_fail:,} wall={time.time()-t_all:.1f}s",
        flush=True,
    )
    print(f"Heuristic failures: {heuristic_fail}", flush=True)
    print(f"Worst best-leaf margin: {worst_best_margin}", flush=True)
    print(f"Worst best-leaf witness: {worst_best_witness}", flush=True)
    if first_exists_fail:
        print(f"First existence failure: {first_exists_fail}", flush=True)
    if first_heuristic_fail:
        print(f"First heuristic failures: {first_heuristic_fail}", flush=True)

    payload = {
        "params": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "all_trees": args.all_trees,
            "tol": args.tol,
            "stop_on_first": args.stop_on_first,
        },
        "summary": {
            "seen": total_seen,
            "considered": total_considered,
            "checked": total_checked,
            "exists_both_nonneg_fail": exists_both_nonneg_fail,
            "all_both_nonneg_fail": all_both_nonneg_fail,
            "heuristic_fail": heuristic_fail,
            "worst_best_leaf_margin": worst_best_margin,
            "worst_best_leaf_witness": worst_best_witness,
            "first_exists_fail": first_exists_fail,
            "first_heuristic_fail": first_heuristic_fail,
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
