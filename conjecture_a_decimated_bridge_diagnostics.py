#!/usr/bin/env python3
"""Decimated bridge diagnostics for Conjecture A.

For d_leaf<=1 trees, this script packages several diagnostics used in the
decimated-peeling proof attempt:

1) Private-neighbor property on every non-empty S subseteq H.
2) Strict edge marginal positivity s(u)-d(h) over heavy-core incidences.
3) Support-outside-N(H) bound max P(u), u in A\\N(H).
4) Bridge lower bound g >= R where
     g = n/3 - mu(T),
     R = F(H) - |A cap N(H)|/6,
   and F(H) = s(N(H)) - d(H) in the decimated core model.
5) Identity decomposition check:
     g = sum_{v in (C\\A)\\(H U N(H))}(1/3-P(v))
       + sum_{u in A\\N(H)}(1/6-P(u)/2)
       + R.
6) Structural stress-test: minimum # non-support neighbors among heavy vertices.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import time
from typing import Any

from conjecture_a_decimation_core_model import decimate_tree, weighted_hard_core_probs
from conjecture_a_hall_subset_scan import hard_core_probs, is_dleaf_le_1, parse_graph6


def decode_subset(mask: int, nodes: list[int]) -> list[int]:
    out: list[int] = []
    rem = mask
    while rem:
        bit = rem & -rem
        i = bit.bit_length() - 1
        rem ^= bit
        out.append(nodes[i])
    return out


def has_private_neighbor_for_all_subsets(h_masks: list[int]) -> tuple[bool, int, int]:
    """Return (ok, subsets_checked, bad_mask)."""
    m = len(h_masks)
    if m == 0:
        return True, 0, 0

    h_pow = 1 << m
    nbr_mask = [0] * h_pow
    checked = 0
    for s in range(1, h_pow):
        lsb = s & -s
        b = lsb.bit_length() - 1
        prev = s ^ lsb
        nbr_mask[s] = nbr_mask[prev] | h_masks[b]
        checked += 1

        rem = s
        has_private = False
        while rem:
            bit = rem & -rem
            i = bit.bit_length() - 1
            rem ^= bit
            if h_masks[i] & ~nbr_mask[s ^ bit]:
                has_private = True
                break
        if not has_private:
            return False, checked, s

    return True, checked, 0


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run decimated bridge diagnostics through a tree range."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=23)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument(
        "--max-heavy",
        type=int,
        default=24,
        help="Skip trees where |H| exceeds this threshold (subset scan guard).",
    )
    parser.add_argument(
        "--verify-original-mu",
        action="store_true",
        help="Also compute mu(T) from original tree marginals and compare.",
    )
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    one_third = 1.0 / 3.0
    one_sixth = 1.0 / 6.0

    total_seen = 0
    total_considered = 0
    total_with_h = 0
    skipped_heavy = 0

    subsets_checked = 0
    private_fail = 0
    first_private_fail: dict[str, Any] | None = None

    heavy_edge_checks = 0
    edge_margin_fail = 0
    min_edge_margin = math.inf
    min_edge_margin_wit: dict[str, Any] | None = None

    support_outside_checks = 0
    support_outside_fail = 0
    max_support_outside_p = -math.inf
    max_support_outside_wit: dict[str, Any] | None = None

    bridge_viol = 0
    min_r_all = math.inf
    min_r_all_wit: dict[str, Any] | None = None
    min_r_nonempty_h = math.inf
    min_r_nonempty_h_wit: dict[str, Any] | None = None

    min_g = math.inf
    min_g_wit: dict[str, Any] | None = None

    decomp_max_err = 0.0
    decomp_max_err_wit: dict[str, Any] | None = None

    mu_core_max_err = 0.0
    mu_core_max_err_wit: dict[str, Any] | None = None
    mu_core_fail = 0

    min_non_support_neighbors = math.inf
    min_non_support_neighbors_wit: dict[str, Any] | None = None

    per_n: dict[str, dict[str, Any]] = {}
    t_all = time.time()

    print("Decimated bridge diagnostics", flush=True)
    print("Checks: private-neighbor, edge margins, support-outside, bridge g>=R", flush=True)
    print("-" * 88, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_with_h = 0
        n_skipped_heavy = 0
        n_subsets_checked = 0
        n_private_fail = 0
        n_edge_checks = 0
        n_edge_margin_fail = 0
        n_support_outside_checks = 0
        n_support_outside_fail = 0
        n_bridge_viol = 0
        n_mu_core_fail = 0

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
            core_non_support = [i for i in range(len(core_adj)) if i not in a_idx]

            # Core mean expression from exact decimation identity.
            mu_core = (
                0.5 * len(a_idx)
                + sum(probs[i] for i in core_non_support)
                + 0.5 * sum(probs[i] for i in a_idx)
            )
            g = n0 * one_third - mu_core
            if g < min_g:
                min_g = g
                min_g_wit = {"n": n0, "g6": g6, "g": g}

            if args.verify_original_mu:
                mu_orig = sum(hard_core_probs(n0, adj))
                mu_err = abs(mu_orig - mu_core)
                if mu_err > mu_core_max_err:
                    mu_core_max_err = mu_err
                    mu_core_max_err_wit = {
                        "n": n0,
                        "g6": g6,
                        "mu_orig": mu_orig,
                        "mu_core": mu_core,
                        "abs_err": mu_err,
                    }
                if mu_err > args.tol:
                    mu_core_fail += 1
                    n_mu_core_fail += 1

            h_nodes = [i for i in core_non_support if probs[i] > one_third + args.tol]
            m = len(h_nodes)
            if m > args.max_heavy:
                skipped_heavy += 1
                n_skipped_heavy += 1
                continue

            if m > 0:
                total_with_h += 1
                n_with_h += 1

            h_set = set(h_nodes)
            nh_nodes = {u for h in h_nodes for u in core_adj[h] if u not in h_set}
            a_inter_nh = sum(1 for u in nh_nodes if u in a_idx)

            # Supply helper.
            def supply(u: int) -> float:
                if u in a_idx:
                    return one_third - 0.5 * probs[u]
                return one_third - probs[u]

            # Edge marginal and structural counts over heavy incidences.
            for h in h_nodes:
                demand_h = probs[h] - one_third
                non_support_nbrs = 0
                for u in core_adj[h]:
                    if u not in a_idx:
                        non_support_nbrs += 1
                    heavy_edge_checks += 1
                    n_edge_checks += 1
                    margin = supply(u) - demand_h
                    if margin < min_edge_margin:
                        min_edge_margin = margin
                        min_edge_margin_wit = {
                            "n": n0,
                            "g6": g6,
                            "h_core_index": h,
                            "u_core_index": u,
                            "u_is_support": u in a_idx,
                            "P_h": probs[h],
                            "P_u": probs[u],
                            "margin": margin,
                        }
                    if margin <= args.tol:
                        edge_margin_fail += 1
                        n_edge_margin_fail += 1

                if non_support_nbrs < min_non_support_neighbors:
                    min_non_support_neighbors = non_support_nbrs
                    min_non_support_neighbors_wit = {
                        "n": n0,
                        "g6": g6,
                        "h_core_index": h,
                        "deg_h": len(core_adj[h]),
                        "non_support_neighbors": non_support_nbrs,
                        "neighbor_is_support": [u in a_idx for u in core_adj[h]],
                        "P_h": probs[h],
                    }

            # Supports outside N(H).
            for u in a_idx:
                if u in nh_nodes:
                    continue
                support_outside_checks += 1
                n_support_outside_checks += 1
                p = probs[u]
                if p > max_support_outside_p:
                    max_support_outside_p = p
                    max_support_outside_wit = {
                        "n": n0,
                        "g6": g6,
                        "u_core_index": u,
                        "P_u": p,
                        "H_size": m,
                        "NH_size": len(nh_nodes),
                    }
                if p > one_third + args.tol:
                    support_outside_fail += 1
                    n_support_outside_fail += 1

            # Decimated Hall slack for full H and bridge term.
            d_h = sum(probs[h] - one_third for h in h_nodes)
            f_h = sum(supply(u) for u in nh_nodes) - d_h
            r = f_h - a_inter_nh * one_sixth

            if r < min_r_all:
                min_r_all = r
                min_r_all_wit = {
                    "n": n0,
                    "g6": g6,
                    "R": r,
                    "F_H": f_h,
                    "A_cap_NH": a_inter_nh,
                    "H_size": m,
                    "NH_size": len(nh_nodes),
                }
            if m > 0 and r < min_r_nonempty_h:
                min_r_nonempty_h = r
                min_r_nonempty_h_wit = {
                    "n": n0,
                    "g6": g6,
                    "R": r,
                    "F_H": f_h,
                    "A_cap_NH": a_inter_nh,
                    "H_size": m,
                    "NH_size": len(nh_nodes),
                }

            if g + args.tol < r:
                bridge_viol += 1
                n_bridge_viol += 1

            # Identity decomposition check for g.
            h_union_nh = h_set | nh_nodes
            term_core = sum(one_third - probs[v] for v in core_non_support if v not in h_union_nh)
            term_support = sum(one_sixth - 0.5 * probs[u] for u in a_idx if u not in nh_nodes)
            rhs = term_core + term_support + r
            err = abs(g - rhs)
            if err > decomp_max_err:
                decomp_max_err = err
                decomp_max_err_wit = {
                    "n": n0,
                    "g6": g6,
                    "g": g,
                    "rhs": rhs,
                    "abs_err": err,
                    "term_core": term_core,
                    "term_support": term_support,
                    "R": r,
                }

            # Private-neighbor scan over all non-empty subsets of H.
            if m > 0:
                u_nodes = sorted(nh_nodes)
                u_idx = {u: i for i, u in enumerate(u_nodes)}
                h_masks: list[int] = []
                for h in h_nodes:
                    mask = 0
                    for u in core_adj[h]:
                        j = u_idx.get(u)
                        if j is not None:
                            mask |= 1 << j
                    h_masks.append(mask)

                ok, checked, bad_mask = has_private_neighbor_for_all_subsets(h_masks)
                subsets_checked += checked
                n_subsets_checked += checked
                if not ok:
                    private_fail += 1
                    n_private_fail += 1
                    if first_private_fail is None:
                        first_private_fail = {
                            "n": n0,
                            "g6": g6,
                            "H_core_indices": h_nodes,
                            "bad_subset_size": bad_mask.bit_count(),
                            "bad_subset": decode_subset(bad_mask, h_nodes),
                        }

        proc.wait()
        per_n[str(n)] = {
            "seen": n_seen,
            "considered": n_considered,
            "with_H": n_with_h,
            "skipped_heavy": n_skipped_heavy,
            "subsets_checked": n_subsets_checked,
            "private_fail": n_private_fail,
            "edge_checks": n_edge_checks,
            "edge_margin_fail": n_edge_margin_fail,
            "support_outside_checks": n_support_outside_checks,
            "support_outside_fail": n_support_outside_fail,
            "bridge_viol": n_bridge_viol,
            "mu_core_fail": n_mu_core_fail,
            "wall_s": time.time() - t0,
        }
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} with_H={n_with_h:7d} "
            f"subsets={n_subsets_checked:10d} priv_fail={n_private_fail:4d} "
            f"edge_fail={n_edge_margin_fail:4d} bridge_viol={n_bridge_viol:4d} "
            f"skip={n_skipped_heavy:4d} ({time.time()-t0:.1f}s)",
            flush=True,
        )

    print("-" * 88, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} with_H={total_with_h:,} "
        f"subsets={subsets_checked:,} private_fail={private_fail:,} edge_checks={heavy_edge_checks:,} "
        f"edge_margin_fail={edge_margin_fail:,} support_outside_checks={support_outside_checks:,} "
        f"support_outside_fail={support_outside_fail:,} bridge_viol={bridge_viol:,} "
        f"skipped_heavy={skipped_heavy:,} wall={time.time()-t_all:.1f}s",
        flush=True,
    )
    print(f"Min edge margin: {min_edge_margin}", flush=True)
    print(f"Max support-outside P: {max_support_outside_p}", flush=True)
    print(f"Min R (all trees): {min_r_all}", flush=True)
    print(f"Min R (H!=empty): {min_r_nonempty_h}", flush=True)
    print(f"Min g=n/3-mu: {min_g}", flush=True)
    print(f"Min non-support neighbors at heavy: {min_non_support_neighbors}", flush=True)
    print(f"Max decomposition error: {decomp_max_err}", flush=True)
    if args.verify_original_mu:
        print(f"Max |mu_orig-mu_core|: {mu_core_max_err}", flush=True)

    payload = {
        "params": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "tol": args.tol,
            "max_heavy": args.max_heavy,
            "verify_original_mu": args.verify_original_mu,
        },
        "summary": {
            "seen": total_seen,
            "considered": total_considered,
            "with_H": total_with_h,
            "skipped_heavy": skipped_heavy,
            "subsets_checked": subsets_checked,
            "private_fail": private_fail,
            "first_private_fail": first_private_fail,
            "edge_checks": heavy_edge_checks,
            "edge_margin_fail": edge_margin_fail,
            "minimum_edge_margin": min_edge_margin,
            "minimum_edge_margin_witness": min_edge_margin_wit,
            "support_outside_checks": support_outside_checks,
            "support_outside_fail": support_outside_fail,
            "maximum_support_outside_p": max_support_outside_p,
            "maximum_support_outside_witness": max_support_outside_wit,
            "bridge_viol": bridge_viol,
            "minimum_r_all": min_r_all,
            "minimum_r_all_witness": min_r_all_wit,
            "minimum_r_nonempty_h": min_r_nonempty_h,
            "minimum_r_nonempty_h_witness": min_r_nonempty_h_wit,
            "minimum_g": min_g,
            "minimum_g_witness": min_g_wit,
            "decomposition_max_abs_err": decomp_max_err,
            "decomposition_max_abs_err_witness": decomp_max_err_wit,
            "mu_core_fail": mu_core_fail,
            "mu_core_max_abs_err": mu_core_max_err,
            "mu_core_max_abs_err_witness": mu_core_max_err_wit,
            "minimum_non_support_neighbors_at_heavy": min_non_support_neighbors,
            "minimum_non_support_neighbors_witness": min_non_support_neighbors_wit,
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
