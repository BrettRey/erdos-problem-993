#!/usr/bin/env python3
"""Generate exact certificates for leaf-heavy WHNC singleton-gap subsets.

For a non-empty heavy subset S:
  F(S) = supply(N(S)) - demand(S),
  supply(u) = 1/3 - P(u),
  demand(h) = P(h) - 1/3.

Let singleton slack be:
  sigma(h) = F({h}).
Define:
  m = min_{h in S} sigma(h),
  credit(h) = sigma(h) - m >= 0,
  overlap_cost(u) = (deg_S(u)-1) * supply(u), for u in N(S).

Canonical identity:
  F(S) - m = (|S|-1)*m + sum_h credit(h) - sum_u overlap_cost(u).

This script emits machine-checkable (exact Fraction) certificates for subsets,
with emphasis on the leaf-heavy obstruction class (all-negative marginals).
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from fractions import Fraction

from attack_conjA_weighted_compensation import (
    hard_core_probs_fraction,
    is_dleaf_le_1,
    parse_graph6,
)


def frac_str(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}"


def frac_float(x: Fraction) -> float:
    return float(x.numerator) / float(x.denominator)


def decode_subset(mask: int, nodes: list[int]) -> list[int]:
    out: list[int] = []
    rem = mask
    while rem:
        lsb = rem & -rem
        i = lsb.bit_length() - 1
        rem ^= lsb
        out.append(nodes[i])
    return out


def analyze_subset(
    adj: list[list[int]],
    probs: list[Fraction],
    heavy: list[int],
    subset: list[int],
) -> dict:
    one_third = Fraction(1, 3)
    h_set = set(heavy)
    s_set = set(subset)

    u_nodes = sorted({u for h in subset for u in adj[h] if u not in h_set})
    supply = {u: one_third - probs[u] for u in u_nodes}
    demand = {h: probs[h] - one_third for h in subset}

    n_s = sorted({u for h in subset for u in adj[h] if u not in h_set})
    f_s = sum(supply[u] for u in n_s) - sum(demand[h] for h in subset)

    singleton = {}
    for h in subset:
        n_h = [u for u in adj[h] if u not in h_set]
        singleton[h] = sum(supply[u] for u in n_h) - demand[h]
    m = min(singleton.values()) if singleton else Fraction(0)

    credit = {h: singleton[h] - m for h in subset}
    deg_s = {u: sum(1 for h in subset if u in adj[h]) for u in u_nodes}
    overlap_cost = {u: (deg_s[u] - 1) * supply[u] for u in u_nodes if deg_s[u] >= 2}

    lhs = f_s - m
    baseline = (len(subset) - 1) * m
    rhs = baseline + sum(credit.values()) - sum(overlap_cost.values())
    identity_ok = lhs == rhs

    # Marginals M(h,S) = F(S)-F(S\{h}) = supply(private_neighbors)-demand(h)
    marginals = {}
    all_negative_marginals = True if subset else False
    private_map = {}
    for h in subset:
        private_u = [u for u in adj[h] if u in u_nodes and deg_s[u] == 1]
        private_map[h] = private_u
        marg = sum(supply[u] for u in private_u) - demand[h]
        marginals[h] = marg
        if marg >= 0:
            all_negative_marginals = False

    leaves = [h for h in subset if len(adj[h]) == 1]
    nonleaves = [h for h in subset if len(adj[h]) > 1]

    return {
        "subset": subset,
        "subset_size": len(subset),
        "leaves_in_subset": leaves,
        "nonleaves_in_subset": nonleaves,
        "leaf_heavy_subset": len(leaves) > 0,
        "all_negative_marginals": all_negative_marginals,
        "F_S": {
            "exact": frac_str(f_s),
            "float": frac_float(f_s),
        },
        "min_singleton": {
            "exact": frac_str(m),
            "float": frac_float(m),
        },
        "gap": {
            "exact": frac_str(lhs),
            "float": frac_float(lhs),
        },
        "singleton_slack": {
            str(h): {"exact": frac_str(v), "float": frac_float(v)}
            for h, v in singleton.items()
        },
        "credit_terms": {
            str(h): {"exact": frac_str(v), "float": frac_float(v)}
            for h, v in credit.items()
        },
        "overlap_cost_terms": {
            str(u): {"exact": frac_str(v), "float": frac_float(v), "deg_S": deg_s[u]}
            for u, v in overlap_cost.items()
        },
        "marginals": {
            str(h): {
                "exact": frac_str(v),
                "float": frac_float(v),
                "private_neighbors": private_map[h],
            }
            for h, v in marginals.items()
        },
        "identity_check": {
            "lhs_gap_exact": frac_str(lhs),
            "rhs_rewrite_exact": frac_str(rhs),
            "baseline_exact": frac_str(baseline),
            "baseline_float": frac_float(baseline),
            "ok": identity_ok,
        },
        "supply": {
            str(u): {"exact": frac_str(v), "float": frac_float(v), "deg_S": deg_s[u]}
            for u, v in supply.items()
        },
        "demand": {
            str(h): {"exact": frac_str(v), "float": frac_float(v)} for h, v in demand.items()
        },
    }


def all_negative_leaf_heavy_subsets(
    adj: list[list[int]],
    probs: list[Fraction],
) -> tuple[list[int], list[list[int]]]:
    one_third = Fraction(1, 3)
    heavy = [v for v, p in enumerate(probs) if p > one_third]
    m = len(heavy)
    if m <= 1:
        return heavy, []

    h_set = set(heavy)
    u_nodes = sorted({u for h in heavy for u in adj[h] if u not in h_set})
    u_idx = {u: i for i, u in enumerate(u_nodes)}

    demand = [probs[h] - one_third for h in heavy]
    masks: list[int] = []
    for h in heavy:
        msk = 0
        for u in adj[h]:
            j = u_idx.get(u)
            if j is not None:
                msk |= 1 << j
        masks.append(msk)

    sup = [one_third - probs[u] for u in u_nodes]
    sup_sum = [Fraction(0)] * (1 << len(u_nodes))
    for s in range(1, 1 << len(u_nodes)):
        lsb = s & -s
        b = lsb.bit_length() - 1
        sup_sum[s] = sup_sum[s ^ lsb] + sup[b]

    union_mask = [0] * (1 << m)
    for s in range(1, 1 << m):
        lsb = s & -s
        b = lsb.bit_length() - 1
        union_mask[s] = union_mask[s ^ lsb] | masks[b]

    out: list[list[int]] = []
    for s in range(1, 1 << m):
        if s & (s - 1) == 0:
            continue
        subset = decode_subset(s, heavy)
        if not any(len(adj[h]) == 1 for h in subset):
            continue

        all_neg = True
        rem = s
        while rem:
            lsb = rem & -rem
            i = lsb.bit_length() - 1
            rem ^= lsb
            priv = masks[i] & ~union_mask[s ^ lsb]
            marg = sup_sum[priv] - demand[i]
            if marg >= 0:
                all_neg = False
                break
        if all_neg:
            out.append(subset)

    return heavy, out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Emit exact certificates for leaf-heavy all-negative subsets."
    )
    parser.add_argument("--min-n", type=int, default=8)
    parser.add_argument("--max-n", type=int, default=20)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument(
        "--max-certs",
        type=int,
        default=25,
        help="Maximum number of subset certificates to emit (smallest gaps kept).",
    )
    parser.add_argument("--out", default="results/whnc_leaf_block_certs.json")
    args = parser.parse_args()

    certs: list[dict] = []
    total_seen = 0
    total_considered = 0
    total_allneg_leaf = 0
    t_all = time.time()

    print("Leaf-block certificate scan", flush=True)
    print("Collecting leaf-heavy all-negative subsets with exact rewrite certificates", flush=True)
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

            probs = hard_core_probs_fraction(n0, adj)
            heavy, subsets = all_negative_leaf_heavy_subsets(adj, probs)
            if not subsets:
                continue

            g6 = line.decode("ascii").strip()
            for subset in subsets:
                n_subsets += 1
                total_allneg_leaf += 1
                cert = analyze_subset(adj, probs, heavy, subset)
                certs.append(
                    {
                        "n": n0,
                        "g6": g6,
                        "heavy": heavy,
                        "certificate": cert,
                    }
                )

        proc.wait()
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
            f"leaf_allneg_subsets={n_subsets:8d} ({time.time() - t0:.1f}s)",
            flush=True,
        )

    # Keep smallest gaps for a compact prototype artifact.
    certs.sort(key=lambda c: c["certificate"]["gap"]["float"])
    kept = certs[: args.max_certs]

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    payload = {
        "params": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "max_certs": args.max_certs,
        },
        "summary": {
            "seen": total_seen,
            "considered": total_considered,
            "leaf_allneg_subsets_total": total_allneg_leaf,
            "certificates_kept": len(kept),
            "wall_s": time.time() - t_all,
        },
        "certificates": kept,
    }
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print("-" * 80, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} "
        f"leaf_allneg_subsets={total_allneg_leaf:,} wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
