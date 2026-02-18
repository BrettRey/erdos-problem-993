#!/usr/bin/env python3
"""Leaf decimation reduction for d_leaf <= 1 trees.

For a tree T at lambda=1 with d_leaf<=1:
  - Let L be leaves.
  - Let A be supports adjacent to leaves (each support has at most one leaf).
  - Let C = V(T) \\ L (leaf-stripped core graph; still a tree when non-empty).

Decimation identity:
  Z_T = 2^|L| * Z_C^(lam),
where Z_C^(lam) is hard-core partition function on C with vertex activities:
  lam(v)=1/2 if v in A, otherwise lam(v)=1.

Marginals:
  P_T(v) = P_C^(lam)(v) for v in C,
  P_T(leaf l with support s) = (1 - P_C^(lam)(s))/2.

Mean identity:
  mu(T) = |A|/2 + sum_{v in C\\A} P_C(v) + (1/2) sum_{v in A} P_C(v).

Hence:
  n/3 - mu(T)
  = (|C|/3 - |A|/6)
    - [sum_{v in C\\A} P_C(v) + (1/2) sum_{v in A} P_C(v)].
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from fractions import Fraction

from conjecture_a_hall_subset_scan import hard_core_probs, is_dleaf_le_1, parse_graph6


def rooted_order(adj: list[list[int]], root: int = 0) -> tuple[list[int], list[int], list[list[int]]]:
    n = len(adj)
    parent = [-1] * n
    children = [[] for _ in range(n)]
    order: list[int] = []
    q = [root]
    seen = [False] * n
    seen[root] = True
    head = 0
    while head < len(q):
        v = q[head]
        head += 1
        order.append(v)
        for u in adj[v]:
            if not seen[u]:
                seen[u] = True
                parent[u] = v
                children[v].append(u)
                q.append(u)
    return order, parent, children


def weighted_partition_tree(adj: list[list[int]], lam: list[Fraction]) -> Fraction:
    n = len(adj)
    if n == 0:
        return Fraction(1, 1)
    order, _parent, children = rooted_order(adj, 0)
    exc = [Fraction(1, 1)] * n
    inc = [Fraction(0, 1)] * n
    for v in reversed(order):
        prod_exc = Fraction(1, 1)
        prod_all = Fraction(1, 1)
        for c in children[v]:
            prod_exc *= exc[c]
            prod_all *= exc[c] + inc[c]
        inc[v] = lam[v] * prod_exc
        exc[v] = prod_all
    return exc[0] + inc[0]


def weighted_hard_core_probs(adj: list[list[int]], lam: list[Fraction]) -> list[float]:
    """Compute occupation probabilities on a tree with vertex-dependent lambda_v."""
    n = len(adj)
    if n == 0:
        return []
    if n == 1:
        lv = float(lam[0])
        return [lv / (1.0 + lv)]

    order, parent, children = rooted_order(adj, 0)

    # R_up[v] = R_{v -> parent(v)}.
    r_up = [0.0] * n
    lam_f = [float(x) for x in lam]
    for v in reversed(order):
        prod = 1.0
        for c in children[v]:
            prod *= 1.0 / (1.0 + r_up[c])
        r_up[v] = lam_f[v] * prod

    # R_down[v] = R_{parent(v) -> v} for non-root v.
    r_down = [0.0] * n
    for v in order:
        m = len(children[v])
        if m == 0:
            continue

        # Prefix/suffix products of 1/(1 + R_{child->v}) for sibling exclusion.
        vals = [1.0 / (1.0 + r_up[c]) for c in children[v]]
        pref = [1.0] * (m + 1)
        suff = [1.0] * (m + 1)
        for i in range(m):
            pref[i + 1] = pref[i] * vals[i]
        for i in range(m - 1, -1, -1):
            suff[i] = suff[i + 1] * vals[i]

        for i, c in enumerate(children[v]):
            sib_prod = pref[i] * suff[i + 1]
            if parent[v] == -1:
                # Root has no parent-side term.
                r_down[c] = lam_f[v] * sib_prod
            else:
                parent_term = 1.0 / (1.0 + r_down[v])
                r_down[c] = lam_f[v] * parent_term * sib_prod

    probs = [0.0] * n
    for v in range(n):
        prod = 1.0
        for c in children[v]:
            prod *= 1.0 + r_up[c]
        if parent[v] != -1:
            prod *= 1.0 + r_down[v]
        r_v = lam_f[v] / prod
        probs[v] = r_v / (1.0 + r_v)
    return probs


def decimate_tree(
    adj: list[list[int]],
) -> tuple[list[list[int]], list[int], set[int], set[int], list[Fraction]]:
    """Return core graph and maps/activity data.

    Returns:
      - core adjacency
      - core_to_orig map
      - leaves set (orig indices)
      - supports set (orig indices in core, adjacent to leaves)
      - lambda list on core (1 or 1/2)
    """
    n = len(adj)
    leaves = {v for v in range(n) if len(adj[v]) == 1}
    supports = {adj[l][0] for l in leaves}

    core_to_orig = [v for v in range(n) if v not in leaves]
    orig_to_core = {v: i for i, v in enumerate(core_to_orig)}
    m = len(core_to_orig)
    core_adj = [[] for _ in range(m)]
    for i, v in enumerate(core_to_orig):
        for u in adj[v]:
            j = orig_to_core.get(u)
            if j is not None:
                core_adj[i].append(j)

    lam = [Fraction(1, 2) if v in supports else Fraction(1, 1) for v in core_to_orig]
    return core_adj, core_to_orig, leaves, supports, lam


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Decimation reduction checks for d_leaf<=1 trees."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=21)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-10)
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    total_seen = 0
    total_considered = 0
    z_fail = 0
    core_prob_fail = 0
    leaf_prob_fail = 0
    gap_identity_fail = 0
    min_gap = None
    min_gap_core = None
    per_n: dict[str, dict] = {}

    t_all = time.time()
    print("Decimation core-model scan", flush=True)
    print("Verifying partition/marginal/gap identities", flush=True)
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_z_fail = 0
        n_core_fail = 0
        n_leaf_fail = 0
        n_gap_fail = 0

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
            probs_t = hard_core_probs(n0, adj)
            mu_t = sum(probs_t)
            gap_t = n0 / 3.0 - mu_t

            core_adj, core_to_orig, leaves, supports, lam = decimate_tree(adj)
            probs_c = weighted_hard_core_probs(core_adj, lam)

            # Partition function identity: Z_T = 2^|L| Z_C^(lam)
            z_t = Fraction(weighted_partition_tree(adj, [Fraction(1, 1)] * n0))
            z_c = weighted_partition_tree(core_adj, lam)
            lhs = z_t
            rhs = (Fraction(2, 1) ** len(leaves)) * z_c
            if lhs != rhs:
                z_fail += 1
                n_z_fail += 1

            # Core marginals identity
            orig_to_core = {v: i for i, v in enumerate(core_to_orig)}
            local_core_fail = False
            for v in core_to_orig:
                i = orig_to_core[v]
                if abs(probs_t[v] - probs_c[i]) > args.tol:
                    core_prob_fail += 1
                    n_core_fail += 1
                    local_core_fail = True
                    break

            # Leaf marginals identity
            local_leaf_fail = False
            if not local_core_fail:
                for l in leaves:
                    s = adj[l][0]
                    i = orig_to_core[s]
                    pred = 0.5 * (1.0 - probs_c[i])
                    if abs(probs_t[l] - pred) > args.tol:
                        leaf_prob_fail += 1
                        n_leaf_fail += 1
                        local_leaf_fail = True
                        break

            # Gap identity
            a_set = set(supports)
            w_occ = 0.0
            for v in core_to_orig:
                i = orig_to_core[v]
                if v in a_set:
                    w_occ += 0.5 * probs_c[i]
                else:
                    w_occ += probs_c[i]
            rhs_gap = (len(core_to_orig) / 3.0 - len(a_set) / 6.0) - w_occ
            if abs(gap_t - rhs_gap) > args.tol:
                gap_identity_fail += 1
                n_gap_fail += 1

            if (min_gap is None) or (gap_t < min_gap["gap_t"]):
                min_gap = {
                    "n": n0,
                    "g6": g6,
                    "gap_t": gap_t,
                    "mu_t": mu_t,
                    "core_size": len(core_to_orig),
                    "leaf_count": len(leaves),
                    "support_count": len(a_set),
                    "rhs_gap": rhs_gap,
                    "w_occ": w_occ,
                }
            if (min_gap_core is None) or (rhs_gap < min_gap_core["rhs_gap"]):
                min_gap_core = {
                    "n": n0,
                    "g6": g6,
                    "gap_t": gap_t,
                    "mu_t": mu_t,
                    "core_size": len(core_to_orig),
                    "leaf_count": len(leaves),
                    "support_count": len(a_set),
                    "rhs_gap": rhs_gap,
                    "w_occ": w_occ,
                }

        proc.wait()
        per_n[str(n)] = {
            "seen": n_seen,
            "considered": n_considered,
            "z_fail": n_z_fail,
            "core_prob_fail": n_core_fail,
            "leaf_prob_fail": n_leaf_fail,
            "gap_identity_fail": n_gap_fail,
            "wall_s": time.time() - t0,
        }
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
            f"z_fail={n_z_fail:5d} core_fail={n_core_fail:5d} "
            f"leaf_fail={n_leaf_fail:5d} gap_fail={n_gap_fail:5d} "
            f"({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} "
        f"z_fail={z_fail:,} core_fail={core_prob_fail:,} leaf_fail={leaf_prob_fail:,} "
        f"gap_fail={gap_identity_fail:,} wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    print(f"Minimum original gap witness: {min_gap}", flush=True)
    print(f"Minimum core-form gap witness: {min_gap_core}", flush=True)

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
                "z_fail": z_fail,
                "core_prob_fail": core_prob_fail,
                "leaf_prob_fail": leaf_prob_fail,
                "gap_identity_fail": gap_identity_fail,
                "min_gap": min_gap,
                "min_gap_core": min_gap_core,
                "wall_s": time.time() - t_all,
            },
            "per_n": per_n,
        }
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
