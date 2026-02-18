#!/usr/bin/env python3
"""Scan private-neighbor properties for heavy subsets in trees.

Definitions (at fugacity lambda=1):
  H = {v : P(v) > 1/3}
  For non-empty S subseteq H, a vertex h in S has a private neighbor if
    exists u ~ h with N(u) cap S = {h}.

Checks:
  strong(S): every h in S has a private neighbor
  weak(S):   at least one h in S has a private neighbor
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time

from conjecture_a_hall_subset_scan import hard_core_probs, is_dleaf_le_1, parse_graph6


def first_private_neighbor_violations(
    adj: list[list[int]],
    probs: list[float],
    tol: float,
) -> tuple[int | None, int | None, list[int], list[bool] | None, list[bool] | None]:
    """Return first strong/weak violating subsets as bitmasks over H indices."""
    one_third = 1.0 / 3.0
    h_nodes = [v for v, pv in enumerate(probs) if pv > one_third + tol]
    m = len(h_nodes)
    if m == 0:
        return None, None, h_nodes, None, None

    h_idx = {h: i for i, h in enumerate(h_nodes)}

    nbr_masks: list[list[int]] = [[] for _ in range(m)]
    for i, h in enumerate(h_nodes):
        for u in adj[h]:
            mask = 0
            for hh in adj[u]:
                j = h_idx.get(hh)
                if j is not None:
                    mask |= 1 << j
            nbr_masks[i].append(mask)

    strong_bad = None
    weak_bad = None
    strong_flags = None
    weak_flags = None

    for s in range(1, 1 << m):
        has_private = [False] * m
        any_private = False
        all_private = True
        rem = s
        while rem:
            lsb = rem & -rem
            i = lsb.bit_length() - 1
            rem ^= lsb
            bit_i = 1 << i

            priv = False
            for mu in nbr_masks[i]:
                if (s & mu) == bit_i:
                    priv = True
                    break
            has_private[i] = priv
            any_private = any_private or priv
            all_private = all_private and priv

        if strong_bad is None and not all_private:
            strong_bad = s
            strong_flags = has_private
        if weak_bad is None and not any_private:
            weak_bad = s
            weak_flags = has_private
        if strong_bad is not None and weak_bad is not None:
            break

    return strong_bad, weak_bad, h_nodes, strong_flags, weak_flags


def decode_subset(mask: int, h_nodes: list[int]) -> list[int]:
    out: list[int] = []
    rem = mask
    while rem:
        lsb = rem & -rem
        i = lsb.bit_length() - 1
        rem ^= lsb
        out.append(h_nodes[i])
    return out


def decode_flags(mask: int, h_nodes: list[int], flags: list[bool] | None) -> dict[int, bool]:
    if flags is None:
        return {}
    out: dict[int, bool] = {}
    rem = mask
    while rem:
        lsb = rem & -rem
        i = lsb.bit_length() - 1
        rem ^= lsb
        out[h_nodes[i]] = flags[i]
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Private-neighbor property scan on heavy subsets."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=23)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument(
        "--all-trees",
        action="store_true",
        help="Do not filter by d_leaf<=1; scan all trees.",
    )
    parser.add_argument(
        "--stop-on-first",
        action="store_true",
        help="Stop once both strong and weak counterexamples are found.",
    )
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    total_seen = 0
    total_considered = 0
    strong_fail_trees = 0
    weak_fail_trees = 0

    first_strong = None
    first_weak = None

    t_all = time.time()
    print("Private-neighbor scan on heavy subsets")
    print("Properties:")
    print("  strong: every h in S has a private neighbor")
    print("  weak:   at least one h in S has a private neighbor")
    print("-" * 80, flush=True)

    done = False
    for n in range(args.min_n, args.max_n + 1):
        if done:
            break
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_strong_fail = 0
        n_weak_fail = 0

        assert proc.stdout is not None
        for line in proc.stdout:
            n0, adj = parse_graph6(line)
            n_seen += 1
            total_seen += 1
            if (not args.all_trees) and (not is_dleaf_le_1(n0, adj)):
                continue

            n_considered += 1
            total_considered += 1

            probs = hard_core_probs(n0, adj)
            strong_bad, weak_bad, h_nodes, strong_flags, weak_flags = (
                first_private_neighbor_violations(adj, probs, args.tol)
            )

            g6 = line.decode("ascii").strip()
            if strong_bad is not None:
                n_strong_fail += 1
                strong_fail_trees += 1
                if first_strong is None:
                    first_strong = {
                        "n": n0,
                        "g6": g6,
                        "H": h_nodes,
                        "S": decode_subset(strong_bad, h_nodes),
                        "private_map_S": decode_flags(strong_bad, h_nodes, strong_flags),
                    }
            if weak_bad is not None:
                n_weak_fail += 1
                weak_fail_trees += 1
                if first_weak is None:
                    first_weak = {
                        "n": n0,
                        "g6": g6,
                        "H": h_nodes,
                        "S": decode_subset(weak_bad, h_nodes),
                        "private_map_S": decode_flags(weak_bad, h_nodes, weak_flags),
                    }

            if args.stop_on_first and first_strong is not None and first_weak is not None:
                done = True
                break

        proc.wait()
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
            f"strong_fail_trees={n_strong_fail:7d} weak_fail_trees={n_weak_fail:7d} "
            f"({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} "
        f"strong_fail_trees={strong_fail_trees:,} weak_fail_trees={weak_fail_trees:,} "
        f"wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    print(f"First strong counterexample: {first_strong}", flush=True)
    print(f"First weak counterexample:   {first_weak}", flush=True)

    if args.out:
        os_dir = args.out.rsplit("/", 1)[0] if "/" in args.out else ""
        if os_dir:
            os.makedirs(os_dir, exist_ok=True)
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
                "strong_fail_trees": strong_fail_trees,
                "weak_fail_trees": weak_fail_trees,
                "first_strong": first_strong,
                "first_weak": first_weak,
            },
        }
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
