#!/usr/bin/env python3
"""Scan private-neighbor property on non-leaf heavy subsets.

Setup:
  H = {v : P(v) > 1/3}
  H_nl = {h in H : deg(h) >= 2}

For non-empty S subseteq H_nl with |S|>=2, check whether some h in S has a
private neighbor relative to S:
  exists u ~ h with N(u) cap S = {h}.

If true for a given h, then by the edge surplus bound
  supply(u) = 1/3 - P(u) > P(h) - 1/3 = demand(h),
so M(h,S) = supply(N(h)\\N(S\\{h})) - demand(h) > 0.
This is the key non-leaf case for singleton-argmin.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time

from conjecture_a_hall_subset_scan import hard_core_probs, is_dleaf_le_1, parse_graph6


def decode_subset(mask: int, nodes: list[int]) -> list[int]:
    out: list[int] = []
    rem = mask
    while rem:
        lsb = rem & -rem
        i = lsb.bit_length() - 1
        rem ^= lsb
        out.append(nodes[i])
    return out


def first_nonleaf_private_failure(
    adj: list[list[int]],
    probs: list[float],
    tol: float,
) -> tuple[bool, list[int] | None, list[int] | None]:
    """Return (ok, H_nonleaf, failing_subset) for one tree."""
    one_third = 1.0 / 3.0
    heavy = [v for v, pv in enumerate(probs) if pv > one_third + tol]
    h_nonleaf = [h for h in heavy if len(adj[h]) > 1]
    k = len(h_nonleaf)
    if k <= 1:
        return True, h_nonleaf, None

    h_set = set(heavy)
    u_nodes = sorted({u for h in heavy for u in adj[h] if u not in h_set})
    u_idx = {u: i for i, u in enumerate(u_nodes)}

    masks: list[int] = []
    for h in h_nonleaf:
        m = 0
        for u in adj[h]:
            j = u_idx.get(u)
            if j is not None:
                m |= 1 << j
        masks.append(m)

    k_pow = 1 << k
    union_mask = [0] * k_pow
    for s in range(1, k_pow):
        lsb = s & -s
        i = lsb.bit_length() - 1
        union_mask[s] = union_mask[s ^ lsb] | masks[i]

    for s in range(1, k_pow):
        if s & (s - 1) == 0:
            continue
        all_no_private = True
        rem = s
        while rem:
            lsb = rem & -rem
            i = lsb.bit_length() - 1
            rem ^= lsb
            priv = masks[i] & ~union_mask[s ^ lsb]
            if priv != 0:
                all_no_private = False
                break
        if all_no_private:
            return False, h_nonleaf, decode_subset(s, h_nonleaf)

    return True, h_nonleaf, None


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Non-leaf heavy private-neighbor scan for singleton-argmin route."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=23)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    total_seen = 0
    total_considered = 0
    bad_trees = 0
    first_bad = None

    per_n: dict[str, dict] = {}

    t_all = time.time()
    print("Non-leaf heavy private-neighbor scan")
    print("Checking: every S subseteq H_nonleaf, |S|>=2 has some private-neighbor vertex")
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_bad = 0

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
            ok, h_nonleaf, failing_s = first_nonleaf_private_failure(adj, probs, args.tol)
            if not ok:
                n_bad += 1
                bad_trees += 1
                if first_bad is None:
                    first_bad = {
                        "n": n0,
                        "g6": line.decode("ascii").strip(),
                        "H_nonleaf": h_nonleaf,
                        "failing_subset": failing_s,
                    }

        proc.wait()
        per_n[str(n)] = {
            "seen": n_seen,
            "considered": n_considered,
            "bad_trees": n_bad,
        }
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
            f"bad_trees={n_bad:7d} ({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} "
        f"bad_trees={bad_trees:,} wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    print(f"First bad witness: {first_bad}", flush=True)

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
                "bad_trees": bad_trees,
                "first_bad": first_bad,
                "wall_s": time.time() - t_all,
            },
            "per_n": per_n,
        }
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
