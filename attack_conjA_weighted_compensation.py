#!/usr/bin/env python3
"""Attack Conjecture A via weighted heavy-neighborhood compensation.

For a tree T with hard-core occupation probabilities P(v) at lambda = 1, define:

  H = {v : P(v) > 1/3}
  N(H) = (union of neighbors of H) \\ H

Candidate inequality (WHNC):

  sum_{h in H} (P(h) - 1/3) <= sum_{u in N(H)} (1/3 - P(u)).

If this holds, then the heavy-vertex excess is fully compensated by deficit on
adjacent non-heavy vertices. Since all vertices in V \\ H have P(v) <= 1/3,
this implies:

  sum_v P(v) <= n/3.

This script verifies WHNC exactly (Fraction arithmetic) on d_leaf <= 1 trees.
"""

from __future__ import annotations

import argparse
import subprocess
import time
from fractions import Fraction


def parse_graph6(line: bytes) -> tuple[int, list[list[int]]]:
    s = line.decode("ascii").strip()
    data = [ord(c) - 63 for c in s]
    n = data[0]
    idx = 1
    adj = [[] for _ in range(n)]
    bit = 5
    word = data[idx] if idx < len(data) else 0
    for j in range(n):
        for i in range(j):
            if word & (1 << bit):
                adj[i].append(j)
                adj[j].append(i)
            if bit == 0:
                bit = 5
                idx += 1
                word = data[idx] if idx < len(data) else 0
            else:
                bit -= 1
    return n, adj


def is_dleaf_le_1(n: int, adj: list[list[int]]) -> bool:
    leaves = {v for v in range(n) if len(adj[v]) == 1}
    for v in range(n):
        leaf_children = sum(1 for u in adj[v] if u in leaves)
        if leaf_children > 1:
            return False
    return True


def hard_core_probs_fraction(n: int, adj: list[list[int]]) -> list[Fraction]:
    if n == 1:
        return [Fraction(1, 2)]

    # Root the tree at 0.
    parent = [-1] * n
    children = [[] for _ in range(n)]
    order: list[int] = []
    queue = [0]
    seen = [False] * n
    seen[0] = True
    q_idx = 0

    while q_idx < len(queue):
        v = queue[q_idx]
        q_idx += 1
        order.append(v)
        for u in adj[v]:
            if not seen[u]:
                seen[u] = True
                parent[u] = v
                children[v].append(u)
                queue.append(u)

    # Bottom-up messages R(child -> parent).
    r_up = [Fraction(0)] * n
    for v in reversed(order):
        prod = Fraction(1, 1)
        for c in children[v]:
            prod *= 1 + r_up[c]
        r_up[v] = Fraction(1, 1) / prod

    # Top-down messages R(parent -> child).
    r_down = [Fraction(0)] * n
    for v in order:
        for c in children[v]:
            prod = Fraction(1, 1)
            if parent[v] != -1:
                prod *= 1 + r_down[v]
            for s in children[v]:
                if s != c:
                    prod *= 1 + r_up[s]
            r_down[c] = Fraction(1, 1) / prod

    probs = [Fraction(0)] * n
    for v in range(n):
        prod = Fraction(1, 1)
        if parent[v] != -1:
            prod *= 1 + r_down[v]
        for c in children[v]:
            prod *= 1 + r_up[c]
        r_full = Fraction(1, 1) / prod
        probs[v] = r_full / (1 + r_full)
    return probs


def whnc_margin(
    probs: list[Fraction], adj: list[list[int]]
) -> tuple[Fraction, Fraction, Fraction, set[int], set[int]]:
    one_third = Fraction(1, 3)
    heavy = {v for v, p in enumerate(probs) if p > one_third}
    nh: set[int] = set()
    for h in heavy:
        nh.update(adj[h])
    nh -= heavy

    excess = sum((probs[h] - one_third) for h in heavy)
    deficit = sum((one_third - probs[u]) for u in nh)
    margin = deficit - excess
    return margin, excess, deficit, heavy, nh


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Verify weighted heavy-neighborhood compensation on d_leaf <= 1 trees."
        )
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=23)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    args = parser.parse_args()

    total_dleaf = 0
    total_fail = 0
    total_equal = 0
    global_min_margin = Fraction(10**9, 1)
    worst_case = None

    t_all = time.time()
    print("WHNC verification on d_leaf <= 1 trees")
    print(
        "Checking: sum_H(P-1/3) <= sum_N(H)(1/3-P), H={v:P(v)>1/3}",
        flush=True,
    )
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        count_dleaf = 0
        fail_n = 0
        equal_n = 0
        min_margin_n: Fraction | None = None

        assert proc.stdout is not None
        for line in proc.stdout:
            nn, adj = parse_graph6(line)
            if not is_dleaf_le_1(nn, adj):
                continue

            count_dleaf += 1
            total_dleaf += 1

            probs = hard_core_probs_fraction(nn, adj)
            margin, excess, deficit, heavy, nh = whnc_margin(probs, adj)

            if min_margin_n is None or margin < min_margin_n:
                min_margin_n = margin
            if margin < global_min_margin:
                global_min_margin = margin
                worst_case = {
                    "n": nn,
                    "g6": line.decode("ascii").strip(),
                    "margin": margin,
                    "excess": excess,
                    "deficit": deficit,
                    "heavy": sorted(heavy),
                    "nh": sorted(nh),
                    "probs": probs,
                }

            if margin < 0:
                total_fail += 1
                fail_n += 1
            elif margin == 0:
                total_equal += 1
                equal_n += 1

        proc.wait()
        if min_margin_n is None:
            margin_str = "N/A"
        else:
            margin_str = f"{float(min_margin_n):.12f}"

        print(
            f"n={n:2d}: d_leaf<=1={count_dleaf:8d}  fail={fail_n:4d}  "
            f"eq={equal_n:4d}  min_margin={margin_str:>14}  "
            f"({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL d_leaf<=1 trees checked: {total_dleaf:,}  "
        f"fail={total_fail}  eq={total_equal}  wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    print(f"Global minimum margin: {float(global_min_margin):.12f}", flush=True)

    if worst_case is not None:
        print("Worst-case example:", flush=True)
        print(f"  n={worst_case['n']} g6={worst_case['g6']}", flush=True)
        print(
            f"  |H|={len(worst_case['heavy'])} |N(H)|={len(worst_case['nh'])}",
            flush=True,
        )
        print(
            "  excess="
            f"{float(worst_case['excess']):.12f} "
            "deficit="
            f"{float(worst_case['deficit']):.12f} "
            "margin="
            f"{float(worst_case['margin']):.12f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
