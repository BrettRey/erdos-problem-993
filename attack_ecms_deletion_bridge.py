#!/usr/bin/env python3
"""Attack ECMS via endpoint-deletion bridge statistics.

For an edge e = uv in tree T, define:
  m_T   = mode(I(T))
  m_C   = mode(I(T/e))
  m_u   = mode(I(T - u))
  m_v   = mode(I(T - v))
  m_*   = max(m_u, m_v)

Candidate bridge inequality:
  |m_C - m_*| <= 1.

This script checks the bridge inequality and reports shift distributions.
It also reports one-sided and full ECMS failures for comparison.
"""

from __future__ import annotations

import argparse
import subprocess
import time
from collections import Counter

from indpoly import independence_poly


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


def mode(poly: list[int]) -> int:
    peak = max(poly)
    for k, coeff in enumerate(poly):
        if coeff == peak:
            return k
    return 0


def induced_without(
    n: int, adj: list[list[int]], removed: int
) -> list[list[int]]:
    keep = [u for u in range(n) if u != removed]
    remap = {u: i for i, u in enumerate(keep)}
    keep_set = set(keep)
    out = [[] for _ in range(n - 1)]
    for u in keep:
        ru = remap[u]
        out[ru] = [remap[w] for w in adj[u] if w in keep_set]
    return out


def contract_edge(
    n: int, adj: list[list[int]], u: int, v: int
) -> list[list[int]]:
    remap = {}
    nxt = 1
    for i in range(n):
        if i == u or i == v:
            remap[i] = 0
        else:
            remap[i] = nxt
            nxt += 1

    out_sets = [set() for _ in range(n - 1)]
    for i in range(n):
        ri = remap[i]
        for j in adj[i]:
            rj = remap[j]
            if ri != rj:
                out_sets[ri].add(rj)
    return [list(s) for s in out_sets]


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Verify endpoint-deletion bridge inequality for contraction modes."
        )
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=18)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    args = parser.parse_args()

    total_edges = 0
    bridge_fail = 0
    one_sided_fail = 0
    ecms_fail = 0
    bridge_dist: Counter[int] = Counter()
    shift_dist: Counter[int] = Counter()
    example_bridge = None
    example_ecms = None

    t_all = time.time()
    print("ECMS deletion-bridge scan")
    print("Checking |mode(T/e) - max(mode(T-u), mode(T-v))| <= 1", flush=True)
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_trees = 0
        n_edges = 0
        n_bridge_fail = 0
        n_ecms_fail = 0

        assert proc.stdout is not None
        for line in proc.stdout:
            n_trees += 1
            nn, adj = parse_graph6(line)

            m_t = mode(independence_poly(nn, adj))

            del_modes = [0] * nn
            for v in range(nn):
                del_adj = induced_without(nn, adj, v)
                del_modes[v] = mode(independence_poly(nn - 1, del_adj))

            for u in range(nn):
                for v in adj[u]:
                    if u >= v:
                        continue
                    n_edges += 1
                    total_edges += 1

                    con_adj = contract_edge(nn, adj, u, v)
                    m_c = mode(independence_poly(nn - 1, con_adj))

                    m_star = max(del_modes[u], del_modes[v])
                    bshift = m_c - m_star
                    tshift = m_c - m_t

                    bridge_dist[bshift] += 1
                    shift_dist[tshift] += 1

                    if abs(bshift) > 1:
                        bridge_fail += 1
                        n_bridge_fail += 1
                        if example_bridge is None:
                            example_bridge = {
                                "n": nn,
                                "g6": line.decode("ascii").strip(),
                                "edge": (u, v),
                                "m_t": m_t,
                                "m_c": m_c,
                                "m_u": del_modes[u],
                                "m_v": del_modes[v],
                            }

                    if m_c < m_t - 1:
                        one_sided_fail += 1

                    if abs(tshift) > 1:
                        ecms_fail += 1
                        n_ecms_fail += 1
                        if example_ecms is None:
                            example_ecms = {
                                "n": nn,
                                "g6": line.decode("ascii").strip(),
                                "edge": (u, v),
                                "m_t": m_t,
                                "m_c": m_c,
                                "m_u": del_modes[u],
                                "m_v": del_modes[v],
                            }

        proc.wait()
        print(
            f"n={n:2d}: trees={n_trees:8d} edges={n_edges:10d} "
            f"bridge_fail={n_bridge_fail:4d} ecms_fail={n_ecms_fail:4d} "
            f"({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"Total edges: {total_edges:,}  bridge_fail={bridge_fail}  "
        f"one_sided_fail={one_sided_fail}  ecms_fail={ecms_fail}  "
        f"wall={time.time() - t_all:.1f}s",
        flush=True,
    )

    print("Bridge shift distribution m_C - m_*:", flush=True)
    for k in sorted(bridge_dist):
        cnt = bridge_dist[k]
        pct = 100.0 * cnt / total_edges if total_edges else 0.0
        print(f"  {k:+d}: {cnt:12d} ({pct:6.3f}%)", flush=True)

    print("ECMS shift distribution m_C - m_T:", flush=True)
    for k in sorted(shift_dist):
        cnt = shift_dist[k]
        pct = 100.0 * cnt / total_edges if total_edges else 0.0
        print(f"  {k:+d}: {cnt:12d} ({pct:6.3f}%)", flush=True)

    if example_bridge is not None:
        print("First bridge violation example:", example_bridge, flush=True)
    if example_ecms is not None:
        print("First ECMS violation example:", example_ecms, flush=True)


if __name__ == "__main__":
    main()
