#!/usr/bin/env python3
"""
Enumerate degree-2-free, leaf-light trees from finite core templates.

Model:
  - choose an unlabeled core tree K on b <= b0 vertices,
  - assign a leaf load l_u per core vertex u with
        max(0, 3 - deg_K(u)) <= l_u <= leaf_cap,
  - attach l_u leaves to u,
  - test independence polynomial unimodality.

Optional: exact dedup via nauty labelg (if installed).
"""

from __future__ import annotations

import argparse
import itertools
import json
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Iterable, List, Tuple

from graph6 import parse_graph6
from indpoly import independence_poly, is_unimodal
from trees import trees


def encode_graph6_small(adj: List[List[int]]) -> str:
    """Encode adjacency list as graph6 for n < 63."""
    n = len(adj)
    if n >= 63:
        raise ValueError(f"n={n} not supported by small graph6 encoder")
    aset = [set(nei) for nei in adj]
    bits: List[int] = []
    for j in range(1, n):
        sj = aset[j]
        for i in range(j):
            bits.append(1 if i in sj else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [chr(n + 63)]
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k : k + 6]:
            v = (v << 1) | b
        out.append(chr(v + 63))
    return "".join(out)


def build_tree_from_core(core_adj: List[List[int]], loads: Tuple[int, ...]) -> List[List[int]]:
    b = len(core_adj)
    n = b + sum(loads)
    adj: List[List[int]] = [nbrs[:] for nbrs in core_adj]
    adj.extend([] for _ in range(n - b))
    nxt = b
    for u, l in enumerate(loads):
        for _ in range(l):
            adj[u].append(nxt)
            adj[nxt].append(u)
            nxt += 1
    return adj


def candidate_g6_stream(b0: int, leaf_cap: int, core_backend: str) -> Iterable[str]:
    for b in range(1, b0 + 1):
        for _, core in trees(b, backend=core_backend):
            ranges: List[range] = []
            ok = True
            for u in range(b):
                lb = max(0, 3 - len(core[u]))
                if lb > leaf_cap:
                    ok = False
                    break
                ranges.append(range(lb, leaf_cap + 1))
            if not ok:
                continue
            for loads in itertools.product(*ranges):
                adj = build_tree_from_core(core, loads)
                yield encode_graph6_small(adj)


def canonicalize_with_labelg(lines: List[str]) -> List[str]:
    if shutil.which("labelg") is None:
        raise RuntimeError("labelg not found; install nauty or use --dedup none")
    with tempfile.TemporaryDirectory() as td:
        inp = Path(td) / "in.g6"
        out = Path(td) / "out.g6"
        inp.write_text("\n".join(lines) + "\n", encoding="ascii")
        subprocess.run(["labelg", "-q", str(inp), str(out)], check=True)
        return out.read_text(encoding="ascii").splitlines()


def evaluate_graphs(g6_iter: Iterable[str], stop_on_fail: bool) -> dict:
    seen = set()
    total = 0
    uniq = 0
    non_uni = 0
    first_bad = None

    for g6 in g6_iter:
        total += 1
        if g6 in seen:
            continue
        seen.add(g6)
        uniq += 1

        n, adj = parse_graph6(g6.encode("ascii"))
        poly = independence_poly(n, adj)
        if not is_unimodal(poly):
            non_uni += 1
            if first_bad is None:
                first_bad = {"n": n, "graph6": g6, "poly": poly}
            if stop_on_fail:
                break

    return {
        "generated": total,
        "unique_checked": uniq,
        "non_unimodal": non_uni,
        "first_counterexample": first_bad,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--b0", type=int, required=True, help="max core size")
    ap.add_argument("--leaf-cap", type=int, required=True, help="max leaves per core vertex")
    ap.add_argument(
        "--core-backend",
        default="networkx",
        choices=["networkx", "geng", "auto"],
        help="backend for core enumeration",
    )
    ap.add_argument("--dedup", default="none", choices=["none", "labelg"])
    ap.add_argument("--stop-on-fail", action="store_true")
    ap.add_argument("--out", default="", help="write JSON results here")
    args = ap.parse_args()

    t0 = time.time()
    g6_lines = list(candidate_g6_stream(args.b0, args.leaf_cap, args.core_backend))
    t1 = time.time()

    if args.dedup == "labelg":
        g6_lines = canonicalize_with_labelg(g6_lines)
    t2 = time.time()

    result = evaluate_graphs(g6_lines, args.stop_on_fail)
    t3 = time.time()

    payload = {
        "b0": args.b0,
        "leaf_cap": args.leaf_cap,
        "core_backend": args.core_backend,
        "dedup": args.dedup,
        **result,
        "time_generate_s": round(t1 - t0, 3),
        "time_canonicalize_s": round(t2 - t1, 3),
        "time_check_s": round(t3 - t2, 3),
        "time_total_s": round(t3 - t0, 3),
    }

    if args.out:
        Path(args.out).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    else:
        print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
