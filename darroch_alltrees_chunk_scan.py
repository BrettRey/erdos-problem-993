#!/usr/bin/env python3
"""Chunked all-tree Darroch scan for a fixed n using geng res/mod splitting.

Checks at lambda=1:
  1) mode(T) <= ceil(mu(T))
  2) mode(T) in {floor(mu(T)), ceil(mu(T))}   (Darroch-style set-membership)

Here mu(T) is the mean independent-set size under uniform weighting:
  mu = I'(1)/I(1) = (sum_k k*i_k)/(sum_k i_k).

Run one chunk via geng's res/mod partition and write a JSON summary.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from typing import Any

from graph6 import parse_graph6
from indpoly import independence_poly


def mode_index_leftmost(poly: list[int]) -> int:
    return max(range(len(poly)), key=lambda i: poly[i])


def main() -> None:
    ap = argparse.ArgumentParser(description="Chunked all-tree Darroch scan.")
    ap.add_argument("--n", type=int, required=True, help="Tree order.")
    ap.add_argument("--res", type=int, default=0, help="Residue class for geng split.")
    ap.add_argument("--mod", type=int, default=1, help="Modulo for geng split.")
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument(
        "--progress-every",
        type=int,
        default=200000,
        help="Print progress every this many trees (0 disables).",
    )
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    if args.n < 1:
        raise ValueError("n must be >= 1")
    if args.mod < 1:
        raise ValueError("mod must be >= 1")
    if not (0 <= args.res < args.mod):
        raise ValueError("res must satisfy 0 <= res < mod")

    cmd = [args.geng, "-q", str(args.n), f"{args.n - 1}:{args.n - 1}", "-c"]
    if args.mod > 1:
        cmd.append(f"{args.res}/{args.mod}")

    t0 = time.time()
    print(
        f"Darroch chunk scan: n={args.n}, chunk={args.res}/{args.mod}",
        flush=True,
    )

    total = 0
    mode_ceil_fail = 0
    darroch_fail = 0

    min_mode_ceil_margin: int | None = None
    min_set_margin: int | None = None

    first_mode_ceil_fail: dict[str, Any] | None = None
    first_darroch_fail: dict[str, Any] | None = None

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    assert proc.stdout is not None

    for line in proc.stdout:
        n0, adj = parse_graph6(line)
        total += 1
        if n0 != args.n:
            raise RuntimeError(f"Unexpected n from geng: got {n0}, expected {args.n}")

        poly = independence_poly(n0, adj)
        m = mode_index_leftmost(poly)

        z = sum(poly)
        s = sum(i * c for i, c in enumerate(poly))
        floor_mu = s // z
        ceil_mu = (s + z - 1) // z

        mode_ceil_margin = ceil_mu - m
        if min_mode_ceil_margin is None or mode_ceil_margin < min_mode_ceil_margin:
            min_mode_ceil_margin = mode_ceil_margin

        set_margin = min(abs(m - floor_mu), abs(m - ceil_mu))
        if min_set_margin is None or set_margin < min_set_margin:
            min_set_margin = set_margin

        if m > ceil_mu:
            mode_ceil_fail += 1
            if first_mode_ceil_fail is None:
                first_mode_ceil_fail = {
                    "g6": line.decode("ascii").strip(),
                    "mode": m,
                    "floor_mu": floor_mu,
                    "ceil_mu": ceil_mu,
                    "mu_float": s / z,
                }

        if m != floor_mu and m != ceil_mu:
            darroch_fail += 1
            if first_darroch_fail is None:
                first_darroch_fail = {
                    "g6": line.decode("ascii").strip(),
                    "mode": m,
                    "floor_mu": floor_mu,
                    "ceil_mu": ceil_mu,
                    "mu_float": s / z,
                }

        if args.progress_every > 0 and total % args.progress_every == 0:
            print(
                f"  chunk {args.res}/{args.mod}: trees={total:,} "
                f"mode<=ceil fails={mode_ceil_fail:,} darroch fails={darroch_fail:,} "
                f"elapsed={time.time()-t0:.1f}s",
                flush=True,
            )

    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"geng exited with code {proc.returncode}")

    summary = {
        "n": args.n,
        "chunk": {"res": args.res, "mod": args.mod},
        "total_trees": total,
        "mode_ceil_fail": mode_ceil_fail,
        "darroch_fail": darroch_fail,
        "min_mode_ceil_margin": min_mode_ceil_margin,
        "min_set_margin": min_set_margin,
        "first_mode_ceil_fail": first_mode_ceil_fail,
        "first_darroch_fail": first_darroch_fail,
        "wall_s": time.time() - t0,
    }

    print(
        f"Done chunk {args.res}/{args.mod}: trees={total:,}, "
        f"mode<=ceil fails={mode_ceil_fail:,}, darroch fails={darroch_fail:,}, "
        f"wall={summary['wall_s']:.1f}s",
        flush=True,
    )

    payload = {"params": vars(args), "summary": summary}
    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
