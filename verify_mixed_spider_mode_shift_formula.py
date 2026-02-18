#!/usr/bin/env python3
"""Check closed-form mode law and parity-step mode shift for mixed spiders.

Family:
  T_{k,j} = S(2^k, 1^j)
  I_{k,j}(x) = (1+2x)^k (1+x)^j + x(1+x)^k.

Claims checked:
  (F) mode(k,j) = floor((4k + 3j + 3)/6)
  (S) mode(k,j+2) = mode(k,j) + 1

Modes are leftmost modes (argmax on exact integer coefficients).
"""

from __future__ import annotations

import argparse
import json
import os
import time
from typing import Any

from conjecture_a_mixed_spider_combined_scan import (
    build_g_coeffs_shift_binom,
    build_h_coeffs_2k,
    mode_lambda_from_fg,
    next_f_times_1px,
)


def maybe_store(dst: list[dict[str, Any]], item: dict[str, Any], cap: int) -> None:
    if len(dst) < cap:
        dst.append(item)


def main() -> None:
    ap = argparse.ArgumentParser(description="Verify mixed-spider mode formula and +2 shift.")
    ap.add_argument("--k-min", type=int, default=1)
    ap.add_argument("--k-max", type=int, default=8000)
    ap.add_argument("--j-max", type=int, default=120)
    ap.add_argument("--store-cap", type=int, default=50)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    if args.k_min < 1:
        raise ValueError("k-min must be >= 1")
    if args.k_max < args.k_min:
        raise ValueError("k-max must be >= k-min")
    if args.j_max < 2:
        raise ValueError("j-max must be >= 2")

    t0 = time.time()
    print(
        f"Mode formula scan: k={args.k_min}..{args.k_max}, j=0..{args.j_max}",
        flush=True,
    )

    mode_formula_fail: list[dict[str, Any]] = []
    mode_shift_fail: list[dict[str, Any]] = []

    # First nontrivial-mode witness (for quick sanity in report).
    first_mode: dict[str, Any] | None = None

    total_mode_checks = 0
    total_shift_checks = 0

    for k in range(args.k_min, args.k_max + 1):
        h = build_h_coeffs_2k(k)
        g = build_g_coeffs_shift_binom(k)
        f = h[:]

        modes: list[int] = []
        for j in range(0, args.j_max + 3):
            m, _, _, _ = mode_lambda_from_fg(f, g)
            modes.append(int(m))

            if first_mode is None and m > 0:
                first_mode = {"k": k, "j": j, "mode": int(m)}

            if j <= args.j_max:
                total_mode_checks += 1
                mf = (4 * k + 3 * j + 3) // 6
                if int(m) != mf:
                    maybe_store(
                        mode_formula_fail,
                        {
                            "k": k,
                            "j": j,
                            "mode": int(m),
                            "mode_formula": mf,
                        },
                        args.store_cap,
                    )

            if j < args.j_max + 2:
                f = next_f_times_1px(f)

        for j in range(0, args.j_max + 1):
            total_shift_checks += 1
            if modes[j + 2] != modes[j] + 1:
                maybe_store(
                    mode_shift_fail,
                    {
                        "k": k,
                        "j": j,
                        "mode_j": modes[j],
                        "mode_j2": modes[j + 2],
                        "mode_step": modes[j + 2] - modes[j],
                    },
                    args.store_cap,
                )

        if (k - args.k_min) % 200 == 0 or k == args.k_max:
            print(
                f"k={k:5d} formula_fail={len(mode_formula_fail)} shift_fail={len(mode_shift_fail)}",
                flush=True,
            )

    payload = {
        "params": vars(args),
        "summary": {
            "mode_formula_checks": total_mode_checks,
            "mode_shift_checks": total_shift_checks,
            "mode_formula_fail_count": len(mode_formula_fail),
            "mode_shift_fail_count": len(mode_shift_fail),
            "first_mode": first_mode,
            "wall_s": time.time() - t0,
        },
        "fails": {
            "mode_formula": mode_formula_fail,
            "mode_shift": mode_shift_fail,
        },
    }

    print("-" * 88, flush=True)
    print(f"mode_formula_checks={total_mode_checks}", flush=True)
    print(f"mode_shift_checks={total_shift_checks}", flush=True)
    print(f"mode_formula_fail_count={len(mode_formula_fail)}", flush=True)
    print(f"mode_shift_fail_count={len(mode_shift_fail)}", flush=True)
    print(f"wall={payload['summary']['wall_s']:.2f}s", flush=True)

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f_out:
            json.dump(payload, f_out, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
