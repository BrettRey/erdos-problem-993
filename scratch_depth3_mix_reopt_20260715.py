#!/usr/bin/env python3
"""Per-size mix re-optimization for the depth-3 champion architecture.

Family shape fixed at (s=2, u=4, t=12): A super-gadgets (center with two
sub-spiders S(2^4), 19 vertices each) + B spiders S(2^12) (25 vertices)
on a common root. At each target size, scan the full mix range from
A = 0 (pure two-phase) to B = 0 (pure depth-3), coarse grid then local
refinement, and report the envelope deficit. Discriminates whether the
observed 0.09/sqrt(n) decay of the proportional ray is a property of the
class envelope or an artifact of mix detuning.
"""

import sys
import time

sys.path.insert(0, ".")

from scratch_depth3_valley_20260715 import tree_poly, size  # noqa: E402
from scripts.valley_search import valley_score  # noqa: E402

S, U, T = 2, 4, 12
SG_SIZE = 1 + S * (1 + 2 * U)   # 19
SP_SIZE = 1 + 2 * T             # 25


def eval_mix(A: int, n_target: int):
    B = max(0, round((n_target - 1 - A * SG_SIZE) / SP_SIZE))
    if A == 0 and B == 0:
        return None
    n = size(A, S, U, B, T)
    poly = tree_poly(A, S, U, B, T)
    vs = valley_score(poly)
    return A, B, n, len(poly) - 1, vs


def optimize(n_target: int, coarse: list[int], refine_steps: int = 2):
    tried = {}
    def run(A):
        if A in tried or A < 0:
            return
        r = eval_mix(A, n_target)
        if r:
            tried[A] = r
    for A in coarse:
        run(A)
    for _ in range(refine_steps):
        bestA = max(tried, key=lambda a: (tried[a][4]["window"],
                                          tried[a][4]["ratio"]))
        step = max(1, bestA // 8)
        for d in (-2 * step, -step, step, 2 * step):
            run(bestA + d)
    bestA = max(tried, key=lambda a: (tried[a][4]["window"],
                                      tried[a][4]["ratio"]))
    return tried[bestA], len(tried)


def main():
    print(f"{'n_tgt':>6} {'A*':>5} {'B*':>5} {'n':>6} {'alpha':>6} "
          f"{'V':>14} {'1-V':>11} {'n(1-V)':>9} {'sq(1-V)':>9} "
          f"{'b/alpha':>8} {'win':>5} {'evals':>5} {'sec':>7}")
    results = []
    plans = [
        (2473, [0, 4, 13, 26, 52, 90, 130]),
        (4945, [0, 8, 26, 52, 104, 180, 260]),
        (9889, [0, 16, 52, 104, 208, 360, 520]),
    ]
    for n_target, coarse in plans:
        t0 = time.time()
        (A, B, n, alpha, vs), evals = optimize(n_target, coarse)
        d = 1.0 - vs["ratio"]
        results.append((n, A, B, d))
        print(f"{n_target:>6} {A:>5} {B:>5} {n:>6} {alpha:>6} "
              f"{vs['ratio']:>14.10f} {d:>11.3e} {n*d:>9.4f} "
              f"{(n**0.5)*d:>9.4f} {vs['pos']/alpha:>8.4f} "
              f"{str(vs['window']):>5s} {evals:>5} "
              f"{time.time()-t0:>7.1f}", flush=True)
        if vs["witness"]:
            print(f"*** WITNESS *** A={A} B={B}")
            return

    # seeded short scan at ~19777 using the optimal-A trend
    if len(results) >= 2 and results[-1][3] > 0:
        ratio = results[-1][1] / max(1, results[-2][1])
        A_seed = max(0, round(results[-1][1] * ratio))
        t0 = time.time()
        tried = {}
        for A in sorted({0, A_seed // 2, A_seed, 3 * A_seed // 2}):
            r = eval_mix(A, 19777)
            if r:
                tried[A] = r
        bestA = max(tried, key=lambda a: (tried[a][4]["window"],
                                          tried[a][4]["ratio"]))
        A, B, n, alpha, vs = tried[bestA]
        d = 1.0 - vs["ratio"]
        print(f"{19777:>6} {A:>5} {B:>5} {n:>6} {alpha:>6} "
              f"{vs['ratio']:>14.10f} {d:>11.3e} {n*d:>9.4f} "
              f"{(n**0.5)*d:>9.4f} {vs['pos']/alpha:>8.4f} "
              f"{str(vs['window']):>5s} {len(tried):>5} "
              f"{time.time()-t0:>7.1f}", flush=True)


if __name__ == "__main__":
    main()
