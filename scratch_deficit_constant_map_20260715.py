#!/usr/bin/env python3
"""Map the valley-deficit constant C(t) = n*(1-V) for two-type hybrids.

Family: a x S(2^2) + b x S(2^t) bouquets. For each phase-gap t and target
size n, optimize the mix fraction b/(a+b) by coarse-to-fine scan and
report the best deficit. If C(t) is bounded in t, the two-phase mechanism
cannot cross V=1 at any size; if C(t) grows, larger t at fixed n could.
"""

import sys
import time

sys.path.insert(0, ".")

from scripts.valley_scaling_probe import bouquet_total  # noqa: E402
from scripts.valley_search import valley_score  # noqa: E402


def eval_mix(a: int, b: int, t: int):
    g = [{2: 2}] * a + [{2: t}] * b
    n = 1 + a * 5 + b * (1 + 2 * t)
    poly = bouquet_total(g)
    vs = valley_score(poly)
    return n, vs


def best_over_mix(t: int, n_target: int):
    """Coarse-to-fine scan over b at roughly constant n."""
    best = None
    tried = set()
    # coarse: b from 1 up to max feasible
    b_max = max(2, (n_target - 6) // (1 + 2 * t))
    grid = sorted(set(max(1, round(b_max * f / 12)) for f in range(1, 13)))
    for rounds in range(3):
        for b in grid:
            if b in tried:
                continue
            tried.add(b)
            a = max(1, (n_target - 1 - b * (1 + 2 * t)) // 5)
            n, vs = eval_mix(a, b, t)
            key = (vs["window"], vs["ratio"])
            rec = (key, b, a, n, vs)
            if best is None or key > best[0]:
                best = rec
        # refine around current best b
        b0 = best[1]
        grid = [max(1, b0 + d) for d in (-3, -2, -1, 1, 2, 3)]
    return best


def main():
    n_targets = [640, 1280, 2560]
    ts = [6, 10, 15, 20, 30, 45, 60, 90]
    print(f"{'t':>4} {'n':>6} {'a':>5} {'b':>4} {'V':>14} {'n(1-V)':>9} "
          f"{'b/alpha':>8} {'win':>5} {'sec':>5}")
    for t in ts:
        for n_target in n_targets:
            t0 = time.time()
            key, b, a, n, vs = best_over_mix(t, n_target)
            alpha_pos = vs["pos"]
            # alpha from poly length is not returned; recompute cheap proxy
            print(f"{t:>4} {n:>6} {a:>5} {b:>4} {vs['ratio']:>14.10f} "
                  f"{n*(1-vs['ratio']):>9.4f} "
                  f"{'-' if vs['pos']<0 else vs['pos']:>8} "
                  f"{str(vs['window']):>5s} {time.time()-t0:>5.1f}",
                  flush=True)
            if vs["witness"]:
                print(f"*** WITNESS at t={t} a={a} b={b} ***")


if __name__ == "__main__":
    main()
