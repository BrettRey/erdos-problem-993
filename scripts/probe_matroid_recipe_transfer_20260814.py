#!/usr/bin/env python3
"""Probe: does the Divoux-Larson-Lowen-Wang Rota-disproof recipe
(arXiv:2608.07342) transfer to tree independence polynomials?

Their construction: (1) severe log-concavity failure in a parametric family
(generalized theta graphs, violation ratio Theta(t)); (2) replication of the
failure by direct sums (convolution); (3) Whittle q-lift, which realizes the
geometric tilt a_k -> q^k a_k inside the matroid category and converts
LC failures into peaks.

This probe measures the tree-side analogue of each step on the
Bautista-Ramos alternating-break family TG_{m,t} (arXiv:2511.00334; builder
validated in probe_absorption_margin_20260811.py against the paper's exact
break sets) and on the two n=26 LC-failure trees (results/analysis_n26.json).

All arithmetic exact (int / Fraction). Emits
results/matroid_recipe_transfer_20260814.json.
"""
import json
import math
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import probe_absorption_margin_20260811 as probe


def lc_break_windows(a):
    """All internal LC breaks k with exact violation ratio and the exact
    tilt window (lo, hi): c in (lo, hi) iff the tilted sequence c^k a_k
    has a strict valley at k."""
    out = []
    for k in range(1, len(a) - 1):
        if a[k] > 0 and a[k - 1] * a[k + 1] > a[k] * a[k]:
            out.append({
                "k": k,
                "ratio": Fraction(a[k - 1] * a[k + 1], a[k] * a[k]),
                "lo": Fraction(a[k], a[k + 1]),
                "hi": Fraction(a[k - 1], a[k]),
            })
    return out


def count_peaks(a):
    """Peak count per arXiv:2608.07342 p.2 (runs collapse; sentinels -inf)."""
    runs = []
    j = 0
    while j < len(a):
        k = j
        while k + 1 < len(a) and a[k + 1] == a[j]:
            k += 1
        runs.append(a[j])
        j = k + 1
    return sum(
        1 for i, v in enumerate(runs)
        if (i == 0 or v > runs[i - 1]) and (i == len(runs) - 1 or v > runs[i + 1])
    )


def is_unimodal(a):
    return count_peaks(a) <= 1


def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return r


def main():
    report = {}

    # --- Step 1 analogue: severity along TG_{m,t} -------------------------
    grid = []
    for m in (2, 4, 8, 12):
        for t in (2, 8, 17, 32):
            tree = probe.TG_mt(m, t)
            F = tree.ipoly()
            alpha = len(F) - 1
            thr = math.ceil((2 * alpha - 1) / 3)  # Levit-Kadrawi prefix bound
            wins = lc_break_windows(F)
            grid.append({
                "m": m, "t": t, "n": tree.n, "alpha": alpha,
                "lk_threshold": thr,
                "breaks": [w["k"] for w in wins],
                "reflected_depths": [alpha - w["k"] for w in wins],
                "max_ratio": str(max((w["ratio"] for w in wins), default=Fraction(0))),
                "max_ratio_float": float(max((w["ratio"] for w in wins), default=Fraction(0))),
                "any_break_below_lk_threshold": any(w["k"] < thr for w in wins),
                "unimodal": is_unimodal(F),
            })
    report["tg_grid"] = grid
    deepest = max((d for row in grid for d in row["reflected_depths"]), default=0)
    report["tg_grid_summary"] = {
        "max_violation_ratio_float": max(r["max_ratio_float"] for r in grid),
        "deepest_reflected_break": deepest,
        "any_break_below_lk_threshold": any(r["any_break_below_lk_threshold"] for r in grid),
        "all_unimodal": all(r["unimodal"] for r in grid),
    }

    # Structural note (checked against the grid): TG break depths are the odd
    # offsets <= 2m-1, while the LK threshold sits at reflected depth
    # ~ alpha/3 = (t+1)m + O(1). 2m-1 < (t+1)m for every t >= 1, so within
    # this family the break band can never reach the prefix.
    report["tg_depth_bound_check"] = all(
        max(row["reflected_depths"], default=0) <= 2 * row["m"] - 1
        for row in grid
    )

    # --- Step 2 analogue: convolution powers of a marginal failure --------
    with open(os.path.join(ROOT, "results", "analysis_n26.json")) as fh:
        d26 = json.load(fh)
    worst = max(d26["lc_failures"], key=lambda f: f["lc_ratio"])
    P = worst["poly"]
    powers = []
    Q = [1]
    for p in range(1, 13):
        Q = pmul(Q, P)
        wins = lc_break_windows(Q)
        powers.append({
            "power": p,
            "break_count": len(wins),
            "max_ratio_float": float(max((w["ratio"] for w in wins), default=Fraction(0))),
            "unimodal": is_unimodal(Q),
        })
    report["n26_worst_powers"] = {
        "base_graph6": worst["graph6"],
        "base_ratio": worst["lc_ratio"],
        "rows": powers,
    }

    # --- Step 3 analogue: how much does a single tilt expose? -------------
    tilt = {}
    for (m, t) in [(8, 8), (12, 8)]:
        F = probe.TG_mt(m, t).ipoly()
        wins = lc_break_windows(F)
        lo_all = max(w["lo"] for w in wins)
        hi_all = min(w["hi"] for w in wins)
        best = (1, None)
        c_min = int(min(w["lo"] for w in wins))
        c_max = int(max(w["hi"] for w in wins)) + 1
        step = max(1, (c_max - c_min) // 600)
        for num in range(c_min, c_max, step):
            tilted = [F[i] * Fraction(num) ** i for i in range(len(F))]
            p = count_peaks(tilted)
            if p > best[0]:
                best = (p, num)
        tilt[f"TG_{m}_{t}"] = {
            "windows": [
                {"k": w["k"], "lo": float(w["lo"]), "hi": float(w["hi"])}
                for w in wins
            ],
            "all_windows_intersect": lo_all < hi_all,
            "best_single_tilt_peaks": best[0],
            "best_c": best[1],
        }
    report["single_tilt"] = tilt

    out = os.path.join(ROOT, "results", "matroid_recipe_transfer_20260814.json")
    with open(out, "w") as fh:
        json.dump(report, fh, indent=1)
    print(json.dumps(report["tg_grid_summary"], indent=1))
    print("tg_depth_bound_check:", report["tg_depth_bound_check"])
    print("n26 worst tree: powers 2..12 all LC:",
          all(r["break_count"] == 0 for r in powers[1:]))
    for k, v in tilt.items():
        print(k, "best single-tilt peaks:", v["best_single_tilt_peaks"],
              "at c =", v["best_c"],
              "| all windows intersect:", v["all_windows_intersect"])
    print("wrote", out)


if __name__ == "__main__":
    main()
