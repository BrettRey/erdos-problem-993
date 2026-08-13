#!/usr/bin/env python3
"""Toy kill-test for the Beurling-majorant difference-test lead.

See notes/beurling_majorant_difference_test_lead_2026-08-12.md, section 4,
step 2. Question: on the known hard trees, is there ANY band width w such
that (a) the angular band |theta| < w around the positive axis on the
saddle circle |z| = x*(k) is free of non-real zeros of I(T;.), and
(b) band-limited Beurling--Selberg majorants of width w resolve the
valley scale (single coefficients)?

Accounting used (recorded in the note):
- Selberg majorant/minorant of an interval at bandwidth w has mass slack
  ~ 2*pi/w in coefficient-count units, so resolving a single coefficient
  needs w > 2*pi.
- Frequencies live on the torus [-pi, pi], so usable w <= pi. Hence the
  PURE lane needs w >= 2*pi > pi: structurally impossible. The empirical
  content of this test is the actual margin: w_avail (widest zero-free
  angular band at each modulus-activity threshold) and the implied
  interval resolution Delta_k = 2*pi / w_avail, which is what a hybrid
  assembly (majorant intervals + a local step) would have to cover.

Certified arithmetic: exact integer coefficients (indpoly DP or the
compressed word recurrence), exact rational saddle bisection (integer
sign tests only), python-flint Arb root isolation with the conservative
non-real test (imaginary ball excludes zero). No float64 root-finding
(DECISIONS.md 2026-07-16).

Run with a python-flint venv:
  <venv>/bin/python scratch_beurling_band_killtest_20260812.py
"""

import json
import math
import sys
from fractions import Fraction

sys.path.insert(0, ".")

from flint import fmpz_poly

from indpoly import independence_poly
from scripts.stress_literature_root_families import (
    balanced_tree_from_word,
    compressed_state_from_word,
    phase_truncated_words,
    poly_add,
)

OUT = "results/beurling_band_killtest_20260812.json"
ETAS = [0.05, 0.1, 0.25, 0.5, 1.0, None]  # None = all non-real zeros
K_FRACS = [Fraction(1, 3), Fraction(1, 2), Fraction(2, 3)]


def edges_to_adj(n, edges):
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    if sum(len(a) for a in adj) != 2 * (n - 1):
        raise ValueError("not a tree: wrong edge count")
    seen = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                stack.append(v)
    if len(seen) != n:
        raise ValueError("not a tree: disconnected")
    return adj


def path_adj(n):
    return edges_to_adj(n, [[i, i + 1] for i in range(n - 1)])


def eval_frac(coeffs, x):
    """Exact Horner evaluation of an integer polynomial at Fraction x."""
    acc = Fraction(0)
    for c in reversed(coeffs):
        acc = acc * x + c
    return acc


def saddle(coeffs, k, iters=70):
    """Exact-bisection solution of x I'(x) = k I(x) on x > 0.

    Sign tests are exact (rational arithmetic on integer coefficients).
    Returns the midpoint Fraction of the final bracket.
    """
    scoef = [(j - k) * c for j, c in enumerate(coeffs)]  # x I' - k I

    def s(x):
        return eval_frac(scoef, x)

    lo, hi = Fraction(0), Fraction(1)
    while s(hi) <= 0:
        hi *= 2
    for _ in range(iters):
        mid = (lo + hi) / 2
        if s(mid) <= 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def tilt_stats(coeffs, x):
    """(mean, variance) of the x-tilted coefficient distribution, exact."""
    a = eval_frac(coeffs, x)
    b = eval_frac([j * c for j, c in enumerate(coeffs)], x)
    c2 = eval_frac([j * j * c for j, c in enumerate(coeffs)], x)
    mu = b / a
    var = c2 / a - mu * mu
    return float(mu), float(var)


def certified_root_data(coeffs):
    """[(modulus, |arg|, certified_nonreal)] from Arb isolation.

    Non-real only when the imaginary ball excludes zero. Roots not
    certified non-real are real (trees have no nonnegative real roots,
    so their angle is pi). Returns also the max ball radius seen.
    """
    roots = fmpz_poly(coeffs).complex_roots()
    data = []
    max_rad = 0.0
    for z, _mult in roots:
        re = float(z.real.mid())
        im = float(z.imag.mid())
        rad = max(float(z.real.rad()), float(z.imag.rad()))
        max_rad = max(max_rad, rad)
        nonreal = float(z.imag.rad()) < abs(float(z.imag.mid()))
        mod = math.hypot(re, im)
        ang = abs(math.atan2(im, re)) if nonreal else math.pi
        data.append((mod, ang, nonreal))
    return data, max_rad


def band_widths(root_data, xstar):
    """w_avail per eta: min |arg| over non-real zeros within the
    log-modulus band |log(r/x*)| <= eta. pi if the band is empty."""
    out = {}
    for eta in ETAS:
        w = math.pi
        for mod, ang, nonreal in root_data:
            if not nonreal or mod <= 0:
                continue
            if eta is None or abs(math.log(mod / xstar)) <= eta:
                w = min(w, ang)
        out["inf" if eta is None else str(eta)] = w
    return out


def analyze(label, coeffs):
    alpha = len(coeffs) - 1
    root_data, max_rad = certified_root_data(coeffs)
    n_nonreal = sum(1 for _, _, nr in root_data if nr)
    ks = sorted({max(1, min(alpha - 1, int(alpha * f))) for f in K_FRACS})
    rows = []
    for k in ks:
        xs = saddle(coeffs, k)
        mu, var = tilt_stats(coeffs, xs)
        xf = float(xs)
        w = band_widths(root_data, xf)
        rows.append({
            "k": k,
            "x_star": xf,
            "tilt_mean": mu,
            "tilt_var": var,
            "w_avail": w,
            "resolution": {e: (2 * math.pi / v if v > 0 else None)
                           for e, v in w.items()},
        })
        print(f"  k={k:4d}  x*={xf:.5f}  var={var:8.2f}  "
              + "  ".join(f"w[{e}]={v:.4f}" for e, v in w.items()))
    return {
        "label": label,
        "alpha": alpha,
        "n_roots": sum(1 for _ in root_data),
        "n_nonreal_certified": n_nonreal,
        "max_ball_radius": max_rad,
        "per_k": rows,
    }


def main():
    cases = []

    # Control: path (claw-free, real-rooted): best possible case.
    print("== control-path-101")
    cases.append(analyze("control-path-101",
                         independence_poly(101, path_adj(101))))

    # Adversarial champions from the July pareto root frontier.
    frontier = json.load(open("results/pareto_root_frontier_20260716.json"))
    for band, rec in sorted(frontier.items()):
        n = rec["n"]
        label = f"pareto-band{band}-n{n}"
        print(f"== {label}")
        adj = edges_to_adj(n, rec["edges"])
        cases.append(analyze(label, independence_poly(n, adj)))

    # Jerrum--Patel phase-truncated family (positive-axis accumulation).
    for label, word in phase_truncated_words(1, 5):
        exc, inc = compressed_state_from_word(word)
        coeffs = poly_add(exc, inc)
        n = len(balanced_tree_from_word(word))
        full = f"{label}-n{n}"
        print(f"== {full}")
        cases.append(analyze(full, coeffs))

    # Unsplit binary control at height 6 (remote positive-axis regime).
    word = (2,) * 6
    exc, inc = compressed_state_from_word(word)
    print("== binary-h6-n127")
    cases.append(analyze("binary-h6-n127", poly_add(exc, inc)))

    summary = {
        "purpose": ("Beurling-majorant band kill-test: w_avail vs the "
                    "w > 2*pi needed for single-coefficient resolution "
                    "(torus cap pi)."),
        "structural_fact": ("usable bandwidth <= pi < 2*pi needed: the "
                            "pure-majorant lane cannot resolve a width-1 "
                            "valley on ANY tree; numbers below give the "
                            "hybrid margin Delta_k = 2*pi/w_avail."),
        "cases": cases,
    }
    with open(OUT, "w") as fh:
        json.dump(summary, fh, indent=1)
    print(f"\nwrote {OUT}")

    # Console verdict table: worst (largest) Delta_k at eta=0.5 per case.
    print("\ncase, worst resolution Delta_k over mid-band k (eta=0.5):")
    for c in cases:
        worst = max(r["resolution"]["0.5"] for r in c["per_k"])
        print(f"  {c['label']:44s}  Delta_k = {worst:8.2f}")


if __name__ == "__main__":
    main()
