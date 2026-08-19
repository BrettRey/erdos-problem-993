#!/usr/bin/env python3
"""Exact counts of unlabeled trees by matching number; size the depth-3 window.

Purpose (lane decision input, 2026-08-19): the dist-3 frontier after the
n <= 32 LC census is depth >= 3, alpha in {17,18,19}, n in [33,38]
(notes/break_depth_lane_2026-08-14.md, consequence 4). Trees are bipartite,
so by Konig alpha(T) = n - mu(T) with mu the matching number. Counting
unlabeled trees by mu therefore sizes the alpha-filtered slice EXACTLY,
without enumerating a single tree, and settles whether the window is a
laptop job, a Modal job, or needs a matching-constrained generator.

Method: Otter-style dissymmetry with a bivariate weight x^n y^mu.

Rooted trees carry a two-value state s in {0,1}: s = 1 iff some maximum
matching of the rooted tree leaves the root exposed. Classical recurrence
(children T_1..T_k of the root):
  - if every child has s_i = 0, the root can stay exposed:
        mu = sum mu_i,      s = 1
  - if some child has s_j = 1, every maximum matching uses the root
    (match it to such a child):
        mu = sum mu_i + 1,  s = 0
Single vertex: mu = 0, s = 1.

Let R1 / R0 be the generating functions of rooted trees with s = 1 / 0.
With MSET the multiset construction,
    R1 = x * MSET(R0)                (all children from R0)
    R0 = x * y * (MSET(R0+R1) - MSET(R0))   (at least one R1 child)

Unrooted trees via Otter t = r - (r^2 - r(x^2))/2, with each pairing term
carrying the edge-join weight: joining rooted trees A, B by a root-root
edge gives mu = mu_A + mu_B + [s_A = 1 and s_B = 1] (use the new edge iff
both roots can be left exposed). So
    P[n][m] (ordered pairs) uses y-bump only on R1 x R1,
    D[n][m] (A paired with itself) puts A at (2 n_A, 2 mu_A + s_A),
    t = r - (P - D)/2.

Self-tests:
  1. Totals sum_mu t[n][mu] against the repo's grounded tree counts
     (CLAUDE.md table = OEIS A000055 for n <= 26; the 2026-08-14 census
     totals for n = 27..32, each of which matched an independent
     `gentreeg -u` count).
  2. Brute force: full (n, mu) histograms from `gentreeg -p` for
     n = 4..BRUTE_MAX, mu computed by greedy leaf matching (exact on
     trees), alpha cross-checked per tree as deg I(T) via the tree DP,
     asserting alpha = n - mu (Konig) tree by tree.

Output: results/depth3_window_sizing_20260819.json
"""

import json
import math
import os
import subprocess
import sys

N_MAX = 44
MU_MAX = 22  # floor(44/2); mu of any tree on n <= 44 vertices is <= 22
BRUTE_MAX = 14

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

# Grounded reference totals (see docstring for provenance).
REF_TOTALS = {
    1: 1, 2: 1, 3: 1, 4: 2, 5: 3, 6: 6, 7: 11, 8: 23, 9: 47, 10: 106,
    15: 7741, 20: 823065, 25: 104636890, 26: 279793450,
    27: 751065460, 28: 2023443032, 29: 5469566585,
    30: 14830871802, 31: 40330829030, 32: 109972410221,
}


def zeros():
    return [[0] * (MU_MAX + 1) for _ in range(N_MAX + 1)]


def mul_item(S, n0, m0, f):
    """Return S * (1 - x^n0 y^m0)^(-f), truncated to (N_MAX, MU_MAX)."""
    R = zeros()
    for n in range(N_MAX + 1):
        jn_max = n // n0
        for m in range(MU_MAX + 1):
            acc = 0
            jm_max = m // m0 if m0 else jn_max
            for j in range(0, min(jn_max, jm_max) + 1):
                s = S[n - j * n0][m - j * m0]
                if s:
                    acc += math.comb(f + j - 1, j) * s
            R[n][m] = acc
    return R


def rooted_counts():
    """Return (r0, r1): counts of rooted trees by (n, mu) and root state."""
    r0, r1 = zeros(), zeros()
    A = zeros()  # MSET over R0 items seen so far
    B = zeros()  # MSET over R0+R1 items seen so far
    A[0][0] = B[0][0] = 1
    for n in range(1, N_MAX + 1):
        for m in range(MU_MAX + 1):
            r1[n][m] = A[n - 1][m]
            if m >= 1:
                r0[n][m] = B[n - 1][m - 1] - A[n - 1][m - 1]
        # make size-n rooted trees available as children of larger trees
        for m in range(MU_MAX + 1):
            if r0[n][m]:
                A = mul_item(A, n, m, r0[n][m])
                B = mul_item(B, n, m, r0[n][m])
            if r1[n][m]:
                B = mul_item(B, n, m, r1[n][m])
    return r0, r1


def unrooted_counts(r0, r1):
    """Otter dissymmetry with edge-join weights; return t[n][mu]."""
    P = zeros()  # ordered pairs, join-weighted
    for na in range(1, N_MAX):
        for ma in range(MU_MAX + 1):
            a0, a1 = r0[na][ma], r1[na][ma]
            if not (a0 or a1):
                continue
            for nb in range(1, N_MAX + 1 - na):
                for mb in range(MU_MAX + 1):
                    b0, b1 = r0[nb][mb], r1[nb][mb]
                    if not (b0 or b1):
                        continue
                    n, m = na + nb, ma + mb
                    if m <= MU_MAX:
                        P[n][m] += a0 * b0 + a0 * b1 + a1 * b0
                    if m + 1 <= MU_MAX:
                        P[n][m + 1] += a1 * b1
    D = zeros()  # A joined with an isomorphic copy of itself
    for na in range(1, N_MAX // 2 + 1):
        for ma in range(MU_MAX + 1):
            n = 2 * na
            if 2 * ma <= MU_MAX:
                D[n][2 * ma] += r0[na][ma]
            if 2 * ma + 1 <= MU_MAX:
                D[n][2 * ma + 1] += r1[na][ma]
    t = zeros()
    for n in range(1, N_MAX + 1):
        for m in range(MU_MAX + 1):
            num = 2 * (r0[n][m] + r1[n][m]) - P[n][m] + D[n][m]
            assert num % 2 == 0, (n, m)
            t[n][m] = num // 2
    return t


# ---------------------------------------------------------------- brute force

def parse_parents(line):
    return [int(x) for x in line.split()]


def matching_number(par, n):
    """Exact mu of a tree given a 1-indexed parent array (greedy from leaves)."""
    children = [[] for _ in range(n + 1)]
    order = []
    for v in range(2, n + 1):
        children[par[v - 1]].append(v)
    # vertices arrive in a root-first order from gentreeg; reverse = leaves first
    matched = [False] * (n + 1)
    mu = 0
    for v in range(n, 1, -1):
        p = par[v - 1]
        if not matched[v] and not matched[p]:
            matched[v] = matched[p] = True
            mu += 1
    return mu


def alpha_via_dp(par, n):
    """deg I(T) via the standard rooted tree DP, exact ints."""
    children = [[] for _ in range(n + 1)]
    for v in range(2, n + 1):
        children[par[v - 1]].append(v)
    # process vertices in reverse (children before parents in gentreeg order)
    excl = [None] * (n + 1)  # max ind set in subtree excluding v
    incl = [None] * (n + 1)  # including v
    for v in range(n, 0, -1):
        e = i = 0
        for c in children[v]:
            e += max(excl[c], incl[c])
            i += excl[c]
        excl[v], incl[v] = e, i + 1
    return max(excl[1], incl[1])


def brute_histograms():
    hists = {}
    for n in range(4, BRUTE_MAX + 1):
        out = subprocess.run(
            ["gentreeg", "-p", "-q", str(n)],
            capture_output=True, text=True, check=True,
        ).stdout
        h = {}
        for line in out.splitlines():
            if not line.strip():
                continue
            par = parse_parents(line)
            assert len(par) == n
            mu = matching_number(par, n)
            alpha = alpha_via_dp(par, n)
            assert alpha == n - mu, (line, alpha, mu)
            h[mu] = h.get(mu, 0) + 1
        hists[n] = h
    return hists


def main():
    r0, r1 = rooted_counts()
    t = unrooted_counts(r0, r1)

    # self-test 1: grounded totals
    for n, ref in sorted(REF_TOTALS.items()):
        tot = sum(t[n])
        assert tot == ref, (n, tot, ref)
    print(f"selftest totals: {len(REF_TOTALS)} grounded orders match "
          f"(incl. census n=27..32)")

    # self-test 2: brute-force (n, mu) histograms
    for n, h in brute_histograms().items():
        for m in range(MU_MAX + 1):
            assert t[n][m] == h.get(m, 0), (n, m, t[n][m], h.get(m, 0))
    print(f"selftest brute force: full (n,mu) histograms match for "
          f"n=4..{BRUTE_MAX} (alpha = n - mu verified per tree)")

    # dist-3 ladder rung d: alpha in {3d+8, 3d+9, 3d+10}, n in [33, 6d+20]
    # (n <= 32 exhausted by the census; alpha >= ceil(n/2) enforced by mu
    # bounds). Rung 3 is the named frontier; rung 4 sizes the ladder's growth.
    def rung(d):
        alphas = (3 * d + 8, 3 * d + 9, 3 * d + 10)
        window = {}
        for n in range(33, min(6 * d + 20, N_MAX) + 1):
            by_alpha = {}
            for alpha in alphas:
                mu = n - alpha
                if 0 <= mu <= MU_MAX:
                    c = t[n][mu]
                    if c:
                        by_alpha[alpha] = c
            window[n] = {
                "total_trees": sum(t[n]),
                "by_alpha": by_alpha,
                "window_trees": sum(by_alpha.values()),
            }
        return alphas, window

    results_by_rung = {}
    for d in (3, 4):
        alphas, window = rung(d)
        grand = sum(w["window_trees"] for w in window.values())
        grand_all = sum(w["total_trees"] for w in window.values())
        results_by_rung[d] = {"alphas": alphas, "window": window,
                              "window_total": grand, "orders_total": grand_all}
        print(f"\nDepth-{d} rung (alpha in {alphas}, "
              f"n in [33,{max(window)}]):")
        for n, w in window.items():
            if not w["window_trees"]:
                continue
            frac = w["window_trees"] / w["total_trees"]
            print(f"  n={n}: {w['window_trees']:>26,} of "
                  f"{w['total_trees']:>26,}  ({frac:.3e})  {w['by_alpha']}")
        print(f"  TOTAL rung {d}: {grand:,} of {grand_all:,} "
              f"({grand / grand_all:.3e})")
    window = results_by_rung[3]["window"]
    grand = results_by_rung[3]["window_total"]
    grand_all = results_by_rung[3]["orders_total"]

    # full mu-spectrum for the window orders, for the note
    spectrum = {n: {str(m): t[n][m] for m in range(MU_MAX + 1) if t[n][m]}
                for n in range(26, N_MAX + 1)}

    out = {
        "date": "2026-08-19",
        "method": "Otter dissymmetry, bivariate x^n y^mu, rooted state "
                  "s=[root exposable in some maximum matching]; "
                  "alpha = n - mu by Konig (trees bipartite)",
        "selftests": {
            "grounded_totals_orders": sorted(REF_TOTALS),
            "brute_force_orders": list(range(4, BRUTE_MAX + 1)),
        },
        "window_definition": "dist-3 ladder rung d: alpha in "
                             "{3d+8,3d+9,3d+10}, n in [33, 6d+20]; rung 3 "
                             "is the frontier after the n<=32 census",
        "rungs": {str(d): r for d, r in results_by_rung.items()},
        "window": window,
        "window_total": grand,
        "orders_total": grand_all,
        "mu_spectrum_n26_44": spectrum,
    }
    path = os.path.join(REPO, "results", "depth3_window_sizing_20260819.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
