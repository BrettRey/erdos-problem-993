#!/usr/bin/env python3
"""Certified root-sector probe for the subcubic-LC lane (2026-08-26).

Sector lemma (elementary, replayed in the lane note): a polynomial with
positive coefficients whose roots all satisfy |arg(-z)| <= pi/3 factors over R
into linear factors (r + x) and quadratic factors (r^2 + 2r cos(phi) x + x^2)
with phi <= pi/3, each log-concave with no internal zeros; LC sequences are
closed under convolution, so the coefficient sequence is log-concave.
Contrapositive: every LC-failing tree has a root with deviation > pi/3 from
the negative real axis.

This probe (all roots certified via python-flint / Arb, per the project's
mandatory root-finding rule; no numpy.roots anywhere):

  1. Self-check: the 21 stored LC-failing trees (2 at n=26 from
     analysis_n26.json, 19 at n=28 from the census shards) must each have a
     certified deviation > pi/3.
  2. Max deviation over ALL subcubic trees at n = 14, 16, 18 (gentreeg -D3),
     and a 1-in-k sample at n = 20.
  3. Structured large subcubic families: binary caterpillars, complete binary
     trees, subdivided binary spiders, up to n ~ 95.

Question: do subcubic trees respect a uniform sector bound (ideally <= pi/3,
which with the lemma would PROVE subcubic LC), and how does the max deviation
move with n and with Delta?

Certified comparisons: an angle claim "dev > pi/3" or "dev < pi/3" is reported
only when the Arb ball for (dev - pi/3) excludes zero. Undecided balls are
counted separately. Midpoints are for reporting only.
"""
import json
import subprocess
import sys
from collections import Counter

sys.setrecursionlimit(100000)

from flint import arb, acb, ctx, fmpz_poly  # noqa: E402

ctx.prec = 128
PI3 = arb.pi() / 3


def ipoly_from_par(par):
    """Exact independence polynomial from a 1-indexed parent array (par[1]=0)."""
    n = len(par) - 1
    children = [[] for _ in range(n + 1)]
    for v in range(2, n + 1):
        children[par[v]].append(v)
    F = {}
    G = {}
    for v in range(n, 0, -1):  # children have larger indices
        f = [1]
        g = [1]
        for c in children[v]:
            fc_gc = [a + b for a, b in
                     zip(F[c] + [0] * (len(G[c]) - len(F[c])),
                         G[c] + [0] * (len(F[c]) - len(G[c])))]
            f = mul(f, fc_gc)
            g = mul(g, F[c])
        F[v] = f
        G[v] = [0] + g  # times x
    top = [a + b for a, b in zip(F[1] + [0] * (len(G[1]) - len(F[1])),
                                 G[1] + [0] * (len(F[1]) - len(G[1])))]
    return top


def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def lc_ok(seq):
    return all(seq[k] * seq[k] >= seq[k - 1] * seq[k + 1]
               for k in range(1, len(seq) - 1))


def max_deviation(coeffs):
    """Certified max deviation |arg(-z)| over roots; returns (arb, verdict).

    verdict in {'gt', 'lt', 'undecided'} comparing the max against pi/3."""
    roots = fmpz_poly(list(coeffs)).complex_roots()
    if sum(m for _z, m in roots) != len(coeffs) - 1:
        raise AssertionError("root multiplicities do not sum to degree")
    devs = []
    for z, _m in roots:
        if z.imag.contains(0):
            # Certified real root: must be negative (positive coefficients).
            if not z.real.upper() < 0:
                raise AssertionError("real root ball not certified negative")
            continue  # deviation 0 from the axis
        devs.append(arb.pi() - abs(z.arg()))
    if not devs:
        return arb(0), "lt"
    m = devs[0]
    for d in devs[1:]:
        m = arb.max(m, d)
    diff = m - PI3
    if diff.lower() > 0:
        return m, "gt"
    if diff.upper() < 0:
        return m, "lt"
    return m, "undecided"


def gentreeg_pars(n, extra=("-D3",)):
    cmd = ["gentreeg", *extra, "-p", "-q", str(n)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        yield [0] + [int(x) for x in line.split()]
    proc.wait()


# ---- structured families (edge lists -> parent arrays via BFS) ----

def par_from_edges(n, edges):
    adj = [[] for _ in range(n + 1)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    order, seen, queue = [], {1}, [1]
    par_map = {1: 0}
    while queue:
        u = queue.pop(0)
        order.append(u)
        for w in adj[u]:
            if w not in seen:
                seen.add(w)
                par_map[w] = u
                queue.append(w)
    idx = {v: i + 1 for i, v in enumerate(order)}
    par = [0] * (n + 1)
    for v in order[1:]:
        par[idx[v]] = idx[par_map[v]]
    return par


def binary_caterpillar(m):
    """Spine P_m, one pendant leaf per spine vertex; n = 2m, Delta = 3."""
    edges = [(i, i + 1) for i in range(1, m)]
    edges += [(i, m + i) for i in range(1, m + 1)]
    return 2 * m, edges


def complete_binary(depth):
    n = 2 ** (depth + 1) - 1
    edges = [(i, 2 * i) for i in range(1, n + 1) if 2 * i <= n]
    edges += [(i, 2 * i + 1) for i in range(1, n + 1) if 2 * i + 1 <= n]
    return n, edges


def binary_spider(arm, k=3):
    """k paths of length `arm` glued at a center; Delta = k."""
    n = 1 + k * arm
    edges = []
    v = 2
    for _ in range(k):
        prev = 1
        for _ in range(arm):
            edges.append((prev, v))
            prev = v
            v += 1
    return n, edges


def main():
    report = {"pi_over_3": float(PI3), "self_check_failures": [],
              "exhaustive": {}, "families": {}}

    # 1. Self-check on the 21 stored LC failures.
    import glob
    import re
    import networkx as nx
    fail_polys = []
    d = json.load(open("results/analysis_n26.json"))
    fail_polys += [f["poly"] for f in d["lc_failures"]]
    for path in sorted(glob.glob("results/lc_census_20260814/n28_s*of64.txt")):
        for line in open(path):
            if line.startswith("FAIL"):
                poly = [int(x) for x in
                        re.search(r"poly=([\d,]+)", line).group(1).split(",")]
                fail_polys.append(poly)
    ok = 0
    for poly in fail_polys:
        m, verdict = max_deviation(poly)
        if verdict != "gt":
            raise AssertionError(
                f"SECTOR LEMMA VIOLATION? failing tree with dev {m} not > pi/3")
        ok += 1
    print(f"[self-check] all {ok} stored LC-failing trees have certified "
          f"deviation > pi/3  (lemma contrapositive holds)")
    report["self_check_count"] = ok

    # 2. Exhaustive subcubic sweeps.
    for n, step in ((14, 1), (16, 1), (18, 1), (20, 8)):
        best = None
        count = verdict_gt = undecided = non_lc = 0
        for i, par in enumerate(gentreeg_pars(n)):
            if i % step:
                continue
            count += 1
            poly = ipoly_from_par(par)
            if not lc_ok(poly):
                non_lc += 1
            m, verdict = max_deviation(poly)
            if verdict == "gt":
                verdict_gt += 1
            elif verdict == "undecided":
                undecided += 1
            key = float(m.mid()) if hasattr(m, "mid") else float(m)
            if best is None or key > best[0]:
                best = (key, par, verdict)
        report["exhaustive"][n] = {
            "trees_checked": count, "sample_step": step,
            "max_deviation_mid": best[0], "argmax_par": best[1],
            "argmax_verdict_vs_pi3": best[2],
            "count_certified_gt_pi3": verdict_gt,
            "count_undecided": undecided, "non_lc": non_lc,
        }
        print(f"[subcubic n={n}] checked {count} (step {step}): "
              f"max dev ~ {best[0]:.5f} rad (pi/3 = {float(PI3):.5f}), "
              f"argmax verdict {best[2]}, certified>pi/3: {verdict_gt}, "
              f"undecided: {undecided}, non-LC: {non_lc}")

    # 3. Structured families.
    fams = []
    for m in (10, 20, 30, 45):
        fams.append((f"binary_caterpillar_m{m}", *binary_caterpillar(m)))
    for depth in (4, 5, 6):
        fams.append((f"complete_binary_d{depth}", *complete_binary(depth)))
    for arm in (10, 20, 30):
        fams.append((f"binary_spider_arm{arm}", *binary_spider(arm)))
    for name, n, edges in fams:
        par = par_from_edges(n, edges)
        poly = ipoly_from_par(par)
        m, verdict = max_deviation(poly)
        report["families"][name] = {
            "n": n, "max_deviation_mid": float(m.mid()),
            "verdict_vs_pi3": verdict, "lc": lc_ok(poly),
        }
        print(f"[{name}] n={n}: max dev ~ {float(m.mid()):.5f} "
              f"({verdict} pi/3), LC={lc_ok(poly)}")

    with open("results/sector_angle_probe_20260826.json", "w") as fh:
        json.dump(report, fh, indent=1)
    print("saved results/sector_angle_probe_20260826.json")


if __name__ == "__main__":
    main()
