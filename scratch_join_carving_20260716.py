#!/usr/bin/env python3
"""Join-vertex engineering: subtractive carving scan.

I(T1 -uv- T2) = I(T1) I(T2) - x^2 I(T1 - N[u]) I(T2 - N[v]).

The subtracted term is the only local mechanism in tree composition that
can CARVE a dip instead of smoothing one (additive mixtures all obey the
C/n barrier measured 2026-07-15/16). The carving profile depends on the
join vertices u, v: high-occupancy vertices maximize the correction.
This scan joins diverse components at structurally distinct vertex
pairs (root / gadget center / mid-leg / leaf) and scores the exact
joined polynomial for witnesses, R_gap, and rise distance c-b.

Exact arithmetic via the generic big-int tree DP.
"""

import sys
import time

sys.path.insert(0, ".")

from indpoly import independence_poly
from scripts.valley_search import bouquet_adj, valley_score
from scratch_rgap_rescore_20260715 import rgap


# ---------------------------------------------------------------------------
# Component library with labeled representative vertices
# ---------------------------------------------------------------------------

def component(spec_gadgets, spec_paths=(), label=""):
    """Return (label, n, adj, reps) with reps = [(name, vertex_id), ...]."""
    spec = (tuple(tuple(sorted(g)) for g in spec_gadgets),
            tuple(spec_paths), 0)
    n, adj = bouquet_adj(*spec)
    reps = [("root", 0)]
    # vertex layout in bouquet_adj: root=0, then per gadget: center,
    # then legs sequentially; then root paths.
    idx = 1
    for gi, legs in enumerate(spec[0]):
        center = idx
        idx += 1
        if gi == 0:
            reps.append(("g0_center", center))
            # first leg's first vertex (adjacent to center)
            if legs:
                reps.append(("g0_leg1", center + 1))
                # deepest vertex of first leg (a leaf)
                reps.append(("g0_leaf", center + legs[0]))
        idx += sum(legs)
    for pi, l in enumerate(spec[1]):
        if pi == 0:
            reps.append(("path1", idx))       # adjacent to root
            reps.append(("path_end", idx + l - 1))
        idx += l
    # dedupe vertex ids
    seen = set()
    out = []
    for name, v in reps:
        if v not in seen and v < n:
            seen.add(v)
            out.append((name, v))
    return label, n, adj, out


def build_components():
    comps = []
    # Galvin bouquets (LC-fragile shapes)
    comps.append(component([(2,) * 6] * 6, (), "T_{6,6,1}"))
    comps.append(component([(2,) * 10] * 10, (), "T_{10,10,1}"))
    # Kadrawi-Levit diagonal
    comps.append(component([(2,) * 3, (2,) * 8, (2,) * 9], (), "T_{3,8,9}"))
    # 3-leg bouquet (dumbbell champion block)
    comps.append(component([(3,) * 6] * 8, (), "8xS(3^6)"))
    # two-type hybrid block
    comps.append(component([(2,) * 2] * 20 + [(2,) * 10] * 4, (),
                           "20xS(2^2)+4xS(2^10)"))
    # subdivided star (dense-low)
    comps.append(component([], [2] * 25 + [3] * 8, "substar(25x2,8x3)"))
    # star-heavy (sharp low mode, high top reach via center-out)
    comps.append(component([(1,) * 30], (), "S(1^30)"))
    comps.append(component([], [1] * 40, "star40"))
    # broom (mid-mode)
    comps.append(component([], [1] * 30 + [12], "broom(12,30)"))
    # path (flat)
    comps.append(component([], [60], "P61"))
    # deep block
    comps.append(component([(4,) * 5] * 6, (), "6xS(4^5)"))
    return comps


def join_adj(nA, adjA, u, nB, adjB, v, ell=1):
    """Join u in A to v in B by a path with ell edges."""
    inter = max(0, ell - 1)
    n = nA + nB + inter
    adj = ([list(x) for x in adjA]
           + [[] for _ in range(inter)]
           + [[x + nA + inter for x in row] for row in adjB])

    def edge(a, b):
        adj[a].append(b)
        adj[b].append(a)

    prev = u
    for i in range(inter):
        edge(prev, nA + i)
        prev = nA + i
    edge(prev, v + nA + inter)
    return n, adj


def main():
    comps = build_components()
    print(f"{len(comps)} components; scanning join-vertex pairs", flush=True)
    t0 = time.time()
    results = []
    count = 0
    for i in range(len(comps)):
        lA, nA, adjA, repsA = comps[i]
        for j in range(i, len(comps)):
            lB, nB, adjB, repsB = comps[j]
            for (nameA, u) in repsA:
                for (nameB, v) in repsB:
                    if j == i and (nameB, v) < (nameA, u):
                        continue
                    for ell in (1, 2):
                        n, adj = join_adj(nA, adjA, u, nB, adjB, v, ell)
                        poly = independence_poly(n, adj)
                        vs = valley_score(poly)
                        count += 1
                        if vs["witness"]:
                            print(f"\n*** WITNESS *** {lA}[{nameA}] --{ell}-- "
                                  f"{lB}[{nameB}] n={n}")
                            print("poly=", poly, flush=True)
                            return
                        r = rgap(poly, thetas=(0.001,))[0.001]
                        if r:
                            results.append(
                                (r["R"], r["c"] - r["b"], n,
                                 f"{lA}[{nameA}] --{ell}-- {lB}[{nameB}]",
                                 r["b"]))
    results.sort(key=lambda t: (-(t[1] > 1), -t[0]))
    print(f"\nscanned {count} joins in {time.time()-t0:.0f}s; no witness")
    multi = [r for r in results if r[1] > 1]
    print(f"joins with rise distance c-b > 1: {len(multi)}")
    print("\ntop by (c-b>1 first, then R):")
    for R, gap, n, label, b in results[:20]:
        print(f"  R={R:.8f} n(1-R)={n*(1-R):8.2f} c-b={gap} b={b} n={n}  {label}")


if __name__ == "__main__":
    main()
