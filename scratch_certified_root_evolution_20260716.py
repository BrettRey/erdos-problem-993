#!/usr/bin/env python3
"""Certified root-herding evolution.

Maximizes the certified near-dominant root angle phi_dom over free tree
mutations, with Arb (python-flint) certified root isolation inside the
scoring loop. No float64 root-finding anywhere. Every child is also
checked exactly for a unimodality witness.

Run with the scratchpad venv python (needs flint + networkx).
"""

import json
import math
import random
import sys
import time

sys.path.insert(0, ".")

import networkx as nx
from flint import fmpz_poly

from indpoly import independence_poly
from scripts.valley_search import bouquet_adj, valley_score
from scratch_valley_freeform_20260715 import mutate_tree
from scratch_rgap_rescore_20260715 import rgap

OUT = "results/certified_root_evolution_20260716.json"


def certified_phi(poly, dom_factor=2.0):
    """Certified (phi_dom, phi_max). Conservative: a root counts as
    non-real only when its imaginary ball excludes zero."""
    P = fmpz_poly(list(poly))
    roots = P.complex_roots()
    mods = [abs(complex(float(z.real), float(z.imag))) for z, m in roots]
    zmin = min(mods)
    bd, ba = 0.0, 0.0
    for (z, m), mod in zip(roots, mods):
        if z.imag.rad() >= abs(z.imag.mid()):
            continue
        phi = math.pi - abs(math.atan2(float(z.imag), float(z.real)))
        ba = max(ba, phi)
        if mod <= dom_factor * zmin:
            bd = max(bd, phi)
    return bd, ba


def seeds():
    out = []

    def bq(g, p, label):
        spec = (tuple(tuple(sorted(x)) for x in g), tuple(p), 0)
        n, adj = bouquet_adj(*spec)
        return [label, n, adj]

    out.append(bq([(2,) * 6] * 6, [], "T_{6,6,1}"))       # phi_dom 0.256
    out.append(bq([(2,) * 8] * 8, [], "T_{8,8,1}"))
    out.append(bq([(2,) * 5] * 7, [], "T_{7,5,1}"))
    out.append(bq([(3,) * 6] * 4, [], "4xS(3^6)"))        # 0.151
    d28 = json.load(open("results/analysis_n28_modal_lc_nm.json"))
    for i, f in enumerate(d28["top_lc_failures"][:3]):
        G = nx.from_graph6_bytes(f["graph6"].encode())
        adj = [sorted(G.neighbors(v)) for v in range(G.number_of_nodes())]
        out.append([f"F28_{i}", len(adj), adj])
    return out


def main():
    budget = float(sys.argv[1]) if len(sys.argv) > 1 else 3000
    rng = random.Random(716)
    pop = []
    for label, n, adj in seeds():
        poly = independence_poly(n, adj)
        pd, pm = certified_phi(poly)
        pop.append([pd, pm, n, adj, label])
        print(f"seed {label}: phi_dom={pd:.5f}", flush=True)
    pop.sort(key=lambda r: -r[0])
    best = pop[0][0]

    t0 = time.time()
    evals = 0
    accepted = 0
    while time.time() - t0 < budget:
        parent = pop[rng.randrange(min(6, len(pop)))]
        m = mutate_tree(parent[2], parent[3], rng)
        if m is None:
            continue
        n2, adj2 = m
        if n2 < 30 or n2 > 150:
            continue
        poly = independence_poly(n2, adj2)
        vs = valley_score(poly)
        if vs["witness"]:
            print("*** WITNESS ***")
            print("poly=", poly)
            json.dump({"witness": True, "n": n2,
                       "edges": [(u, v) for u in range(n2)
                                 for v in adj2[u] if u < v],
                       "poly": [str(c) for c in poly]},
                      open(OUT, "w"))
            return
        try:
            pd, pm = certified_phi(poly)
        except Exception as e:
            print("root isolation failed:", e, flush=True)
            continue
        evals += 1
        if pd > pop[-1][0]:
            pop.append([pd, pm, n2, adj2, "evolved"])
            pop.sort(key=lambda r: -r[0])
            pop = pop[:10]
            accepted += 1
            if pd > best:
                best = pd
                r = rgap(poly, thetas=(0.001,))[0.001]
                rr = (f"R={r['R']:.6f} c-b={r['c']-r['b']}"
                      if r else "no descent")
                print(f"[{time.time()-t0:7.0f}s] new best phi_dom={pd:.5f} "
                      f"phi_max={pm:.5f} n={n2} {rr} "
                      f"(evals={evals})", flush=True)
                json.dump({"witness": False, "phi_dom": pd, "phi_max": pm,
                           "n": n2,
                           "edges": [(u, v) for u in range(n2)
                                     for v in adj2[u] if u < v]},
                          open(OUT, "w"))
    print(f"\ndone: evals={evals} accepted={accepted} "
          f"best certified phi_dom={best:.5f} "
          f"({time.time()-t0:.0f}s)", flush=True)
    for pd, pm, n, adj, label in pop[:5]:
        print(f"  phi_dom={pd:.5f} phi_max={pm:.5f} n={n} {label}")


if __name__ == "__main__":
    main()
