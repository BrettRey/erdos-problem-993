#!/usr/bin/env python3
"""Attack: aperiodic composition families.

Trees built by iterating rooted growth operations under a control word.
Regular families iterate one map, so their root geometry is confined to
the regular-family region of the (dominance ratio, angle) plane;
aperiodic control (Fibonacci, Thue-Morse) composes different maps in a
non-periodic order, the standard source of exotic spectra in analogous
operator settings. Certified Arb root isolation throughout; exact
witness and R_gap checks on every constructed tree.

Ops on rooted (E, I) with parallel adjacency construction:
  P      new root, old root as only child
  S_w    new root: old root + w fresh leaves as children
  G_t    new root: old root + fresh spider S(2^t) center as children
  D_j    decorate current root with j leaves (no new root)

Run with the scratchpad flint venv.
"""

import json
import math
import sys
import time

sys.path.insert(0, ".")

from flint import fmpz_poly

from indpoly import independence_poly
from scripts.valley_search import valley_score
from scripts.valley_scaling_probe import kmul, kadd, shift, kpow, path_poly
from scratch_rgap_rescore_20260715 import rgap

OUT = "results/aperiodic_composition_20260716.json"


class T:
    """Rooted tree with exact (E, I) polys and explicit adjacency."""

    def __init__(self):
        self.E = [1]
        self.I = [0, 1]
        self.adj = [[]]
        self.root = 0

    def n(self):
        return len(self.adj)

    def _new(self, attach):
        self.adj.append([])
        v = len(self.adj) - 1
        self.adj[attach].append(v)
        self.adj[v].append(attach)
        return v

    def op_P(self):
        E, I = self.E, self.I
        self.E, self.I = kadd(E, I), shift(E)
        old = self.root
        self.adj.append([])
        v = len(self.adj) - 1
        self.adj[old].append(v)
        self.adj[v].append(old)
        self.root = v

    def op_S(self, w):
        E, I = self.E, self.I
        self.E = kmul(kadd(E, I), kpow([1, 1], w))
        self.I = shift(E)
        old = self.root
        self.adj.append([])
        v = len(self.adj) - 1
        self.adj[old].append(v)
        self.adj[v].append(old)
        self.root = v
        for _ in range(w):
            self._new(v)

    def op_G(self, t):
        Es = kpow(path_poly(2), t)
        Is = shift(kpow([1, 1], t))
        E, I = self.E, self.I
        self.E = kmul(kadd(E, I), kadd(Es, Is))
        self.I = kmul(shift(E), Es)
        old = self.root
        self.adj.append([])
        v = len(self.adj) - 1
        self.adj[old].append(v)
        self.adj[v].append(old)
        self.root = v
        c = self._new(v)
        for _ in range(t):
            x1 = self._new(c)
            self._new(x1)

    def op_D(self, j):
        self.E = kmul(self.E, kpow([1, 1], j))
        for _ in range(j):
            self._new(self.root)

    def poly(self):
        return kadd(self.E, self.I)


OPS = {
    "P": lambda tr: tr.op_P(),
    "S2": lambda tr: tr.op_S(2),
    "S5": lambda tr: tr.op_S(5),
    "G2": lambda tr: tr.op_G(2),
    "G3": lambda tr: tr.op_G(3),
    "G5": lambda tr: tr.op_G(5),
    "PD2": lambda tr: (tr.op_P(), tr.op_D(2)),
}


def fib_word(m):
    a, b = "a", "ab"
    while len(b) < m:
        a, b = b, b + a
    return b[:m]


def tm_word(m):
    s = "a"
    while len(s) < m:
        s = "".join("ab" if c == "a" else "ba" for c in s)
    return s[:m]


def words(m):
    import random
    rng = random.Random(31)
    return {
        "fib": fib_word(m),
        "thue-morse": tm_word(m),
        "period-ab": ("ab" * m)[:m],
        "period-aabb": ("aabb" * m)[:m],
        "pure-a": "a" * m,
        "pure-b": "b" * m,
        "random": "".join(rng.choice("ab") for _ in range(m)),
    }


def certified_points(poly):
    """All certified non-real roots as (modulus_ratio, angle) points."""
    P = fmpz_poly(list(poly))
    roots = P.complex_roots()
    mods = [abs(complex(float(z.real), float(z.imag))) for z, m in roots]
    zmin = min(mods)
    pts = []
    for (z, m), mod in zip(roots, mods):
        if z.imag.rad() >= abs(z.imag.mid()):
            continue
        phi = math.pi - abs(math.atan2(float(z.imag), float(z.real)))
        pts.append((mod / zmin, phi))
    return pts


def main():
    pairs = [("S2", "G3"), ("P", "G5"), ("P", "S5"), ("G2", "G5"),
             ("PD2", "G3")]
    sizes = (60, 130, 260)
    best = {}       # control -> best phi at ratio thresholds
    records = []
    t0 = time.time()
    for opa, opb in pairs:
        for wname, word in words(200).items():
            tr = T()
            checkpoints = set()
            for ch in word:
                (OPS[opa] if ch == "a" else OPS[opb])(tr)
                if tr.n() > max(sizes) + 30:
                    break
                for s in sizes:
                    if tr.n() >= s and s not in checkpoints:
                        checkpoints.add(s)
                        poly = tr.poly()
                        assert poly == independence_poly(
                            tr.n(), tr.adj), "map/adjacency mismatch"
                        vs = valley_score(poly)
                        if vs["witness"]:
                            print(f"*** WITNESS *** {opa}/{opb} {wname} "
                                  f"n={tr.n()}")
                            print("poly=", poly, flush=True)
                            json.dump({"witness": True,
                                       "edges": [(u, v) for u in
                                                 range(tr.n())
                                                 for v in tr.adj[u]
                                                 if u < v]},
                                      open(OUT, "w"))
                            return
                        pts = certified_points(poly)
                        r = rgap(poly, thetas=(0.001,))[0.001]
                        rec = {
                            "ops": f"{opa}/{opb}", "word": wname,
                            "n": tr.n(),
                            "points": [(round(a, 4), round(b, 5))
                                       for a, b in pts],
                            "R": r["R"] if r else None,
                            "cb": (r["c"] - r["b"]) if r else None,
                            "edges": [(u, v) for u in range(tr.n())
                                      for v in tr.adj[u] if u < v],
                        }
                        records.append(rec)
                        for thr in (1.1, 1.2, 1.5, 2.0):
                            cand = max((phi for ra, phi in pts
                                        if ra <= thr), default=0.0)
                            key = (wname, thr)
                            if cand > best.get(key, 0.0):
                                best[key] = cand

    print(f"{len(records)} constructions scored in {time.time()-t0:.0f}s\n")
    print(f"{'control':>12} " + " ".join(f"phi@r<={t:<4}" for t in
                                         (1.1, 1.2, 1.5, 2.0)))
    for wname in words(10):
        row = [best.get((wname, t), 0.0) for t in (1.1, 1.2, 1.5, 2.0)]
        print(f"{wname:>12} " + " ".join(f"{v:9.5f}" for v in row))
    mx = max((r["R"] or 0) for r in records)
    anyc = max((r["cb"] or 0) for r in records)
    print(f"\nbest R_gap={mx:.6f}; max rise distance c-b={anyc}")
    json.dump(records, open(OUT, "w"))
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
