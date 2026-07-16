#!/usr/bin/env python3
"""Scaling study of valley deficit 1 - V(n) for champion architectures.

For each champion family from the 2026-07-15 valley campaign, scale the
size parameter geometrically and measure the exact valley ratio V(n).
Decay diagnostics: n*(1-V) and sqrt(n)*(1-V). If 1-V ~ C/n^theta with
theta <= 1 and C > 0 stable, the family saturates below 1 (broom-like)
and cannot cross. Acceleration of the deficit toward 0 faster than any
stable power law would mark counterexample territory.

Exact arithmetic throughout; Kronecker substitution (packing coefficients
into one big integer) makes degree-2000 exact convolutions fast.
"""

import argparse
import json
import os
import sys
import time

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO)

from scripts.valley_search import valley_score  # noqa: E402


# ---------------------------------------------------------------------------
# Fast exact polynomial arithmetic via Kronecker substitution
# ---------------------------------------------------------------------------

def kmul(a: list[int], b: list[int]) -> list[int]:
    """Exact convolution of nonnegative-int polys via one bigint multiply."""
    if not a or not b:
        return []
    amax = max(a)
    bmax = max(b)
    if amax == 0 or bmax == 0:
        return [0] * (len(a) + len(b) - 1)
    # each output coeff < min(la,lb) * amax * bmax
    bound = min(len(a), len(b)) * amax * bmax
    B = bound.bit_length() + 1
    A = 0
    for c in reversed(a):
        A = (A << B) | c
    Bv = 0
    for c in reversed(b):
        Bv = (Bv << B) | c
    P = A * Bv
    mask = (1 << B) - 1
    out = []
    for _ in range(len(a) + len(b) - 1):
        out.append(P & mask)
        P >>= B
    return out


def kpow(a: list[int], e: int) -> list[int]:
    result = [1]
    base = a
    while e:
        if e & 1:
            result = kmul(result, base)
        base = kmul(base, base)
        e >>= 1
    return result


def kadd(a, b):
    if len(a) < len(b):
        a, b = b, a
    out = list(a)
    for i, c in enumerate(b):
        out[i] += c
    return out


def shift(a, k=1):
    return [0] * k + a


# Path polys (exact, small)
_P: list[list[int]] = [[1], [1, 1]]


def path_poly(l: int) -> list[int]:
    while len(_P) <= l:
        _P.append(kadd(_P[-1], shift(_P[-2])))
    return _P[l]


def gadget_EI(legs: dict[int, int]) -> tuple[list[int], list[int]]:
    """(center-excluded, center-included) for a spider gadget.

    legs: {leg_length: count}
    """
    E = [1]
    I = [1]
    for l, cnt in legs.items():
        E = kmul(E, kpow(path_poly(l), cnt))
        I = kmul(I, kpow(path_poly(l - 1), cnt))
    return E, shift(I)


def bouquet_from_gadgets(gadget_list: list[dict[int, int]]) -> tuple[list[int], list[int]]:
    """(root-excluded, root-included) polys for a bouquet of gadgets."""
    ex = [1]
    inc = [1]
    for legs in gadget_list:
        E, I = gadget_EI(legs)
        ex = kmul(ex, kadd(E, I))
        inc = kmul(inc, E)
    return ex, shift(inc)


def bouquet_total(gadget_list) -> list[int]:
    ex, inc = bouquet_from_gadgets(gadget_list)
    return kadd(ex, inc)


def dumbbell_total(gA: list[dict[int, int]], gB: list[dict[int, int]],
                   ell: int) -> list[int]:
    """Two bouquets joined root-to-root by a path with ell edges."""
    e, i = bouquet_from_gadgets(gA)
    # walk ell-1 intermediate path vertices
    for _ in range(ell - 1):
        e, i = kadd(e, i), shift(e)
    # B root: children = its gadgets + the chain end
    exB = [1]
    incB = [1]
    for legs in gB:
        E, I = gadget_EI(legs)
        exB = kmul(exB, kadd(E, I))
        incB = kmul(incB, E)
    total_ex = kmul(exB, kadd(e, i))
    total_inc = shift(kmul(incB, e))
    return kadd(total_ex, total_inc)


def gsize(legs: dict[int, int]) -> int:
    return 1 + sum(l * c for l, c in legs.items())


# ---------------------------------------------------------------------------
# Families
# ---------------------------------------------------------------------------

def family_point(name: str, scale: int, extra: dict) -> tuple[str, int, list[int]]:
    """Return (label, n, poly) for the given family at the given scale."""
    if name == "T3MN_h":
        # [S(2^3), S(2^M), S(2^M + 3^h)] champion, h fixed
        h = extra.get("h", 3)
        M = scale
        g = [{2: 3}, {2: M}, {2: M, 3: h}]
        n = 1 + sum(gsize(x) for x in g)
        return f"S(2^3)+S(2^{M})+S(2^{M}+3^{h})", n, bouquet_total(g)
    if name == "T3MN_pure":
        M = scale
        g = [{2: 3}, {2: M}, {2: M}]
        n = 1 + sum(gsize(x) for x in g)
        return f"T_(3,{M},{M})", n, bouquet_total(g)
    if name == "hybrid_2t":
        # a x S(2^2) + b x S(2^10), a:b ~ 34:7 champion ratio
        a = scale
        b = max(1, round(scale * 7 / 34))
        g = [{2: 2}] * a + [{2: 10}] * b
        n = 1 + sum(gsize(x) for x in g)
        return f"{a}xS(2^2)+{b}xS(2^10)", n, bouquet_total(g)
    if name == "dumbbell_3legs":
        # c x S(3^6) -- ell -- c x S(3^6)
        c = scale
        ell = extra.get("ell", 2)
        gA = [{3: 6}] * c
        n = 2 * (1 + c * gsize({3: 6}) - c * 0) + (ell - 1)
        # gsize({3:6}) = 19; bouquet size = 1 + 19c
        n = 2 * (1 + 19 * c) + (ell - 1)
        return f"{c}xS(3^6)--{ell}--{c}xS(3^6)", n, dumbbell_total(gA, gA, ell)
    if name == "dumbbell_mixed":
        # c x S(2^7) -- ell -- b x S(3^8): two different densities
        c = scale
        b = max(1, round(scale * 8 / 14))
        ell = extra.get("ell", 2)
        gA = [{2: 7}] * c
        gB = [{3: 8}] * b
        n = (1 + 15 * c) + (1 + 25 * b) + (ell - 1)
        return f"{c}xS(2^7)--{ell}--{b}xS(3^8)", n, dumbbell_total(gA, gB, ell)
    raise ValueError(name)


FAMILY_SCALES = {
    "T3MN_h": [40, 75, 150, 300, 600, 1200],
    "T3MN_pure": [40, 75, 150, 300, 600, 1200],
    "hybrid_2t": [34, 68, 136, 272, 544],
    "dumbbell_3legs": [12, 24, 48, 96, 192],
    "dumbbell_mixed": [14, 28, 56, 112, 224],
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--families", default="all")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    names = (list(FAMILY_SCALES) if args.families == "all"
             else args.families.split(","))
    rows = []
    for name in names:
        print(f"\n=== family {name} ===", flush=True)
        print(f"{'label':>42s} {'n':>6} {'alpha':>6} {'V':>14} "
              f"{'1-V':>12} {'n(1-V)':>10} {'sq(n)(1-V)':>10} "
              f"{'b/alpha':>8} {'win':>5} {'sec':>6}")
        for scale in FAMILY_SCALES[name]:
            t0 = time.time()
            label, n, poly = family_point(name, scale, {})
            vs = valley_score(poly)
            dt = time.time() - t0
            alpha = len(poly) - 1
            deficit = 1.0 - vs["ratio"]
            row = {
                "family": name, "label": label, "n": n, "alpha": alpha,
                "V": vs["ratio"], "deficit": deficit,
                "pos": vs["pos"], "rise_pos": vs["rise_pos"],
                "window": vs["window"], "witness": vs["witness"],
                "sec": round(dt, 2),
            }
            rows.append(row)
            print(f"{label:>42s} {n:>6} {alpha:>6} {vs['ratio']:>14.10f} "
                  f"{deficit:>12.3e} {n*deficit:>10.4f} "
                  f"{(n**0.5)*deficit:>10.4f} "
                  f"{vs['pos']/alpha if alpha else 0:>8.4f} "
                  f"{str(vs['window']):>5s} {dt:>6.1f}", flush=True)
            if vs["witness"]:
                print(f"*** WITNESS *** {label} poly written to out file")
                row["poly"] = poly

    if args.out:
        with open(args.out, "w") as f:
            json.dump(rows, f, indent=2)
        print(f"\nsaved {args.out}")


if __name__ == "__main__":
    main()
