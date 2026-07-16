#!/usr/bin/env python3
"""Direct valley (non-unimodality) search over the spider-bouquet grammar.

Disproof-oriented search for Erdos Problem #993. Unlike targeted.py and
nm_optimizer.py, which rank by log-concavity defect or tail near-miss,
this scores every tree by a direct valley ratio:

    V(T) = max_b  min( max_{a<b} i_a , max_{c>b} i_c ) / i_b

V(T) > 1 is exactly a unimodality counterexample (a dip at b with higher
coefficients on both sides). The comparison is done in exact integer
arithmetic; the float ratio is only used for ranking.

Grammar: a root joined to k "spider gadgets" (a center vertex carrying
paths of prescribed lengths), plus optional path arms and pendant leaves
on the root itself. This one grammar realizes:

  - Galvin T_{m,t,1}        = root + m copies of gadget S(2^t)
  - Kadrawi-Levit T_{3,m,n} = root + gadgets S(2^3), S(2^m), S(2^n)
  - T*_{3,m,n}              = same with one leg of the first gadget
                              lengthened from 2 to 4 (Li 2026, Fig. 1.2)
  - multi-arm stars, spiders, brooms, subdivided stars as degenerate cases

Gadget polynomials are cached by leg multiset, so large sweeps that reuse
gadgets are much faster than generic tree DP.

Usage:
  python3 scripts/valley_search.py --selftest
  python3 scripts/valley_search.py --sweep A --max-n 320 --out results.json
  python3 scripts/valley_search.py --sweep all --max-n 320 --hillclimb 40
"""

import argparse
import heapq
import json
import os
import random
import sys
import time

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO)

from indpoly import independence_poly, is_unimodal  # noqa: E402


# ---------------------------------------------------------------------------
# Exact polynomial helpers (arbitrary precision, plain lists)
# ---------------------------------------------------------------------------

def _polymul(a: list[int], b: list[int]) -> list[int]:
    la, lb = len(a), len(b)
    out = [0] * (la + lb - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                out[i + j] += ai * bj
    return out


def _polyadd(a: list[int], b: list[int]) -> list[int]:
    if len(a) < len(b):
        a, b = b, a
    out = list(a)
    for i, bi in enumerate(b):
        out[i] += bi
    return out


# Path polynomials. _PATH[l] = I(P_l; x); P_0 is the empty graph.
_PATH: list[list[int]] = [[1], [1, 1]]


def _path_poly(length: int) -> list[int]:
    while len(_PATH) <= length:
        # I(P_l) = I(P_{l-1}) + x I(P_{l-2})  (delete/include an end vertex)
        nxt = _polyadd(_PATH[-1], [0] + _PATH[-2])
        _PATH.append(nxt)
    return _PATH[length]


# Gadget cache: leg multiset -> (E, I) with E = center excluded, I = included.
_GADGET: dict[tuple[int, ...], tuple[list[int], list[int]]] = {}


def _gadget_polys(legs: tuple[int, ...]) -> tuple[list[int], list[int]]:
    key = tuple(sorted(legs))
    hit = _GADGET.get(key)
    if hit is not None:
        return hit
    excl = [1]
    incl = [1]
    for l in key:
        excl = _polymul(excl, _path_poly(l))
        incl = _polymul(incl, _path_poly(l - 1))
    incl = [0] + incl  # factor x for the center itself
    _GADGET[key] = (excl, incl)
    return excl, incl


def bouquet_poly(gadgets: list[tuple[int, ...]],
                 root_paths: tuple[int, ...] = (),
                 root_leaves: int = 0) -> list[int]:
    """Independence polynomial of the bouquet tree, exactly."""
    root_excl = [1]
    root_incl = [1]
    for legs in gadgets:
        e, i = _gadget_polys(legs)
        root_excl = _polymul(root_excl, _polyadd(e, i))
        root_incl = _polymul(root_incl, e)
    for l in root_paths:
        root_excl = _polymul(root_excl, _path_poly(l))
        root_incl = _polymul(root_incl, _path_poly(l - 1))
    if root_leaves:
        onepx = [1, 1]
        for _ in range(root_leaves):
            root_excl = _polymul(root_excl, onepx)
        # leaves must be excluded when the root is in the set: factor 1
    root_incl = [0] + root_incl  # factor x for the root
    return _polyadd(root_excl, root_incl)


def bouquet_size(gadgets, root_paths=(), root_leaves=0) -> int:
    return (1 + sum(1 + sum(legs) for legs in gadgets)
            + sum(root_paths) + root_leaves)


def bouquet_adj(gadgets, root_paths=(), root_leaves=0):
    """Explicit adjacency list (for cross-checking and export)."""
    n = bouquet_size(gadgets, root_paths, root_leaves)
    adj: list[list[int]] = [[] for _ in range(n)]

    def edge(u, v):
        adj[u].append(v)
        adj[v].append(u)

    nxt = 1
    for legs in gadgets:
        center = nxt
        nxt += 1
        edge(0, center)
        for l in legs:
            prev = center
            for _ in range(l):
                edge(prev, nxt)
                prev = nxt
                nxt += 1
    for l in root_paths:
        prev = 0
        for _ in range(l):
            edge(prev, nxt)
            prev = nxt
            nxt += 1
    for _ in range(root_leaves):
        edge(0, nxt)
        nxt += 1
    assert nxt == n
    return n, adj


# ---------------------------------------------------------------------------
# Valley score
# ---------------------------------------------------------------------------

def valley_score(seq: list[int]) -> dict:
    """Direct non-unimodality score, exact comparison, float ranking.

    Returns dict with:
      ratio      max_b min(prefixmax, suffixmax)/i_b  (float)
      pos        argmax b
      rise_pos   first c > pos with i_c = suffix max at pos
      window     True if rise_pos < ceil((2*alpha - 1)/3), the region where
                 an ascent is not already excluded by the decreasing-tail
                 theorem (Levit-Mandrescu; see Basit-Galvin restatement)
      witness    True iff min(prefixmax, suffixmax) > i_b exactly
    """
    m = len(seq)
    if m < 3:
        return {"ratio": 0.0, "pos": -1, "rise_pos": -1,
                "window": False, "witness": False}
    alpha = m - 1
    pref = [0] * m
    best_prefix = seq[0]
    for b in range(1, m):
        pref[b] = best_prefix
        if seq[b] > best_prefix:
            best_prefix = seq[b]
    suff = [0] * m
    best_suffix = seq[-1]
    for b in range(m - 2, -1, -1):
        suff[b] = best_suffix
        if seq[b] > best_suffix:
            best_suffix = seq[b]

    best_ratio = 0.0
    best_b = -1
    witness = False
    for b in range(1, m - 1):
        lo = pref[b] if pref[b] < suff[b] else suff[b]
        if seq[b] == 0:
            continue
        if lo > seq[b]:
            witness = True
        r = lo / seq[b]
        if r > best_ratio:
            best_ratio = r
            best_b = b

    rise_pos = -1
    if best_b >= 0:
        target = suff[best_b]
        for c in range(best_b + 1, m):
            if seq[c] == target:
                rise_pos = c
                break
    tail_start = -((-(2 * alpha - 1)) // 3)  # ceil((2 alpha - 1)/3)
    window = 0 <= rise_pos < tail_start
    return {"ratio": best_ratio, "pos": best_b, "rise_pos": rise_pos,
            "window": window, "witness": witness}


# ---------------------------------------------------------------------------
# Spec handling
# ---------------------------------------------------------------------------

def canon(spec) -> tuple:
    gadgets, root_paths, root_leaves = spec
    g = tuple(sorted(tuple(sorted(l)) for l in gadgets))
    return (g, tuple(sorted(root_paths)), root_leaves)


def spec_label(spec) -> str:
    gadgets, root_paths, root_leaves = spec
    parts = []
    from collections import Counter
    cnt = Counter(tuple(sorted(l)) for l in gadgets)
    for legs, k in sorted(cnt.items()):
        legc = Counter(legs)
        legstr = "+".join(f"{v}^{c}" if c > 1 else f"{v}"
                          for v, c in sorted(legc.items()))
        parts.append(f"{k}xS({legstr})")
    if root_paths:
        parts.append("paths" + str(tuple(sorted(root_paths))))
    if root_leaves:
        parts.append(f"{root_leaves}leaves")
    return " ".join(parts) if parts else "empty"


def evaluate_spec(spec) -> dict:
    gadgets, root_paths, root_leaves = spec
    poly = bouquet_poly(list(gadgets), tuple(root_paths), root_leaves)
    vs = valley_score(poly)
    n = bouquet_size(gadgets, root_paths, root_leaves)
    out = {
        "label": spec_label(spec),
        "spec": [ [list(l) for l in spec[0]], list(spec[1]), spec[2] ],
        "n": n,
        "alpha": len(poly) - 1,
        **vs,
    }
    if vs["witness"]:
        out["poly"] = poly
    return out


# ---------------------------------------------------------------------------
# Sweeps
# ---------------------------------------------------------------------------

def sweep_A(max_n: int):
    """Perturbed T_{3,m,n} / T*_{3,m,n}."""
    for m in range(1, 90):
        for n_ in range(m, 90):
            base = [(2,) * 3, (2,) * m, (2,) * n_]
            if bouquet_size(base) > max_n:
                break
            # unperturbed + one-leg extensions of the small gadget (2->3..6)
            yield (tuple(map(tuple, base)), (), 0)
            for ext in (3, 4, 5, 6):
                g0 = tuple([ext] + [2] * 2)
                yield ((g0, (2,) * m, (2,) * n_), (), 0)
            # one extra leg on the small gadget of length 1..3
            for extra in (1, 2, 3):
                g0 = tuple([2, 2, 2, extra])
                yield ((g0, (2,) * m, (2,) * n_), (), 0)
            # extra root decorations
            for p in (1, 2, 3):
                yield (tuple(map(tuple, base)), (p,), 0)
            # fourth small gadget
            for g in (2, 3, 4):
                spec = (tuple(map(tuple, base)) + ((2,) * g,), (), 0)
                if bouquet_size(spec[0]) <= max_n:
                    yield spec
            # mixed leg lengths in one big gadget: a legs of 2, b of 3
            for b3 in (1, 2, 3):
                if n_ - b3 >= 1:
                    g2 = tuple([2] * (n_ - b3) + [3] * b3)
                    yield (((2,) * 3, (2,) * m, g2), (), 0)


def sweep_B(max_n: int):
    """Hub-bouquets: s copies of S(2^t) + exceptional gadget(s)."""
    menu = []
    for g in range(1, 9):
        menu.append((2,) * g)
        menu.append((1,) * g)
        menu.append(tuple([2] * g + [3]))
        menu.append(tuple([3] * g))
    for L in (4, 6, 8, 12, 16):
        menu.append((L,))
        menu.append((L, 2, 2))
    menu = sorted(set(menu))
    for t in range(2, 13):
        for s in range(3, 41):
            base = [(2,) * t] * s
            if bouquet_size(base) > max_n:
                break
            yield (tuple(base), (), 0)
            for exc in menu:
                spec_g = tuple(base) + (exc,)
                if bouquet_size(spec_g) <= max_n:
                    yield (spec_g, (), 0)
            for p in (1, 2, 3, 5):
                yield (tuple(base), (p,), 0)


def sweep_C(max_n: int):
    """Two-type Galvin hybrids: a copies of S(2^t1) + b copies of S(2^t2)."""
    for t1 in range(2, 12):
        for t2 in range(t1 + 1, 14):
            for a in range(1, 40):
                if bouquet_size([(2,) * t1] * a) > max_n:
                    break
                for b in range(1, 40):
                    spec_g = tuple([(2,) * t1] * a + [(2,) * t2] * b)
                    if bouquet_size(spec_g) > max_n:
                        break
                    yield (spec_g, (), 0)


def sweep_D(max_n: int):
    """Mixed subdivided stars at the root + one long arm."""
    for a2 in range(0, 60):
        for a3 in range(0, 40):
            if a2 + a3 < 3:
                continue
            base_paths = [2] * a2 + [3] * a3
            if 1 + sum(base_paths) > max_n:
                break
            for L in (0, 1, 4, 7, 10, 15):
                paths = base_paths + ([L] if L else [])
                if 1 + sum(paths) > max_n:
                    continue
                yield ((), tuple(paths), 0)


SWEEPS = {"A": sweep_A, "B": sweep_B, "C": sweep_C, "D": sweep_D}


# ---------------------------------------------------------------------------
# Hill climbing on champion specs
# ---------------------------------------------------------------------------

def mutate_spec(spec, rng: random.Random, max_n: int):
    gadgets, root_paths, root_leaves = spec
    gadgets = [list(l) for l in gadgets]
    root_paths = list(root_paths)
    for _ in range(30):
        g2 = [list(l) for l in gadgets]
        p2 = list(root_paths)
        l2 = root_leaves
        op = rng.randrange(8)
        if op == 0 and g2:  # +-1 a leg length
            gi = rng.randrange(len(g2))
            li = rng.randrange(len(g2[gi]))
            g2[gi][li] = max(1, g2[gi][li] + rng.choice((-1, 1)))
        elif op == 1 and g2:  # add a leg
            gi = rng.randrange(len(g2))
            g2[gi].append(rng.choice((1, 2, 2, 2, 3)))
        elif op == 2 and g2:  # remove a leg
            gi = rng.randrange(len(g2))
            if len(g2[gi]) > 1:
                g2[gi].pop(rng.randrange(len(g2[gi])))
            else:
                g2.pop(gi)
        elif op == 3 and g2:  # duplicate a gadget
            g2.append(list(rng.choice(g2)))
        elif op == 4 and len(g2) > 1:  # drop a gadget
            g2.pop(rng.randrange(len(g2)))
        elif op == 5:  # root path tweak
            if p2 and rng.random() < 0.5:
                pi = rng.randrange(len(p2))
                p2[pi] = max(1, p2[pi] + rng.choice((-1, 1)))
            else:
                p2.append(rng.choice((1, 2, 3, 4)))
        elif op == 6:  # root leaves tweak
            l2 = max(0, l2 + rng.choice((-1, 1)))
        elif op == 7 and g2:  # clone-with-perturb: copy gadget, change one leg
            src = list(rng.choice(g2))
            if src:
                src[rng.randrange(len(src))] = rng.choice((1, 2, 3, 4))
                g2.append(src)
        cand = (tuple(tuple(l) for l in g2 if l), tuple(p2), l2)
        if not cand[0] and not cand[1]:
            continue
        if bouquet_size(cand[0], cand[1], cand[2]) <= max_n:
            return cand
    return spec


def hill_climb(seeds: list[dict], max_n: int, generations: int,
               pop: int, rng: random.Random, seen: set):
    """Evolutionary refinement of the best specs by valley ratio."""
    def key(r):
        return (r["window"], r["ratio"])

    population = sorted(seeds, key=key, reverse=True)[:pop]
    best = dict(population[0]) if population else None
    for gen in range(generations):
        children = []
        for r in population:
            spec = (tuple(tuple(l) for l in r["spec"][0]),
                    tuple(r["spec"][1]), r["spec"][2])
            for _ in range(6):
                child = mutate_spec(spec, rng, max_n)
                c = canon(child)
                if c in seen:
                    continue
                seen.add(c)
                res = evaluate_spec(child)
                children.append(res)
                if res["witness"]:
                    return res, population
        population = sorted(population + children, key=key, reverse=True)[:pop]
        if population and key(population[0]) > key(best):
            best = dict(population[0])
        if gen % 10 == 0:
            b = population[0]
            print(f"  gen {gen:3d}  best ratio={b['ratio']:.10f} "
                  f"window={b['window']} n={b['n']} {b['label']}",
                  flush=True)
    return None, population


# ---------------------------------------------------------------------------
# Self test
# ---------------------------------------------------------------------------

def selftest():
    rng = random.Random(7)
    # 1. bouquet_poly vs generic tree DP on random specs
    for trial in range(60):
        k = rng.randint(0, 4)
        gadgets = tuple(tuple(rng.randint(1, 4)
                              for _ in range(rng.randint(1, 5)))
                        for _ in range(k))
        paths = tuple(rng.randint(1, 5) for _ in range(rng.randint(0, 3)))
        leaves = rng.randint(0, 3)
        if not gadgets and not paths and not leaves:
            continue
        fast = bouquet_poly(list(gadgets), paths, leaves)
        n, adj = bouquet_adj(gadgets, paths, leaves)
        slow = independence_poly(n, adj)
        assert fast == slow, (gadgets, paths, leaves, fast, slow)
    # 2. valley score sanity: known unimodal and a synthetic valley
    assert not valley_score([1, 5, 9, 7, 3])["witness"]
    v = valley_score([1, 5, 4, 6, 2])
    assert v["witness"] and v["pos"] == 2
    assert valley_score([1, 5, 5, 5, 2])["witness"] is False
    # 3. T_{3,m,n} matches Li's structure: n = 1 + 3 + 3+m+n + 2*(3+m+n)?
    #    size = 1 root + 3 centers + (3+m+n) legs * 2 vertices each
    for m, n_ in ((4, 4), (3, 4), (5, 7)):
        spec = [(2,) * 3, (2,) * m, (2,) * n_]
        assert bouquet_size(spec) == 4 + 2 * (3 + m + n_)
        p = bouquet_poly(spec)
        assert is_unimodal(p), (m, n_)  # Li 2026 Theorems 1.4-1.5
    # 4. Galvin T_{m,t,1} reproduces the targeted.py generator
    from targeted import make_T_m_t_1
    for m, t in ((4, 3), (6, 6), (5, 2)):
        n, adj = make_T_m_t_1(m, t)
        slow = independence_poly(n, adj)
        fast = bouquet_poly([(2,) * t] * m)
        assert fast == slow, (m, t)
    print("selftest: all checks passed")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sweep", default="all",
                    help="A|B|C|D|all (comma-separated ok)")
    ap.add_argument("--max-n", type=int, default=320)
    ap.add_argument("--hillclimb", type=int, default=0,
                    help="generations of refinement on sweep champions")
    ap.add_argument("--pop", type=int, default=24)
    ap.add_argument("--top", type=int, default=40)
    ap.add_argument("--out", default="")
    ap.add_argument("--seed", type=int, default=993)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        selftest()
        return

    names = list("ABCD") if args.sweep == "all" else args.sweep.split(",")
    rng = random.Random(args.seed)
    seen: set = set()
    heap: list = []   # (ratio_key, counter, result) min-heap of top results
    counter = 0
    counterexamples = []
    t0 = time.time()
    tested = 0

    for name in names:
        gen = SWEEPS[name](args.max_n)
        t1 = time.time()
        count = 0
        for spec in gen:
            c = canon(spec)
            if c in seen:
                continue
            seen.add(c)
            res = evaluate_spec(spec)
            res["sweep"] = name
            tested += 1
            count += 1
            if res["witness"]:
                counterexamples.append(res)
                print(f"\n*** VALLEY WITNESS *** {res['label']} n={res['n']}")
                print(f"    poly={res['poly']}\n", flush=True)
            k = (res["window"], res["ratio"])
            counter += 1
            if len(heap) < args.top:
                heapq.heappush(heap, (k, counter, res))
            elif k > heap[0][0]:
                heapq.heapreplace(heap, (k, counter, res))
        print(f"sweep {name}: {count} specs, {time.time()-t1:.1f}s",
              flush=True)

    top = sorted((h[2] for h in heap),
                 key=lambda r: (r["window"], r["ratio"]), reverse=True)

    if args.hillclimb and top:
        print(f"\nhill climb: {args.hillclimb} generations, pop {args.pop}",
              flush=True)
        witness, population = hill_climb(top, args.max_n, args.hillclimb,
                                         args.pop, rng, seen)
        if witness:
            counterexamples.append(witness)
            print(f"\n*** VALLEY WITNESS (hill climb) *** {witness['label']}")
            print(f"    poly={witness['poly']}", flush=True)
        merged = {id(r): r for r in top + population}
        top = sorted(merged.values(),
                     key=lambda r: (r["window"], r["ratio"]),
                     reverse=True)[:args.top]

    print(f"\n{'='*74}")
    print(f"tested {tested} specs in {time.time()-t0:.1f}s; "
          f"counterexamples: {len(counterexamples)}")
    print(f"top valley ratios (window=True means rise precedes the "
          f"decreasing-tail threshold):")
    for r in top[:args.top]:
        print(f"  V={r['ratio']:.10f} win={str(r['window']):5s} "
              f"b={r['pos']:>3} rise={r['rise_pos']:>3} "
              f"n={r['n']:>4} a={r['alpha']:>4} [{r.get('sweep','hc')}] "
              f"{r['label']}")

    if args.out:
        payload = {
            "max_n": args.max_n,
            "sweeps": names,
            "tested": tested,
            "counterexamples": counterexamples,
            "top": [{k: v for k, v in r.items() if k != "poly"}
                    for r in top],
            "elapsed_s": round(time.time() - t0, 1),
        }
        with open(args.out, "w") as f:
            json.dump(payload, f, indent=2)
        print(f"saved {args.out}")


if __name__ == "__main__":
    main()
