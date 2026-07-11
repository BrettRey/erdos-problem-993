#!/usr/bin/env python3
"""Adversarial exact search for prefix extension-count variance bounds.

This is deliberately a scratch driver.  It scores a tree at ranks r satisfying

    r + 5 <= 2 (alpha-r)

by the exact ratio Var(e_r)/E[e_r].  No floating-point comparison is used for
selection.  A second, slower audit computes the worst (minimum-denominator)
maximum matching for the matching-block ABV bound.
"""

from __future__ import annotations

import argparse
import heapq
import json
import random
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

from graph6 import parse_graph6
from indpoly import _polymul, independence_poly
from scratch_extension_variance_dp_20260711 import extension_moment_jets


Adj = list[list[int]]


@dataclass(frozen=True)
class Score:
    num: int
    den: int
    rank: int
    alpha: int

    def ratio(self) -> Fraction:
        return Fraction(self.num, self.den)


def edges(adj: Adj) -> list[tuple[int, int]]:
    return [(u, v) for u, ns in enumerate(adj) for v in ns if u < v]


def from_edges(n: int, es: list[tuple[int, int]]) -> Adj:
    out = [[] for _ in range(n)]
    for u, v in es:
        out[u].append(v)
        out[v].append(u)
    return out


def prufer_decode(seq: list[int]) -> Adj:
    n = len(seq) + 2
    degree = [1] * n
    for v in seq:
        degree[v] += 1
    leaves = [v for v, d in enumerate(degree) if d == 1]
    heapq.heapify(leaves)
    es: list[tuple[int, int]] = []
    for v in seq:
        u = heapq.heappop(leaves)
        es.append((u, v))
        degree[u] -= 1
        degree[v] -= 1
        if degree[v] == 1:
            heapq.heappush(leaves, v)
    u = heapq.heappop(leaves)
    v = heapq.heappop(leaves)
    es.append((u, v))
    return from_edges(n, es)


def prufer_encode(adj: Adj) -> list[int]:
    n = len(adj)
    if n <= 2:
        return []
    degree = [len(ns) for ns in adj]
    leaves = [v for v, d in enumerate(degree) if d == 1]
    heapq.heapify(leaves)
    out: list[int] = []
    for _ in range(n - 2):
        u = heapq.heappop(leaves)
        v = next(v for v in adj[u] if degree[v] > 0)
        out.append(v)
        degree[u] = 0
        degree[v] -= 1
        if degree[v] == 1:
            heapq.heappush(leaves, v)
    return out


def best_score(adj: Adj, min_rank: int = 1) -> Score | None:
    jet = extension_moment_jets(adj)
    alpha = len(jet.value) - 1
    best: Score | None = None
    for r in range(min_rank, alpha + 1):
        if r + 5 > 2 * (alpha - r):
            break
        z, d1, d2 = jet.value[r], jet.first[r], jet.second[r]
        if d1 == 0:
            continue
        # Var(e)/E[e] = (z(d2+d1)-d1^2)/(z d1).
        score = Score(z * (d2 + d1) - d1 * d1, z * d1, r, alpha)
        if best is None or score.num * best.den > best.num * score.den:
            best = score
    return best


def score_at(adj: Adj, rank: int) -> Score:
    jet = extension_moment_jets(adj)
    alpha = len(jet.value) - 1
    z, d1, d2 = jet.value[rank], jet.first[rank], jet.second[rank]
    return Score(z * (d2 + d1) - d1 * d1, z * d1, rank, alpha)


def best_global_slope_score(
    adj: Adj, include_eta: bool, prefix_only: bool = False,
    min_rank: int = 0,
) -> Score | None:
    """Maximize Var(e)/(2 mu + 2 eta) or Var(e)/(2 mu), all ranks."""
    jet = extension_moment_jets(adj)
    alpha = len(jet.value) - 1
    best: Score | None = None
    for r, z in enumerate(jet.value):
        if r < min_rank:
            continue
        if prefix_only and r + 5 > 2 * (alpha - r):
            continue
        if r >= len(jet.first) or jet.first[r] == 0:
            continue
        d1 = jet.first[r]
        d2 = jet.second[r] if r < len(jet.second) else 0
        z2 = jet.value[r + 2] if r + 2 < len(jet.value) else 0
        hsum2 = d2 - (r + 2) * (r + 1) * z2
        assert hsum2 >= 0 and hsum2 % 2 == 0
        hsum = hsum2 // 2
        num = z * (d2 + d1) - d1 * d1
        den = 2 * z * (d1 + hsum if include_eta else d1)
        score = Score(num, den, r, alpha)
        if best is None or score.num * best.den > best.num * score.den:
            best = score
    return best


def best_mean_increment_score(
    adj: Adj, prefix_only: bool = False, min_rank: int = 0
) -> Score | None:
    """Maximize mu_{r+1}-mu_r, whose target upper bound is one."""
    p = independence_poly(len(adj), adj)
    alpha = len(p) - 1
    best: Score | None = None
    for r in range(min_rank, alpha - 1):
        if prefix_only and r + 5 > 2 * (alpha - r):
            continue
        if p[r] == 0 or p[r + 1] == 0:
            continue
        num = (r + 2) * p[r + 2] * p[r] - (r + 1) * p[r + 1] * p[r + 1]
        den = p[r + 1] * p[r]
        score = Score(num, den, r, alpha)
        if best is None or score.num * best.den > best.num * score.den:
            best = score
    return best


def component_poly(adj: Adj, keep: list[bool]) -> list[int]:
    """Independence polynomial of the induced forest on keep vertices."""
    n = len(adj)
    seen = [False] * n
    total = [1]
    for root in range(n):
        if not keep[root] or seen[root]:
            continue
        parent = {root: -1}
        order = [root]
        seen[root] = True
        for u in order:
            for v in adj[u]:
                if keep[v] and not seen[v]:
                    seen[v] = True
                    parent[v] = u
                    order.append(v)
        e0: dict[int, list[int]] = {}
        e1: dict[int, list[int]] = {}
        for u in reversed(order):
            a, b = [1], [0, 1]
            for v in adj[u]:
                if parent.get(v) != u:
                    continue
                child_total = [0] * max(len(e0[v]), len(e1[v]))
                for i, x in enumerate(e0[v]):
                    child_total[i] += x
                for i, x in enumerate(e1[v]):
                    child_total[i] += x
                while len(child_total) > 1 and child_total[-1] == 0:
                    child_total.pop()
                a = _polymul(a, child_total)
                b = _polymul(b, e0[v])
            e0[u], e1[u] = a, b
        comp = [0] * max(len(e0[root]), len(e1[root]))
        for i, x in enumerate(e0[root]):
            comp[i] += x
        for i, x in enumerate(e1[root]):
            comp[i] += x
        total = _polymul(total, comp)
    return total


def edge_free_weights(adj: Adj, rank: int) -> dict[tuple[int, int], int]:
    """w_uv=# rank-r independent sets leaving both endpoints addable."""
    n = len(adj)
    out: dict[tuple[int, int], int] = {}
    for u, v in edges(adj):
        banned = {u, v, *adj[u], *adj[v]}
        keep = [x not in banned for x in range(n)]
        p = component_poly(adj, keep)
        out[(u, v)] = p[rank] if rank < len(p) else 0
    return out


def extremal_maximum_matching(
    adj: Adj, weights: dict[tuple[int, int], int], minimize: bool
) -> tuple[int, int, list[tuple[int, int]]]:
    """Return (cardinality, weight, edges) for extremal max matching."""
    n = len(adj)
    parent = [-2] * n
    parent[0] = -1
    order = [0]
    for u in order:
        for v in adj[u]:
            if parent[v] == -2:
                parent[v] = u
                order.append(v)

    # state 0: u not matched to parent; state 1: u is matched to parent.
    dp: list[list[tuple[int, int, list[tuple[int, int]]]]] = [
        [(-10**9, 0, []), (-10**9, 0, [])] for _ in range(n)
    ]

    def better(a: tuple[int, int, list], b: tuple[int, int, list]) -> bool:
        if a[0] != b[0]:
            return a[0] > b[0]
        return a[1] < b[1] if minimize else a[1] > b[1]

    for u in reversed(order):
        children = [v for v in adj[u] if parent[v] == u]
        # If u is matched upward, no child may match u.
        card = sum(dp[v][0][0] for v in children)
        weight = sum(dp[v][0][1] for v in children)
        cert = sum((dp[v][0][2] for v in children), [])
        dp[u][1] = (card, weight, cert)
        # State 0: either unmatched, or matched to exactly one child.
        best = (card, weight, list(cert))
        for chosen in children:
            base_card = 1
            base_weight = weights[tuple(sorted((u, chosen)))]
            base_cert = [(u, chosen)]
            for v in children:
                state = 1 if v == chosen else 0
                base_card += dp[v][state][0]
                base_weight += dp[v][state][1]
                base_cert += dp[v][state][2]
            cand = (base_card, base_weight, base_cert)
            if better(cand, best):
                best = cand
        dp[u][0] = best
    return dp[0][0]


def abv_audit(adj: Adj, rank: int) -> dict:
    jet = extension_moment_jets(adj)
    z, d1, d2 = jet.value[rank], jet.first[rank], jet.second[rank]
    varnum = z * (d2 + d1) - d1 * d1
    weights = edge_free_weights(adj, rank)
    lo = extremal_maximum_matching(adj, weights, True)
    hi = extremal_maximum_matching(adj, weights, False)
    assert lo[0] == hi[0]
    # E sum Y_j^2 = (d1+2 sum_{matched edge} free_count)/z.
    lo_den = z * (d1 + 2 * lo[1])
    hi_den = z * (d1 + 2 * hi[1])
    return {
        "rank": rank,
        "z": z,
        "d1": d1,
        "d2": d2,
        "var_num": varnum,
        "matching_cardinality": lo[0],
        "min_bsum": lo[1],
        "max_bsum": hi[1],
        "worst_ratio": [varnum, lo_den],
        "best_ratio": [varnum, hi_den],
        "min_matching": lo[2],
        "max_matching": hi[2],
    }


def random_prufer(n: int, rng: random.Random, hubs: int | None = None) -> list[int]:
    if hubs is None:
        return [rng.randrange(n) for _ in range(n - 2)]
    centers = rng.sample(range(n), hubs)
    p = []
    for _ in range(n - 2):
        p.append(rng.choice(centers) if rng.random() < 0.82 else rng.randrange(n))
    return p


def mutate(seq: list[int], rng: random.Random, strength: int = 1) -> list[int]:
    out = list(seq)
    n = len(seq) + 2
    for _ in range(strength):
        mode = rng.randrange(4)
        if mode == 0:
            out[rng.randrange(len(out))] = rng.randrange(n)
        elif mode == 1:
            # Copy a code symbol to amplify/erode a hub.
            out[rng.randrange(len(out))] = out[rng.randrange(len(out))]
        elif mode == 2:
            a, b = rng.sample(range(len(out)), 2)
            out[a], out[b] = out[b], out[a]
        else:
            center = rng.randrange(n)
            for _ in range(rng.randint(1, max(1, len(out) // 20))):
                out[rng.randrange(len(out))] = center
    return out


def evolutionary(
    n: int, generations: int, popsize: int, seed: int, min_rank: int,
    global_target: str | None = None, global_prefix: bool = False,
) -> list[dict]:
    rng = random.Random(seed)
    pop: list[list[int]] = []
    # Exact star and a portfolio of uniform / hub-biased Prüfer seeds.
    pop.append([0] * (n - 2))
    # Seed the sharp perfect-matching lifts from D13.  Starting only from a
    # star makes the prefix-GSB objective artificially easy (its extension
    # count is deterministic beyond rank one), while these n=2*alpha trees
    # sit close to the actual three-coefficient boundary.
    if n % 2 == 0 and n >= 8:
        alpha = n // 2
        corona_edges = [(0, u) for u in range(1, alpha)]
        corona_edges += [(u, alpha + u) for u in range(alpha)]
        pop.append(prufer_encode(from_edges(n, corona_edges)))

        # Misaligned lift of a nearly balanced double star.  Hub blocks are
        # (0,1),(2,3); the bridge deliberately hits endpoint 3.
        p = (alpha - 2) // 2
        q = alpha - 2 - p
        lift_edges: list[tuple[int, int]] = [(0, 1), (2, 3), (0, 3)]
        next_vertex = 4
        for hub, count in ((0, p), (2, q)):
            for _ in range(count):
                base, leaf = next_vertex, next_vertex + 1
                next_vertex += 2
                lift_edges.extend(((base, leaf), (hub, base)))
        assert next_vertex == n
        pop.append(prufer_encode(from_edges(n, lift_edges)))
    # Seed the global search with padded exact ordered-LC witnesses so that
    # Prüfer mutation tests whether their deep-tail variance mechanism can be
    # amplified rather than relying on random discovery of that morphology.
    if global_target and n >= 26:
        try:
            rows = json.loads(Path("results/analysis_n26.json").read_text())["lc_failures"]
            for row in rows:
                _, base = parse_graph6(row["graph6"].encode())
                for mode in range(3):
                    es = edges(base)
                    nn = len(base)
                    last = max(range(nn), key=lambda u: len(base[u]))
                    while nn < n:
                        if mode == 0:
                            anchor = last
                        elif mode == 1:
                            anchor = rng.randrange(nn)
                        else:
                            anchor = max(range(nn), key=lambda u: sum(1 for x, _ in es if x == u) + sum(1 for _, x in es if x == u))
                        es.append((anchor, nn))
                        last = nn
                        nn += 1
                    pop.append(prufer_encode(from_edges(n, es)))
        except (OSError, KeyError, ValueError):
            pass
    for i in range(popsize - len(pop)):
        hubs = None if i % 3 == 0 else rng.randint(2, min(12, n))
        pop.append(random_prufer(n, rng, hubs))

    hall: list[tuple[Score, list[int]]] = []
    for gen in range(generations):
        scored: list[tuple[Score, list[int]]] = []
        for seq in pop:
            adj = prufer_decode(seq)
            if global_target == "D":
                score = best_mean_increment_score(adj, global_prefix, min_rank)
            elif global_target:
                score = best_global_slope_score(
                    adj, global_target == "A", global_prefix, min_rank
                )
            else:
                score = best_score(adj, min_rank)
            if score is not None:
                scored.append((score, seq))
        scored.sort(key=lambda x: x[0].ratio(), reverse=True)
        hall.extend(scored[:8])
        hall.sort(key=lambda x: x[0].ratio(), reverse=True)
        hall = hall[:40]
        best = scored[0][0]
        print(json.dumps({
            "generation": gen, "n": n, "rank": best.rank,
            "alpha": best.alpha, "ratio": float(best.ratio()),
            "fraction": [best.num, best.den],
        }), flush=True)
        elites = [seq for _, seq in scored[: max(4, popsize // 6)]]
        pop = list(elites)
        while len(pop) < popsize:
            parent = rng.choice(elites)
            pop.append(mutate(parent, rng, 1 if rng.random() < 0.75 else rng.randint(2, 6)))

    # Deduplicate by edge certificate (label-dependent is okay for audit).
    seen: set[tuple[tuple[int, int], ...]] = set()
    out = []
    for score, seq in hall:
        adj = prufer_decode(seq)
        cert = tuple(edges(adj))
        if cert in seen:
            continue
        seen.add(cert)
        out.append({
            "n": n, "rank": score.rank, "alpha": score.alpha,
            "ratio": [score.num, score.den], "edges": list(cert),
            "degree_sequence": sorted((len(ns) for ns in adj), reverse=True),
        })
    return out


def pad_tree(adj: Adj, kind: str, amount: int, anchor: int = 0) -> Adj:
    es = edges(adj)
    n = len(adj)
    if kind == "leaves":
        for _ in range(amount):
            es.append((anchor, n))
            n += 1
    elif kind == "path":
        last = anchor
        for _ in range(amount):
            es.append((last, n))
            last = n
            n += 1
    elif kind == "corona":
        old_n = n
        for u in range(old_n):
            for _ in range(amount):
                es.append((u, n))
                n += 1
    else:
        raise ValueError(kind)
    return from_edges(n, es)


def add_path_arm(es: list[tuple[int, int]], n: int, root: int, length: int) -> int:
    last = root
    for _ in range(length):
        es.append((last, n))
        last = n
        n += 1
    return n


def decorated_grammar_tree(rng: random.Random, max_n: int) -> Adj:
    """Random multi-hub core with arms and repeated rooted substitutions."""
    h = rng.randint(1, min(14, max_n))
    if h == 1:
        core: list[tuple[int, int]] = []
    elif rng.random() < 0.4:
        core = [(i - 1, i) for i in range(1, h)]
    elif rng.random() < 0.5:
        core = [(0, i) for i in range(1, h)]
    else:
        core = edges(prufer_decode(random_prufer(h, rng, rng.randint(2, min(h, 5)))))
    es = list(core)
    n = h
    # Heavy-tailed arm loads produce stars, subdivided stars, and multi-hubs.
    for u in range(h):
        load = rng.choice((0, 0, 1, 2, 3, 5, 8, 13, 21, 34))
        for _ in range(load):
            if n >= max_n:
                break
            length = rng.choices((1, 2, 3, 4), weights=(7, 6, 3, 1))[0]
            length = min(length, max_n - n)
            n = add_path_arm(es, n, u, length)
        # Rooted substitution: a child carrying t pendant P2 paths.
        copies = rng.choice((0, 0, 0, 1, 1, 2, 3))
        for _ in range(copies):
            t = rng.randint(1, 8)
            need = 1 + 2 * t
            if n + need > max_n:
                break
            branch = n
            es.append((u, branch))
            n += 1
            for _ in range(t):
                n = add_path_arm(es, n, branch, 2)
    return from_edges(n, es)


def grammar_search(samples: int, max_n: int, seed: int) -> list[dict]:
    rng = random.Random(seed)
    bands = (1, 2, 3, 5, 10, 20, 40)
    hall: dict[int, list[tuple[Score, Adj]]] = {b: [] for b in bands}
    failures: list[dict] = []
    for sample in range(samples):
        adj = decorated_grammar_tree(rng, max_n)
        jet = extension_moment_jets(adj)
        alpha = len(jet.value) - 1
        for r in range(1, alpha + 1):
            if r + 5 > 2 * (alpha - r):
                break
            if r >= len(jet.first) or jet.first[r] == 0:
                continue
            z, d1, d2 = jet.value[r], jet.first[r], jet.second[r]
            score = Score(z * (d2 + d1) - d1 * d1, z * d1, r, alpha)
            if score.num > score.den:
                failures.append({
                    "sample": sample, "n": len(adj), "rank": r, "alpha": alpha,
                    "ratio": [score.num, score.den], "edges": edges(adj),
                })
                return failures
            for band in bands:
                if r < band:
                    continue
                rows = hall[band]
                rows.append((score, adj))
                rows.sort(key=lambda x: x[0].ratio(), reverse=True)
                del rows[8:]
        if (sample + 1) % 500 == 0:
            print(json.dumps({
                "sample": sample + 1,
                "best": {str(b): float(hall[b][0][0].ratio()) if hall[b] else None for b in bands},
            }), flush=True)
    out = []
    for band in bands:
        for score, adj in hall[band]:
            out.append({
                "band": band, "n": len(adj), "rank": score.rank,
                "alpha": score.alpha, "ratio": [score.num, score.den],
                "degree_sequence": sorted((len(ns) for ns in adj), reverse=True),
                "edges": edges(adj),
            })
    out.sort(key=lambda x: (x["band"], -float(Fraction(*x["ratio"]))))
    return out


def audit_known_padding(max_amount: int) -> list[dict]:
    data = json.loads(Path("results/analysis_n26.json").read_text())
    out = []
    for row in data["lc_failures"]:
        n, base = parse_graph6(row["graph6"].encode())
        assert n == 26
        anchors = sorted(range(n), key=lambda u: len(base[u]), reverse=True)[:5]
        for kind in ("leaves", "path"):
            for amount in range(max_amount + 1):
                for anchor in anchors:
                    adj = pad_tree(base, kind, amount, anchor)
                    s = best_score(adj)
                    if s is None:
                        continue
                    out.append({
                        "source": row["graph6"], "kind": kind,
                        "amount": amount, "anchor": anchor, "n": len(adj),
                        "rank": s.rank, "alpha": s.alpha,
                        "ratio": [s.num, s.den], "edges": edges(adj),
                    })
    out.sort(key=lambda x: Fraction(*x["ratio"]), reverse=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=80)
    ap.add_argument("--generations", type=int, default=30)
    ap.add_argument("--popsize", type=int, default=100)
    ap.add_argument("--seed", type=int, default=993)
    ap.add_argument("--min-rank", type=int, default=2)
    ap.add_argument("--padding", type=int, default=-1)
    ap.add_argument("--grammar-samples", type=int, default=-1)
    ap.add_argument("--max-n", type=int, default=300)
    ap.add_argument("--global-target", choices=("A", "B", "D"))
    ap.add_argument("--global-prefix", action="store_true")
    ap.add_argument("--out", default="/tmp/variance_search.json")
    args = ap.parse_args()
    if args.padding >= 0:
        rows = audit_known_padding(args.padding)
    elif args.grammar_samples >= 0:
        rows = grammar_search(args.grammar_samples, args.max_n, args.seed)
    else:
        rows = evolutionary(
            args.n, args.generations, args.popsize, args.seed, args.min_rank,
            args.global_target, args.global_prefix,
        )
    Path(args.out).write_text(json.dumps(rows, indent=2) + "\n")
    if rows:
        # The grammar output is grouped by rank band, so select the global
        # exact maximum for the final expensive ABV audit.
        best = max(rows, key=lambda x: Fraction(*x["ratio"]))
        print(json.dumps({k: v for k, v in best.items() if k != "edges"}, indent=2))
        adj = from_edges(best["n"], [tuple(e) for e in best["edges"]])
        print(json.dumps({"abv": abv_audit(adj, best["rank"])}))


if __name__ == "__main__":
    main()
