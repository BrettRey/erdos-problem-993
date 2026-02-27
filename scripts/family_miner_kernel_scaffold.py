#!/usr/bin/env python3
"""family_miner_kernel_scaffold.py

Mining pipeline for kernel+scaffold families at a canonical root u.

Model:
  P(x) = \\prod_i F_i(x)^{k_i}
  Q(x) = x \\prod_i G_i(x)^{k_i}
  I(x) = (1+2x)P(x) + (1+x)Q(x)

We mine:
  - kernel pairs: two rooted types (F,G) with identical ratio r(x)=G(x)/F(x)
    (symbolic equality; optional sampled-lambda equality)
    but different [x]F (= component size).
  - scaffold knobs: two additional rooted types.

For each candidate family (kernel pair + 2 scaffolds), scan a bounded parameter box
for same-(m,lambda) and same-(m,lambda,rho) collisions, where
  m = leftmost mode index of I,
  lambda = i_{m-1}/i_m,
  rho = Q(lambda)/P(lambda).

Outputs ranked JSON diagnostics.

This script is designed to be drop-in for repo usage but is self-contained.
It does NOT depend on project-specific enumerators.

Exact rationals only: all lambda/rho are Fractions; no floats.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, List, Tuple, Optional, Iterable, Any

import networkx as nx
from networkx.algorithms.isomorphism import is_isomorphic

# No floating computations are used anywhere.
# We avoid heavy symbolic normalisation; kernel equality is certified by
# exact polynomial cross-multiplication over Z.


def poly_add(a: List[int], b: List[int]) -> List[int]:
    n = max(len(a), len(b))
    out = [0]*n
    for i in range(n):
        if i < len(a):
            out[i] += a[i]
        if i < len(b):
            out[i] += b[i]
    # trim
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_mul(a: List[int], b: List[int]) -> List[int]:
    if a == [0] or b == [0]:
        return [0]
    out = [0]*(len(a)+len(b)-1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj:
                out[i+j] += ai*bj
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_pow(a: List[int], k: int) -> List[int]:
    # exponentiation by squaring
    if k == 0:
        return [1]
    if k == 1:
        return a
    res = [1]
    base = a
    e = k
    while e > 0:
        if e & 1:
            res = poly_mul(res, base)
        e >>= 1
        if e:
            base = poly_mul(base, base)
    return res


def poly_eval_frac(a: List[int], lam: Fraction) -> Fraction:
    # Horner
    res = Fraction(0, 1)
    for coeff in reversed(a):
        res = res*lam + Fraction(coeff, 1)
    return res


def leftmost_mode(coeffs: List[int]) -> int:
    mx = max(coeffs)
    for i, v in enumerate(coeffs):
        if v == mx:
            return i
    raise RuntimeError("empty coeff list")


def frac_to_numden(q: Fraction) -> Dict[str, int]:
    return {"num": int(q.numerator), "den": int(q.denominator)}


def sign_mix_score(int_list: List[int]) -> Tuple[bool, int]:
    # returns (mixed?, sign_changes ignoring zeros)
    s = [0 if v == 0 else (1 if v > 0 else -1) for v in int_list]
    nonzero = [v for v in s if v != 0]
    if not nonzero:
        return (False, 0)
    mixed = (min(nonzero) < 0) and (max(nonzero) > 0)
    changes = 0
    for i in range(1, len(nonzero)):
        if nonzero[i] != nonzero[i-1]:
            changes += 1
    return mixed, changes


@dataclass(frozen=True)
class RootedType:
    # rooted component type
    n: int
    # root degree in component
    deg_root: int
    # independence polynomial of component
    F: Tuple[int, ...]
    # independence polynomial of component with root removed (= dp0[root])
    G: Tuple[int, ...]
    # sampled ratio signature r(λ)=G(λ)/F(λ) at fixed rational λ samples
    r_sig: Tuple[Tuple[int, int], ...]
    # diagnostic hash label
    key: str


def independence_poly_tree(G: nx.Graph, root: int) -> Tuple[List[int], List[int]]:
    """Return (dp0[root], dp1[root]) as coefficient lists for a rooted tree."""
    parent = {root: -1}
    order = [root]
    for v in order:
        for w in G.neighbors(v):
            if w == parent[v]:
                continue
            parent[w] = v
            order.append(w)
    # postorder
    order.reverse()

    dp0: Dict[int, List[int]] = {}
    dp1: Dict[int, List[int]] = {}

    for v in order:
        children = [w for w in G.neighbors(v) if parent.get(w, -2) == v]
        # dp0[v] = prod (dp0[c]+dp1[c])
        p0 = [1]
        for c in children:
            p0 = poly_mul(p0, poly_add(dp0[c], dp1[c]))
        # dp1[v] = x * prod dp0[c]
        p1 = [0, 1]  # x
        for c in children:
            p1 = poly_mul(p1, dp0[c])
        dp0[v] = p0
        dp1[v] = p1

    return dp0[root], dp1[root]


def attachable_dleaf_le_1(G: nx.Graph, root: int) -> bool:
    """Check whether attaching an extra neighbor to root (making it non-leaf unless n=1)
    yields a tree in which every vertex in the component has <=1 leaf neighbor.

    This is a sufficient per-component check for membership in the global gated class.
    """
    n = G.number_of_nodes()
    # leaf status after attachment
    deg_after = {v: G.degree(v) for v in G.nodes()}
    deg_after[root] += 1

    leaf_after = {v: (deg_after[v] == 1) for v in G.nodes()}

    # count leaf neighbors within component
    for v in G.nodes():
        cnt = 0
        for w in G.neighbors(v):
            if leaf_after[w]:
                cnt += 1
                if cnt > 1:
                    return False
    return True


def ratio_signature(F: Tuple[int, ...], G: Tuple[int, ...], samples: List[Fraction]) -> Tuple[Tuple[int, int], ...]:
    """Return a hashable sampled signature for r(λ)=G(λ)/F(λ) over a fixed list of rational samples."""
    out: List[Tuple[int, int]] = []
    for lam in samples:
        Fv = poly_eval_frac(list(F), lam)
        Gv = poly_eval_frac(list(G), lam)
        # Fv>0 for λ>0 on independence polynomials.
        r = Gv / Fv
        out.append((int(r.numerator), int(r.denominator)))
    return tuple(out)


def equal_ratio_symbolic(F1: Tuple[int, ...], G1: Tuple[int, ...], F2: Tuple[int, ...], G2: Tuple[int, ...]) -> bool:
    """Certify equality G1/F1 == G2/F2 as rational functions by cross-multiplication over Z[x]."""
    return poly_mul(list(G1), list(F2)) == poly_mul(list(G2), list(F1))


DEFAULT_RATIO_SAMPLES = [Fraction(1, 2), Fraction(2, 3), Fraction(3, 4), Fraction(4, 5)]


def rooted_type_library(max_n: int, ratio_samples: Optional[List[Fraction]] = None) -> List[RootedType]:
    """Generate attachable rooted tree types up to size max_n."""
    types: List[RootedType] = []
    seen_hash = set()

    samples = ratio_samples if ratio_samples is not None else DEFAULT_RATIO_SAMPLES

    for n in range(1, max_n+1):
        # networkx.nonisomorphic_trees is defined for n>=2.
        if n == 1:
            Tlist = [nx.Graph()]
            Tlist[0].add_node(0)
        else:
            Tlist = list(nx.generators.nonisomorphic_trees(n))

        for T in Tlist:
            # relabel to 0..n-1 (already)
            for root in range(n):
                if not attachable_dleaf_le_1(T, root):
                    continue

                # rooted WL hash to deduplicate
                H = T.copy()
                nx.set_node_attributes(H, {v: (1 if v == root else 0) for v in H.nodes()}, "root")
                h = nx.weisfeiler_lehman_graph_hash(H, node_attr="root")
                sig = (n, h)
                if sig in seen_hash:
                    continue
                seen_hash.add(sig)

                dp0, dp1 = independence_poly_tree(T, root)
                F = tuple(poly_add(dp0, dp1))
                G = tuple(dp0)  # dp0[root] = I(T - root)
                r_sig = ratio_signature(F, G, samples)

                rt = RootedType(
                    n=n,
                    deg_root=int(T.degree(root)),
                    F=F,
                    G=G,
                    r_sig=r_sig,
                    key=f"n{n}_deg{int(T.degree(root))}_wl{h[:10]}"
                )
                types.append(rt)

    return types


@dataclass
class CandidateFamily:
    kernel_ratio_key: str
    K1: RootedType
    K2: RootedType
    S1: RootedType
    S2: RootedType


def build_polys_for_params(fam: CandidateFamily, a: int, b: int, s1: int, s2: int) -> Tuple[List[int], List[int], List[int]]:
    # P = F1^a F2^b Fs1^s1 Fs2^s2
    P = [1]
    for (typ, k) in [(fam.K1, a), (fam.K2, b), (fam.S1, s1), (fam.S2, s2)]:
        if k:
            P = poly_mul(P, list(poly_pow(list(typ.F), k)))
    # Q = x * G1^a G2^b Gs1^s1 Gs2^s2
    Q = [0, 1]
    for (typ, k) in [(fam.K1, a), (fam.K2, b), (fam.S1, s1), (fam.S2, s2)]:
        if k:
            Q = poly_mul(Q, list(poly_pow(list(typ.G), k)))

    # I = (1+2x)P + (1+x)Q
    one_2x = [1, 2]
    one_x = [1, 1]
    I = poly_add(poly_mul(one_2x, P), poly_mul(one_x, Q))
    return P, Q, I


def invariants_from_PQI(P: List[int], Q: List[int], I: List[int], m_min: int) -> Optional[Dict[str, Any]]:
    m = leftmost_mode(I)
    if m < m_min or m == 0:
        return None
    lam = Fraction(I[m-1], I[m])
    # rho = Q(lam)/P(lam)
    Pval = poly_eval_frac(P, lam)
    Qval = poly_eval_frac(Q, lam)
    if Pval == 0:
        return None
    rho = Qval / Pval
    return {
        "m": m,
        "lambda": lam,
        "rho": rho,
    }


def mine_families(types: List[RootedType],
                  m_min: int,
                  max_kernel_total: int,
                  max_scaffold: int,
                  max_families: int,
                  require_adjacent_kernel: bool) -> List[Dict[str, Any]]:

    # group by sampled ratio signature
    groups: Dict[Tuple[Tuple[int, int], ...], List[RootedType]] = {}
    for t in types:
        groups.setdefault(t.r_sig, []).append(t)

    # kernel candidate pairs
    kernel_pairs: List[Tuple[Tuple[Tuple[int, int], ...], RootedType, RootedType]] = []
    for sig, lst in groups.items():
        # need >=2 and different sizes, and certify true ratio equality by cross-multiplication
        lst2 = sorted(lst, key=lambda z: z.n)
        for i in range(len(lst2)):
            for j in range(i+1, len(lst2)):
                if lst2[i].n == lst2[j].n:
                    continue
                if require_adjacent_kernel and abs(lst2[i].n - lst2[j].n) != 1:
                    continue
                if not equal_ratio_symbolic(lst2[i].F, lst2[i].G, lst2[j].F, lst2[j].G):
                    continue
                kernel_pairs.append((sig, lst2[i], lst2[j]))

    # scaffold pool: prefer small n
    scaffold_pool = sorted(types, key=lambda t: (t.n, t.deg_root))

    candidates: List[CandidateFamily] = []
    for rkey, k1, k2 in kernel_pairs:
        # pick two scaffold types with different ratio keys (to give knobs)
        # avoid using kernel-ratio scaffolds
        for i in range(min(len(scaffold_pool), 60)):
            s1 = scaffold_pool[i]
            if s1.r_sig == rkey:
                continue
            for j in range(i+1, min(len(scaffold_pool), 80)):
                s2 = scaffold_pool[j]
                if s2.r_sig == rkey:
                    continue
                if s2.r_sig == s1.r_sig:
                    continue
                candidates.append(CandidateFamily(
                    kernel_ratio_key=f"sig={rkey}",
                    K1=k1,
                    K2=k2,
                    S1=s1,
                    S2=s2,
                ))
                if len(candidates) >= max_families:
                    return evaluate_candidates(candidates, m_min, max_kernel_total, max_scaffold)

    return evaluate_candidates(candidates, m_min, max_kernel_total, max_scaffold)


def evaluate_candidates(cands: List[CandidateFamily], m_min: int, max_kernel_total: int, max_scaffold: int) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []

    for fam in cands:
        # scan parameter box
        # kernel total tK in 1..max_kernel_total
        # distribution a=0..tK, b=tK-a
        # scaffolds s1,s2 in 0..max_scaffold
        hit_mlr = None
        hit_ml = None
        hit_ml_kernel_line = None
        hit_mlr_kernel_line = None

        # store maps for collision detection
        # value is a compact witness record
        seen_ml: Dict[Tuple[int, Fraction], Dict[str, Any]] = {}
        seen_mlr: Dict[Tuple[int, Fraction, Fraction], Dict[str, Any]] = {}

        # per-line: detect collisions specifically along the kernel direction
        # (fixed tK, fixed scaffold counts, varying mixture a/b)

        # diagnostics
        min_lambda_gap: Optional[Fraction] = None
        best_signmix = (False, 0)

        for tK in range(1, max_kernel_total+1):
            for s1c in range(0, max_scaffold+1):
                for s2c in range(0, max_scaffold+1):
                    line_seen_ml: Dict[Tuple[int, Fraction], Dict[str, Any]] = {}
                    line_seen_mlr: Dict[Tuple[int, Fraction, Fraction], Dict[str, Any]] = {}
                    # along kernel line
                    lambdas_line: List[Fraction] = []
                    invs_line: List[Tuple[int, Fraction, List[int]]] = []
                    for a in range(0, tK+1):
                        b = tK - a
                        P, Q, I = build_polys_for_params(fam, a, b, s1c, s2c)
                        inv = invariants_from_PQI(P, Q, I, m_min)
                        if inv is None:
                            continue
                        m = inv["m"]
                        lam = inv["lambda"]
                        rho = inv["rho"]
                        N = a*fam.K1.n + b*fam.K2.n + s1c*fam.S1.n + s2c*fam.S2.n

                        witness = {
                            "params": {"a": a, "b": b, "s1": s1c, "s2": s2c, "tK": tK},
                            "N": N,
                            "m": m,
                            "lambda": frac_to_numden(lam),
                            "rho": frac_to_numden(rho),
                            "P": P,
                            "Q": Q,
                            "I": I,
                        }

                        # global collision checks
                        key_ml = (m, lam)
                        key_mlr = (m, lam, rho)

                        # kernel-line collision checks (same fixed tK,s1c,s2c)
                        if hit_ml_kernel_line is None and key_ml in line_seen_ml and line_seen_ml[key_ml]["N"] != N:
                            hit_ml_kernel_line = {
                                "line": {"tK": tK, "s1": s1c, "s2": s2c},
                                "key": (m, frac_to_numden(lam)),
                                "A": line_seen_ml[key_ml],
                                "B": witness,
                            }
                        else:
                            line_seen_ml.setdefault(key_ml, witness)

                        if hit_mlr_kernel_line is None and key_mlr in line_seen_mlr and line_seen_mlr[key_mlr]["N"] != N:
                            hit_mlr_kernel_line = {
                                "line": {"tK": tK, "s1": s1c, "s2": s2c},
                                "key": (m, frac_to_numden(lam), frac_to_numden(rho)),
                                "A": line_seen_mlr[key_mlr],
                                "B": witness,
                            }
                        else:
                            line_seen_mlr.setdefault(key_mlr, witness)

                        if key_ml in seen_ml and seen_ml[key_ml]["N"] != N and hit_ml is None:
                            hit_ml = {
                                "key": (m, frac_to_numden(lam)),
                                "A": seen_ml[key_ml],
                                "B": witness,
                            }
                        else:
                            seen_ml.setdefault(key_ml, witness)

                        if key_mlr in seen_mlr and seen_mlr[key_mlr]["N"] != N and hit_mlr is None:
                            hit_mlr = {
                                "key": (m, frac_to_numden(lam), frac_to_numden(rho)),
                                "A": seen_mlr[key_mlr],
                                "B": witness,
                            }
                        else:
                            seen_mlr.setdefault(key_mlr, witness)

                        lambdas_line.append(lam)
                        invs_line.append((m, lam, I))

                    # line diagnostics
                    if len(lambdas_line) >= 2:
                        # minimal gap between consecutive lambdas (in sorted order)
                        srt = sorted(lambdas_line)
                        for u, v in zip(srt, srt[1:]):
                            gap = v - u
                            if gap > 0 and (min_lambda_gap is None or gap < min_lambda_gap):
                                min_lambda_gap = gap

                        # signmix between two extremes if both exist
                        # choose first two invs on line with same m if possible
                        invs_by_m: Dict[int, List[Tuple[Fraction, List[int]]]] = {}
                        for m, lam, I in invs_line:
                            invs_by_m.setdefault(m, []).append((lam, I))
                        for m0, lst in invs_by_m.items():
                            if len(lst) >= 2:
                                # take two with farthest lambdas
                                lst2 = sorted(lst, key=lambda t: t[0])
                                I1 = lst2[0][1]
                                I2 = lst2[-1][1]
                                # pad
                                L = max(len(I1), len(I2))
                                D = [ (I2[i] if i < len(I2) else 0) - (I1[i] if i < len(I1) else 0) for i in range(L) ]
                                sm = sign_mix_score(D)
                                if sm[0] and (not best_signmix[0] or sm[1] > best_signmix[1]):
                                    best_signmix = sm

        # score / pack
        score = 0
        if hit_mlr is not None:
            score += 10_000
        if hit_ml is not None:
            score += 1_000
        if hit_mlr_kernel_line is not None:
            score += 5_000
        if hit_ml_kernel_line is not None:
            score += 500
        if min_lambda_gap is not None:
            # smaller gap => higher score
            # compare by denominator size: use 1/gap as heuristic integer
            # since Fractions, compute floor(1/gap)
            invgap = (Fraction(1, 1) / min_lambda_gap)
            score += min(int(invgap), 500)

        if best_signmix[0]:
            score += 50 + best_signmix[1]

        out.append({
            "score": score,
            "kernel": {
                "ratio_sig": [list(p) for p in fam.K1.r_sig],
                "ratio_repr": {"num": list(fam.K1.G), "den": list(fam.K1.F)},
                "K1": {"key": fam.K1.key, "n": fam.K1.n, "deg_root": fam.K1.deg_root, "F": list(fam.K1.F), "G": list(fam.K1.G)},
                "K2": {"key": fam.K2.key, "n": fam.K2.n, "deg_root": fam.K2.deg_root, "F": list(fam.K2.F), "G": list(fam.K2.G)},
                "delta_n": fam.K2.n - fam.K1.n,
            },
            "scaffold": {
                "S1": {"key": fam.S1.key, "n": fam.S1.n, "deg_root": fam.S1.deg_root, "F": list(fam.S1.F), "G": list(fam.S1.G)},
                "S2": {"key": fam.S2.key, "n": fam.S2.n, "deg_root": fam.S2.deg_root, "F": list(fam.S2.F), "G": list(fam.S2.G)},
            },
            "scan": {
                "m_min": m_min,
                "max_kernel_total": max_kernel_total,
                "max_scaffold": max_scaffold,
            },
            "hits": {
                "ml_collision": hit_ml,
                "mlr_collision": hit_mlr,
                "ml_collision_kernel_line": hit_ml_kernel_line,
                "mlr_collision_kernel_line": hit_mlr_kernel_line,
            },
            "diagnostics": {
                "min_lambda_gap": None if min_lambda_gap is None else frac_to_numden(min_lambda_gap),
                "sign_mixed": bool(best_signmix[0]),
                "sign_changes": int(best_signmix[1]),
            }
        })

    # sort by score desc
    out.sort(key=lambda d: d["score"], reverse=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-comp-n", type=int, default=10)
    ap.add_argument("--m-min", type=int, default=4)
    ap.add_argument("--max-kernel-total", type=int, default=6)
    ap.add_argument("--max-scaffold", type=int, default=3)
    ap.add_argument("--max-families", type=int, default=200)
    ap.add_argument("--require-adjacent-kernel", action="store_true")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    t0 = time.perf_counter()

    types = rooted_type_library(args.max_comp_n)

    ranked = mine_families(
        types=types,
        m_min=args.m_min,
        max_kernel_total=args.max_kernel_total,
        max_scaffold=args.max_scaffold,
        max_families=args.max_families,
        require_adjacent_kernel=args.require_adjacent_kernel,
    )

    runtime_s = time.perf_counter() - t0

    out = {
        "params": {
            "max_comp_n": args.max_comp_n,
            "m_min": args.m_min,
            "max_kernel_total": args.max_kernel_total,
            "max_scaffold": args.max_scaffold,
            "max_families": args.max_families,
            "require_adjacent_kernel": bool(args.require_adjacent_kernel),
        },
        "library": {
            "rooted_types": len(types),
        },
        "results": {
            "candidates": len(ranked),
            "any_mlr_split_in_box": any(c["hits"]["mlr_collision"] is not None for c in ranked),
        },
        "top": ranked[:10],
        "all": ranked,
        "runtime_s": runtime_s,
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
