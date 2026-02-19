#!/usr/bin/env python3
"""Diagnostic decomposition of the degree-2 bridge for focused tie-fugacity.

For each d_leaf<=1 tree T, choose a leaf l whose support s has degree 2.
Let u be the other neighbour of s (so s is adjacent to l and u).

Define:
  A = T - l
  B = T - {l, s}

When deg(s) = 2, the IS polynomials factor through the hub polynomials
P = dp_B[u][0] (u excluded) and Q = dp_B[u][1] (u included), giving:

  I(T) = (1+2x)P + (1+x)Q
  I(A) = (1+x)P + Q
  I(B) = P + Q

The Phi functional decomposes as:

  Phi_m(T; lam) = (1+lam)*Phi_m(B; lam) + lam*Z_B(lam) + lam*Phi_{m-1}(P; lam)

where the middle term lam*Z_B is the "pendant bonus" (always positive).

This script records, for each tree:
1. P, Q polynomial metadata (degree, LC status, mode)
2. Three-term decomposition values at lambda = lambda_m^T
3. STRONG C2 margin: lambda_m^T - lambda_{m-1}^B
4. Cross-tree comparison: lambda_m^A vs lambda_{m-1}^B
5. Per-term signs and ratios
6. Tightest witnesses with full structural data
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from collections import Counter
from typing import Any

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from graph6 import parse_graph6
from indpoly import independence_poly, is_log_concave, is_unimodal, _polyadd, _polymul


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def mode_index_leftmost(poly: list[int]) -> int:
    return max(range(len(poly)), key=lambda i: poly[i])


def remove_vertices(adj: list[list[int]], remove_set: set[int]) -> list[list[int]]:
    keep = [v for v in range(len(adj)) if v not in remove_set]
    idx = {v: i for i, v in enumerate(keep)}
    out = [[] for _ in keep]
    for v in keep:
        vv = idx[v]
        for u in adj[v]:
            if u in idx:
                out[vv].append(idx[u])
    return out


def phi_q(poly: list[int], q: int, lam: float) -> float:
    """Phi_q(poly; lam) = Z(lam) * (mu(lam) - (q-1))."""
    z = 0.0
    mu_num = 0.0
    p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p
        z += w
        mu_num += k * w
        p *= lam
    if z == 0:
        return 0.0
    mu = mu_num / z
    return z * (mu - (q - 1))


def eval_poly(poly: list[int], lam: float) -> float:
    """Evaluate polynomial at lam."""
    val = 0.0
    p = 1.0
    for ck in poly:
        val += ck * p
        p *= lam
    return val


def weighted_mean(poly: list[int], lam: float) -> float:
    """Weighted mean mu(lam) = sum(k * c_k * lam^k) / sum(c_k * lam^k)."""
    z = 0.0
    mu_num = 0.0
    p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p
        z += w
        mu_num += k * w
        p *= lam
    return mu_num / z if z != 0 else 0.0


def compute_hub_polys(adj_B: list[list[int]], u_in_B: int) -> tuple[list[int], list[int]]:
    """Compute P = dp[u][0] and Q = dp[u][1] by rooting B at u.

    Returns (P, Q) as coefficient lists.
    """
    n = len(adj_B)
    if n == 0:
        return [1], []
    if n == 1:
        return [1], [0, 1]

    # BFS to build tree rooted at u_in_B
    parent = [-1] * n
    children = [[] for _ in range(n)]
    visited = [False] * n
    visited[u_in_B] = True
    bfs_queue = [u_in_B]
    head = 0
    while head < len(bfs_queue):
        v = bfs_queue[head]
        head += 1
        for w in adj_B[v]:
            if not visited[w]:
                visited[w] = True
                parent[w] = v
                children[v].append(w)
                bfs_queue.append(w)

    # Post-order traversal
    order = []
    stack = [(u_in_B, False)]
    while stack:
        v, processed = stack.pop()
        if processed:
            order.append(v)
            continue
        stack.append((v, True))
        for c in children[v]:
            stack.append((c, False))

    dp0 = [[] for _ in range(n)]
    dp1 = [[] for _ in range(n)]

    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
        else:
            prod = [1]
            for c in children[v]:
                summand = _polyadd(dp0[c], dp1[c])
                prod = _polymul(prod, summand)
            dp0[v] = prod

            prod = [1]
            for c in children[v]:
                prod = _polymul(prod, dp0[c])
            dp1[v] = [0] + prod

    return dp0[u_in_B], dp1[u_in_B]


def verify_PQ_identity(poly_T: list[int], poly_A: list[int], poly_B: list[int],
                        P: list[int], Q: list[int]) -> tuple[bool, bool, bool]:
    """Verify the three identities: I(T)=(1+2x)P+(1+x)Q, I(A)=(1+x)P+Q, I(B)=P+Q."""
    # I(B) = P + Q
    check_B = _polyadd(P, Q)
    ok_B = check_B == poly_B or (len(check_B) >= len(poly_B) and
           all(check_B[i] == (poly_B[i] if i < len(poly_B) else 0) for i in range(len(check_B))))

    # I(A) = (1+x)P + Q = P + xP + Q
    xP = [0] + P  # multiply by x
    term1 = _polyadd(P, xP)  # (1+x)P
    check_A = _polyadd(term1, Q)
    ok_A = all(check_A[i] == (poly_A[i] if i < len(poly_A) else 0)
               for i in range(max(len(check_A), len(poly_A))))

    # I(T) = (1+2x)P + (1+x)Q = P + 2xP + Q + xQ
    twoxP = [0] + [2 * c for c in P]
    xQ = [0] + Q
    check_T = _polyadd(_polyadd(P, twoxP), _polyadd(Q, xQ))
    ok_T = all(check_T[i] == (poly_T[i] if i < len(poly_T) else 0)
               for i in range(max(len(check_T), len(poly_T))))

    return ok_T, ok_A, ok_B


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def fresh_stats() -> dict[str, Any]:
    return {
        "seen": 0,
        "considered": 0,
        "checked": 0,
        "no_deg2_support": 0,
        "identity_fail": 0,
        # Three-term decomposition
        "term1_neg": 0,    # (1+lam)*Phi_m(B) < 0
        "term2_pos": 0,    # lam*Z_B always positive (sanity)
        "term3_neg": 0,    # lam*Phi_{m-1}(P) < 0
        "both_neg": 0,     # term1 and term3 both negative
        # B properties
        "B_not_dleaf1": 0,
        "A_not_dleaf1": 0,
        # STRONG C2
        "strong_c2_fail": 0,
        "min_strong_c2_margin": None,
        "min_strong_c2_witness": None,
        # Cross-tree
        "cross_tree_fail": 0,
        "min_cross_margin": None,
        "min_cross_witness": None,
        # Three-term minima
        "min_term1": None,
        "min_term1_witness": None,
        "min_term3": None,
        "min_term3_witness": None,
        "min_pendant_ratio": None,  # min(pendant_bonus / |negative_sum|)
        "min_pendant_ratio_witness": None,
        # P, Q properties
        "P_not_lc": 0,
        "Q_not_lc": 0,
        "P_not_unimodal": 0,
        "phi_m1_P_neg": 0,  # Phi_{m-1}(P; lam) < 0
        # Overall failure
        "fail_phi_total": 0,
        "wall_s": 0.0,
    }


def maybe_update_min(stats: dict[str, Any], key: str, witness_key: str,
                     value: float, witness: dict[str, Any]) -> None:
    cur = stats[key]
    if cur is None or value < cur:
        stats[key] = value
        stats[witness_key] = witness


def inc(d: dict[str, int], key: str, amount: int = 1) -> None:
    d[key] = d.get(key, 0) + amount


def merge_stats(dst: dict[str, Any], src: dict[str, Any]) -> None:
    for k in ["seen", "considered", "checked", "no_deg2_support", "identity_fail",
              "term1_neg", "term2_pos", "term3_neg", "both_neg",
              "B_not_dleaf1", "A_not_dleaf1",
              "strong_c2_fail", "cross_tree_fail",
              "P_not_lc", "Q_not_lc", "P_not_unimodal", "phi_m1_P_neg",
              "fail_phi_total"]:
        dst[k] += int(src.get(k, 0))

    for key, witness_key in [
        ("min_strong_c2_margin", "min_strong_c2_witness"),
        ("min_cross_margin", "min_cross_witness"),
        ("min_term1", "min_term1_witness"),
        ("min_term3", "min_term3_witness"),
        ("min_pendant_ratio", "min_pendant_ratio_witness"),
    ]:
        if src.get(key) is not None:
            maybe_update_min(dst, key, witness_key, float(src[key]), src[witness_key])


def aggregate(per_n: dict[str, dict[str, Any]]) -> dict[str, Any]:
    out = fresh_stats()
    for key in sorted(per_n, key=lambda s: int(s)):
        merge_stats(out, per_n[key])
    out.pop("wall_s", None)
    return out


def write_payload(out_path: str, params: dict[str, Any],
                  per_n: dict[str, dict[str, Any]]) -> None:
    payload = {
        "params": params,
        "summary": aggregate(per_n),
        "per_n": {k: per_n[k] for k in sorted(per_n, key=lambda s: int(s))},
    }
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def fmt_opt(val: float | None) -> str:
    return "None" if val is None else f"{val:.6g}"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Diagnose degree-2 bridge decomposition for tie-fugacity."
    )
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--all-trees", action="store_true")
    ap.add_argument("--out", default="results/bridge_decomposition_diagnostic_n20.json")
    ap.add_argument("--no-resume", action="store_true")
    args = ap.parse_args()

    params = {
        "min_n": args.min_n,
        "max_n": args.max_n,
        "tol": args.tol,
        "all_trees": args.all_trees,
    }

    per_n: dict[str, dict[str, Any]] = {}
    if (not args.no_resume) and args.out and os.path.exists(args.out):
        with open(args.out, "r", encoding="utf-8") as f:
            old = json.load(f)
        per_n = old.get("per_n", {})
        done = ", ".join(sorted(per_n, key=lambda s: int(s)))
        print(f"Resuming from {args.out}; completed n: {done}", flush=True)

    scope = "all trees" if args.all_trees else "d_leaf<=1"
    print(f"Bridge decomposition diagnostic on {scope}, n={args.min_n}..{args.max_n}",
          flush=True)
    print(f"Output: {args.out}", flush=True)
    print("-" * 100, flush=True)

    t_all = time.time()
    for n in range(args.min_n, args.max_n + 1):
        n_key = str(n)
        if n_key in per_n:
            print(f"n={n:2d}: already complete, skipping", flush=True)
            continue

        t0 = time.time()
        stats = fresh_stats()

        proc = subprocess.Popen(
            [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"],
            stdout=subprocess.PIPE,
        )
        assert proc.stdout is not None

        for line in proc.stdout:
            nn, adj = parse_graph6(line)
            stats["seen"] += 1

            if (not args.all_trees) and (not is_dleaf_le_1(nn, adj)):
                continue
            stats["considered"] += 1

            poly_T = independence_poly(nn, adj)
            m = mode_index_leftmost(poly_T)
            if m == 0 or poly_T[m - 1] <= 0 or poly_T[m] <= 0:
                continue

            lam = poly_T[m - 1] / poly_T[m]
            deg = [len(nb) for nb in adj]
            g6 = line.decode("ascii").strip()

            # Find a leaf with degree-2 support (prefer smallest leaf index)
            best_leaf = -1
            best_support = -1
            for v in range(nn):
                if deg[v] == 1:
                    s = adj[v][0]
                    if deg[s] == 2:
                        best_leaf = v
                        best_support = s
                        break

            if best_leaf < 0:
                stats["no_deg2_support"] += 1
                continue

            stats["checked"] += 1
            leaf = best_leaf
            support = best_support

            # u = the other neighbour of s (not l)
            u_orig = [w for w in adj[support] if w != leaf][0]

            # Build A = T - l, B = T - {l, s}
            a_adj = remove_vertices(adj, {leaf})
            b_adj = remove_vertices(adj, {leaf, support})
            poly_A = independence_poly(len(a_adj), a_adj)
            poly_B = independence_poly(len(b_adj), b_adj)

            # Find u's index in B
            keep_B = [v for v in range(nn) if v not in {leaf, support}]
            idx_B = {v: i for i, v in enumerate(keep_B)}
            u_in_B = idx_B[u_orig]

            # Check if A and B remain d_leaf ≤ 1
            if not is_dleaf_le_1(len(b_adj), b_adj):
                stats["B_not_dleaf1"] += 1
            if not is_dleaf_le_1(len(a_adj), a_adj):
                stats["A_not_dleaf1"] += 1

            # Compute P, Q (hub-off, hub-on at u in B)
            P, Q = compute_hub_polys(b_adj, u_in_B)

            # Verify identities
            ok_T, ok_A, ok_B = verify_PQ_identity(poly_T, poly_A, poly_B, P, Q)
            if not (ok_T and ok_A and ok_B):
                stats["identity_fail"] += 1

            # P, Q properties
            if len(P) > 2 and not is_log_concave(P):
                stats["P_not_lc"] += 1
            if len(Q) > 2 and not is_log_concave(Q):
                stats["Q_not_lc"] += 1

            # Check P unimodality
            if not is_unimodal(P):
                stats["P_not_unimodal"] += 1

            mode_P = mode_index_leftmost(P) if P else 0
            mode_Q = mode_index_leftmost(Q) if Q and any(c > 0 for c in Q) else 0
            mode_B = mode_index_leftmost(poly_B) if poly_B else 0
            mode_A = mode_index_leftmost(poly_A) if poly_A else 0

            # --- Three-term decomposition ---
            # Phi_m(T; lam) = (1+lam)*Phi_m(B; lam) + lam*Z_B(lam) + lam*Phi_{m-1}(P; lam)
            phi_m_B = phi_q(poly_B, m, lam)
            z_B = eval_poly(poly_B, lam)
            phi_m1_P = phi_q(P, m - 1, lam)

            term1 = (1 + lam) * phi_m_B       # (1+lam)*Phi_m(B; lam)
            term2 = lam * z_B                  # lam*Z_B(lam) -- pendant bonus
            term3 = lam * phi_m1_P             # lam*Phi_{m-1}(P; lam)
            total = term1 + term2 + term3

            # Verify total matches Phi_m(T; lam)
            phi_m_T = phi_q(poly_T, m, lam)

            if term1 < -args.tol:
                stats["term1_neg"] += 1
            if term2 > args.tol:
                stats["term2_pos"] += 1
            if term3 < -args.tol:
                stats["term3_neg"] += 1
            if term1 < -args.tol and term3 < -args.tol:
                stats["both_neg"] += 1
            if phi_m1_P < -args.tol:
                stats["phi_m1_P_neg"] += 1
            if total < -args.tol:
                stats["fail_phi_total"] += 1

            # Pendant ratio: term2 / |negative part|
            neg_part = min(term1, 0) + min(term3, 0)
            if neg_part < -args.tol:
                pendant_ratio = term2 / abs(neg_part)
            else:
                pendant_ratio = float("inf")

            base_witness = {
                "n": nn, "g6": g6, "mode_T": m, "lam": lam,
                "leaf": leaf, "support": support, "u": u_orig,
                "deg_sig": dict(sorted(Counter(deg).items())),
                "mode_A": mode_A, "mode_B": mode_B,
                "mode_P": mode_P, "mode_Q": mode_Q,
                "P_lc": is_log_concave(P) if len(P) > 2 else True,
                "Q_lc": is_log_concave(Q) if len(Q) > 2 else True,
                "term1": term1, "term2": term2, "term3": term3,
                "total": total, "phi_m_T": phi_m_T,
            }

            maybe_update_min(stats, "min_term1", "min_term1_witness", term1, base_witness)
            maybe_update_min(stats, "min_term3", "min_term3_witness", term3, base_witness)
            if pendant_ratio != float("inf"):
                maybe_update_min(stats, "min_pendant_ratio", "min_pendant_ratio_witness",
                                 pendant_ratio, base_witness)

            # --- STRONG C2: lambda_m^T >= lambda_{m-1}^B ---
            if m - 2 >= 0 and m - 1 < len(poly_B) and poly_B[m - 1] > 0:
                lam_B = poly_B[m - 2] / poly_B[m - 1]
                c2_margin = lam - lam_B
                c2_witness = dict(base_witness)
                c2_witness["lam_B"] = lam_B
                c2_witness["c2_margin"] = c2_margin
                if c2_margin < -args.tol:
                    stats["strong_c2_fail"] += 1
                maybe_update_min(stats, "min_strong_c2_margin",
                                 "min_strong_c2_witness", c2_margin, c2_witness)

            # --- Cross-tree: lambda_m^A vs lambda_{m-1}^B ---
            if (m - 1 < len(poly_A) and m < len(poly_A) and poly_A[m] > 0
                    and m - 2 >= 0 and m - 1 < len(poly_B) and poly_B[m - 1] > 0):
                lam_A = poly_A[m - 1] / poly_A[m]
                lam_B2 = poly_B[m - 2] / poly_B[m - 1]
                cross_margin = lam_A - lam_B2
                cross_witness = dict(base_witness)
                cross_witness["lam_A"] = lam_A
                cross_witness["lam_B"] = lam_B2
                cross_witness["cross_margin"] = cross_margin
                if cross_margin < -args.tol:
                    stats["cross_tree_fail"] += 1
                maybe_update_min(stats, "min_cross_margin",
                                 "min_cross_witness", cross_margin, cross_witness)

        proc.wait()
        stats["wall_s"] = time.time() - t0
        per_n[n_key] = stats
        write_payload(args.out, params, per_n)

        print(
            f"n={n:2d}: seen={stats['seen']:9d} ck={stats['checked']:8d} "
            f"t1-={stats['term1_neg']:5d} t3-={stats['term3_neg']:5d} "
            f"both-={stats['both_neg']:5d} P!lc={stats['P_not_lc']:4d} "
            f"c2F={stats['strong_c2_fail']:3d} xF={stats['cross_tree_fail']:3d} "
            f"totF={stats['fail_phi_total']:3d} "
            f"min_c2={fmt_opt(stats['min_strong_c2_margin'])} "
            f"min_pr={fmt_opt(stats['min_pendant_ratio'])} "
            f"({stats['wall_s']:.1f}s)",
            flush=True,
        )

    summary = aggregate(per_n)
    wall = time.time() - t_all
    print("=" * 100, flush=True)
    print(f"TOTAL checked={summary['checked']:,}  wall={wall:.1f}s", flush=True)
    print(f"  identity_fail={summary['identity_fail']}", flush=True)
    print(f"  term1_neg={summary['term1_neg']}  term3_neg={summary['term3_neg']}  "
          f"both_neg={summary['both_neg']}", flush=True)
    print(f"  B_not_dleaf1={summary['B_not_dleaf1']}  A_not_dleaf1={summary['A_not_dleaf1']}",
          flush=True)
    print(f"  P_not_lc={summary['P_not_lc']}  Q_not_lc={summary['Q_not_lc']}  "
          f"P_not_unimodal={summary['P_not_unimodal']}", flush=True)
    print(f"  phi_m1_P_neg={summary['phi_m1_P_neg']}", flush=True)
    print(f"  strong_c2_fail={summary['strong_c2_fail']}  "
          f"cross_tree_fail={summary['cross_tree_fail']}", flush=True)
    print(f"  fail_phi_total={summary['fail_phi_total']}", flush=True)
    print(f"  min_strong_c2_margin={summary['min_strong_c2_margin']}", flush=True)
    print(f"  min_cross_margin={summary['min_cross_margin']}", flush=True)
    print(f"  min_term1={summary['min_term1']}", flush=True)
    print(f"  min_term3={summary['min_term3']}", flush=True)
    print(f"  min_pendant_ratio={summary['min_pendant_ratio']}", flush=True)
    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
