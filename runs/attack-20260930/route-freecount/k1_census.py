"""K1: exhaustive exact kill-test of Lemma PV on all trees of order n.

Usage: python3 k1_census.py N [res mod] > out.json
Trees from `gentreeg -p N [res/mod]`.

For every tree, vertex v and level k:
  PV    (slack):  k i_{k-1}(T) j_k <= (k+1) i_k(T) j_{k-1}
  PV0   (no slack, Mason-type diagnostic only): i_{k-1} j_k <= i_k j_{k-1}
  RDs   (ratio dominance with slack at v):   k j_k e_{k-1} <= (k+1) j_{k-1} e_k
  WLC   (weak LC of T-N[v] at k-1):          k j_k j_{k-2} <= (k+1) j_{k-1}^2
where j = I(T-N[v]), e = I(T-v). PV is checked on the window [ceil(n/4), q]
and, for the "where does PV first fail" diagnostic, at every k in [1, alpha].
All comparisons are exact integer cross-multiplications.
"""

import json
import subprocess
import sys
import time
from fractions import Fraction

from pv_lib import co, lc_ok, parse_parent_array, pv_lhs_rhs, tree_data, window


def main():
    N = int(sys.argv[1])
    cmd = ["gentreeg", "-q", "-p", str(N)]
    if len(sys.argv) >= 4:
        cmd.append(f"{sys.argv[2]}/{sys.argv[3]}")
    t0 = time.time()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    stats = {
        "n": N,
        "shard": sys.argv[2:4],
        "trees": 0,
        "window_checks_vk": 0,
        "window_levels": 0,
        "lc_fail_window": 0,
        "pv_fail_window": 0,
        "pv0_fail_window": 0,
        "rds_fail_window": 0,
        "wlc_fail_window": 0,
        "mediant_implication_violations": 0,
        "pv_fail_anywhere": 0,
        "pv_fail_anywhere_trees": 0,
        "pv_fail_offsets_k_minus_alpha": {},
        "pv_fail_offsets_k_minus_q": {},
        "lc_fail_anywhere": 0,
        "max_pv_ratio_window": None,
        "max_pv_ratio_window_float": None,
        "max_pv_ratio_witness": None,
        "max_pv0_ratio_window_float": None,
        "min_window_k_minus_q_of_pv_fail": None,
        "empty_window_trees": 0,
    }
    best = Fraction(-1)
    best0 = Fraction(-1)
    for line in proc.stdout:
        toks = line.split()
        if not toks:
            continue
        adj = parse_parent_array(toks)
        n = len(adj)
        I, J, E, alpha = tree_data(adj)
        stats["trees"] += 1
        lo, q = window(n, alpha)
        lo = max(lo, 1)
        if lo > q:
            stats["empty_window_trees"] += 1
        tree_fail_any = False
        for k in range(1, alpha + 1):
            in_win = lo <= k <= q
            if not lc_ok(I, k):
                stats["lc_fail_anywhere"] += 1
                if in_win:
                    stats["lc_fail_window"] += 1
            if in_win:
                stats["window_levels"] += 1
            for v in range(n):
                Jv, Ev = J[v], E[v]
                lhs, rhs = pv_lhs_rhs(I, Jv, k)
                pv_ok = lhs <= rhs
                if not pv_ok:
                    stats["pv_fail_anywhere"] += 1
                    tree_fail_any = True
                    d = str(k - alpha)
                    stats["pv_fail_offsets_k_minus_alpha"][d] = (
                        stats["pv_fail_offsets_k_minus_alpha"].get(d, 0) + 1)
                    dq = k - q
                    stats["pv_fail_offsets_k_minus_q"][str(dq)] = (
                        stats["pv_fail_offsets_k_minus_q"].get(str(dq), 0) + 1)
                    m = stats["min_window_k_minus_q_of_pv_fail"]
                    if m is None or dq < m:
                        stats["min_window_k_minus_q_of_pv_fail"] = dq
                if not in_win:
                    continue
                stats["window_checks_vk"] += 1
                if not pv_ok:
                    stats["pv_fail_window"] += 1
                if rhs > 0:
                    r = Fraction(lhs, rhs)
                    if r > best:
                        best = r
                        stats["max_pv_ratio_witness"] = {
                            "parent_array": toks, "v": v, "deg_v": len(adj[v]),
                            "k": k, "alpha": alpha, "q": q, "lo": lo,
                            "I": I, "J_v": Jv}
                l0, r0 = pv_lhs_rhs(I, Jv, k, slack=False)
                if l0 > r0:
                    stats["pv0_fail_window"] += 1
                if r0 > 0:
                    r = Fraction(l0, r0)
                    if r > best0:
                        best0 = r
                # ratio dominance with slack at v
                rds = k * co(Jv, k) * co(Ev, k - 1) <= (k + 1) * co(Jv, k - 1) * co(Ev, k)
                wlc = (k < 2) or (k * co(Jv, k) * co(Jv, k - 2) <= (k + 1) * co(Jv, k - 1) ** 2)
                if not rds:
                    stats["rds_fail_window"] += 1
                if not wlc:
                    stats["wlc_fail_window"] += 1
                if rds and wlc and not pv_ok:
                    stats["mediant_implication_violations"] += 1
        if tree_fail_any:
            stats["pv_fail_anywhere_trees"] += 1
    proc.wait()
    if best >= 0:
        stats["max_pv_ratio_window"] = f"{best.numerator}/{best.denominator}"
        stats["max_pv_ratio_window_float"] = float(best)  # diagnostic only
    if best0 >= 0:
        stats["max_pv0_ratio_window_float"] = float(best0)  # diagnostic only
    stats["seconds"] = round(time.time() - t0, 2)
    json.dump(stats, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
