"""Merge r1_census shard outputs, verify tree counts against OEIS A000055,
recompute every reported tightest case EXACTLY (Fraction) with pv_lib and
with an independent non-DP recursion (vertex deletion recurrence
I(G) = I(G - v) + x I(G - N[v]) on vertex sets, memoised), and recheck any
R1/PV violation the same way.
Usage: python3 merge.py N0 N1 SHARDS > results.json
"""
import glob
import json
import sys
from fractions import Fraction
from functools import lru_cache

sys.path.insert(0, "/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/runs/attack-20260930/route-freecount")
from pv_lib import co, parse_parent_array, tree_data, window  # noqa: E402

HERE = "/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/runs/attack-20260930/r1-exhaustive"
OEIS = {19: 317955, 20: 823065, 21: 2144505, 22: 5623756, 23: 14828074, 24: 39299897,
        25: 104636890, 26: 279793450, 27: 751065460, 28: 2023443032}
NREPORT = 10


def indep_by_recursion(adj, verts):
    """Independence polynomial of the induced subgraph on `verts`, by the
    deletion recurrence (no tree DP)."""
    @lru_cache(maxsize=None)
    def rec(S):
        if not S:
            return (1,)
        # choose a vertex of max degree inside S
        v = max(S, key=lambda x: sum(1 for w in adj[x] if w in S))
        a = rec(S - {v})
        b = rec(S - {v} - frozenset(w for w in adj[v] if w in S))
        out = [0] * max(len(a), len(b) + 1)
        for i, c in enumerate(a):
            out[i] += c
        for i, c in enumerate(b):
            out[i + 1] += c
        return tuple(out)
    return list(rec(frozenset(verts)))


def exact_case(par_tokens, k, u_c=None, recursion_check=True):
    """Exact L_u(k), normalisations, for C vertex u_c (1-indexed)."""
    adj = parse_parent_array(par_tokens)
    n = len(adj)
    I, J, E, alpha = tree_data(adj)
    lo, q = window(n, alpha)
    lo = max(lo, 1)
    D = [k * co(I, k - 1) * co(J[v], k) - (k + 1) * co(I, k) * co(J[v], k - 1) for v in range(n)]
    L = []
    for u in range(n):
        s = Fraction(D[u], len(adj[u]) + 1)
        for v in adj[u]:
            s += Fraction(D[v], len(adj[v]) + 1)
        L.append(s)
    absS = sum(abs(d) for d in D)
    res = {"n": n, "alpha": alpha, "window": [lo, q], "k": k, "I": I}
    if u_c is not None:
        u = u_c - 1
        res.update({
            "u_c_1indexed": u_c, "deg_u": len(adj[u]),
            "L_u": str(L[u]),
            "norm1_exact": str(L[u] / (k * co(I, k - 1) * co(I, k))),
            "norm2_exact": str(L[u] * n / absS) if absS else None,
            "norm1_float_diag": float(L[u] / (k * co(I, k - 1) * co(I, k))),
            "norm2_float_diag": float(L[u] * n / absS) if absS else None,
            "max_L_over_u": str(max(L)),
            "D_u": D[u],
        })
    if recursion_check:
        V = set(range(n))
        I2 = indep_by_recursion(adj, V)
        assert I2 == I, ("recursion I mismatch", I2, I)
        for v in range(n):
            Jv = indep_by_recursion(adj, V - {v} - set(adj[v]))
            b = list(J[v])
            while len(b) > 1 and b[-1] == 0:
                b.pop()
            assert Jv == b, ("recursion J mismatch", v)
        res["independent_recursion_check"] = "I and every J_v agree"
    return res


def parse_kv(line):
    return dict(x.split("=", 1) for x in line.split()[1:])


def shard_files(n, m, d):
    fs = sorted(glob.glob(f"{HERE}/{d}/n{n}_s*of{m}.txt"))
    return fs if len(fs) == m else None


def main():
    n0, n1, m = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    out = {"per_n": {}, "violations": []}
    for n in range(n0, n1 + 1):
        f_raw, f_ext = shard_files(n, m, "raw"), shard_files(n, m, "raw_ext")
        files = f_raw or f_ext
        if files is None:
            out["per_n"][n] = {"status": "incomplete: no complete shard set"}
            continue
        agg = {}
        ext = {}
        ext_hist = {}
        ext_min = None
        tops = {"TOP_NORM1": [], "TOP_NORM2": []}
        pvmax = None
        viol_lines = []
        # identity check between the two binaries where both ran
        same = None
        if f_raw and f_ext:
            same = True
            for a, b in zip(f_raw, f_ext):
                la = [x for x in open(a) if not x.startswith("EXT")]
                lb = [x for x in open(b) if not x.startswith("EXT")]
                if la != lb:
                    same = False
        if f_ext:
            for f in f_ext:
                for line in open(f):
                    if line.startswith("EXTSTATS"):
                        for key, val in parse_kv(line).items():
                            if key != "n":
                                ext[key] = ext.get(key, 0) + int(val)
                    elif line.startswith("EXT_R1_MIN_OFFSET"):
                        kv = parse_kv(line)
                        if ext_min is None or int(kv["k_minus_q"]) < int(ext_min["k_minus_q"]):
                            ext_min = kv
                    elif line.startswith("EXT_"):
                        tag = line.split()[0]
                        kv = parse_kv(line)
                        key = [k for k in kv if k not in ("n", "trees")][0]
                        ext_hist.setdefault(tag, {})
                        ext_hist[tag][kv[key]] = ext_hist[tag].get(kv[key], 0) + int(kv["trees"])
        for f in files:
            for line in open(f):
                if line.startswith("STATS"):
                    kv = parse_kv(line)
                    for key, val in kv.items():
                        if key in ("n", "LAM"):
                            agg[key] = val
                        else:
                            agg[key] = agg.get(key, 0) + int(val)
                elif line.startswith("TOP_NORM"):
                    tag = line.split()[0]
                    tops[tag].append(parse_kv(line))
                elif line.startswith("PVMAX"):
                    kv = parse_kv(line)
                    r = Fraction(int(kv["lhs"]), int(kv["rhs"]))
                    if pvmax is None or r > pvmax[0]:
                        pvmax = (r, kv)
                elif line.startswith("R1VIOL ") or line.startswith("PVVIOL "):
                    viol_lines.append(line.strip())
        rec = {"source_dir": "raw (r1_census)" if f_raw else "raw_ext (r1_ext)",
               "r1_census_vs_r1_ext_inwindow_identical": same,
               "stats": agg, "oeis_A000055": OEIS.get(n), "count_ok": int(agg["trees"]) == OEIS.get(n)}
        if f_ext:
            rec["outside_window_ext"] = {"stats": ext, "hist": ext_hist, "r1_min_offset_above_q": ext_min}
        if pvmax:
            rec["pv_max_ratio_lhs_over_rhs"] = {"exact": f"{pvmax[0].numerator}/{pvmax[0].denominator}",
                                                "float_diag": float(pvmax[0]), **pvmax[1]}
        for tag, lst in tops.items():
            lst.sort(key=lambda kv: -float(kv["float"]))
            rep = []
            for i, kv in enumerate(lst[:NREPORT]):
                ex = exact_case(kv["par"].split(","), int(kv["k"]), int(kv["u"]), recursion_check=(i < 3))
                # consistency: C's LL / LAM equals exact L_u
                assert Fraction(int(kv["LL"]), int(agg["LAM"])) == Fraction(ex["L_u"]), (tag, kv)
                ex["par"] = kv["par"]
                ex["float_rank_value_C"] = float(kv["float"])
                rep.append(ex)
            rec[tag] = rep
        rec["violation_lines_first"] = viol_lines[:20]
        rec["violation_lines_total_printed"] = len(viol_lines)
        out["per_n"][n] = rec
        for line in viol_lines[:10]:
            kv = parse_kv(line)
            u = int(kv.get("u", kv.get("v")))
            ex = exact_case(kv["par"].split(","), int(kv["k"]), u, recursion_check=True)
            ex["line"] = line
            out["violations"].append(ex)
        print(f"merged n={n}", file=sys.stderr, flush=True)
    json.dump(out, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
