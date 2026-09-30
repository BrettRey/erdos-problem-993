#!/usr/bin/env python3
"""Part 2: tilted scale V_k * delta_k on the central window (FLOAT DIAGNOSTIC).

Independent Python re-implementation (bisection on log lambda, treepoly.tilt)
of the C program's safeguarded-Newton diagnostic.

(a) For every census n: recompute V_k delta_k for the extremal trees recorded
    by the C census (exact-min star, window-top minimiser, Newton-ratio
    minimiser, the C VDMIN minimiser and the per-V-bin minimisers) and check
    the C VDMIN value.
(b) For the Bernoulli samples (hash of the parent array, rate ~2000/A000055(n))
    at n in {16, 20, 24}: per-tree min over W5 of V_k delta_k, its location,
    and V delta at the k minimising delta; summary statistics.
"""
import json
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from treepoly import adj_from_parent, ipoly_from_adj, margins, q_of, tilt, tilted_Q  # noqa: E402


def profile(poly, n):
    rows = []
    ms = margins(poly, n, "n5")
    for k, d in ms.items():
        lam, V = tilt(poly, k, iters=120)
        rows.append(dict(k=k, delta=float(d), lam=lam, V=V, Vd=V * float(d), Q=tilted_Q(poly, k, lam, V)))
    return rows


def poly_of_par(par_str):
    par = [int(x) for x in par_str.split(",")]
    return ipoly_from_adj(adj_from_parent(par))


def extremal_checks(ns):
    out = {}
    for n in ns:
        path = os.path.join(HERE, "data", f"census_n{n}.json")
        if not os.path.exists(path):
            continue
        c = json.load(open(path))
        rec = {}
        cands = {
            "min_delta_star": c["TOP5"][0]["par"],
            "window_top_min": c["window_top_min_delta_q"]["par"],
            "newton_min": c["newton_ratio_min_exact"]["par"],
            "VDMIN_C": c["VDMIN_FLOAT_diagnostic"]["par"],
        }
        for b, r in c.get("VD_by_V_bin_FLOAT_diagnostic", {}).items():
            cands[f"VDbin_{b}"] = r["par"]
        for name, par in cands.items():
            poly = poly_of_par(par)
            prof = profile(poly, n)
            best = min(prof, key=lambda r: r["Vd"])
            md = min(prof, key=lambda r: r["delta"])
            qm = min(prof, key=lambda r: r["Q"])
            rec[name] = dict(par=par, alpha=len(poly) - 1, q=q_of(len(poly) - 1),
                             Q_min=qm["Q"], Q_min_k=qm["k"], Q_min_V=qm["V"], Q_at_min_delta=md["Q"],
                             Vd_min=best["Vd"], Vd_min_k=best["k"], Vd_min_V=best["V"],
                             Vd_min_lambda=best["lam"],
                             Vd_at_min_delta=md["Vd"], k_min_delta=md["k"], V_at_min_delta=md["V"])
        cv = c["VDMIN_FLOAT_diagnostic"]["value"]
        pv = rec["VDMIN_C"]["Vd_min"]
        rec["C_vs_python_VDMIN"] = dict(C=cv, python=pv, absdiff=abs(cv - pv))
        out[n] = rec
        print(f"n={n}: C VDMIN={cv:.9f} python={pv:.9f} | star Vd@argmin={rec['min_delta_star']['Vd_at_min_delta']:.4f} "
              f"(V={rec['min_delta_star']['V_at_min_delta']:.2f}) | newton-min tree Vd_min={rec['newton_min']['Vd_min']:.4f}",
              flush=True)
    return out


def sample_stats(n):
    path = os.path.join(HERE, "data", f"sample_n{n}.jsonl")
    trees = [json.loads(l) for l in open(path)]
    per = []
    for t in trees:
        poly = [int(x) for x in t["poly"].split(",")]
        prof = profile(poly, n)
        best = min(prof, key=lambda r: r["Vd"])
        md = min(prof, key=lambda r: r["delta"])
        alpha = len(poly) - 1
        per.append(dict(par=t["par"], alpha=alpha, q=q_of(alpha), Vd_min=best["Vd"], k=best["k"],
                        q_minus_k=q_of(alpha) - best["k"], V=best["V"], lam=best["lam"],
                        Vd_at_min_delta=md["Vd"], n_min_delta=n * md["delta"], V_at_min_delta=md["V"],
                        maxV=max(r["V"] for r in prof),
                        Q_min=min(r["Q"] for r in prof),
                        Vd_minV4=min((r["Vd"] for r in prof if r["V"] >= 4), default=None)))
    vds = sorted(r["Vd_min"] for r in per)
    worst = min(per, key=lambda r: r["Vd_min"])
    frac_at_top = sum(r["q_minus_k"] == 0 for r in per) / len(per)
    v4 = [r["Vd_minV4"] for r in per if r["Vd_minV4"] is not None]
    summ = dict(n=n, sampled=len(per), Vd_min=vds[0], Vd_p01=vds[len(vds) // 100], Vd_median=statistics.median(vds),
                Vd_max=vds[-1], worst=worst, frac_min_at_window_top=frac_at_top,
                Vd_at_min_delta_min=min(r["Vd_at_min_delta"] for r in per),
                Vd_at_min_delta_median=statistics.median(r["Vd_at_min_delta"] for r in per),
                n_min_delta_min=min(r["n_min_delta"] for r in per),
                n_min_delta_median=statistics.median(r["n_min_delta"] for r in per),
                Vd_min_restricted_V_ge_4=(min(v4) if v4 else None), trees_with_V_ge_4=len(v4),
                max_V_seen=max(r["maxV"] for r in per),
                below_quarter=sum(v < 0.25 for v in vds),
                Q_min=min(r["Q_min"] for r in per), Q_median=statistics.median(r["Q_min"] for r in per),
                Q_nonpositive=sum(r["Q_min"] <= 0 for r in per))
    json.dump(per, open(os.path.join(HERE, "data", f"tilted_sample_n{n}.json"), "w"), indent=0)
    print(f"sample n={n}: {len(per)} trees  Vd_min={summ['Vd_min']:.4f} p01={summ['Vd_p01']:.4f} "
          f"median={summ['Vd_median']:.4f} | at-top frac={frac_at_top:.2f} | Vd@argmin-delta min={summ['Vd_at_min_delta_min']:.4f} "
          f"| n*minDelta min={summ['n_min_delta_min']:.3f} median={summ['n_min_delta_median']:.3f} "
          f"| V>=4: {len(v4)} trees, min Vd={summ['Vd_min_restricted_V_ge_4']} | maxV={summ['max_V_seen']:.2f} "
          f"| Qmin={summ['Q_min']:.4f} Qmed={summ['Q_median']:.4f} Q<=0: {summ['Q_nonpositive']}",
          flush=True)
    return summ


def main():
    ns = [int(x) for x in sys.argv[1:]] or list(range(10, 28))
    ext = extremal_checks(ns)
    samp = {n: sample_stats(n) for n in (16, 20, 24)
            if os.path.exists(os.path.join(HERE, "data", f"sample_n{n}.jsonl"))}
    json.dump(dict(extremal=ext, samples=samp), open(os.path.join(HERE, "data", "tilted_summary.json"), "w"),
              indent=1)


if __name__ == "__main__":
    main()
