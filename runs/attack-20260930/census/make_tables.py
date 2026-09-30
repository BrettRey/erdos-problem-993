#!/usr/bin/env python3
"""Print markdown tables (for the returned report) and write data/per_n_summary.csv from data/*.json,
and verify the star closed form against every exhaustive minimum."""
import glob
import json
import os
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")


def star_min(n):
    # K_{1,n-1}: delta_k = n/((k+1)(n-k)) for k >= 3; minimum at the mode
    return Fraction(4, n + 2) if n % 2 == 0 else Fraction(4 * n, (n + 1) ** 2)


def census_rows():
    rows = []
    for p in sorted(glob.glob(os.path.join(D, "census_n*.json")), key=lambda s: int(s.split("_n")[1].split(".")[0])):
        rows.append(json.load(open(p)))
    return rows


def main():
    out = []
    cr = census_rows()
    out.append("### Table 1. Exhaustive census (exact)\n")
    out.append("| n | trees (= A000055) | central LC violations W5 / W4 | min δ over W5 (exact) | n·min δ | = star closed form | argmin | runner-up n·m(T) | window-top min n·δ_q | Newton ratio min (exact, float) | trees failing Newton in W5 |")
    out.append("|---|---|---|---|---|---|---|---|---|---|---|")
    allstar = True
    for c in cr:
        n = c["n"]
        t = c["TOP5"][0]
        m = Fraction(t["margin"])
        ok = (m == star_min(n))
        allstar &= ok
        is_star = t["par"] == ",".join(["0"] + ["1"] * (n - 1))
        ru = c["TOP5"][1]
        tq = c["window_top_min_delta_q"]
        nw = c["newton_ratio_min_exact"]
        out.append(f"| {n} | {c['trees']:,} ({'ok' if c['count_ok'] else 'MISMATCH'}) | {c['viol5']} / {c['viol4']} | {t['margin']} | "
                   f"{t['n_times_margin_float']:.4f} | {'yes' if ok else 'NO'} | {'star K_{1,n-1}' if is_star else t['par']}, α={t['alpha']}, k={t['k']} | "
                   f"{ru['n_times_margin_float']:.4f} (α={ru['alpha']}, k={ru['k']}) | {tq['n_times_value_float']:.4f} ({'star' if tq['is_star'] else 'α=' + str(tq['alpha'])}) | "
                   f"{nw['value_float']:.4f} (α={nw['alpha']}, k={nw['k']}, q={nw['q']}) | {c.get('newton_fail_trees', 'n/a'):,} |")
    out.append(f"\nStar closed form matches the exhaustive minimum at every n above: **{allstar}**.\n")

    out.append("### Table 2. Tilted diagnostics from the exhaustive runs (FLOAT)\n")
    out.append("| n | min V·δ over all trees, k∈W5 (V, α, k, q) | min V·δ with V∈[2,4) | min V·δ with V≥4 | star V·δ at its argmin |")
    out.append("|---|---|---|---|---|")
    ts = json.load(open(os.path.join(D, "tilted_summary.json"))) if os.path.exists(os.path.join(D, "tilted_summary.json")) else {"extremal": {}}
    for c in cr:
        n = c["n"]
        v = c["VDMIN_FLOAT_diagnostic"]
        b = c.get("VD_by_V_bin_FLOAT_diagnostic", {})
        b2 = b.get("2.0")
        b4 = [r for k, r in b.items() if float(k) >= 4]
        b4m = min(b4, key=lambda r: r["value"]) if b4 else None
        st = ts["extremal"].get(str(n), {}).get("min_delta_star", {})
        out.append(f"| {n} | {v['value']:.4f} (V={v['V']:.2f}, α={v['alpha']}, k={v['k']}, q={v['q']}) | "
                   f"{format(b2['value'], '.4f') if b2 else '-'} | "
                   f"{(format(b4m['value'], '.4f') + ' (V=' + format(b4m['V'], '.2f') + ')') if b4m else '-'} | "
                   f"{(format(st['Vd_at_min_delta'], '.4f') + ' (V=' + format(st['V_at_min_delta'], '.2f') + ')') if st else '-'} |")

    qs = sorted(glob.glob(os.path.join(D, "qstat_n*.json")), key=lambda s: int(s.split("_n")[1].split(".")[0]))
    if qs:
        out.append("\n### Table 3. Fang-type statistic Q_k = V^{3/2}(2p_k − p_{k−1} − p_{k+1}) at the tilt with mean k (FLOAT)\n")
        out.append("| n | trees | trees with some Q_k ≤ 0 in W5 | min Q (V, α, k, q) | min Q with V≥2 | min Q with V≥4 |")
        out.append("|---|---|---|---|---|---|")
        for p in qs:
            q = json.load(open(p))
            bins = q["Q_by_V_bin_FLOAT"]
            allb = min(bins.values(), key=lambda r: r["Q"])
            ge2 = [r for k, r in bins.items() if float(k) >= 2]
            ge4 = [r for k, r in bins.items() if float(k) >= 4]
            g2 = min(ge2, key=lambda r: r["Q"]) if ge2 else None
            g4 = min(ge4, key=lambda r: r["Q"]) if ge4 else None
            out.append(f"| {q['n']} | {q['trees']:,} | {q['q_nonpositive_trees']} | {allb['Q']:.4f} (V={allb['V']:.2f}, α={allb['alpha']}, k={allb['k']}, q={allb['q']}) | "
                       f"{(format(g2['Q'], '.4f')) if g2 else '-'} | {(format(g4['Q'], '.4f') + ' (V=' + format(g4['V'], '.2f') + ')') if g4 else '-'} |")

    fp = os.path.join(D, "families.json")
    if os.path.exists(fp):
        fam = json.load(open(fp))["rows"]
        by = {}
        for r in fam:
            by.setdefault(r["family"], []).append(r)
        out.append("\n### Table 4. Families (exact δ; FLOAT tilted columns)\n")
        out.append("| family | members | max n | n·min δ at max n | worst n·min δ (n) | n·δ_q at max n | min Newton ratio | min V·δ (V there) | V·δ at argmin, max n | min Q | Q at argmin, max n | central viol. | non-unimodal |")
        out.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for f, rs in by.items():
            rs.sort(key=lambda r: r["n"])
            last = rs[-1]
            w = min(rs, key=lambda r: r["n_min_delta_W5"])
            vd = min(rs, key=lambda r: r["Vdelta_min"])
            out.append(f"| {f} | {len(rs)} | {last['n']} | {last['n_min_delta_W5']:.4f} | {w['n_min_delta_W5']:.4f} ({w['n']}) | {last['n_delta_q']:.3f} | "
                       f"{min(r['newton_ratio_min'] for r in rs):.4f} | {vd['Vdelta_min']:.4f} ({vd['Vdelta_min_V']:.2f}) | {last['Vdelta_at_argmin']:.4f} | "
                       f"{min(r.get('Q_min', 9) for r in rs):.4f} | {last.get('Q_at_argmin', float('nan')):.4f} | "
                       f"{sum(r['central_violation'] for r in rs)} | {sum(not r['unimodal'] for r in rs)} |")
        # large-V regime
        big = [(b, v) for r in fam for b, v in r.get("Vdelta_min_by_Vbin", {}).items()]
        out.append("\nFamily minima of V·δ by V bin (all members, all window k): " +
                   ", ".join(f"V≥{b}: {min(v for bb, v in big if bb == b):.4f}" for b in sorted({b for b, _ in big}, key=float)))
        # star comparison at equal n
        stars = {r["n"]: r for r in by.get("star", [])}
        worse = []
        for r in fam:
            n = r["n"]
            if r["family"] != "star" and r["min_delta_W5_float"] < float(star_min(n)) - 1e-15:
                worse.append((r["family"], n, r["n_min_delta_W5"], float(n * star_min(n))))
        out.append(f"\nFamily members (any n) whose central minimum lies strictly below the star's value at the same n: {len(worse)} {worse[:10]}")
    txt = "\n".join(out)
    print(txt)
    # machine-readable per-n summary
    import csv
    with open(os.path.join(D, "per_n_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["n", "trees", "count_ok", "viol5", "viol4", "min_delta_W5_exact", "n_times_min_delta",
                    "equals_star_closed_form", "argmin_par", "argmin_alpha", "argmin_k", "runner_up_n_times",
                    "window_top_min_n_delta_q", "newton_min_exact", "newton_min_float", "newton_fail_trees",
                    "VD_min_float", "VD_min_V", "VD_min_alpha", "VD_min_k", "VD_min_q"])
        for c in cr:
            n = c["n"]; t = c["TOP5"][0]; v = c["VDMIN_FLOAT_diagnostic"]; nw = c["newton_ratio_min_exact"]
            w.writerow([n, c["trees"], c["count_ok"], c["viol5"], c["viol4"], t["margin"], t["n_times_margin_float"],
                        Fraction(t["margin"]) == star_min(n), t["par"], t["alpha"], t["k"],
                        c["TOP5"][1]["n_times_margin_float"], c["window_top_min_delta_q"]["n_times_value_float"],
                        nw["value"], nw["value_float"], c.get("newton_fail_trees"),
                        v["value"], v["V"], v["alpha"], v["k"], v["q"]])


if __name__ == "__main__":
    main()
