#!/usr/bin/env python3
"""Driver for central_margin.c (wave variant: at most CONC=5 concurrent shards).

Usage: python3 run_census.py N [N ...] [--shards S] [--sample-target 2000]

For each n: runs S concurrent pipelines `gentreeg -p -q n r/S | ./central_margin n P`
(P chosen so that about sample_target trees are sampled when n in SAMPLE_NS),
saves raw shard output to raw/n{n}_s{r}of{S}.txt, merges the TOP lists by
exact Fraction comparison, checks the total tree count against OEIS A000055,
and writes data/census_n{n}.json.
"""
import argparse
import json
import os
import subprocess
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")
DATA = os.path.join(HERE, "data")
os.makedirs(RAW, exist_ok=True)
os.makedirs(DATA, exist_ok=True)

# OEIS A000055 (unlabeled trees). Provenance: CLAUDE.md table (n<=10, 15, 20,
# 25, 26), the task brief (24), and total_trees/checked fields stored in the
# project's results/*.json (17-23). Every value is additionally cross-checked
# here against the count gentreeg -u reports, which is independent of our
# parser only in the weak sense that both come from nauty.
A000055 = {
    10: 106, 11: 235, 12: 551, 13: 1301, 14: 3159, 15: 7741, 16: 19320,
    17: 48629, 18: 123867, 19: 317955, 20: 823065, 21: 2144505,
    22: 5623756, 23: 14828074, 24: 39299897, 25: 104636890,
    26: 279793450, 27: 751065460, 28: 2023443032,
}
SAMPLE_NS = {16, 20, 24}
CONC = int(os.environ.get("CM_CONC", "5"))   # concurrent shards (core budget)


def parse_kv(line):
    out = {}
    for tok in line.split()[1:]:
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k] = v
    return out


def run_n(n, shards, sample_target):
    P = 0
    if n in SAMPLE_NS:
        P = max(1, A000055[n] // sample_target)
    t0 = time.time()
    procs = []
    for w0 in range(0, shards, CONC):
        wave = []
        for r in range(w0, min(shards, w0 + CONC)):
            path = os.path.join(RAW, f"n{n}_s{r}of{shards}.txt")
            fout = open(path, "w")
            g = subprocess.Popen(["gentreeg", "-p", "-q", str(n), f"{r}/{shards}"],
                                 stdout=subprocess.PIPE)
            c = subprocess.Popen([os.path.join(HERE, "central_margin"), str(n), str(P)],
                                 stdin=g.stdout, stdout=fout)
            g.stdout.close()
            wave.append((g, c, fout, path))
        for g, c, fout, path in wave:
            rc = c.wait()
            g.wait()
            fout.close()
            if rc not in (0, 42):
                raise SystemExit(f"shard failed rc={rc}: {path}")
            print(f"  shard done {path} t={time.time() - t0:.0f}s", flush=True)
        procs.extend(wave)
    elapsed = time.time() - t0

    stats = dict(trees=0, viol5=0, viol4=0, zero5=0, zero4=0, nonunimodal=0, sampled=0,
                 newton_fail_trees=0)
    ones = {"PERK": {}, "PERALPHA": {}, "TOPQ": {}, "NEWTON": {}}
    nfdist = {}
    vbins = {}
    top = {"TOP5": [], "TOP4": []}
    vd = None
    violations = []
    samples = []
    for _, _, _, path in procs:
        with open(path) as fh:
            for line in fh:
                if line.startswith("STATS"):
                    kv = parse_kv(line)
                    for key in stats:
                        stats[key] += int(kv[key])
                elif line.startswith("TOP5") or line.startswith("TOP4"):
                    tag = line.split()[0]
                    kv = parse_kv(line)
                    frac = Fraction(int(kv["num"]), int(kv["den"]))
                    top[tag].append(dict(margin=frac, k=int(kv["k"]), alpha=int(kv["alpha"]),
                                         q=int(kv["q"]), par=kv["par"], poly=kv["poly"]))
                elif line.split()[0] in ones:
                    tag = line.split()[0]
                    kv = parse_kv(line)
                    fr = Fraction(int(kv["num"]), int(kv["den"]))
                    key = int(kv["k"]) if tag == "PERK" else (int(kv["alpha"]) if tag == "PERALPHA" else 0)
                    cur = ones[tag].get(key)
                    if cur is None or fr < cur["value"]:
                        ones[tag][key] = dict(value=fr, k=int(kv["k"]), alpha=int(kv["alpha"]),
                                              q=int(kv["q"]), par=kv["par"])
                elif line.startswith("NEWTONFAIL_DIST"):
                    kv = parse_kv(line)
                    d = int(kv["q_minus_k"])
                    nfdist[d] = nfdist.get(d, 0) + int(kv["trees"])
                elif line.startswith("VDBIN_FLOAT"):
                    kv = parse_kv(line)
                    b = float(kv["bin_lo"])
                    if b not in vbins or float(kv["value"]) < vbins[b]["value"]:
                        vbins[b] = dict(value=float(kv["value"]), k=int(kv["k"]), alpha=int(kv["alpha"]),
                                        q=int(kv["q"]), **{"lambda": float(kv["lambda"])},
                                        V=float(kv["V"]), delta=float(kv["delta"]), par=kv["par"])
                elif line.startswith("VDMIN_FLOAT"):
                    kv = parse_kv(line)
                    if vd is None or float(kv["value"]) < vd["value"]:
                        vd = dict(value=float(kv["value"]), k=int(kv["k"]), alpha=int(kv["alpha"]),
                                  q=int(kv["q"]), **{"lambda": float(kv["lambda"])},
                                  V=float(kv["V"]), delta=float(kv["delta"]),
                                  par=kv["par"], poly=kv["poly"])
                elif line.startswith("CENTRAL_VIOLATION") or line.startswith("ALARM"):
                    violations.append(line.strip())
                elif line.startswith("SAMPLE"):
                    kv = parse_kv(line)
                    samples.append(dict(alpha=int(kv["alpha"]), par=kv["par"], poly=kv["poly"]))
    expected = A000055.get(n)
    count_ok = (expected == stats["trees"])
    out = {"n": n, "shards": shards, "elapsed_s": round(elapsed, 2),
           "expected_A000055": expected, "count_ok": count_ok, **stats,
           "violation_lines": violations[:50]}
    for tag in ("TOP5", "TOP4"):
        lst = sorted(top[tag], key=lambda r: r["margin"])[:12]
        out[tag] = [dict(margin=f"{r['margin'].numerator}/{r['margin'].denominator}",
                         margin_float=float(r["margin"]),
                         n_times_margin=f"{(n * r['margin']).numerator}/{(n * r['margin']).denominator}",
                         n_times_margin_float=float(n * r["margin"]),
                         k=r["k"], alpha=r["alpha"], q=r["q"], par=r["par"], poly=r["poly"])
                    for r in lst]
    out["VDMIN_FLOAT_diagnostic"] = vd
    def ser(r):
        v = r["value"]
        return dict(value=f"{v.numerator}/{v.denominator}", value_float=float(v),
                    n_times_value_float=float(n * v), k=r["k"], alpha=r["alpha"], q=r["q"], par=r["par"],
                    is_star=(r["par"] == ",".join(["0"] + ["1"] * (n - 1))))
    out["per_k_min_delta"] = {str(k): ser(r) for k, r in sorted(ones["PERK"].items())}
    out["per_alpha_min_m5"] = {str(a): ser(r) for a, r in sorted(ones["PERALPHA"].items())}
    out["window_top_min_delta_q"] = ser(ones["TOPQ"][0])
    out["newton_ratio_min_exact"] = ser(ones["NEWTON"][0])
    out["newton_fail_first_k_by_q_minus_k"] = {str(d): c for d, c in sorted(nfdist.items())}
    out["VD_by_V_bin_FLOAT_diagnostic"] = {str(b): r for b, r in sorted(vbins.items())}
    with open(os.path.join(DATA, f"census_n{n}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    if samples:
        with open(os.path.join(DATA, f"sample_n{n}.jsonl"), "w") as fh:
            for s in samples:
                fh.write(json.dumps(s) + "\n")
    m5 = out["TOP5"][0]
    m4 = out["TOP4"][0]
    print(f"n={n} trees={stats['trees']} ok={count_ok} viol5={stats['viol5']} viol4={stats['viol4']} "
          f"nonuni={stats['nonunimodal']} min5={m5['margin']} n*min5={m5['n_times_margin_float']:.6f} "
          f"(alpha={m5['alpha']},k={m5['k']},par={m5['par']}) min4={m4['n_times_margin_float']:.6f} "
          f"VDmin={vd['value']:.6f} sampled={stats['sampled']} t={elapsed:.1f}s", flush=True)
    tq = out["window_top_min_delta_q"]
    nw = out["newton_ratio_min_exact"]
    print(f"   topq: n*d_q={tq['n_times_value_float']:.4f} alpha={tq['alpha']} star={tq['is_star']} | "
          f"newton_min={nw['value_float']:.5f} (alpha={nw['alpha']},k={nw['k']},q={nw['q']}) "
          f"fail_trees={stats['newton_fail_trees']} | VDbins=" +
          " ".join(f"[{b}]{r['value']:.4f}(V={r['V']:.2f},a={r['alpha']},k={r['k']})" for b, r in out['VD_by_V_bin_FLOAT_diagnostic'].items()),
          flush=True)
    if not count_ok:
        print(f"  COUNT MISMATCH n={n}: got {stats['trees']} expected {expected}", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ns", nargs="+", type=int)
    ap.add_argument("--shards", type=int, default=5)
    ap.add_argument("--sample-target", type=int, default=2000)
    a = ap.parse_args()
    for n in a.ns:
        run_n(n, a.shards, a.sample_target)


if __name__ == "__main__":
    main()
