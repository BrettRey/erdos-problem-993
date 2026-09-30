#!/usr/bin/env python3
"""Driver for central_margin_q.c: exhaustive Fang-style concavity-at-mean statistic
Q_k = V_k^{3/2}(2p_k - p_{k-1} - p_{k+1}) (FLOAT DIAGNOSTIC) over W5, by V bin.

Usage: python3 run_q.py N [N ...] [--shards S] [--conc C]
Writes raw_q/n{n}_s{r}of{S}.txt and data/qstat_n{n}.json.
"""
import argparse
import json
import os
import subprocess
import time

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw_q")
DATA = os.path.join(HERE, "data")
os.makedirs(RAW, exist_ok=True)
A000055 = {10: 106, 11: 235, 12: 551, 13: 1301, 14: 3159, 15: 7741, 16: 19320, 17: 48629,
           18: 123867, 19: 317955, 20: 823065, 21: 2144505, 22: 5623756, 23: 14828074,
           24: 39299897, 25: 104636890, 26: 279793450}


def kv(line):
    return dict(t.split("=", 1) for t in line.split()[1:] if "=" in t)


def run(n, shards, conc):
    t0 = time.time()
    paths = []
    for w0 in range(0, shards, conc):
        wave = []
        for r in range(w0, min(shards, w0 + conc)):
            path = os.path.join(RAW, f"n{n}_s{r}of{shards}.txt")
            fout = open(path, "w")
            g = subprocess.Popen(["gentreeg", "-p", "-q", str(n), f"{r}/{shards}"], stdout=subprocess.PIPE)
            c = subprocess.Popen([os.path.join(HERE, "central_margin_q"), str(n), "0"], stdin=g.stdout, stdout=fout)
            g.stdout.close()
            wave.append((g, c, fout, path))
        for g, c, fout, path in wave:
            rc = c.wait(); g.wait(); fout.close()
            if rc not in (0, 42):
                raise SystemExit(f"fail {path} rc={rc}")
            paths.append(path)
    trees = 0
    qneg = 0
    bins = {}
    for p in paths:
        for line in open(p):
            if line.startswith("STATS"):
                trees += int(kv(line)["trees"])
            elif line.startswith("QNONPOS_TREES"):
                qneg += int(kv(line)["count"])
            elif line.startswith("QBIN_FLOAT"):
                d = kv(line)
                b = d["bin_lo"]
                if b not in bins or float(d["Q"]) < bins[b]["Q"]:
                    bins[b] = dict(Q=float(d["Q"]), V=float(d["V"]), k=int(d["k"]), alpha=int(d["alpha"]),
                                   q=int(d["q"]), delta=float(d["delta"]), par=d["par"], **{"lambda": float(d["lambda"])})
    out = dict(n=n, trees=trees, count_ok=(trees == A000055[n]), q_nonpositive_trees=qneg,
               Q_by_V_bin_FLOAT=bins, Q_min_overall=min(b["Q"] for b in bins.values()),
               fang_limit=1 / (2 * 3.141592653589793) ** 0.5, elapsed_s=round(time.time() - t0, 1))
    json.dump(out, open(os.path.join(DATA, f"qstat_n{n}.json"), "w"), indent=1)
    print(f"n={n} trees={trees} ok={out['count_ok']} Qnonpos_trees={qneg} Qmin={out['Q_min_overall']:.4f} | " +
          " ".join(f"[{b}]{r['Q']:.4f}(V={r['V']:.2f},a={r['alpha']},k={r['k']},q={r['q']})" for b, r in bins.items()) +
          f" t={out['elapsed_s']}s", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("ns", nargs="+", type=int)
    ap.add_argument("--shards", type=int, default=1)
    ap.add_argument("--conc", type=int, default=1)
    a = ap.parse_args()
    for n in a.ns:
        run(n, a.shards, a.conc)
