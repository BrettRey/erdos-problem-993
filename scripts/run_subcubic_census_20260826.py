#!/usr/bin/env python3
"""Subcubic (max degree <= 3) LC census, n = 33..36 (2026-08-26).

Motivation: the matroid-intersection reading of Schweitzer arXiv:2608.23262.
A tree's independence system is an intersection of exactly Delta(T) matroids
and no fewer (notes/matroid_intersection_number_lane_2026-08-26.md), so
subcubic trees are the tree realizations of the three-matroid rung, where
Schweitzer's general LC failure lives. The complete n <= 32 census has NO
failure with Delta <= 3 (minimum over all 1,230 stored failures is Delta = 4).
This sweep pushes the subcubic frontier past the general census.

Engine: scripts/lc_census_alpha (n <= 38, 128-bit Turan products), run with
no alpha filter. Producer: gentreeg -D3. Expected totals are OEIS A000672
(b-file fetched 2026-08-26; offset validated against locally generated counts
at n = 10, 15, 20, 24 exactly).

Gates before the ladder:
  A (positive control): the two stored n=26 LC failures, re-rooted to
    parent-array form, must each produce a FAIL line from lc_census_alpha.
  B/C: full subcubic runs at n=26 and n=28 must match A000672 exactly and
    produce zero failures (every stored failure at those orders has
    Delta >= 4, so a subcubic failure would be an engine/pipe bug).

Restartable: shards whose file ends in a STATS line are skipped.

Usage:
  python3 scripts/run_subcubic_census_20260826.py            # gates + ladder
  python3 scripts/run_subcubic_census_20260826.py --merge-only
"""
import argparse
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CENSUS = os.path.join(HERE, "lc_census_alpha")
OUTDIR = os.path.join(ROOT, "results", "subcubic_census_20260826")

# OEIS A000672 (number of trees with all degrees <= 3), b-file fetched
# 2026-08-26. Offset validated: gentreeg -D3 -u gives 37, 1132, 52233,
# 1265579 at n = 10, 15, 20, 24, matching b(10), b(15), b(20), b(24).
EXPECTED = {
    26: 6408674,
    28: 32935002,
    33: 2076217086,
    34: 4790669518,
    35: 11077270335,
    36: 25664800714,
    37: 59574334817,
}

LADDER = [33, 34, 35, 36]
MOD = 32
WORKERS = 7


def shard_path(n, s):
    return os.path.join(OUTDIR, f"n{n}_s{s}of{MOD}.txt")


def shard_done(path):
    if not os.path.exists(path):
        return False
    with open(path, "rb") as fh:
        try:
            fh.seek(-512, 2)
        except OSError:
            fh.seek(0)
        return b"\nSTATS " in b"\n" + fh.read()


def run_shard(n, s):
    out = shard_path(n, s)
    if shard_done(out):
        return (s, "cached")
    cmd = (f"gentreeg -D3 -p -q {n} {s}/{MOD} 2>/dev/null | "
           f"'{CENSUS}' {n} > '{out}.part' && mv '{out}.part' '{out}'")
    r = subprocess.run(["bash", "-c", cmd])
    return (s, "done" if r.returncode == 0 else f"EXIT{r.returncode}")


def merge(n, mod=MOD, prefix="n"):
    trees = failures = alarms = 0
    fails, missing = [], []
    for s in range(mod):
        p = os.path.join(OUTDIR, f"{prefix}{n}_s{s}of{mod}.txt")
        if not shard_done(p):
            missing.append(s)
            continue
        for line in open(p):
            if line.startswith("STATS"):
                d = dict(f.split("=") for f in line.split()[1:])
                trees += int(d["trees"])
                failures += int(d["failures"])
                alarms += int(d["alarms"])
            elif line.startswith(("FAIL", "ALARM")):
                fails.append(line.strip())
    exp = EXPECTED[n]
    ok = not missing and trees == exp
    summary = {
        "n": n, "mod": mod, "max_degree": 3, "shards_missing": missing,
        "trees": trees, "expected_trees_A000672": exp, "complete": ok,
        "failures": failures, "alarms": alarms, "fail_lines": fails,
    }
    with open(os.path.join(OUTDIR, f"summary_n{n}.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(f"[subcubic] n={n}: trees={trees} expected={exp} complete={ok} "
          f"failures={failures} alarms={alarms}", flush=True)
    if alarms:
        print("!" * 70 + "\nALARM: NON-UNIMODAL SUBCUBIC TREE — inspect NOW\n" + "!" * 70)
    return summary


def gate_positive_control():
    """The two stored n=26 failures must FAIL under lc_census_alpha."""
    import networkx as nx
    stored = json.load(open(os.path.join(ROOT, "results", "analysis_n26.json")))
    lines = []
    for f in stored["lc_failures"]:
        G = nx.from_graph6_bytes(f["graph6"].encode())
        assert nx.is_tree(G) and G.number_of_nodes() == 26
        # Relabel by BFS so parent[i] < i, 1-indexed, parent[1] = 0.
        order = list(nx.bfs_tree(G, min(G.nodes())).nodes())
        idx = {v: i + 1 for i, v in enumerate(order)}
        par = [0] * 27
        for u, v in nx.bfs_edges(G, min(G.nodes())):
            par[idx[v]] = idx[u]
        lines.append(" ".join(str(par[i]) for i in range(1, 27)))
    inp = ("\n".join(lines) + "\n").encode()
    r = subprocess.run([CENSUS, "26"], input=inp, capture_output=True)
    out = r.stdout.decode()
    n_fail = out.count("FAIL ")
    if n_fail != 2:
        print(out)
        raise SystemExit(f"gate A FAILED: expected 2 FAIL lines, got {n_fail}")
    for f, line in zip(stored["lc_failures"], out.splitlines()):
        want = ",".join(str(c) for c in f["poly"])
        if f"poly={want}" not in out:
            raise SystemExit(f"gate A FAILED: stored poly {want} not reproduced")
    print("[gate A] positive control passed: both stored n=26 failures reproduced "
          "with matching polynomials", flush=True)


def run_order(n):
    todo = [s for s in range(MOD) if not shard_done(shard_path(n, s))]
    print(f"[subcubic] n={n}: {len(todo)}/{MOD} shards to run", flush=True)
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(run_shard, n, s): s for s in todo}
        for fut in as_completed(futs):
            s, status = fut.result()
            if status.startswith("EXIT"):
                raise SystemExit(f"n={n} shard {s}: {status}")
    return merge(n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--merge-only", action="store_true")
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    if args.merge_only:
        for n in LADDER:
            merge(n)
        return

    gate_positive_control()
    for n, gate in ((26, "B"), (28, "C")):
        s = run_order(n)
        if not s["complete"] or s["failures"] != 0 or s["alarms"] != 0:
            raise SystemExit(f"gate {gate} FAILED at n={n}: {s}")
        print(f"[gate {gate}] n={n} subcubic: complete, zero failures, "
              f"total matches A000672", flush=True)

    for n in LADDER:
        s = run_order(n)
        if not s["complete"]:
            raise SystemExit(f"n={n} incomplete: {s['shards_missing']}")
    print("[subcubic] ladder complete", flush=True)


if __name__ == "__main__":
    main()
