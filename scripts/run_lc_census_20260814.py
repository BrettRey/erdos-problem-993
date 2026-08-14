#!/usr/bin/env python3
"""Sharded driver for the n<=32 LC-failure census (lc_census.c).

Runs gentreeg res/mod shards through the C census engine with a worker
pool. Restartable: a shard whose output file ends in a STATS line is
skipped. Completeness is verified two ways: every shard must carry its
STATS trailer, and the summed tree count must equal an independent
`gentreeg -u <n>` count (cached in the output directory, so the expected
total never comes from anyone's memory).

Usage:
  python3 scripts/run_lc_census_20260814.py --n 29 --mod 64 --workers 9
  python3 scripts/run_lc_census_20260814.py --n 29 --merge-only
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
CENSUS = os.path.join(HERE, "lc_census")
OUTDIR = os.path.join(ROOT, "results", "lc_census_20260814")


def expected_count(n):
    cache = os.path.join(OUTDIR, f"expected_n{n}.txt")
    if os.path.exists(cache):
        return int(open(cache).read().strip())
    print(f"[driver] counting trees at n={n} via gentreeg -u ...", flush=True)
    r = subprocess.run(["gentreeg", "-u", str(n)], capture_output=True,
                       text=True, check=True)
    m = re.search(r">Z\s+(\d+)\s+trees generated", r.stderr + r.stdout)
    if not m:
        raise RuntimeError(f"cannot parse gentreeg -u output: {r.stderr!r}")
    cnt = int(m.group(1))
    with open(cache, "w") as fh:
        fh.write(str(cnt))
    return cnt


def shard_path(n, s, mod):
    return os.path.join(OUTDIR, f"n{n}_s{s}of{mod}.txt")


def shard_done(path):
    if not os.path.exists(path):
        return False
    with open(path, "rb") as fh:
        try:
            fh.seek(-256, 2)
        except OSError:
            fh.seek(0)
        return b"\nSTATS " in b"\n" + fh.read()


def run_shard(n, s, mod):
    out = shard_path(n, s, mod)
    if shard_done(out):
        return (s, "cached")
    cmd = (f"gentreeg -p -q {n} {s}/{mod} 2>/dev/null | "
           f"'{CENSUS}' {n} > '{out}.part' && mv '{out}.part' '{out}'")
    r = subprocess.run(["bash", "-c", cmd])
    if r.returncode != 0:
        return (s, f"EXIT{r.returncode}")
    return (s, "done")


def merge(n, mod):
    trees = failures = alarms = 0
    fails = []
    missing = []
    for s in range(mod):
        p = shard_path(n, s, mod)
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
    exp = expected_count(n)
    ok = not missing and trees == exp
    summary = {
        "n": n, "mod": mod, "shards_missing": missing,
        "trees": trees, "expected_trees": exp, "complete": ok,
        "failures": failures, "alarms": alarms,
        "fail_lines": fails,
    }
    out = os.path.join(OUTDIR, f"summary_n{n}.json")
    with open(out, "w") as fh:
        json.dump(summary, fh, indent=1)
    print(f"[driver] n={n}: trees={trees} expected={exp} complete={ok} "
          f"failures={failures} alarms={alarms} -> {out}", flush=True)
    if alarms:
        print("!" * 70 + "\nALARM: NON-UNIMODAL TREE(S) FOUND — see summary\n"
              + "!" * 70, flush=True)
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--mod", type=int, default=64)
    ap.add_argument("--workers", type=int, default=9)
    ap.add_argument("--merge-only", action="store_true")
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    if not args.merge_only:
        if not os.path.exists(CENSUS):
            sys.exit("build first: cc -O3 -march=native -o scripts/lc_census "
                     "scripts/lc_census.c")
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(run_shard, args.n, s, args.mod): s
                    for s in range(args.mod)}
            done = 0
            for fut in as_completed(futs):
                s, status = fut.result()
                done += 1
                if status.startswith("EXIT"):
                    print(f"[driver] shard {s}: {status}", flush=True)
                if done % 8 == 0 or done == args.mod:
                    print(f"[driver] n={args.n}: {done}/{args.mod} shards",
                          flush=True)
    merge(args.n, args.mod)


if __name__ == "__main__":
    main()
