"""Run r1_census over all trees for each n in [N0, N1], SHARDS parallel shards
(gentreeg res/mod), writing raw/n{n}_s{r}of{m}.txt per shard as it finishes.
Usage: python3 run_census.py N0 N1 SHARDS [BINARY] [RAWDIR]
(BINARY default r1_census; r1_ext adds outside-window EXT lines, identical in-window code.)
"""
import subprocess
import sys
import time

HERE = "/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/runs/attack-20260930/r1-exhaustive"


def main():
    n0, n1, m = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    binary = sys.argv[4] if len(sys.argv) > 4 else "r1_census"
    rawdir = sys.argv[5] if len(sys.argv) > 5 else "raw"
    for n in range(n0, n1 + 1):
        t0 = time.time()
        procs = []
        for r in range(m):
            out = open(f"{HERE}/{rawdir}/n{n}_s{r}of{m}.txt.part", "w")
            g = subprocess.Popen(["gentreeg", "-p", "-q", str(n), f"{r}/{m}"], stdout=subprocess.PIPE)
            c = subprocess.Popen([f"{HERE}/{binary}", str(n)], stdin=g.stdout, stdout=out)
            g.stdout.close()
            procs.append((r, g, c, out))
        for r, g, c, out in procs:
            rc = c.wait(); g.wait(); out.close()
            import os
            os.rename(f"{HERE}/{rawdir}/n{n}_s{r}of{m}.txt.part", f"{HERE}/{rawdir}/n{n}_s{r}of{m}.txt")
            print(f"n={n} shard={r}/{m} rc={rc} gen_rc={g.returncode}", flush=True)
        print(f"n={n} done t={time.time() - t0:.1f}s at {time.strftime('%H:%M:%S')}", flush=True)


if __name__ == "__main__":
    main()
