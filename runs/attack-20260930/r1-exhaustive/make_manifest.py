"""Write manifest.yaml for the r1-exhaustive lane: model, git HEAD, SHA-256 of
inputs read, scripts, binaries and every raw/result file."""
import glob
import hashlib
import os
import subprocess
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = "/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993"
ATT = f"{PROJ}/runs/attack-20260930"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    head = subprocess.check_output(["git", "-C", PROJ, "rev-parse", "HEAD"], text=True).strip()
    inputs = [f"{ATT}/route-freecount/pv_lib.py", f"{ATT}/route-freecount/RETURN.json",
              f"{ATT}/route-freecount/k1_census.py", f"{ATT}/route-freecount/k2_aggregate.py",
              f"{ATT}/route-freecount/k4_hillclimb.py", f"{ATT}/route-freecount/k3_families.py",
              f"{ATT}/census/central_margin.c", "/opt/homebrew/bin/gentreeg"]
    own = sorted(glob.glob(f"{HERE}/*.c") + glob.glob(f"{HERE}/*.py") + [f"{HERE}/r1_census", f"{HERE}/r1_ext"])
    outs = sorted(set(glob.glob(f"{HERE}/*.json") + glob.glob(f"{HERE}/*.log") + glob.glob(f"{HERE}/*.txt")
                      + glob.glob(f"{HERE}/raw/*.txt") + glob.glob(f"{HERE}/raw_ext/*.txt")))
    L = ["# manifest for runs/attack-20260930/r1-exhaustive",
         f"written: {time.strftime('%Y-%m-%dT%H:%M:%S%z')}",
         "model: claude-opus-5-5 (Claude Code subagent)",
         f"git_head: {head}",
         "compiler: " + subprocess.check_output(["cc", "--version"], text=True).splitlines()[0],
         "build: cc -O3 -o r1_census r1_census.c ; cc -O3 -o r1_ext r1_ext.c",
         "inputs:"]
    for p in inputs:
        L.append(f"  - path: {p}\n    sha256: {sha(p)}")
    L.append("scripts_and_binaries:")
    for p in own:
        L.append(f"  - path: {os.path.relpath(p, HERE)}\n    sha256: {sha(p)}")
    L.append("outputs:")
    for p in outs:
        if p.endswith(".part"):
            continue
        L.append(f"  - path: {os.path.relpath(p, HERE)}\n    sha256: {sha(p)}")
    with open(f"{HERE}/manifest.yaml", "w") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
