#!/usr/bin/env python3
"""Write manifest.yaml: model, project HEAD, SHA-256 of inputs and outputs."""
import hashlib, subprocess, os, glob
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = "/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993"
RUN = PROJ + "/runs/attack-20260930"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
head = subprocess.run(["git", "-C", PROJ, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
inputs = [RUN + "/route-lclt/RETURN.json", RUN + "/route-lclt/scripts/kt4_mixing_clt.py", RUN + "/route-lclt/scripts/lclt_lib.py",
          RUN + "/route-lclt/data/kt5_census_n8_16.jsonl", RUN + "/audit-a/RETURN.json", RUN + "/audit-b/RETURN.json",
          RUN + "/audit-b/n0_chains.py", RUN + "/audit-b/data/lemma42_certificates.json", RUN + "/audit-b/data/n0_chains_stdout.txt",
          RUN + "/audit-a/fang.txt", RUN + "/route-freecount/RETURN.json"]
outs = sorted(f for f in glob.glob(HERE + "/*") if os.path.isfile(f) and not f.endswith("manifest.yaml"))
with open(HERE + "/manifest.yaml", "w") as fh:
    fh.write("lane: lclt-check (independent check of route-lclt binomial-smoothing certificate; hybrid N1)\n")
    fh.write("model: claude-opus-5-5 (Claude Code subagent)\n")
    fh.write(f"project_git_head: {head}\n")
    fh.write("date: 2026-09-30\n")
    fh.write("arithmetic: exact_AC.py / exact_large.py exact (Fraction, int); kernel_constants.py mpmath 40-60 digits (diagnostic) + python-flint arb enclosure of int|phi'''| and c_*; n1_hybrid.py mpmath order-of-magnitude (diagnostic, imports audit-b/n0_chains.py constants unchanged); kt5_sharp.py float64 diagnostic reusing route-lclt joint_JB/cert_one\n")
    fh.write("report_note: full report is the structured return (RETURN.json written by the harness); no .md report written\n")
    fh.write("inputs:\n")
    for p in inputs:
        if os.path.exists(p): fh.write(f"  - path: {p}\n    sha256: {sha(p)}\n")
    fh.write("outputs:\n")
    for p in outs:
        fh.write(f"  - path: {p}\n    sha256: {sha(p)}\n")
print(open(HERE + "/manifest.yaml").read())
