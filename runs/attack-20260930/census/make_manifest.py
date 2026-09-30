#!/usr/bin/env python3
"""Write manifest.yaml: model, git HEAD, SHA-256 of inputs read and of outputs."""
import glob
import hashlib
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
LIT = "/Users/brettreynolds/projects/LLM-CLI-projects/literature"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


head = subprocess.check_output(["git", "-C", PROJ, "rev-parse", "HEAD"], text=True).strip()
inputs = [
    os.path.join(PROJ, "scripts/lc_census.c"),
    os.path.join(PROJ, "scripts/probe_absorption_margin_20260811.py"),
    os.path.join(PROJ, "scripts/search_break_depth_20260814.py"),
    os.path.join(PROJ, "notes/mason_free_count_reformulation_2026-09-02.md"),
    os.path.join(PROJ, "notes/literature/arxiv_2603_17114.txt"),
    os.path.join(PROJ, "notes/literature/arxiv_2511_00334.txt"),
    os.path.join(PROJ, "paper/poisson_binomial/main.tex"),
    os.path.join(LIT, "fang_2026_unimodality_large_forests.pdf"),
]
outputs = sorted(
    [p for p in glob.glob(os.path.join(HERE, "*")) if os.path.isfile(p) and not p.endswith("manifest.yaml")]
    + glob.glob(os.path.join(HERE, "data", "*.json*"))
    + glob.glob(os.path.join(HERE, "data", "hillclimb", "*.json"))
    + glob.glob(os.path.join(HERE, "sources", "*"))
)
lines = [
    "# manifest for the central-window margin census (attack-20260930/census)",
    "model: claude-opus-5-5 (Claude Code subagent)",
    "date: 2026-09-30",
    f"project_git_head: {head}",
    "tools: [gentreeg (nauty, /opt/homebrew/bin/gentreeg), cc -O3 (central_margin.c, central_margin_q.c), python3 (Fractions; networkx for the independent cross-check)]",
    "dependence_cluster: single agent; exhaustive results cross-checked by an independent Python/networkx enumeration for n<=14; tilted diagnostics cross-checked C (Newton) vs Python (bisection)",
    "inputs_read_sha256:",
]
for p in inputs:
    lines.append(f"  - path: {p}\n    sha256: {sha(p)}")
lines.append("outputs_sha256:")
for p in outputs:
    lines.append(f"  - path: {os.path.relpath(p, HERE)}\n    sha256: {sha(p)}")
lines.append("raw_shard_outputs: raw/ and raw_q/ (not hashed individually; merged results in data/)")
open(os.path.join(HERE, "manifest.yaml"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines[:12]))
print(f"... {len(outputs)} outputs hashed")
