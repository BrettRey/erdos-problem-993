#!/usr/bin/env python3
"""Write manifest.yaml: model, repo HEAD, SHA-256 of inputs read, scripts and data."""
import hashlib, os, subprocess, glob, datetime
ROOT = "/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993"
OUT = os.path.join(ROOT, "runs/attack-20260930/route-lclt")
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()
inputs = [
    "/Users/brettreynolds/projects/LLM-CLI-projects/literature/fang_2026_unimodality_large_forests.pdf",
    ROOT + "/notes/hard_core_attack.md",
    ROOT + "/notes/tie_fugacity_leaf_decomp_analysis.md",
    ROOT + "/notes/mode_mean_tiepoint_2026-02-18.md",
    ROOT + "/notes/real_collar_conjecture_2026-07-16.md",
    ROOT + "/notes/mason_free_count_reformulation_2026-09-02.md",
    ROOT + "/notes/fang_lu_nevo_yao_zheng_2026_large_forests_2026-09-25.md",
    ROOT + "/notes/why_trees_resist_2026-07-16.md",
    ROOT + "/notes/beurling_majorant_difference_test_lead_2026-08-12.md",
    ROOT + "/gpt_attack/bridge_window_unimodality/outcome_2026-07-16/bridge_lemma_report.md",
    ROOT + "/scripts/mason_b_census_gentreeg_20260902.py",
]
head = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
lines = ["# route-lclt run manifest",
         f"generated: {datetime.datetime.now().isoformat(timespec='seconds')}",
         "model: claude-opus-5-5 (Claude Opus 5.5, workflow subagent)",
         f"repo_head: {head}",
         "report_md: NOT WRITTEN (harness blocked report-file writes; full report returned inline to the parent)",
         "tools: python3 (system, 3.14) with mpmath 1.3.0, numpy 2.4.2, scipy 1.17.0; nauty gentreeg",
         "exactness: kt12_census, check_identity_exact are exact (int/Fraction); kt12_families LC sign exact, Gamma/rho 60-digit mpmath diagnostics; kt3 80-digit mpmath diagnostic; kt4/kt5 float64 diagnostics",
         "inputs_read:"]
for p in inputs:
    lines.append(f"  - path: {p}\n    sha256: {sha(p)}")
lines.append("outputs:")
for p in sorted(glob.glob(OUT + "/scripts/*.py") + glob.glob(OUT + "/data/*") + [OUT + "/fang.txt"]):
    if os.path.isfile(p):
        lines.append(f"  - path: {os.path.relpath(p, OUT)}\n    sha256: {sha(p)}")
open(os.path.join(OUT, "manifest.yaml"), "w").write("\n".join(lines) + "\n")
print("manifest written", len(lines))
