"""Write manifest.yaml: model, git HEAD, SHA-256 of inputs read, scripts, data."""
import glob
import hashlib
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = "/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993"
LIT = "/Users/brettreynolds/projects/LLM-CLI-projects/literature"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def main():
    head = subprocess.check_output(["git", "-C", PROJ, "rev-parse", "HEAD"], text=True).strip()
    inputs = [f"{PROJ}/notes/{n}" for n in [
        "mason_free_count_reformulation_2026-09-02.md",
        "attack3_mean_lift_2026-02-19.md",
        "decimation_weighted_whnc_2026-02-18.md",
        "fang_lu_nevo_yao_zheng_2026_large_forests_2026-09-25.md",
        "ratio_dominance_discovery_2026-03-01.md",
        "condition_C_proof_structure_2026-02-28.md",
        "scc_false_n28_2026-03-01.md",
        "d22_global_switching_obstruction_packet_2026-07-11.md",
        "joint_density_analytic_sprint_2026-09-04.md",
        "arxiv_astra_transfer_2026-09-04.md",
    ]] + [f"{LIT}/fang_2026_unimodality_large_forests.pdf", f"{HERE}/fang.txt"]
    shards = sorted(glob.glob(f"{PROJ}/results/lc_census_20260814/n*_s*of64.txt"))
    h = hashlib.sha256()
    for s in shards:
        with open(s, "rb") as f:
            h.update(f.read())
    scripts = sorted(glob.glob(f"{HERE}/*.py"))
    data = sorted(glob.glob(f"{HERE}/data/*"))
    lines = [
        "# Manifest for route-freecount (attack-20260930)",
        "model: claude-opus-5-5",
        f"git_head: {head}",
        "report: REPORT.md not written (the harness blocks report files from subagents); the full report is the structured return of this run",
        "compute: python3 (exact integer / Fraction arithmetic), at most 2 cores, each run under 20 minutes",
        "inputs_read:",
    ]
    for p in inputs:
        lines.append(f"  - path: {p}\n    sha256: {sha(p)}")
    lines.append("census_shards_concatenated:")
    lines.append(f"  glob: {PROJ}/results/lc_census_20260814/n*_s*of64.txt")
    lines.append(f"  count: {len(shards)}")
    lines.append(f"  sha256_of_sorted_concatenation: {h.hexdigest()}")
    lines.append("scripts:")
    for p in scripts:
        lines.append(f"  - path: {os.path.relpath(p, HERE)}\n    sha256: {sha(p)}")
    lines.append("data:")
    for p in data:
        lines.append(f"  - path: {os.path.relpath(p, HERE)}\n    sha256: {sha(p)}")
    with open(f"{HERE}/manifest.yaml", "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
