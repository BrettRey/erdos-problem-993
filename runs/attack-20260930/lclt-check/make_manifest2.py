#!/usr/bin/env python3
"""Write manifest.yaml for the resumed lclt-check run: model, project HEAD, SHA-256 of
key inputs and every output file in this directory (including rerun/)."""
import hashlib, subprocess, os, glob
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = "/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993"
RUN = PROJ + "/runs/attack-20260930"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
head = subprocess.run(["git", "-C", PROJ, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
inputs = [RUN + "/route-lclt/RETURN.json", RUN + "/route-lclt/fang.txt", RUN + "/route-lclt/scripts/kt4_mixing_clt.py",
          RUN + "/route-lclt/scripts/lclt_lib.py", RUN + "/audit-a/RETURN.json", RUN + "/audit-b/RETURN.json",
          RUN + "/audit-b/n0_chains.py", RUN + "/SUMMARY.md",
          "/Users/brettreynolds/projects/LLM-CLI-projects/literature/fang_2026_unimodality_large_forests.pdf"]
outs = sorted(f for f in glob.glob(HERE + "/*") + glob.glob(HERE + "/rerun/*")
              if os.path.isfile(f) and not f.endswith("manifest.yaml") and "__pycache__" not in f)
with open(HERE + "/manifest.yaml", "w") as fh:
    fh.write("lane: lclt-check (independent check of route-lclt claims A, B, C; hybrid N1; Lemma L gap)\n")
    fh.write("model: claude-opus-5-5 (Claude Code subagent; first attempt stopped after 19 min, resumed 2026-09-30 15:50 EDT)\n")
    fh.write(f"project_git_head: {head}\n")
    fh.write("date: 2026-09-30\n")
    fh.write("provenance: files at top level were written by the stopped first attempt and re-verified (rerun/ holds the reruns; "
             "census and exact_large outputs byte-identical). New in the resumed run: kernel_check_exact.py, tv_exact_point.py, "
             "tv_sup_scan.py, n1_check.py, exact_large2.py, census_range.py (n=13,14), make_manifest2.py, rerun/*. kt5_sharp.py (first attempt) reuses route-lclt joint_JB/cert_one, so its n=16 rerun is NOT independent of route-lclt code. n1_hybrid.py/json are the first attempt's; "
             "n1_check.py recomputes them without importing audit-b code. make_manifest.py/return_draft.json are first-attempt leftovers.\n")
    fh.write("arithmetic: EXACT (Fraction/int): exact_AC.py, exact_large.py, exact_large2.py, kernel_check_exact.py TV rows, tv_exact_point.py; "
             "CERTIFIED (python-flint 0.9.0 Arb, 200 bits): I3 and c_* in kernel_check_exact.py; "
             "DIAGNOSTIC (floats): kernel_constants.py (mpmath), tv_sup_scan.py, n1_check.py, n1_hybrid.py, kt5_sharp.py\n")
    fh.write("report_note: full report is the subagent's final message (harness blocks .md report files)\n")
    fh.write("inputs:\n")
    for p in inputs:
        if os.path.exists(p): fh.write(f"  - path: {p}\n    sha256: {sha(p)}\n")
    fh.write("outputs:\n")
    for p in outs:
        fh.write(f"  - path: {p}\n    sha256: {sha(p)}\n")
print(open(HERE + "/manifest.yaml").read()[:1500])
