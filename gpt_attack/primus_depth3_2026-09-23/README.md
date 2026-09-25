# Erdős #993: a bounded depth-three proof task

Prepared for Brett Reynolds on 23 September 2026. This is a focused research
handoff, not a claimed solution of Erdős Problem #993.

## Start here

1. Read `PROMPT.md`; it is also the text to paste into the Primus project prompt.
2. Read `MATHEMATICS.md` for the exact target and available results.
3. Run the checks below, then read `FAILED_ROUTES.md`.
4. Use `lean/` only when its formal interfaces or verification are useful.

Suggested project title: **Tree independence sequences: close the residual
depth-three window**. Select the Math mode if available. No service-specific
upload schema or file limit is assumed. If ZIP extraction is not supported,
extract this archive locally and attach its component files as the interface
permits; start with the prompt and mathematical brief.

## Quick checks

Python 3.11 or later and NetworkX 3.6.1 are required for the mathematical
replays. The integrity check uses the Python standard library only.

```sh
python3 check_manifest.py
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -B -m unittest -v test_packet
.venv/bin/python -B replay.py
```

`replay.py` rechecks three fixed graph witnesses; it does not launch a search.
It recomputes the independence polynomial and the extendable/blocked counts,
then checks the stored exact values. Its output distinguishes a failed
auxiliary inequality, the depth-three target, and unimodality. An additional
graph can be checked using `replay.py --graph6 'GRAPH6_STRING'`.

The unit tests include brute-force set counting on small fixed regression
cases. These test the packaged implementation, not the universal target.

## Contents

| Path | Purpose |
|---|---|
| `PROMPT.md` | Paste-ready assignment and required result grading |
| `MATHEMATICS.md` | Definitions, target, remaining regime, and theorem interfaces |
| `FAILED_ROUTES.md` | Known false shortcuts and carefully qualified optional leads |
| `replay.py`, `test_packet.py` | Exact witness replay and regression checks |
| `indpoly.py`, `scripts/` | Unmodified relevant numerical helpers from the project |
| `evidence/` | Two original JSON certificates; only three selected witnesses are replayed |
| `lean/` | Established Lean source, pinned dependencies, and a new specification only |
| `PROVENANCE.md`, `CHECKS.md` | Origin, historical verification, and packaging checks |
| `MANIFEST.sha256`, `check_manifest.py` | File integrity verification |
| `LICENSE` | Original project software licence |

The files in `evidence/` preserve historical search context. Their aggregate
search counts are not newly replayed by the three-witness check. Several
historical fields contain approximate display ratios; acceptance decisions in
`replay.py` use integers only.

## Lean, if needed

The included project pins Lean/mathlib 4.28.0 and the dependency commits in
`lean/lake-manifest.json`. Dependencies and build caches are deliberately absent.
From `lean/`, with the pinned Lean toolchain installed:

```sh
lake exe cache get
lake build
lake env lean HistoricalAudit.lean
lake env lean PrimusSpec.lean
```

These commands may download dependencies and take time. `PrimusSpec.lean`
defines the open propositions; it does not prove them. The historical audit
checks the already established subcases, not the open target. See `CHECKS.md`
for precisely what was run while making this ZIP.

Do not upload `.venv/`, `.lake/`, caches, or a repository history. The supplied
archive contains no confidential referee material or account credentials.
