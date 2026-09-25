# Provenance and scope

Package date: 2026-09-23. Prepared by Codex for Brett Reynolds after approval of
the focused handoff brief. No Primus job was created or uploaded by this step.

Source repository revision: `b85d5dc01bdcaf3d6853147b28c350d74d3ede21`.
The working tree had unrelated changes; those were preserved and excluded.
The package is a curated snapshot, not a repository export.

## Unmodified files copied from the project

- `lean/RequestProject/`, `lean/lakefile.toml`, `lean/lake-manifest.json`, and
  `lean/lean-toolchain`: from
  `formalization/depth3_b1_zero_20260906/checkpoint-timeout-1/`.
- `lean/HistoricalAudit.lean`: from
  `formalization/depth3_b1_zero_20260906/LocalAudit.lean` (renamed only).
- `indpoly.py`: the project-root numerical helper.
- `scripts/audit_b1_zero_forbidden_core_20260905.py`,
  `scripts/audit_pendant_p2_boundary_20260904.py`, and
  `scripts/probe_blocked_profile_depth3_20260828.py`: corresponding project
  scripts, retained unchanged for the functions imported by the replay.
- `evidence/b1_positive_sign_obstructions_20260905.json` and
  `evidence/cross_reserve_witness_lifts_20260904.json`: corresponding original
  files under `results/`.
- `LICENSE`: project-root MIT licence.

The wrapper, tests, integrity checker, brief, and `lean/PrimusSpec.lean` are new
packaging files. The new specification does not change the historical
definitions or add a proved theorem.

## Historical verification, not a new claim of replay

The project's 6 September records report a fresh Lean build and exact-target
axiom audit of the b_1=0, low-density, and blocked-shadow graph theorems. Their
audited axioms were `propext`, `Classical.choice`, and `Quot.sound`. The source
was also reviewed for counting semantics and coverage. These records concern
the supplied established subcases, not W, R, or S in the new brief.

The frozen historical `DepthThreeSpec.lean` hash is
`aaf283ca986d1119706e68169cff73fab3c050f99b66ce68736b353e6a6da090`.
The historical `RequestProject/B1Zero/Target.lean` hash is
`aa7be0df446ceea33b4ca6f670d3f899bf2fd6ec8a4ef8dba1324820f75252fe`.
The full historical returned-source hash list is retained in the originating
repository; packaging compares the copied Lean sources with that record.

`CHECKS.md` reports what was actually checked during packaging. Do not turn a
hash comparison into a claim that a new full Lean build was performed.

## Exclusions

No confidential referee reports, correspondence, submission history, personal
account identifiers, credentials, API task metadata, private project-management
records, unrelated manuscripts, `.git` history, downloaded dependencies, Lean
build caches, Python environments, or third-party paper PDFs are included.
The original status journal and provider-generated narrative summaries are
excluded to avoid stale status or overbroad claims.

## Returning results

Return the exact revised/additional source files, a short result report, and
reproduction commands. Identify the exact statement proved or refuted.
If using Lean, report exact theorem ascriptions and `#print axioms` output;
do not silently alter supplied definitions, dependency pins, or hypotheses.
Any claimed result will require independent replay before adoption.
