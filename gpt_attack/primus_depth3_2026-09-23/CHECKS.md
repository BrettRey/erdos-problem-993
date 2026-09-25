# Packaging checks — 23 September 2026

These checks validate the handoff, not the open mathematical target.

## Exact numerical checks run for this package

Environment: Python 3.14.7, NetworkX 3.6.1. The replay disables NumPy and uses
Python integers for every decisive arithmetic comparison.

| Check | Result |
|---|---|
| `python3 -B -m unittest -v test_packet` | PASS: all six tests |
| Direct set-count comparisons on the 48 nonisomorphic trees of orders 1–8 | PASS: polynomial and extendable profiles agree |
| Empty/disconnected forest regression cases | PASS |
| Relabelling, invalid-input rejection, and unimodality plateau handling | PASS |
| `python3 -B replay.py` | PASS: all three saved witnesses recomputed |

The recomputed witness margins are:

| Witness | Combined correction D | Full margin s_3^2-s_2 s_4 |
|---|---:|---:|
| Pair-only defect-one adverse example | -56,227,048 | 96,167,097,495 |
| Unary-only defect-one adverse example | -178,212,783 | 94,423,068,753 |
| Blocked-profile log-concavity failure | -1,336,360,568 | 462,192,427,400 |

All three are trees in the target window, and all three are in the already
covered low-density region. They refute tempting universal shortcuts, but
are not witnesses in the remaining high-density regime R. Their full
independence sequences are unimodal. The historical aggregate lift search
was not rerun.

## Lean checks run for this package

- All 52 copied historical Lean modules and the three pinned configuration
  files match the original 6 September SHA-256 record.
- The copied `RequestProject/` source tree is byte-for-byte identical to the
  historical returned source tree.
- `PrimusSpec.lean` was elaborated successfully with the existing pinned
  Lean 4.28.0 environment. Its printed definitions were checked against the
  brief, including integer-valued signed correction and the non-strict W target.
  It states open propositions; compilation is not a proof of them.
- `HistoricalAudit.lean` was rerun successfully in that environment. Exact
  ascriptions of the established three graph theorems passed. They and the
  supporting declarations printed only `propext`, `Classical.choice`, and
  `Quot.sound` as axioms.
- A source scan found no `sorry`, `admit`, added `axiom`, `native_decide`,
  `implemented_by`, `unsafe`, or `extern` tokens in the included Lean files.

The two Lean commands used the existing compiled dependencies in the original
pinned project. This is a new elaboration/axiom check, **not** a new clean build
of all 52 modules or of mathlib. The previous fresh build is historical evidence,
described in `PROVENANCE.md`. Instructions for a fresh rebuild are in `README.md`.

## Content and integrity checks

The package is assembled from an explicit selection of files, not from a
recursive repository archive. The selection excludes confidential reviews,
correspondence, account metadata, credentials, caches, and repository history.
A path/credential-pattern scan found no such data; the only matches for
"referee" were the documents' statements that referee material is excluded.
There are no symbolic links.

`MANIFEST.sha256` covers all payload files except itself. Verify it using
`python3 check_manifest.py`. This checks integrity, not mathematical correctness.
