---
schema_version: 1
project: papers/queue/erdos-problem-993
intended_use: journal-submission
target: Electronic Communications in Probability (EJMS), research article
record_status: unresolved
canonical_source: paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.tex
canonical_pdf: paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.pdf
source_commit: 1b66c60c694e293418dc2faa14f605fdf1d435b3
source_sha256: 1e950299b8c1838d391d84f6249d706b6e1cc9a83592c14bc096728be549f642
pdf_sha256: 3c4bd14e3f71d0e2a8cecc12dc543ac102a4f4072a43d7b92acff617bd0e7d24
checked_at: '2026-10-03T21:51:57-04:00'
dirty_paths: []
decision_owner: Brett Reynolds
decision_status: pending
accepted_risks: []
evidence_artifacts:
  - path: submission/pre-submission-checklist-2026-10-03.md
    sha256: 15667feca845dd58fbb79b990b2dc29fadfdbcc066b9f34afb153b5ebff56c68
    supports: The exact source builds cleanly to the pinned 13-page PDF. Its abstract, keywords, links, metadata and declarations match the portal record and the ECP instructions of 3 October 2026, and the supplement DOI is live with a checksum equal to the local archive.
    relation: composition
    inspected_source_sha256: 1e950299b8c1838d391d84f6249d706b6e1cc9a83592c14bc096728be549f642
  - path: notes/passes/2026-10-03-preflight-reverification.md
    sha256: 11baf7ddc5c78061e8ffbae1bc6f782314dd645285d197544c03553dce67b5d4
    supports: Every registered manuscript pass is current on the exact source. The last two commits and the final sentence reorder were checked against each pass's dimension, and one editorial scar was found and fixed.
    relation: composition
    inspected_source_sha256: 1e950299b8c1838d391d84f6249d706b6e1cc9a83592c14bc096728be549f642
  - path: paper/poisson_binomial/CERTIFICATE.md
    sha256: e421c7974138ec03b8cb7221769b13163c52daaa0d3480d7c26b14acb7ada352
    supports: The computations the proof needs and the Lean formalizations of Proposition 3.1 and Theorem 1.1 are listed with their replay commands and records. The copy in the published Zenodo archive is byte-identical. Every numbered cross-reference resolves to the same item in the exact source.
    relation: convergence
    inspected_source_sha256: 1e950299b8c1838d391d84f6249d706b6e1cc9a83592c14bc096728be549f642
  - path: runs/pb-novelty-search-20261003/REPORT.md
    sha256: 8bff08f8f0abd358a9e152d54db44739fa1853cf80259c8554704077f828194e
    supports: In 37 listed zbMATH Open and arXiv queries, including a delta since 2026-07-01, no document states or implies Theorem 1.1 as stated in the exact source. This is a bounded negative and not a substitute for the MathSciNet pass.
    relation: convergence
    inspected_source_sha256: 1e950299b8c1838d391d84f6249d706b6e1cc9a83592c14bc096728be549f642
review_artifacts: []
unresolved_items:
  - The venue record requires an authenticated MathSciNet pass, which has not been done. A free zbMATH and arXiv sweep found nothing. Brett decides whether to run the pass or accept the risk.
  - The venue record requires an expert novelty / proof-interface check, which has not been done. Brett must confirm whether his acceptance of no human audit covers it.
  - The venue record asks for Brett's own forecast before submission. It is not yet given.
  - Portal file roles and the EJMS-generated PDF have not yet been previewed. Check them after upload and before the final submit.
  - No human mathematician has audited the paper. Brett accepted this risk on 2026-10-03 ("Happy to submit"), but before the final sentence reorder existed. It becomes an accepted_risk with decision_status risk-accepted when he gives the final go on this exact object.
invalidation_triggers:
  - Canonical source or PDF changes, including a rebuild (pdfTeX embeds the build time).
  - A load-bearing source, analysis, computation or Lean statement is corrected.
  - The intended use, target, or submitted package changes.
  - The Zenodo record or the ECP instructions change before submission.
---
# Paper Assurance Record
<!-- SUMMARY: Version-pinned assurance record for the ECP submission of the Poisson–binomial paper; evidence current on source 1e950299…f642 / PDF 3c4bd14e…7d24; venue-record novelty conditions open; awaiting Brett's decisions, portal preview and final go · status: unresolved · updated: 2026-10-03 -->

## Declared use

This record covers the first submission of the exact PDF pinned above to
*Electronic Communications in Probability* through EJMS, as a research article
by a single author. The supplement is cited by its published Zenodo DOI. This
record does not itself authorize the portal submit action. Brett's explicit go
in the session does.

## Bounded readiness claim

For the exact objects pinned above, the recorded checks support three claims:
the build and package are ready, they comply with ECP's public instructions,
and every registered manuscript pass is current. They also show that the
archived computations and Lean proofs the paper relies on are identified and
their replay is recorded. They do not establish that the paper is correct
beyond what the Lean kernel and the exact computations check, that a human
expert would agree, or that ECP will accept it.

## Object

Commit `1b66c60c694e293418dc2faa14f605fdf1d435b3` has a clean tree for the
canonical source. The PDF is a gitignored build product, pinned by its
SHA-256. The source bundle (SHA-256 `e1abd93e…07bd16`) is held for acceptance
and is not part of this submission.

## Evidence ledger

| Check | Bounded conclusion | Artifact | Inspected version | Status |
| --- | --- | --- | --- | --- |
| Build, metadata, venue | Clean build of the exact source to the pinned PDF. The abstract and declarations match the portal record and the ECP instructions. The DOI is live and its checksum matches | `submission/pre-submission-checklist-2026-10-03.md` | source `1e950299…f642` | supported |
| Manuscript passes | All registered passes are current on the exact source, and the last diffs were checked against each pass's dimension | `notes/passes/2026-10-03-preflight-reverification.md` | source `1e950299…f642` | supported |
| Mathematical certificate | The exact computations and the Lean proofs (PBScalar, PBDeduction; Lean 4.28.0, Mathlib v4.28.0, standard axioms only) are identified with replay records. The cross-references match the exact source | `paper/poisson_binomial/CERTIFICATE.md` | source `1e950299…f642` (cross-references) | supported |
| Free novelty sweep | 37 zbMATH Open and arXiv queries (delta since 2026-07-01 included); nothing states or implies Theorem 1.1. Bounded negative; the MathSciNet pass is still open | `runs/pb-novelty-search-20261003/REPORT.md` | source `1e950299…f642` (theorem statement) | supported, MathSciNet open |

The local replay records named in CERTIFICATE.md are
`formalization/pb_scalar_inequality_aristotle_result/LOCAL_REPLAY.md`
(SHA-256 `4200aa58…94c0`) and
`formalization/pb_deduction_aristotle_result/LOCAL_REPLAY.md`
(`52eabf36…d6bf`). They inspected the Lean projects, not the manuscript
source, so they support this record through CERTIFICATE.md and are not listed
as separate evidence.

## Interface audit

| Inputs | Relation | Justification | Result |
| --- | --- | --- | --- |
| Build/package checklist + pass re-check | composition | Both inspected the same source hash. The checklist certifies the object the passes read, so their endpoints align | supported |
| Manuscript checks + mathematical certificate | convergence | Prose and package checks probe presentation defects. Exact computation and the Lean kernel probe mathematical defects. Neither stands in for the other | supported |

## Review lineage

No model review is used as assurance evidence. ChatGPT (three reviews),
Elicit, and Codex (GPT-6.1-sol) reviews were triaged in DECISIONS.md as
prompts for checks, and their raw text is kept under `runs/`. Their
agreement contributes nothing to this record.

## Unresolved and defeated conditions

Five items are unresolved, listed in the front matter:
- the MathSciNet pass;
- the expert check;
- Brett's forecast;
- the portal preview;
- recording Brett's risk acceptance on this exact object.

The first three are conditions in the venue decision record. The first
version of this record left them out, because the pre-flight had not read
the record's unchecked boxes. Nothing is defeated: the free novelty sweep
turned up no conflicting result.

## Invalidation triggers

Any rebuild of the PDF invalidates the record, because the hash changes even
when the text doesn't. So does any change to the source, the Zenodo record
or the ECP instructions.

## Decision

Brett Reynolds is the human decision owner. The decision is pending for this
exact object. After the portal preview, and when Brett gives the final go, the
record moves to `record_status: supported` and `decision_status:
risk-accepted`, with the no-human-audit risk in `accepted_risks`.
`paper_assurance.py check --gate` must then report CURRENT and gate PASS
before the submit button is pressed.
