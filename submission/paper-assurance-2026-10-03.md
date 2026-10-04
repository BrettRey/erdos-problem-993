---
schema_version: 1
project: papers/queue/erdos-problem-993
intended_use: journal-submission
target: Electronic Communications in Probability (EJMS), research article
record_status: unresolved
canonical_source: paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.tex
canonical_pdf: paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.pdf
source_commit: df330141a619260aa53f3e3c06994f8ed4b2d66a
source_sha256: 0687e054618274f119133d23f9c430ab5ce1a8c21bc38cc53768a1c49d4ffc4c
pdf_sha256: aa6a5f4f59bd119239081e5f99b95c32d53363fd836b48d4ef59f2ff4927359e
checked_at: '2026-10-04T08:58:44-04:00'
dirty_paths: []
decision_owner: Brett Reynolds
decision_status: pending
accepted_risks: []
evidence_artifacts:
- path: submission/pre-submission-checklist-2026-10-03.md
  sha256: 66faa249eb058a74e660d510a597ea5456a977aaa0a50d1c9527e6f100e3952b
  supports: The exact source builds cleanly to the pinned 13-page PDF. Its abstract,
    keywords, links, metadata and declarations match the portal record and the ECP
    instructions of 3 October 2026, the supplement DOI is live with a checksum equal
    to the local archive, and all ten EJMS source checkboxes hold literally.
  relation: composition
  inspected_source_sha256: 0687e054618274f119133d23f9c430ab5ce1a8c21bc38cc53768a1c49d4ffc4c
- path: notes/passes/2026-10-04-ejms-compliance-reverification.md
  sha256: 9496ac7eb0f19cad4b1b66b3c8adc547913a040dd62d2cc09eae105bd1b14f9c
  supports: Every registered manuscript pass is current on the exact source. The EJMS
    compliance diff (journal mode, template metadata, Table 1 moved) changed no wording,
    and each pass was judged against it. This continues the 2026-10-03 pre-flight re-check
    of the earlier diffs.
  relation: composition
  inspected_source_sha256: 0687e054618274f119133d23f9c430ab5ce1a8c21bc38cc53768a1c49d4ffc4c
- path: paper/poisson_binomial/CERTIFICATE.md
  sha256: e421c7974138ec03b8cb7221769b13163c52daaa0d3480d7c26b14acb7ada352
  supports: The computations the proof needs and the Lean formalizations of Proposition
    3.1 and Theorem 1.1 are listed with their replay commands and records. The copy
    in the published Zenodo archive is byte-identical. Every numbered cross-reference
    resolves to the same item in the exact source.
  relation: convergence
  inspected_source_sha256: 0687e054618274f119133d23f9c430ab5ce1a8c21bc38cc53768a1c49d4ffc4c
- path: runs/pb-novelty-search-20261003/REPORT.md
  sha256: 8bff08f8f0abd358a9e152d54db44739fa1853cf80259c8554704077f828194e
  supports: In 37 listed zbMATH Open and arXiv queries, including a delta since 2026-07-01,
    no document states or implies Theorem 1.1 as stated in the exact source. This
    is a bounded negative and not a substitute for the MathSciNet pass.
  relation: convergence
  inspected_source_sha256: 0687e054618274f119133d23f9c430ab5ce1a8c21bc38cc53768a1c49d4ffc4c
- path: runs/pb-novelty-search-20261003/MATHSCINET.md
  sha256: 50268d990ad2e1c6c35124fd87846bd993e35cd8d7a6e1460e1ba27d8458b5c8
  supports: An authenticated MathSciNet pass found nothing that states or implies
    Theorem 1.1 as stated in the exact source. It covered keyword searches and the
    forward citations of Hillion–Johnson 2016, Johnson 2017, Pitman 1997 and Darroch
    1964. This meets the venue record's MathSciNet condition; it is a bounded negative.
  relation: convergence
  inspected_source_sha256: 0687e054618274f119133d23f9c430ab5ce1a8c21bc38cc53768a1c49d4ffc4c
review_artifacts: []
unresolved_items:
- Portal file roles and the EJMS-generated PDF have not yet been previewed. Check
  them after upload and before the final submit.
- Two risks are to be recorded as accepted_risks, with decision_status risk-accepted,
  when Brett gives the final go on this exact object. (1) No human mathematician has
  audited the paper; Brett accepted this on 2026-10-03 ("Happy to submit"), before
  the sentence reorder (which he approved on 2026-10-04) and the EJMS compliance
  changes (which he chose on 2026-10-04) existed. (2) No expert novelty / proof-interface check
  has been done; Brett accepted this on 2026-10-04 ("I have no way to get an expert
  check, so that's not in the cards").
invalidation_triggers:
- Canonical source or PDF changes, including a rebuild (pdfTeX embeds the build time).
- A load-bearing source, analysis, computation or Lean statement is corrected.
- The intended use, target, or submitted package changes.
- The Zenodo record or the ECP instructions change before submission.
---
# Paper Assurance Record
<!-- SUMMARY: Version-pinned assurance record for the ECP submission of the Poisson–binomial paper; evidence current on source 0687e054…ffc4c / PDF aa6a5f4f…359e (journal mode, EJMS checkboxes verified); venue-record conditions met or risk-accepted; awaiting portal preview and Brett's final go · status: unresolved · updated: 2026-10-04 -->

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

Commit `df330141a619260aa53f3e3c06994f8ed4b2d66a` has a clean tree for the
canonical source. The PDF is a gitignored build product, pinned by its
SHA-256. The source bundle (SHA-256 `926d432a…b2b9`) is held for acceptance
and is not part of this submission.

The object changed on 2026-10-04 from source `1e950299…f642` / PDF
`3c4bd14e…7d24`. The change put the paper in journal mode, restored the
template metadata and moved Table 1, to satisfy the EJMS source checkboxes.
No wording, result number or theorem statement changed: a word-multiset
comparison of the PDF text shows only journal furniture added. So the
certificate's cross-references and the two novelty searches carry over to the
exact source unchanged.

## Evidence ledger

| Check | Bounded conclusion | Artifact | Inspected version | Status |
| --- | --- | --- | --- | --- |
| Build, metadata, venue | Clean build of the exact source to the pinned PDF. The abstract and declarations match the portal record and the ECP instructions. The DOI is live and its checksum matches | `submission/pre-submission-checklist-2026-10-03.md` | source `0687e054…ffc4c` | supported |
| Manuscript passes | All registered passes are current on the exact source, and the last diffs were checked against each pass's dimension | `notes/passes/2026-10-04-ejms-compliance-reverification.md` (continuing `2026-10-03-preflight-reverification.md`) | source `0687e054…ffc4c` | supported |
| Mathematical certificate | The exact computations and the Lean proofs (PBScalar, PBDeduction; Lean 4.28.0, Mathlib v4.28.0, standard axioms only) are identified with replay records. The cross-references match the exact source | `paper/poisson_binomial/CERTIFICATE.md` | source `0687e054…ffc4c` (cross-references) | supported |
| Free novelty sweep | 37 zbMATH Open and arXiv queries (delta since 2026-07-01 included); nothing states or implies Theorem 1.1 | `runs/pb-novelty-search-20261003/REPORT.md` | source `0687e054…ffc4c` (theorem statement) | supported |
| MathSciNet pass | Keyword searches plus four forward-citation lists through Brett's U of T session; nothing states or implies Theorem 1.1 | `runs/pb-novelty-search-20261003/MATHSCINET.md` | source `0687e054…ffc4c` (theorem statement) | supported |

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

Two items are unresolved, listed in the front matter: the portal preview, and
recording Brett's two risk acceptances (no human audit, no expert check) on
this exact object.

Every venue-record condition is now met or accepted:
- the MathSciNet pass was met on 2026-10-04;
- Brett gave his forecast on 2026-10-04;
- he accepted the missing expert check as a risk on 2026-10-04. The first
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
