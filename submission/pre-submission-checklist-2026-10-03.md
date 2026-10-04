# Pre-submission checklist: ECP, Poisson–binomial paper
<!-- SUMMARY: final pre-flight for the ECP submission via EJMS; all gates pass on source 1e950299…f642 / PDF 3c4bd14e…7d24; risk of no human audit accepted by Brett; MathSciNet pass done 2026-10-04 (nothing found); expert check and Brett's forecast open · status: blocked on Brett's decisions · updated: 2026-10-03 -->

Completed from `Project-Management/templates/paper-submission-checklist.md`.
Items that don't apply are marked N/A with the reason.

## Submission record

- Project: `papers/queue/erdos-problem-993`. This is the Poisson–binomial spin-off; the #993 line itself was dropped on 2026-10-03.
- Submission type: journal, first submission.
- Target: Electronic Communications in Probability (ECP), via EJMS, <https://www.e-publications.org/ims/submission/>. Brett's account `brettrey` already holds the ECP Author role (checked in the browser on 2026-10-03).
- Canonical source: `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.tex`, SHA-256 `1e950299b8c1838d391d84f6249d706b6e1cc9a83592c14bc096728be549f642`.
- Canonical PDF: `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.pdf` (gitignored build product), SHA-256 `3c4bd14e3f71d0e2a8cecc12dc543ac102a4f4072a43d7b92acff617bd0e7d24`, 13 pages.
  - Built 2026-10-03 at 21:44 EDT with `pdflatex`, `bibtex`, `pdflatex`, `pdflatex`.
  - **This is the file to upload.** Rebuilding changes the hash, because pdfTeX embeds the build time.
- Source package: `paper/poisson_binomial/ecp_source_bundle.zip`, SHA-256 `e1abd93ea80d5c83b04a55a61deee0b9a9e499b46025c877df7ab0883307bd16`. It is held for acceptance, since ECP asks only for a PDF at submission.
- Supplement: Zenodo, DOI 10.5281/zenodo.23129480 (published). It holds the same file as `paper/poisson_binomial/poisson_binomial_certificate_supplement.zip`: MD5 `5c670f4c405ed79eacb7372e938b9e41`, 567,017 bytes, SHA-256 `b1e6462d…309d`.
- Venue decision record: `submission/venue-decision-2026-07-16.md`.
- Portal-fields record: `submission/portal-fields-ECP-2026-10-03.md`.
- Paper assurance record: `submission/paper-assurance-2026-10-03.md`.
- Decision owner: Brett Reynolds. Assisting agent: Claude Code (Opus 5.5).
- Date checked: 2026-10-03, 21:30–22:00 EDT.

## 1. Submission decision

- [x] The venue decision record exists, with a journal-reader contract, alternatives and a risk test.
- [ ] **The venue record's own pre-submission conditions.** The first version of this checklist did not read them; a side agent caught the gap.
  - *Independent audit of the full proof:* open. This is covered by Brett's acceptance of no human audit. Since July, Theorem 1.1 has been fully kernel-checked in Lean.
  - *Authenticated MathSciNet pass and submission-date delta search:* **done 2026-10-04**, nothing found. The MathSciNet pass ran through Brett's U of T session at his instruction and covered keyword searches plus four forward-citation lists (`runs/pb-novelty-search-20261003/MATHSCINET.md`). The delta since 2026-07-01 is in the zbMATH and arXiv sweep (`REPORT.md`, 37 queries listed). Neither found anything that states or implies Theorem 1.1.
  - *Expert novelty / proof-interface check:* **open**. Brett must confirm whether his acceptance of no human audit covers it.
  - *Durable DOI archive:* done (Zenodo 10.5281/zenodo.23129480).
  - *Brett's own forecast:* the record asks for it before submission. **Pending.**
- [x] Desk and reviewer risks. The strongest open risk is that no human mathematician has audited the paper. **Brett accepted it, 2026-10-03: "Happy to submit."** The cross-family model checks are recorded in DECISIONS.md: ChatGPT Pro reconstructed all 275 Bernstein coefficients, and Codex (GPT-6.1-sol) redid the computations of Section 4. Theorem 1.1 is fully kernel-checked in Lean (PBScalar and PBDeduction, replayed locally, standard axioms only).
- [x] Current venue instructions. ECP journal page and Manuscript Submission page, pasted by Brett on 2026-10-03 and saved verbatim in `submission/ecp-author-instructions-2026-10-03.md`.
- [x] Not under consideration elsewhere, and no prior submission of this paper.
- [x] Public history to disclose if asked: the source has been in the public GitHub repository since commit `8435919`, and the supplement is on Zenodo. There is no preprint. ECP permits preprints.
- [x] Coauthors: N/A (single author).
- [x] Constraints. No collaborators, participants or restricted data. Brett can't post to math arXiv, so Zenodo only.

## 2. Objection and source gate

- [x] External reviews were triaged in DECISIONS.md: two ChatGPT reviews, a third ChatGPT review, Elicit, and two Codex reviews. Their raw text is in `runs/`. Each suggestion was adopted with credit or set aside with a reason. None is treated as evidence.
- [x] Sources reread before citing. The three references added today (Saumard–Wellner, Brändén, Liggett) are in `notes/passes/2026-10-03-review4-reverification.md`. Liggett is behind a CAPTCHA, so it is cited only for the definition, which it is the standard citation for.
- [x] Numbers. The `numbers-audit` pass is current. The printed margins of Section 4.1 come from `scripts/verify_pb_compact_monotone.py`.
- [x] Projectibility audit: N/A (mathematics paper).

## 3. Manuscript content gate

- [x] Title, abstract and keywords match the contribution. The abstract matches the portal record word for word (182 words), checked by script.
- [x] Abstract style: one paragraph, plain prose, ECP length.
- [x] Presentation passes are all current on the pass ledger after today's diff-scoped re-check (`notes/passes/2026-10-03-preflight-reverification.md`).
  - That re-check found one editorial scar, the sentence order left by demoting the window remark, and fixed it by reordering existing sentences.
  - `review-board` has never been run and is not required by the gate. Brett asked for no new cold reads.
- [x] No placeholders in the source or the PDF text (TODO, PENDING, XXX, ??, TBD).
- [x] Figures: none. Table 1 is cited in the text.
- [x] Declarations:
  - competing interests: none (Brett, 2026-10-03);
  - funding: none (Brett);
  - acknowledgements: no people to thank (Brett);
  - ethics: N/A;
  - data and code: Zenodo DOI, printed in the supplement box.
- [x] AI disclosure. A page-1 footnote, plus the acknowledgement after the references, which names GPT-5.6, ChatGPT, Elicit, Codex (GPT-6.1-sol), Aristotle and Claude Opus 5.5 with their specific contributions. Neither ECP page has an AI-specific rule.

## 4. Blind review

- N/A. ECP shows no double-anonymous requirement; the venue record checked the public instructions. Recheck if the portal asks.

## 5. Build and file gate

- [x] House-style scan: N/A. The paper uses the journal's `ejpecp` class and conventions, which Brett accepted on 2026-10-03.
- [x] `git diff --check` is clean.
- [x] Clean rebuild: 0 lines beginning `!`, 0 undefined references, 0 overfull, 0 underfull, no rerun request. The log's other "warning" hits are a hyperref option note from the class and a bibtex function count.
- [x] `pdftotext` check. Page 1 (title, author, abstract, keywords, MSC) and the supplement box (archive DOI) render correctly, and all 14 references are present.
- [x] `pdfinfo`: 13 pages. The title, author, subject and keywords metadata are set by the class.
- [x] PDF hyperlinks, read from the decompressed PDF: the DOI (twice), ORCID `0000-0003-0073-7195`, the email, the MSC page, and 14 MathSciNet links, one per reference. The DOI resolves (doi.org 302 to Zenodo).
- [x] Visual check: rendered pages 1, 3, 12 and 13 were inspected.
- [x] Source bundle. It was rebuilt after today's fix and compiles alone (pdflatex ×3): 13 pp, 0 errors, and text identical to the manuscript PDF.
- [x] Supplement. Zenodo API: status published, creator ORCID correct, affiliations Humber Polytechnic and University of Toronto, file checksum equal to the local build. The archive contains no manuscript file, so today's sentence reorder doesn't affect it.
- [x] Mechanical gates.
  - `check-terms.py --gate`: its one acronym flag, "ST", is the product $S\cdot T$ in a display in Section 4.2, a false positive.
  - `check-quotes.py --gate`: the paper contains no quotations.
  - The skill's `submission_gate_check.py` assumes a `main.tex` layout, which the project's no-`main.*` rule forbids, so its file-layout failures don't apply.

## 6. Preprint gate

- N/A. There is no preprint. The supplement deposit was handled in its own decision entries.

## 7. Journal-specific gate

- [x] The portal-fields record is complete. Every remaining open item is answered or left to the portal by Brett's decision. Brett's home postal code is entered only in the portal and is not recorded in this public repository.
- [x] Instructions checked on submission day.
- [x] Length: 13 pp, the stated maximum per the 2025 IMS editorial report. Brett accepted it.
- [x] Cover letter `submission/cover-letter-ECP-2026-10-03.md` has no placeholders. It gives the DOI, affiliations and ORCID.
- [x] The data availability statement points to a live DOI.
- [x] Reviewer suggestions: none. That field is left to the editor.
- [ ] Portal file roles and the portal-generated PDF preview are checked before the final submit. **Do in session.**
- [ ] Confirmation, manuscript ID and acknowledgement email are saved under `submission/`. **After submit.**

## 8. Assurance composition gate

- [x] `submission/paper-assurance-2026-10-03.md` pins the hashes above. `paper_assurance.py check --gate` must report CURRENT and gate PASS; the result is recorded in the record's commit.
- [x] Accepted residual risk: no independent human audit. It is named in `accepted_risks` and was accepted by Brett.

## 9. Final submit gate

- [ ] The uploaded PDF hash is confirmed equal to `3c4bd14e…7d24`.
- [ ] **Brett gives the final go.** The agent fills the form and stops at the final submit button.

## Stop conditions

**Active:** the venue record's expert novelty / proof-interface check is open, and Brett has not accepted it as a risk. The MathSciNet pass was done on 2026-10-04 and found nothing. Brett's forecast is pending. The stop condition is "Risks accepted explicitly by Brett only if a listed gate remains open." The no-human-audit risk is accepted.
