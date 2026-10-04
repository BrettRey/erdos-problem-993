# Portal Fields Record: ECP, Poisson–binomial paper
<!-- SUMMARY: Answer sheet for the EJMS/ECP portal for the Poisson–binomial first-descent paper; filled 2026-10-03 from variance-local-log-concavity-poisson-binomial.tex and the venue record; section 11 resolved 2026-10-04 (MathSciNet done; expert check risk-accepted; forecast given); upload PDF aa6a5f4f…359e · status: ready for portal · updated: 2026-10-04 -->

## Record

- [x] Project: `papers/queue/erdos-problem-993` (spin-off manuscript; the #993 line itself was dropped 2026-10-03)
- [x] Manuscript title (short form): Variance and local log-concavity of Bernoulli sums
- [x] Venue: Electronic Communications in Probability (ECP)
- [x] Portal platform: EJMS (IMS Electronic Journal Management System)
- [x] Portal URL: <https://www.e-publications.org/ims/submission/> (named on the ECP Manuscript Submission page; Brett logged in there on 2026-10-03)
- [x] Author instructions URL: <https://imstat.org/journals-and-publications/electronic-communications-in-probability/> and its Manuscript Submission page. Brett pasted both texts on 2026-10-03; they are saved verbatim in `submission/ecp-author-instructions-2026-10-03.md`.
- [x] Date instructions checked: 2026-10-03 (Brett, browser). Recheck only if submission slips by weeks.
- [x] Venue decision record path: `submission/venue-decision-2026-07-16.md`
- [x] Pre-submission checklist path: `runs/pb-submission-gate-20261003/` (referee read, revision verifier, cold reads) and DECISIONS 2026-10-03
- [x] Paper assurance record path: `paper/poisson_binomial/CERTIFICATE.md` (inside the supplement)
- [x] Canonical source file: `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.tex`
- [x] Canonical PDF: `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.pdf` (gitignored; built with `pdflatex`, `bibtex`, `pdflatex` ×2).
  - The final build is 2026-10-04, 08:54 EDT, from source SHA-256 `0687e054…ffc4c`. Its SHA-256 is `aa6a5f4f59bd119239081e5f99b95c32d53363fd836b48d4ef59f2ff4927359e`.
  - It has 13 pages. Brett accepted 13, the journal's stated maximum, on 2026-10-03.
  - It is in journal mode (`\documentclass[ECP]{ejpecp}`, with the template's placeholder volume, paper number and DOI), and Table 1 sits beside its first callout. Both changes make the EJMS source checkboxes literally true; Brett's call, 2026-10-04.
  - It prints the published supplement DOI.
  - **Upload this file without rebuilding**: pdfTeX embeds the build time, so a rebuild changes the hash.
- [x] Account the submission is made under: Brett's existing EJMS account `brettrey`, which holds the ECP Author role (Roles page, checked 2026-10-03)
- [x] Decision owner: Brett Reynolds
- [x] Assisting agent/model: Claude Code (Opus 5.5)

## 1. Routing

| Field | Value | Limit | Source |
|---|---|---|---|
| Journal | Electronic Communications in Probability | | venue record |
| Article type | Research article | ECP "normally 12 pages or less"; the 2025 IMS editorial report starts review only at ≤12, at most 13 | venue record, Fit Evidence |
| Section / category | none known | | |
| Review model | single-anonymized at the public-instruction level (no double-anonymous requirement found); recheck on login | | venue record |

Page count of the canonical PDF: 13 (`pdfinfo`, final build, 2026-10-03).

## 2. Title, abstract, keywords

**Title**

> Variance and local log-concavity of Poisson–binomial laws

- Source: `variance-local-log-concavity-poisson-binomial.tex:5–6`

**Short title**

> Variance and local log-concavity of Bernoulli sums

- Source: `variance-local-log-concavity-poisson-binomial.tex:3`

**Abstract.** Plain text with LaTeX math, as EJMS accepts TeX in the abstract box; confirm on login.

> Let $W$ be the number of successes in independent trials with unequal success probabilities, with probability mass function $f$ and variance $V\geq1$, and let $D$ be the first index at which $f$ decreases. The normalized Turán deficit $\delta_k=1-f_{k-1}f_{k+1}/f_k^2$, a measure of local log-concavity, is a bounded transform of the discrete curvature of $\log f$ at $k$; for a normal density of variance $V$, that curvature is $1/V$. Quantitative log-concavity is used to derive concentration and log-Sobolev inequalities for discrete laws, but for Bernoulli sums the known lower bounds on $\delta_k$ depend on the number of summands or on the success probabilities. We prove $V\delta_D\geq1/4$, and binomial laws show that no constant above $1/3$ is possible. With cubic inequalities of Hillion and Johnson, this gives $\delta_k\geq1/(4V+|k-D|)$ on the whole support. Ultra-log-concavity alone gives no variance-scaled lower bound on $\delta_D$. For the number of eigenvalues of a Haar unitary matrix in a half circle, our bound at $D$ is of order $1/\log N$, against order $1/N$ from ultra-log-concavity. The proof of the main theorem is computer-assisted and is verified in the Lean proof assistant.

- Source: `variance-local-log-concavity-poisson-binomial.tex` `\ABSTRACT{...}` (no displayed formula in the current abstract).
- [x] Rechecked word for word at the final pre-flight (2026-10-03). A script compared the `.tex` `\ABSTRACT` with this block and found them identical (182 words). The page-1 text of the final PDF matches.

**Keywords**

> Poisson–binomial law; log-concavity; Turán inequality; ultra-log-concavity; computer-assisted proof; formal verification

- Source: `variance-local-log-concavity-poisson-binomial.tex:15–17`, semicolon-separated as in the class

**MSC2020**

> Primary 60E15; Secondary 60C05, 05A20

- Source: `variance-local-log-concavity-poisson-binomial.tex:19–20`

## 3. Authors

| # | Name | Email | Institution | Department | Postal code | Country | ORCID | Corresponding |
|---|---|---|---|---|---|---|---|---|
| 1 | Brett Reynolds | brett.reynolds@humber.ca | Humber Polytechnic, Toronto; University of Toronto | Faculty of Liberal Arts and Sciences (Humber); Department of Linguistics (U of T) | supplied by Brett 2026-10-03 (home address); entered in the portal only, not recorded in this public repository | Canada | 0000-0003-0073-7195 | yes |

- Source: `variance-local-log-concavity-poisson-binomial.tex:10–13`
- Single author; no coauthor approvals needed.

## 4. Files and portal item types

| Local path | Portal item type | Reviewer sees it | Notes |
|---|---|---|---|
| `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.pdf` | manuscript PDF | yes | final build, SHA-256 `aa6a5f4f…359e`, 13 pp |
| `paper/poisson_binomial/ecp_source_bundle.zip` (the `.tex` with the bibliography pasted in, plus `ejpecp.cls` v1.12.3, 2026/03/04) | **not needed at submission**: the ECP Manuscript Submission page asks only for a PDF; hold it for acceptance | no | Rebuilt 2026-10-04 after the EJMS compliance changes by `scripts/build_pb_ecp_source_bundle.py`, SHA-256 `926d432a…b2b9`. Compiles alone (pdflatex ×3): 13 pp, 0 errors, text identical to the manuscript PDF. It is in journal mode, with the template's placeholder volume, paper number and DOI, which the managing editor sets |
| `paper/poisson_binomial/poisson_binomial_certificate_supplement.zip` | supplementary material | yes | SHA-256 `b1e6462d…309d` (2026-10-03 build with the Section 4.1 check program and the Lean proofs of Proposition 3.1 and Theorem 1.1); same file as the Zenodo deposit |
| `submission/cover-letter-ECP-2026-10-03.md` | cover letter / comments to editor | no | |

- [ ] The uploaded PDF hash matches the final build recorded above.
- [ ] Portal proof PDF previewed before final submit.

## 5. Declarations

**Competing interests**

> The author declares no competing interests.

- Source: composed here; confirmed by Brett 2026-10-03 ("no competing interests")

**Funding**

> No specific funding was received for this work.

- Source: Brett, 2026-10-03: "no funding". Confirmed. The manuscript carries no funding line, which ECP does not require.

**Acknowledgement of people**

- Brett, 2026-10-03: "nobody to thank". The acknowledgement credits only the AI systems and tools, as it does now.

**Data and code availability**

> The programs and exact certificate data supporting Proposition 3.1 are in the supplementary archive, deposited at https://doi.org/10.5281/zenodo.23129480. The development repository is https://github.com/BrettRey/erdos-problem-993.

- [x] Do not submit until the Zenodo record is live and the DOI resolves. (Live and resolving, verified 2026-10-03.)

**Ethics.** No human participants or data.

**AI-use disclosure.** On page 1 (title footnote) and in the acknowledgements ("AI use and author responsibility"), which come after the references. The acknowledgements credit specifically: GPT-5.6 via Codex (proof search for Theorem 1.1); an AI-generated ChatGPT review (proposed the framing, Corollary 1.3, Proposition 1.4 with proof, Proposition 1.5, Example 1.7); a second ChatGPT review (the Pitman bound (5.1); independent reconstruction of the Bernstein computations); Aristotle (the Lean proofs, and the monotonicity argument of Section 4.1); further AI-generated reviews (ChatGPT, Elicit, Codex with GPT-6.1-sol) that checked computations and Lean statements, found a boundary case, suggested references, and prompted corrections and the reorganization of Section 1; Claude Opus 5.5 via Claude Code (revision, proof checking, verification code, source checking); Claude, ChatGPT, Codex and Gemini in earlier sessions. Dates and model labels ("Latest" at Pro effort for the two ChatGPT reviews, 3 October 2026) are in the archive (CERTIFICATE.md, "Provenance of AI contributions"). No IMS/ECP-specific AI rule was found on 2026-07-16; **recheck on submission day** and move or extend the disclosure if a portal field asks for it.

## 6. Reviewers

None suggested. Leave the field to the editor unless the portal requires names. If it does, Brett picks names that have been verified as real, currently placeable people; nobody is invented here.

## 7. History and overlap

- Preprint: none on arXiv or any preprint server.
- **Public source:** the manuscript source has been visible in the public GitHub repository since commit `8435919` ("Add ECP Poisson-binomial paper"). If the portal asks about prior public versions, disclose this. ECP permits preprints.
- Prior submission of this paper: none.
- Related work by the author: *Mean bounds, structural reductions, and exhaustive verification for tree independence polynomial unimodality* (Zenodo 19100781; declined by E-JC 2026-09-13). It does not contain the theorem of this paper. Its only Poisson–binomial content is Darroch's mode theorem, cited (`paper/main_v2.tex` L320, 326, 655, 659, 797).
- Not under consideration elsewhere: confirmed.

## 8. Publishing options and cost

- Open access: ECP is free to publish and free to read, CC BY 4.0 (venue record, Fit Evidence).
- APC: none.

## 11. Open before login

| Field | What is missing | Who decides | Resolved |
|---|---|---|---|
| Supplement DOI | **Done 2026-10-03/04:** 10.5281/zenodo.23129480 published (record 23129480); DOI resolves; public file MD5 5c670f4c405ed79eacb7372e938b9e41, 567,017 bytes, matching the local build; creator ORCID 0000-0003-0073-7195 with affiliations Humber Polytechnic and University of Toronto | Brett | **yes** |
| EJMS account | **Done:** existing account `brettrey` holds the ECP Author role (2026-10-03) | Brett | **yes** |
| Author instructions | **Done 2026-10-03:** PDF submission via EJMS (add the "ECP author" role if already registered in EJMS); ECP template; no AI-specific rule on either page; CC BY for articles; no fees | Brett | **yes** |
| Postal code | **Done:** Brett supplied his home postal code on 2026-10-03. It is entered in the portal only, because this repository is public | Brett | **yes** |
| Competing interests | **Done:** "The author declares no competing interests." Confirmed by Brett on 2026-10-03 | Brett | **yes** |
| Independent human audit | no human mathematician has reviewed the paper; a cross-family model check (ChatGPT Pro, 2026-10-03) found no error and reconstructed all 275 coefficients of the Bernstein certificate (now the supplement's second check of 3 < H ≤ 16; the paper's Section 4.1 uses the monotonicity route, kernel-checked in Lean but not reviewed by ChatGPT Pro); proceeding without a human audit is Brett's call and should be recorded in DECISIONS. **Risk accepted by Brett, 2026-10-03 ("Happy to submit")**; recorded in DECISIONS and in `accepted_risks` of `submission/paper-assurance-2026-10-03.md` | Brett | **yes (risk accepted)** |
| MathSciNet novelty search | **Done 2026-10-04.** Run through Brett's U of T session at his instruction ("Drive it"). It covered keyword searches plus the forward citations of Hillion–Johnson 2016, Johnson 2017, Pitman 1997 and Darroch 1964. Nothing states or implies Theorem 1.1 (`runs/pb-novelty-search-20261003/MATHSCINET.md`). The delta since July is in the arXiv sweep in `REPORT.md`. An earlier row on 2026-10-03 had wrongly called this resolved without Brett's decision | Brett | **yes** |
| Expert novelty / proof-interface check | **Risk accepted by Brett, 2026-10-04:** "I have no way to get an expert check, so that's not in the cards." | Brett | **yes (risk accepted)** |
| Brett's forecast | **Given 2026-10-04**; the numbers are kept in the private prediction ledger. To be logged with the event at submission | Brett | **yes** |

- [x] Every row is resolved. Two risks are accepted by Brett: no human audit, and no expert check.
