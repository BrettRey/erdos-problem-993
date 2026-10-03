# Portal Fields Record: ECP, Poisson–binomial paper
<!-- SUMMARY: Answer sheet for the EJMS/ECP portal for the Poisson–binomial first-descent paper; filled 2026-10-03 from variance-local-log-concavity-poisson-binomial.tex and the venue record; section 11 lists what Brett must supply before login · status: open (section 11 not empty) · updated: 2026-10-03 -->

## Record

- [x] Project: `papers/queue/erdos-problem-993` (spin-off manuscript; the #993 line itself was dropped 2026-10-03)
- [x] Manuscript title (short form): Variance and local log-concavity of Bernoulli sums
- [x] Venue: Electronic Communications in Probability (ECP)
- [x] Portal platform: EJMS (IMS Electronic Journal Management System)
- [ ] Portal URL: <https://www.e-publications.org/ims/submission/> (from a search snippet of the imstat author page, 2026-10-03; confirm on login)
- [ ] Author instructions URL: imstat ECP author page. Automated access is blocked by an anti-bot challenge, which was not circumvented. Brett must open it in a normal browser on submission day.
- [ ] Date instructions checked: 2026-07-16 (venue record); **recheck on submission day**
- [x] Venue decision record path: `submission/venue-decision-2026-07-16.md`
- [x] Pre-submission checklist path: `runs/pb-submission-gate-20261003/` (referee read, revision verifier, cold reads) and DECISIONS 2026-10-03
- [x] Paper assurance record path: `paper/poisson_binomial/CERTIFICATE.md` (inside the supplement)
- [x] Canonical source file: `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.tex`
- [ ] Canonical PDF: `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.pdf` (gitignored; rebuild with `pdflatex`, `bibtex`, `pdflatex` ×2). The pre-DOI build of 2026-10-03 (after the metaphor, paragraph-opening and redundancy passes) has SHA-256 `aceeb103b81235c59a74e29625babd417f96d63510cd279eb80d7f1ff1e53983`, 12 pages, and prints `ZENODO-DOI-PENDING`. Record the final hash after the DOI is inserted. Source bundle: `python3 scripts/build_pb_ecp_source_bundle.py`, which refuses to run while the placeholder remains.
- [ ] Account the submission is made under: Brett's EJMS account (not yet created)
- [x] Decision owner: Brett Reynolds
- [x] Assisting agent/model: Claude Code (Opus 5.5)

## 1. Routing

| Field | Value | Limit | Source |
|---|---|---|---|
| Journal | Electronic Communications in Probability | | venue record |
| Article type | Research article | ECP "normally 12 pages or less"; the 2025 IMS editorial report starts review only at ≤12, at most 13 | venue record, Fit Evidence |
| Section / category | none known | | |
| Review model | single-anonymized at the public-instruction level (no double-anonymous requirement found); recheck on login | | venue record |

Page count of the canonical PDF: 12, the last page references only (`pdfinfo`, 2026-10-03).

## 2. Title, abstract, keywords

**Title**

> Variance and local log-concavity of Poisson–binomial laws

- Source: `variance-local-log-concavity-poisson-binomial.tex:5–6`

**Short title**

> Variance and local log-concavity of Bernoulli sums

- Source: `variance-local-log-concavity-poisson-binomial.tex:3`

**Abstract.** Plain text with LaTeX math, as EJMS accepts TeX in the abstract box; confirm on login.

> Let $W$ be a finite sum of independent Bernoulli variables with probability mass function $f$ and variance $V\geq1$, and let $D$ be the first index at which $f$ decreases. The normalized Turán deficit $\delta_k=1-f_{k-1}f_{k+1}/f_k^2$ is a bounded transform of the discrete curvature of $\log f$ at $k$, and for a normal density of variance $V$ that curvature is $1/V$. We prove $V\delta_D\geq1/4$, and an explicit family of binomial laws shows that no constant above $1/3$ is possible, so the best constant lies between $1/4$ and $1/3$. Combined with cubic inequalities of Hillion and Johnson, this gives $\delta_k\geq1/(4V+|k-D|)$ at every $k$ in the support, and at the rightmost mode $c$ we have $1/(4V+1)\leq\delta_c<2/V$. Ultra-log-concavity alone gives no variance-scaled lower bound on $\delta_D$. For the number of eigenvalues of an $N\times N$ Haar unitary matrix in a half circle, the bound at $D$ is of order $1/\log N$, against order $1/N$ from ultra-log-concavity. The proof combines the cubic inequalities with a maximal-mass bound of Bobkov, Marsiglietti, and Melbourne to reduce the main inequality to a one-variable inequality, which is proved by Bernstein expansions, with a computer-assisted step in exact rational arithmetic.

- Source: `variance-local-log-concavity-poisson-binomial.tex` `\ABSTRACT{...}` (no displayed formula in the current abstract).
- [ ] Recheck word for word against `pdftotext -f 1 -l 1 variance-local-log-concavity-poisson-binomial.pdf` at the final build.

**Keywords**

> Poisson–binomial law; Bernoulli sum; log-concavity; Turán inequality; ultra-log-concavity; modal index; random unitary matrix; computer-assisted proof

- Source: `variance-local-log-concavity-poisson-binomial.tex:15–17`, semicolon-separated as in the class

**MSC2020**

> Primary 60E15; Secondary 60C05, 05A20

- Source: `variance-local-log-concavity-poisson-binomial.tex:19–20`

## 3. Authors

| # | Name | Email | Institution | Department | Postal code | Country | ORCID | Corresponding |
|---|---|---|---|---|---|---|---|---|
| 1 | Brett Reynolds | brett.reynolds@humber.ca | Humber College, Toronto | **open** | **open** | Canada | 0000-0003-2407-9448 | yes |

- Source: `variance-local-log-concavity-poisson-binomial.tex:10–13`
- Single author; no coauthor approvals needed.

## 4. Files and portal item types

| Local path | Portal item type | Reviewer sees it | Notes |
|---|---|---|---|
| `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.pdf` | manuscript PDF | yes | final build after the DOI insert |
| `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.tex` + `ejpecp.cls` + bbl inlined | source files (supplementary or at acceptance) | per portal | the ECP sample asks for the bibliography inside the document (`sample.tex` L439–442); inline `variance-local-log-concavity-poisson-binomial.bbl` in the source bundle |
| `paper/poisson_binomial/poisson_binomial_certificate_supplement.zip` | supplementary material | yes | SHA-256 `097652b0…ee0e` (2026-10-03 build with the CUE threshold script); same file as the Zenodo deposit |
| `submission/cover-letter-ECP-2026-10-03.md` | cover letter / comments to editor | no | |

- [ ] The uploaded PDF hash matches the final build recorded above.
- [ ] Portal proof PDF previewed before final submit.

## 5. Declarations

**Competing interests**

> The author declares no competing interests.

- Source: composed here; **Brett to confirm**

**Funding**

> No specific funding was received for this work.

- Source: composed here; **Brett to confirm**

**Data and code availability**

> The programs and exact certificate data supporting Proposition 3.1 are in the supplementary archive, deposited at [Zenodo DOI]. The development repository is https://github.com/BrettRey/erdos-problem-993.

- [ ] Do not submit until the Zenodo record is live and the DOI resolves.

**Ethics.** No human participants or data.

**AI-use disclosure.** On page 1 (title footnote) and in the acknowledgements ("AI use and author responsibility"). The acknowledgements credit specifically: GPT-5.6 via Codex (proof search for Theorem 1.1); two ChatGPT sessions (model labelled "Latest", 3 October 2026, Pro effort): a review (proposed the framing, Corollary 1.3, Proposition 1.4 with proof, Proposition 1.5, Example 1.6) and a referee report (the Pitman (20) bound and Table 1, independent reconstruction of the computations); an Elicit referee report (clarifications); Claude Opus 5.5 via Claude Code (revision, verification code, source checking); Claude, ChatGPT, Codex and Gemini in earlier sessions; Aristotle (Lean proofs). No IMS/ECP-specific AI rule was found on 2026-07-16; **recheck on submission day** and move or extend the disclosure if a portal field asks for it.

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
| Supplement DOI | Zenodo upload of `poisson_binomial_certificate_supplement.zip` (reserve the DOI first, then replace `ZENODO-DOI-PENDING` in `variance-local-log-concavity-poisson-binomial.tex` and rebuild) | Brett (his account) | no |
| EJMS account | not yet created | Brett | no |
| Author instructions | live check of the imstat author page, including field limits, source-file rules, and any AI-use rule | Brett (browser) | no |
| Department and postal code | not in the manuscript | Brett | no |
| Competing interests and funding | composed defaults above need confirmation | Brett | no |
| Independent human audit | no human mathematician has reviewed the paper; a cross-family model check (ChatGPT Pro, 2026-10-03) found no error and reconstructed all 275 coefficients; proceeding without a human audit is Brett's call and should be recorded in DECISIONS | Brett | no |
| MathSciNet novelty search | optional; needs an institutional login | Brett | no |

- [ ] This table is empty, or every remaining row is explicitly accepted by Brett as answerable live in the portal.
