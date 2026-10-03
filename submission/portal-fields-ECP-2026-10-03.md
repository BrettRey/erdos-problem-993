# Portal Fields Record: ECP, Poisson–binomial paper
<!-- SUMMARY: Answer sheet for the EJMS/ECP portal for the Poisson–binomial first-descent paper; filled 2026-10-03 from variance-scaled-turan-first-descent.tex and the venue record; section 11 lists what Brett must supply before login · status: open (section 11 not empty) · updated: 2026-10-03 -->

## Record

- [x] Project: `papers/queue/erdos-problem-993` (spin-off manuscript; the #993 line itself was dropped 2026-10-03)
- [x] Manuscript title (short form): Variance-scaled Turán inequality at first descent
- [x] Venue: Electronic Communications in Probability (ECP)
- [x] Portal platform: EJMS (IMS Electronic Journal Management System)
- [ ] Portal URL: <https://www.e-publications.org/ims/submission/> (from a search snippet of the imstat author page, 2026-10-03; confirm on login)
- [ ] Author instructions URL: imstat ECP author page. Automated access is blocked by an anti-bot challenge, which was not circumvented. Brett must open it in a normal browser on submission day.
- [ ] Date instructions checked: 2026-07-16 (venue record); **recheck on submission day**
- [x] Venue decision record path: `submission/venue-decision-2026-07-16.md`
- [x] Pre-submission checklist path: `runs/pb-submission-gate-20261003/` (referee read, revision verifier, cold reads) and DECISIONS 2026-10-03
- [x] Paper assurance record path: `paper/poisson_binomial/CERTIFICATE.md` (inside the supplement)
- [x] Canonical source file: `paper/poisson_binomial/variance-scaled-turan-first-descent.tex`
- [ ] Canonical PDF: `paper/poisson_binomial/variance-scaled-turan-first-descent.pdf` (gitignored; rebuild with `pdflatex`, `bibtex`, `pdflatex` ×2). The pre-DOI build of 2026-10-03 has SHA-256 `0b207438fc2767f1b91e3efc084966e29099d7d10d5380d3e55e7dbe543e535d`, 11 pages, and prints `ZENODO-DOI-PENDING`. Record the final hash after the DOI is inserted. Source bundle: `python3 scripts/build_pb_ecp_source_bundle.py`, which refuses to run while the placeholder remains.
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

Page count of the canonical PDF: 11 (`pdfinfo`, 2026-10-03).

## 2. Title, abstract, keywords

**Title**

> A variance-scaled Turán inequality at the first descent of a Poisson–binomial mass function

- Source: `variance-scaled-turan-first-descent.tex:5–6`

**Short title**

> A variance-scaled Turán inequality at the first descent

- Source: `variance-scaled-turan-first-descent.tex:3`

**Abstract.** Plain text with LaTeX math, as EJMS accepts TeX in the abstract box; confirm on login.

> Let $W$ be a finite sum of independent Bernoulli summands with probability mass function (pmf) $f=(f_k)$, where $f_k=\mathbb{P}(W=k)$, and variance $V=\operatorname{Var}(W)\geq1$. If $D$ is the first-descent index of this pmf (the least $k$ with $f_k<f_{k-1}$), we prove $V(1-f_{D-1}f_{D+1}/f_D^2)\geq1/4$. The bracket is the normalized slack in the log-concavity (Turán) inequality $f_D^2\geq f_{D-1}f_{D+1}$. Consequently $f_{D+r}/f_D\leq\exp(-r/(4V))$ for every support index $D+r$ with $r\geq1$. An explicit family of binomial laws shows that no universal constant larger than $1/3$ is possible. The proof uses cubic inequalities of Hillion and Johnson to bound the masses near the mode from below, turns these bounds into a lower bound on $V$, and closes with a maximal-mass bound of Bobkov, Marsiglietti, and Melbourne. The resulting one-variable inequality is proved in exact arithmetic by Bernstein expansions.

- Source: `variance-scaled-turan-first-descent.tex` `\ABSTRACT{...}`; the displayed formula is inlined here.
- [ ] Recheck word for word against `pdftotext -f 1 -l 1 variance-scaled-turan-first-descent.pdf` at the final build.

**Keywords**

> Poisson–binomial law; Bernoulli sum; probability mass function; log-concavity; Turán inequality; modal index; computer-assisted proof

- Source: `variance-scaled-turan-first-descent.tex:15–16`, semicolon-separated as in the class

**MSC2020**

> Primary 60E15; Secondary 60C05, 05A20

- Source: `variance-scaled-turan-first-descent.tex:18–19`

## 3. Authors

| # | Name | Email | Institution | Department | Postal code | Country | ORCID | Corresponding |
|---|---|---|---|---|---|---|---|---|
| 1 | Brett Reynolds | brett.reynolds@humber.ca | Humber College, Toronto | **open** | **open** | Canada | 0000-0003-2407-9448 | yes |

- Source: `variance-scaled-turan-first-descent.tex:10–13`
- Single author; no coauthor approvals needed.

## 4. Files and portal item types

| Local path | Portal item type | Reviewer sees it | Notes |
|---|---|---|---|
| `paper/poisson_binomial/variance-scaled-turan-first-descent.pdf` | manuscript PDF | yes | final build after the DOI insert |
| `paper/poisson_binomial/variance-scaled-turan-first-descent.tex` + `ejpecp.cls` + bbl inlined | source files (supplementary or at acceptance) | per portal | the ECP sample asks for the bibliography inside the document (`sample.tex` L439–442); inline `variance-scaled-turan-first-descent.bbl` in the source bundle |
| `paper/poisson_binomial/poisson_binomial_certificate_supplement.zip` | supplementary material | yes | SHA-256 `69dc523b…17af` (2026-10-03 build); same file as the Zenodo deposit |
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

**AI-use disclosure.** On page 1 (title footnote) and in the acknowledgements ("AI use and author responsibility"), naming GPT-5.6 via Codex, Claude, ChatGPT, Codex, Gemini, and Aristotle (Harmonic). No IMS/ECP-specific AI rule was found on 2026-07-16; **recheck on submission day** and move or extend the disclosure if a portal field asks for it.

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
| Supplement DOI | Zenodo upload of `poisson_binomial_certificate_supplement.zip` (reserve the DOI first, then replace `ZENODO-DOI-PENDING` in `variance-scaled-turan-first-descent.tex` and rebuild) | Brett (his account) | no |
| EJMS account | not yet created | Brett | no |
| Author instructions | live check of the imstat author page, including field limits, source-file rules, and any AI-use rule | Brett (browser) | no |
| Department and postal code | not in the manuscript | Brett | no |
| Competing interests and funding | composed defaults above need confirmation | Brett | no |
| Independent human audit | the paper has had no human mathematician's review; proceeding without one is Brett's call and should be recorded in DECISIONS | Brett | no |
| MathSciNet novelty search | optional; needs an institutional login | Brett | no |

- [ ] This table is empty, or every remaining row is explicitly accepted by Brett as answerable live in the portal.
