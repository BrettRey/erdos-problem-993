# Zenodo preprint record: field values
<!-- SUMMARY: values for the Zenodo preprint record of the ECP-submitted Poisson–binomial paper, each with its source; published by Brett 2026-10-04 as 10.5281/zenodo.23137930, verified by API · status: published · updated: 2026-10-04 -->

Brett asked on 2026-10-04 for the preprint PDF to be built and the Zenodo
record drafted. Zenodo is the route because he cannot post to math arXiv. The
draft stops before **Publish**, as the supplement record did.

## File

- `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial-preprint.pdf`
  - Built by `scripts/build_pb_preprint.py` from the submitted source (commit 66990f9), with the class's `preprint` option on. The canonical `.tex` is unchanged.
  - 13 pages, 293,168 bytes. SHA-256 `425c6bc45a1fac4e87f52f3e367f7c5d36cd863caee306f0cfb5b00c8a5142bd`.
  - Its text equals the submitted PDF's text minus the journal furniture (header, placeholder DOI, ISSN, footers, "Submitted to ECP" line), checked by word-multiset comparison.
  - Gitignored as a build product; rebuild with the script.

## Metadata

| Field | Value | Source |
|---|---|---|
| Resource type | Publication / Preprint | Brett: preprint |
| Title | Variance and local log-concavity of Poisson–binomial laws | `.tex` `\TITLE` |
| Publication date | 2026-10-04 | today; the day the paper was submitted |
| Creator | Reynolds, Brett; ORCID 0000-0003-0073-7195; affiliations Humber Polytechnic, University of Toronto | as on the supplement record (verified by API 2026-10-03) |
| Description | the text below | the `.tex` abstract in Unicode, plus a status line |
| License | Creative Commons Attribution 4.0 International | the licence ECP publishes under; Brett can change it before publishing |
| Keywords | Poisson–binomial law; log-concavity; Turán inequality; ultra-log-concavity; computer-assisted proof; formal verification | `.tex` `\KEYWORDS`, as on the supplement record |
| Version | v1 | |
| Language | English | |
| Related work | Is supplemented by 10.5281/zenodo.23129480 (DOI, Software) | the supplement record |
| DOI | reserve a new Zenodo DOI | |
| Funding | none | Brett, 2026-10-03 |

### Description text

> Let W be the number of successes in independent trials with unequal success probabilities, with probability mass function f and variance V ≥ 1, and let D be the first index at which f decreases. The normalized Turán deficit δ_k = 1 − f_{k−1}f_{k+1}/f_k², a measure of local log-concavity, is a bounded transform of the discrete curvature of log f at k; for a normal density of variance V, that curvature is 1/V. Quantitative log-concavity is used to derive concentration and log-Sobolev inequalities for discrete laws, but for Bernoulli sums the known lower bounds on δ_k depend on the number of summands or on the success probabilities. We prove V·δ_D ≥ 1/4, and binomial laws show that no constant above 1/3 is possible. With cubic inequalities of Hillion and Johnson, this gives δ_k ≥ 1/(4V + |k − D|) on the whole support. Ultra-log-concavity alone gives no variance-scaled lower bound on δ_D. For the number of eigenvalues of a Haar unitary matrix in a half circle, our bound at D is of order 1/log N, against order 1/N from ultra-log-concavity. The proof of the main theorem is computer-assisted and is verified in the Lean proof assistant.
>
> Preprint, submitted to Electronic Communications in Probability (October 2026). MSC2020: 60E15 (primary); 60C05, 05A20 (secondary). Code, exact certificate data and the Lean formalizations: https://doi.org/10.5281/zenodo.23129480.

## Draft state (2026-10-04)

- **Zenodo draft 23137930**, reserved DOI **10.5281/zenodo.23137930**. Saved and **not published**.
- **File:** the uploaded PDF shows md5 `01c88afb9cab41c5c9a220cb74c0996b` and 293.17 KB on Zenodo, equal to the local build.
- **Fields:** entered as in the table above. Two details:
  - The keyword is written "Poisson-binomial law" with a hyphen, to match the supplement record.
  - The author was found by ORCID. Zenodo's ORCID entry had "Humber College", which was replaced by Humber Polytechnic and University of Toronto (ROR).
- **Preview** (`/records/23137930?preview=1`) shows: Preprint, Open, v1; the author with two affiliations; the description; the file; and "Is supplemented by: Software 10.5281/zenodo.23129480 (DOI)".
- **Licence and language:** CC BY 4.0 (Zenodo's default, kept) and English.
- **Not done; Brett decides:**
  - press **Publish**, after which the files are fixed;
  - add the reverse link on the supplement record ("Is supplement to" 10.5281/zenodo.23137930);
  - replace the July draft on brettreynolds.ca.

## Published (2026-10-04)

Brett published the record. It was verified through the Zenodo API on
2026-10-04:
- status published; DOI 10.5281/zenodo.23137930, which resolves (doi.org 302 to Zenodo);
- Preprint, `cc-by-4.0`, v1, eng;
- creator ORCID 0000-0003-0073-7195, affiliations Humber Polytechnic and University of Toronto;
- the six keywords; `isSupplementedBy` 10.5281/zenodo.23129480;
- file of 293,168 bytes, md5 `01c88afb9cab41c5c9a220cb74c0996b`, equal to the local build.

## Versions (2026-10-04)

From the Zenodo API:
- **Preprint:** version DOI 10.5281/zenodo.23137930; concept ("all versions") DOI **10.5281/zenodo.23137929**.
- **Supplement:** version DOI 10.5281/zenodo.23129480; concept DOI 10.5281/zenodo.23129479.

Rules:
- A revised preprint (after referee changes) goes in as **New version** on record 23137930, never as a fresh upload, which would split links and citations. Build it with `scripts/build_pb_preprint.py`.
- Website and CV links use the concept DOI 10.5281/zenodo.23137929, so they always resolve to the newest version.
- The paper keeps citing the supplement's version DOI (…23129480), because the proof rests on that exact certificate.
