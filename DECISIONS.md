# Decisions Log

Append-only record of project decisions. Agents: add an entry whenever a non-trivial decision is made during a session (structural changes, venue choices, theoretical commitments, scope changes, reviewer feedback acted on). Keep entries short.

Format: `## YYYY-MM-DD` then bullet points with **bold topic** and brief rationale.

---

## 2026-04-25

- **Submit to E-JC after review-board cleanup** Round-3 review shifted the manuscript to submission-ready. Chose the Electronic Journal of Combinatorics as next venue after the Experimental Mathematics desk rejection. Submitted single-PDF initial submission `60492-1`; no source files uploaded at this stage per E-JC guidelines. Final PDF removes the ORCID icon to avoid logo-like title-block material and builds without overfull boxes.

## 2026-04-23

- **Desk-rejected from Experimental Mathematics** Editor-in-Chief Alexander Kasprzyk returned a standard "not suitable for publication" notice on 2026-04-23, no reviews, no reason given. Submission ID 264515082, lodged 2026-03-29 (roughly four-week turnaround). Candidate next venues remain European Journal of Combinatorics and Discrete Mathematics, as flagged in the 2026-03-17/18 notes; venue choice deferred pending Brett's call. PORTFOLIO.md, CV, and publications.html updated to drop the "under review" framing.

## 2026-03-17/18

- **Aristotle (harmonic.fun) for Lean 4 formalization** Used Aristotle to formalize proved results in Lean 4. Four runs: P3 leaf-swap injection (8 sorries filled), J ≤ E subgraph monotonicity (4 sorries), star+star w_2 + binomial LC (8 sorries), STP2 closure (the open problem — unsolved, as expected). All auxiliary proofs verified, main conjecture remains open.
- **Consolidate into Formal/ not lean/** The existing `Formal/` directory (pinned Mathlib v4.28.0, lakefile.toml) is the canonical formalization home. Aristotle's `lean/` directory was a scratch workspace. New files (JleE.lean, Algebra.lean, STP2Closure.lean) ported into `Formal/` with correct namespaces.
- **isLC definition has ℕ subtraction bug** Aristotle identified that `isLC` using ℕ subtraction at k=0 creates the constraint f(0)·f(1) ≤ f(0)², which is more restrictive than standard mathematical LC. The classical result "LC preserved under convolution" (`lc_conv`) is FALSE under this definition. Must use k ≥ 1 guard or work over ℤ.
- **Aristotle can't solve open problems** Three hours on the STP2 closure theorem produced useful base cases (leaf×leaf, degree-0) but not the general proof. Tool is effective for formalizing known proofs, not discovering new ones.
- **Use official Aristotle CLI/SDK instead of the web app for routine runs** Installed `aristotlelib` via `uv tool install` and added a repo-local wrapper script so Lean-project submissions, listings, and result downloads can be driven from the terminal with repo defaults.
- **Merge Aristotle proof of `tree_has_pendant` into canonical `Formal/` tree** Submitted the repo-root Lean v4.28.0 project with a narrow prompt for `Formal/P3.lean`; Aristotle discharged the standard finite-tree leaf existence lemma using `SimpleGraph.IsTree.exists_vert_degree_one_of_nontrivial`, closing the last real `sorry` in `Formal/P3.lean`.
- **Direct Zenodo versioning for paper snapshots** The GitHub release `paper-v2-2026-03-18-doi` was published correctly, but Zenodo ingest failed again. The paper DOI refresh was completed manually on Zenodo instead, minting `10.5281/zenodo.19100781` under concept DOI `10.5281/zenodo.18745546`. Future paper-only DOI refreshes should use Zenodo's direct `New version` flow rather than GitHub retry releases.
- **Submit without waiting for private feedback** Emailed Ohr Kadrawi with the current manuscript, but decided not to hold journal submission on outreach replies. The practical next submission targets are `European Journal of Combinatorics` and `Discrete Mathematics`.
- **Shutdown posture documented** End-of-session state was normalized into the repo docs: current paper DOI, current GitHub release tag, outreach status, and the preferred Zenodo workflow are now recorded so the project can be safely parked.

## 2026-03-12

- **Forum update positioned as follow-up** Because the `erdosproblems.com` thread for `#993` already contained a January 7, 2026 comment announcing `n <= 29` verification, the new public comment was framed as a follow-up emphasizing the public repo, exact `n=28` LC / near-miss audit, and the current structural manuscript status rather than re-announcing the frontier.
- **Conservative wiki ask** Chose to propose an AI-contributions wiki entry in section `2(e)` ("numerical exploration") rather than a stronger `1(d)` classification, to avoid overclaiming while the manuscript is not yet publicly posted.
- **No YAML PR yet** Decided not to open a `data/problems.yaml` PR at this stage; the metadata edit is low-value compared to the wiki issue and can wait for maintainer feedback.
- **Endorsement outreach kept private** Kept the arXiv endorsement request off the public forum/wiki channels and prepared it only for direct one-to-one outreach.

## 2026-03-11

- **Prospect-aware search scoring** Added archive-plus-prospect scoring to `nm_optimizer.py` and `scripts/lc_breaker_optimizer.py` rather than trying to reconstruct unavailable AlphaEvolve code; the transferable idea is to reward lineages that empirically produce stronger children while retaining best distinct exact trees.
- **Ramsey paper citation placement** Added the Nagda--Raghavan--Thakurta Ramsey/AlphaEvolve paper to `paper/references.bib` and cited it in `paper/main_v2.tex` as broader AI-assisted extremal-combinatorics context, not as direct tree-unimodality prior art.
- **No `n=29` LC/NM burn** Estimated Modal cost for a full `n=29` LC + near-miss run was on the order of $1k before credits, which is outside budget; an initial dispatch was stopped immediately and no sustained large compute will be pursued.
- **Project shelved** With exhaustive unimodality already closed through `n=29` and no affordable high-value computation left, the project is being put to bed; only low-cost paper/admin cleanup would justify reopening.
