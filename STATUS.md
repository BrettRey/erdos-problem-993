# Erdos Problem #993 -- Independent Set Sequence Unimodality for Trees

**Problem:** Is the independent set sequence of every tree unimodal?

**Source:** Alavi, Malde, Schwenk, and Erdos (1987)
**Reference:** [erdosproblems.com #993](https://erdosproblems.com/993)

## Source of truth

The current manuscript is `paper/main_v2.tex` (XeLaTeX + biber). Numeric snapshots live in `results/*.json` where available. The main proof-status references are `notes/one_private_status.md` and `notes/conjecture_A_analysis.md`; subdivision identity details live in `subdivision_correct.py` and `verify_subdivision_formula.py`.

## Current state (2026-07-11)

### Session notes (broad prefix-GSB variance push, numbered D1--D27 round)

- **Lane.** This is the variance/moment attack (D6/D13), separate from the
  2026-07-10 signed-reserve Poisson-binomial lane (D4). The numbered search
  through D27 was run against one strict resolution gate; the full index is
  `notes/solve_today_registry_2026-07-11.md`. No resolution: no counterexample
  and no full proof. The manuscript was not changed.
- **Central reduction (D13).** Prefix GSB, `Var(e_r) <= 2 mu_r + 2 eta_r` for
  `r <= L-2` with `mu_r=(r+1)i_(r+1)/i_r`, implies full unimodality: it makes
  `d_r=mu_r-(r+1)` nonincreasing through the uncertified prefix, and the known
  decreasing tail from `L` absorbs the rest, with plateaus automatic. In
  coefficients this is `(r+2)i_r i_(r+2) <= (r+1)i_(r+1)^2 + i_r i_(r+1)`. It is
  strictly weaker than ordered log-concavity and survives the 26-vertex
  deep-tail LC failure.
- **Base cases proved (r <= 3).** D6 proves `Var_r(e) <= E_r(e)`
  unconditionally for `r=1` and `r=2` (exact discriminant argument), so any
  minimal failure has `r>=3`. D12 proves the rank-three prefix-GSB theorem
  `alpha>=7 ==> 5 i_3 i_5 <= 4 i_4^2 + i_3 i_4` (equivalently `mu_4<=mu_3+1`)
  with an exact symbolic certificate. Packet
  `notes/rank3_prefix_gsb_theorem_packet_2026-07-11.md`, replay
  `scratch_rank3_gsb_certificate_20260711.py` (independently audited; local
  replay prints `symbolic_certificates: passed` with positive finite-base
  minima through order 14). A minimal failure therefore has `r>=4`, `k>=4`.
- **Rank-four method limit (D12).** The rank-three proof cannot be generalized
  from the same recorded radius-two statistics alone. Two order-13 trees (`L??????_B_H__y`,
  `L??????o@_DA@{`) share identical radius-two moments but have `i_6=170` and
  `168`. This blocks extension of that certificate from those statistics; it
  does not rule out other rank-four mechanisms. Replay
  `scratch_rank4_radius2_obstruction_20260711.py`.
- **Spectral route blocked (D14).** The down--up chain reduces prefix GSB to a
  gap bound plus an aggregate energy inequality, but the separated energy
  estimate is false at an exact, audited `n=210` witness (join two rooted
  spherical trees with branching `(2,3,2,1,2,1,1)`, `alpha=124`, last prefix
  rank `r=81`). That tree is still unimodal (true GSB ratio approximately
  0.214). Frozen packet `notes/d14_ed_obstruction_packet_2026-07-11.md`,
  certificate `scratch_d14_ed_obstruction_certificate_20260711.py` passed.
- **Earlier families, mostly negative.** D1 shifted-ratio domination is exactly
  false at rooted `n=29`. D2 minimal-counterexample corridor is
  theorem-strength blocked (subsumes open vertex-deletion mode stability;
  adjacent summand modes are not universal). D3 connector engineering
  converged with no counterexample (best order 171, no later ascent). D7
  pair-switching Hall/min-cut models each fail at small order and are blocked
  pending a new block-tree cut inequality. The unedged `Var<=2 mu` is false at
  Galvin's `T_(6,6,1)`, and unrestricted GSB is false at Galvin's `T_(14,8,1)`,
  so a global monotone-`d_r` proof is blocked and the prefix window is
  essential. D10 and D16 are blocked.
- **D19--D24 localization.** Decorated paths compress back to one
  product-minus-corner bridge (D19); direct leaves suppress the correcting
  branch and nested stars binomially smooth hard bumps (D20); block drift is
  algebraically equivalent to GSB (D21); bare nonedge signs and componentwise
  collision allocation fail (D22); stable or marginal periodic phase grammars
  cannot maintain standardized phase separation (D23); and separated or
  termwise-positive Bencs deficit extractions fail (D24).
- **D25--D26 local-allocation limits.** The sign `q_nonedge<=0` fails exactly on
  `T_(40,20,1)` at `n=1641,r=20`. Parent-oriented and symmetric inverse-degree
  fork allocations fail on an exact 9,418-vertex prefix composition at
  `r=3279`, while global ordered LC and GSB remain positive. The subsequent
  bounded `AF_lambda` cycle found one exact all-parameter obstruction: a center
  joined to five copies of `G(60,18)`, with `n=11106,alpha=5701,r=3799`, has a
  negative gap already at `lambda=0` and a negative slope, so every
  `lambda in [0,1]` fails. Bare support exchange fails at `n=20,r=7`; its
  one-fallback correction has extensive finite support but no proof. Full
  marked networks saturated through `n=11`, but their total-capacity cut is
  exactly GSB, so this is diagnostic rather than a reduction.
- **Remaining candidate (D27).** Fixed-orbit and fixed-common-set proofs of
  `q_far<=0` are false, while the aggregate prefix sign remains compatible with
  finite evidence through every tree of order 18. The path-Rayleigh and Bencs
  identities do not control the required fixed-rank slice, so this is still
  theorem-strength rather than a demonstrated simpler route.
- **Computational falsification (evidence only).** The exact bivariate
  extension-moment DP finds zero prefix failures of `Var_r(e)<=E_r(e)` across
  all 3,490,529 trees through order 21, including all 2,144,505 at order 21.
  This is falsification force only; it is not a proof and does not reduce a
  minimal counterexample into the computed range.
- **Posture.** These are private research results. No manuscript change, no
  release, no Erdős Problems post, no solution claim. Do not open a broad D28
  round. The bounded augmented-fork cycle is complete and refuted. Retain
  aggregate `q_far<=0` only as a separate parked proposition unless a genuinely
  global rank-window compensation mechanism appears. Radius-three moments
  remain a possible later experiment, not a demonstrated necessity; D14 and
  the full marked-flow route stay frozen.
- **Shipped.** Commit `b5cf7ed` ("Add 2026-07-11 obstruction search") pushed to
  `origin/master` with all 78 session artifacts. No build ran (no changes touch
  `paper/main_v2.tex`).
- **Follow-up package.** These session notes include the D25 all-parameter
  obstruction packet, exact replay, compact JSON certificate, and the D19--D27
  status corrections. The exact replay and full test suite pass; the paper is
  still untouched.

## Current state (2026-07-10)

### Session notes (GPT-5.6 Ultra signed-reserve breakthrough)

- **Universal finite PB theorem proved.** For every finite Poisson-binomial
  law with variance `V>=1` and a supported first strict descent `D`,
  `V * Delta_eff(D) >= 1/4`; consequently the raw post-descent reserve also
  satisfies `V * reserve >= 1/4`. Shifting `X-Y` to an ordinary PB law closes
  the full nonterminal finite signed bridge. Terminal raw reserve is immediate;
  its effective quotient is undefined.
- **Proof mechanism.** Hillion--Johnson's cubic inequalities propagate local
  Turán curvature into explicit two-sided mass windows. A pairwise variance
  bound and a sharp max-atom variance bound reduce the theorem to one scalar
  inequality. The original exact Sturm replay closes `3<H<=16`; a second,
  Lean-friendly certificate now replaces all 13 finite cells by 275/275
  strictly positive Bernstein coefficients. An analytic
  Bonferroni--Bernstein argument closes `H>=16`.
- **Exact replay.** `scripts/verify_universal_pb_effective_drop.py` certifies
  the endpoint recurrence, mass-window algebra, asymmetric degree-10 cell,
  twelve compact Sturm cells, and the full Bernstein tail. Its separate exact
  11,320-vector scan found zero failures. Certificate:
  `results/universal_pb_effective_drop_certificate_2026-07-10.json`.
- **Lean-friendly finite replay.**
  `scripts/verify_universal_pb_finite_bernstein.py` exactly reconstructs all
  13 finite-cell numerators in the full Bernstein basis, verifies all 275
  coefficients are positive, and emits Lean 4.28 identities that compile.
  Certificate:
  `results/universal_pb_finite_bernstein_certificate_2026-07-10.json`.
- **Aristotle core formalization complete.** Project
  `8d59a353-5c7b-4071-b837-9ab7bf561be3`, task
  `eb147738-3627-4d2b-96a5-6922b28b92e6`, filled all nine proof holes in the
  minimal Lean 4.28 packet at `formalization/pb_effective_drop_aristotle/`.
  Local replay builds with no `sorry` or new axioms. This verifies the
  endpoint-aware curvature propagation, endpoint exclusion, crossing ratio,
  and raw-from-effective layer conditional on the normalized
  Hillion--Johnson recurrence; it is not yet an end-to-end formalization of
  the universal PB theorem.
- **Independent special-case proofs.** The full Skellam limit and finite laws
  with one and two reflected Bernoulli factors were proved and independently
  audited before the universal argument was found. Their replay harnesses
  passed 2,691 high-precision Skellam rows, 47,850 exact one-factor rows, and
  97,488 exact two-factor rows. A redundant three-factor direct proof was also
  found, but the universal theorem supersedes it as a dependency.
- **Product term closed, issue #5 still open.** For
  `A=(1+x)^sQ`, `V_A=s/4+V_Q`, and `s>=4`, every post-descent ratio is at most
  `1-1/(s+4V_Q)`. The remaining hub-bouquet problem is perturbation by
  `B=xR`: prove separate pre-descent and post-descent increment budgets,
  especially when broom arms grow with `s`. This does not close Case B, PNP,
  or Erdős #993.
- **Audit corrections preserved.** The old conditional side-bound-only target
  is exactly false, and strict-first-descent continuity fails at plateaus.
  Exact regression tests now cover both failures; the full suite passes
  (`55` tests).
- **Manuscript posture.** `paper/main_v2.tex` was not changed. The universal PB
  theorem is a new private result and would require a deliberate standalone or
  revision decision, a Hillion--Johnson citation, and the exact certificate
  material before manuscript migration.
- **Shipped and externally tracked.** Commits `6b078d2` and `011c65b` are on
  `origin/master`. Issue #5 was updated but deliberately left open for the
  `A+xR` perturbation problem; issue #7 now tracks end-to-end Lean
  formalization. No release, manuscript update, Erdős Problems post, or
  solution announcement was made.
- **Clean stopping point.** The next bounded mathematical target is the
  fixed-arm `A+xR` perturbation theorem, with separate pre-descent and
  post-descent increment budgets. Growing-broom asymptotics and full Lean
  formalization should not be taken first. Stale Aristotle/Lean build children
  were terminated; no project, Aristotle, harness, or agent process remains.
- **Verification:** `python3 -m unittest test_all.py -v` passed; all five new
  theorem/certificate harnesses passed; the emitted finite-cell Lean file and
  completed Aristotle core build under Lean 4.28; `python3 -m py_compile` on
  the new scripts passed; `git diff --check` passed.

## Current state (2026-07-04)

### Session notes (2026-07-04, morning, signed-reserve audit hardening)
- **Mode switch:** issue #5 signed-reserve work is now in audit-hardening mode, not forward proof-push mode. Do not close issue #5 or claim signed reserve, hub-bouquet reserve, or Erdos #993 from the current notes.
- **Claim-status ledger added:** `notes/literature/audit_hardening_claims_ledger_2026-07-04.md` is the current reliability surface. It separates theorem-level claims, finite algebra, computational evidence, superseded targets, regression coverage, and permitted next work.
- **Prompt 4 patched:** `notes/literature/corrected_side_bound_audit_2026-07-04.md` now records the code/formula audit outcome and fixes the overclaim findings: missing `V=(m+n)/4 >= 1` scope, finite-search caveats, terminal-descent caveat, max-identity-error caveat, and float-vs-exact certification boundaries.
- **Numeric guardrails added:** `first_descent` now uses a relative-drop tolerance to avoid float-resolved plateau ties; signed analysis tracks no-descent, terminal-descent, nonpositive-mass, and non-finite diagnostic cases explicitly.
- **Regression coverage expanded:** `test_all.py` covers the corrected X-side `beta_X` boundary, exact fair-binomial constants (`5/8`, `3/4` at total count 4), plateau tolerance, terminal-descent counting, half-heavy+dust near misses, and an exact rational half-heavy+dust row below a `0.8` fallback constant.
- **Shipped:** commit `a2c1abd` (`Harden signed reserve audit surface`) pushed to `origin/master`. Issue #5 updated at `https://github.com/BrettRey/erdos-problem-993/issues/5#issuecomment-4882264719`.
- **Verification:** `python3 -m unittest test_all.py -v` passed (`53` tests); `python3 -m py_compile scripts/analyze_signed_conditionals.py scripts/probe_signed_pb_reserve.py scripts/signed_ratio_drop_breaker.py test_all.py` passed; hardened analyzer rerun on `results/signed_ratio_drop_breaker_extended_2026-07-04.json` processed `68` rows with `max_identity_error=6.661e-16`, `terminal=0`, and `nonfinite=0`.
- **Next action:** exact-certify or explicitly quarantine the large float side-bound breaker, then externally audit the one-sided effective-drop route before any theorem-level manuscript migration.

## Current state (2026-07-01)

### Session notes (2026-07-01, STP2 formal hygiene and public targets)
- **Lean LC/STP2 boundary fixed.** `Formal/STP2Closure.lean` now guards both `isLC` and `isSTP2` by `1 <= k`, removing the artificial `k = 0` condition caused by truncated subtraction on `Nat`.
- **Boundary regression added.** The file records `leafI_self_conv_stp2`, checking that `leafI * leafI = [1, 2, 1]` satisfies guarded STP2 rather than failing at the old artificial boundary.
- **Abstract STP2 closure quarantined.** The contiguous-support toy pair `toyI = [1,4,1]`, `toyE = [1,1,1]` is formalized: it satisfies the abstract shape hypotheses, including `noGaps`, but its self-convolution violates guarded STP2 at `k = 3`.
- **Missing invariant clarified.** The old broad closure shells are renamed as `stp2_conv_closure_tree_realizable_conjecture` and `stp2_multi_child_closure_tree_realizable_conjecture`, guarded by the deliberately empty placeholder `treeDPPair`. The next STP2 target is a genuine tree-DP realizability invariant, not another coefficient-shape hypothesis.
- **Public coordination updated.** `README.md` now lists the July 2026 public targets: tree-DP realizability, fixed-r certificate emission, forest/product valley search, and the remaining `n = 29` LC/near-miss audit.
- **Verification:** `lake env lean Formal/STP2Closure.lean`, `lake env lean gpt_attack/axiom_fixed_r_certificate/problem.lean`, `lake build`, and `git diff --check` passed. `Formal/STP2Closure.lean` now contains no `sorry`, `admit`, `axiom`, or `unsafe`.

## Current state (2026-06-12)

### Session notes (2026-06-12, afternoon, fixed-`r` Lean certificate bridge)
- **Fable packet replayed with local verification.** Processed `/Users/brettreynolds/Downloads/files(4).zip` through `files(7).zip`; retained useful patches only after `lake env lean gpt_attack/axiom_fixed_r_certificate/problem.lean` passed locally.
- **Gibbs bridge closed.** `gibbsProb` and `gibbsMean` are marked `noncomputable`; the concrete derivative identity `deriv_gibbsMean_eq_variance_div`, positivity lemma `gibbsZ_pos`, and `mean_shift_bound_for_independence_polynomial` all check.
- **Abstract Route-2 bridge packaged.** Added checked composition and certificate wrappers through `fixed_r_certificate_composition`, `Route2Certificate`, `fixed_r_certificate_composition_split`, `Route2SplitCertificate`, `Route2SplitCertificateFor`, `route2_family_from_finite_and_tail`, `Route2FamilyCertificate`, and `route2_of_family_certificate`.
- **Tex-gap audit outcome.** The small Lean gaps identified in `fixed_r_certificate_target.tex` are closed at the abstract level: the `B` versus `F^-` split, hub-on mode preservation, static comparison, widened budget, and finite-prefix/tail threshold split are all represented in Lean.
- **Literature update.** Downloaded and read Li 2026 plus Li--Yang--Zhang--Zhang 2025 under `notes/literature/`. Useful takeaway: partial 2-Schur positivity plus tail monotonicity may be a research lead, but global 2-Schur positivity is not a new invariant between log-concavity and unimodality.
- **Remaining work.** The abstract bridge is now closed. Remaining Lean work is the larger spider-polynomial model (`P_r`, coefficient extraction, `I(T_{a,r}) = F + G`) or adapting certificate scripts to emit `Route2FamilyCertificate` / `Route2SplitCertificateFor` instances. Broad #993 proof search remains parked.
- **Verification:** `lake env lean gpt_attack/axiom_fixed_r_certificate/problem.lean` passed; `git diff --check` passed; no `sorry`, `admit`, `axiom`, or `unsafe` occurs in the checked Lean file.

## Current state (2026-06-05)

### Session notes (2026-06-05, morning, LEAP/fixed-`r` formalization)
- **Reviewed arXiv:2606.03303v2 (LEAP).** The useful impact is workflow-level: bounded Lean-agent decomposition with verified parent sketches, not a new mathematical ingredient for Erdős #993.
- **Reopen decision:** do not reopen broad #993 proof search. Reopen only the fixed-`r` certificate bridge lane if work remains narrowly aimed at formalization or refutation of exact bridge lemmas.
- **Roughdraft checkpoint:** created and reviewed `notes/leap_reopen_assessment_2026-06-05.md`; Roughdraft returned no CriticMarkup comments.
- **Lean-facing progress:** updated `gpt_attack/axiom_fixed_r_certificate/problem.lean` to replace the opaque mean-shift Lipschitz certificate with proved finite-support variance lemmas, a derivative-to-Lipschitz bridge, and `mean_shift_bound_from_finite_gibbs_distribution`.
- **Remaining explicit subgoal:** prove the concrete Gibbs derivative identity `deriv mu t = finiteVariance N (p t) / t` for the polynomial-weight distribution.
- **Verification:** `lake env lean gpt_attack/axiom_fixed_r_certificate/problem.lean` passed; `git diff --check` passed. No broad `lake build` was run.

## Current state (2026-05-27)

### Session notes (2026-05-27, AxiomProver scan)
- **Scanned recent AxiomProver/Axiom Math paper pattern.** Recorded the scan in `notes/axiomprover_scan_2026-05-27.md`. The takeaway is workflow-level: the transferable pattern is narrow theorem/certificate targets with Lean verification, not broad informal proof discovery.
- **Best next proof-search target:** fixed-`r` Route-2 certificate machinery, especially the adjacent-to-global `F` mode margin, fugacity perturbation, and mean perturbation bridge lemmas.
- **New proof-agent packet:** added `gpt_attack/axiom_fixed_r_certificate/` with a README, natural-language target, task prompt, and Lean-facing arithmetic bridge lemmas. `problem.lean` checks with no `sorry`s; the next target is the finite-support Gibbs mean-shift lemma and full certificate criterion.
- **Parking decision:** mothball the main project again. Reopen only for E-JC reviews, an AxiomProver/Lean-agent run against `gpt_attack/axiom_fixed_r_certificate/`, or a tightly bounded proof session aimed at the finite-support Gibbs mean-shift lemma.
- **Lean check note:** a `lake build` attempt reached the expected `Formal/STP2Closure.lean` `sorry` warnings and then hung during `Formal/P3.lean` replay; the build was terminated. Do not treat it as a successful verification run.

## Current state (2026-05-24)

### Session notes (2026-05-24, afternoon, AlphaProof Nexus scan)
- **Reviewed `/Users/brettreynolds/Downloads/2605.22763v1.md`.** The DeepMind AlphaProof Nexus paper reports 9 solved Erdős problems, but not #993; no new independence-polynomial or tree-unimodality theorem was found to import.
- **Project impact:** no manuscript edit is indicated. The useful lesson is methodological: formal proof search should be aimed at narrow, frozen Lean/certificate subgoals, especially the existing `Formal/STP2Closure.lean` sorries or fixed-`r` Route-2 certificate inequalities.
- **Graph-theory side result:** the paper's unrelated spanning-tree/leaves/local-independence result was checked as nearby context but does not bear directly on the mode-mean or Conjecture A bottlenecks.
- **Working tree:** unchanged apart from shutdown tracking files. Pre-existing untracked artifacts remain: `erdos_993_gpt_attack_2026-04-24.zip`, `gpt_attack/main_v2.pdf`, and `results/route2_spider_lane_j0_1_a200_r80_float.json`.

## Current state (2026-05-21)

### Session notes (2026-05-21, fixed-`r` Route-2 workspace)
- **Fixed-`r` proof workspace shipped.** Commit `6d2c29c` (`Add fixed-r route2 proof workspace`) pushed to `origin/master`, adding `gpt_attack/` scripts, fixed-`r`/Route-2 notes, and replayable JSON certificate snapshots.
- **Mathematical posture:** the clean working note is `notes/fixed_r_proof_note_clean_2026-05-21.md`. It states a fixed-`r` computable-threshold criterion for structured spider lanes `S(2^a,r)`, not a manuscript-ready unconditional theorem for all fixed `r` and not a proof for all `d_leaf<=1` trees.
- **Submission posture:** do not submit the fixed-`r` material now. The E-JC submission remains the public manuscript track; the fixed-`r` work is a private technical reserve unless the fragile proof points are audited and the certificate machinery is made manuscript-grade.
- **Verification before shipping:** `python3 test_all.py` passed; `paper/main_v2.tex` rebuilt with `xelatex + biber + xelatex + xelatex` with only existing underfull-box warnings; `python3 -m py_compile gpt_attack/*.py` passed.
- **Untracked artifacts intentionally left out:** `erdos_993_gpt_attack_2026-04-24.zip`, `gpt_attack/main_v2.pdf`, and `results/route2_spider_lane_j0_1_a200_r80_float.json`.
- **Next action:** park the project unless E-JC reviews arrive or a proof-audit session specifically targets the global `F` margin lemma, the shifted-positive termination criterion, and the `r=2,3` boundary cases.

## Current state (2026-04-25)

### Session notes (2026-04-25, E-JC submission)
- **Submitted to the Electronic Journal of Combinatorics.** Submission number `60492-1`; initial file uploaded was the single PDF `paper/main_v2.pdf`. Acknowledgement received same day; OJS submission ID `15526`. Author dashboard: https://www.combinatorics.org/ojs/index.php/eljc/authorDashboard/submission/15526. E-JC requests no status enquiries before six months elapse.
- **Submission build:** final PDF built from `paper/main_v2.tex` after reviewer-driven revisions; no overfull boxes, no unresolved references/citations. ORCID icon removed from the title block to avoid logo-like material in the E-JC submission PDF.
- **Submission metadata:** OJS abstract saved in `paper/abstract_submission.txt` and `paper/abstract_submission.md`; comments for editor noted public GitHub/Zenodo artifacts, AI disclosure, manual checking, and that the paper is not under consideration elsewhere.

## Current state (2026-04-23)

### Session notes (2026-04-23, midday)
- **Desk-rejected from Experimental Mathematics.** Editor-in-Chief Alexander Kasprzyk emailed a standard "not suitable" notice on 2026-04-23, no reviews, no specific reason. Roughly four-week turnaround from the 2026-03-29 submission (ID 264515082).
- **Tracking surfaces updated:** `PORTFOLIO.md` row flipped to "Next venue TBD"; `personal/CV/main.tex` and `personal/personal-website/publications.html` both drop the "Under review at Experimental Mathematics" framing and present the paper as a Zenodo preprint again; `cv.pdf` rebuilt and synced into the website copy.
- **Repo manuscript unchanged** — the paper itself hasn't been touched this session; `paper/main_v2.tex` is still the live draft matching Zenodo version DOI `10.5281/zenodo.19100781`.
- **Next action:** Brett's call on the next venue. Candidates from the 2026-03-18 notes: European Journal of Combinatorics, Discrete Mathematics. Worth a short look at scope fit and typical desk-rejection rates before resubmitting.

## Current state (2026-03-18)

### Session notes (2026-03-18, shutdown prep)
- **Manuscript state:** `paper/main_v2.tex` is the live draft and incorporates the final minor-revision polish pass from the review-board synthesis. The paper build is clean.
- **Paper snapshot release:** published GitHub release `paper-v2-2026-03-18-doi` from commit `487096954fef98f7ff93a352d241f3f8d62ef0e5`.
- **Zenodo paper record:** the GitHub -> Zenodo ingest failed again, so the paper-only DOI refresh was completed manually on Zenodo instead:
  - concept DOI `10.5281/zenodo.18745546`
  - version DOI `10.5281/zenodo.19100781`
  - record `https://zenodo.org/records/19100781`
- **Zenodo workflow note:** if another paper-only DOI refresh is needed, use Zenodo's direct `New version` flow rather than retrying the GitHub release bridge.
- **Outreach:** emailed Ohr Kadrawi (`orka@ariel.ac.il`) with the current paper on `2026-03-18`.
- **Submission posture:** do not wait on outreach replies before submitting. Best-fit next venues are `European Journal of Combinatorics` and `Discrete Mathematics`.
- **Formalization state:** `Formal/P3.lean`, `Formal/JleE.lean`, and `Formal/Algebra.lean` are fully proved; `Formal/STP2Closure.lean` still isolates the two genuine open closure theorems.
- **Shutdown posture:** no local background jobs intentionally left running; the repo is now at a documented stopping point.

### Session notes (2026-03-17/18, evening–morning)
- **Aristotle formalization tool** (aristotle.harmonic.fun): first access. Lean 4 proof assistant that fills `sorry`s.
- **Four Aristotle runs:**
  1. P3 (leaf-swap injection): 8 sorries → 0. Built and verified.
  2. J ≤ E (subgraph monotonicity): 4 sorries → 0. Built and verified.
  3. Algebra (star+star w_2 + binomial LC): 8 sorries → 0. Built and verified.
  4. STP2 closure (THE OPEN PROBLEM): ~3 hours. Proved base cases (leaf×leaf, degree-0). Main theorem still `sorry`. Identified false lemmas in our formalization.
- **Consolidation:** Ported JleE.lean, Algebra.lean, STP2Closure.lean into `Formal/` (the canonical Lean project, pinned Mathlib v4.28.0). The `lean/` directory is a scratch copy.
- **Key finding:** `isLC` definition using ℕ subtraction at k=0 makes `lc_conv` FALSE. Must guard with k ≥ 1.
- **Formal/ directory now contains:**
  - `Basic.lean` (610 lines, algebraic infrastructure, 0 sorries)
  - `P3.lean` (201 lines, leaf-swap injection, 1 sorry: `tree_has_pendant`)
  - `JleE.lean` (58 lines, J ≤ E, 0 sorries)
  - `Algebra.lean` (88 lines, star+star w_2 + binomial LC, 0 sorries)
  - `STP2Closure.lean` (245 lines, open problem, 2 sorries: the main theorem + multi-child corollary)
- Project disposition unchanged: shelved pending proof breakthrough or paper submission.

## Current state (2026-03-14)

### Update (2026-03-14)
- **AI-contributions wiki entry accepted.** GitHub issue #257 on `teorth/erdosproblems` closed as completed. Entry added to Section 2(e) ("AI tools used to perform numerical exploration") with two thumbs up. Conservative placement was the right call.

## Current state (2026-03-12)

### Session notes (2026-03-12)
- Prepared copy-ready outreach texts under:
  - `outreach/forum_comment_993.txt`
  - `outreach/erdosproblems_wiki_993.txt`
  - `outreach/arxiv_endorsement_request.txt`
- Confirmed that public comments for problem `#993` are posted on the forum thread:
  - `https://www.erdosproblems.com/forum/thread/993`
  - posting requires forum login; the main problem page itself does not expose the comment box directly
- Checked the existing thread state:
  - Jake Mallen had already posted a January 7, 2026 note claiming `n <= 29` verification
  - accordingly, the new forum message was framed as a follow-up centered on the public repo, exact `n=28` LC / near-miss figures, and the current structural manuscript status
- Posted the forum follow-up on the `#993` thread.
- Opened a GitHub issue on `teorth/erdosproblems` proposing a conservative AI-contributions wiki addition for `#993` in section `2(e)`; no `data/problems.yaml` PR was opened.
- Prepared private arXiv endorsement-request text using code `GTGLTK`; David Galvin had already been contacted previously and had not replied by this session.
- Immediate next step is passive:
  - wait for responses on the forum issue / GitHub issue / private endorsement outreach rather than opening more public threads

## Current state (2026-03-11)

### Session notes (2026-03-11)
- Added prospect-aware archive scoring to:
  - `nm_optimizer.py`
  - `scripts/lc_breaker_optimizer.py`
- Added the March 2026 Ramsey/AlphaEvolve paper to:
  - `paper/references.bib`
  - `paper/main_v2.tex`
- Rebuilt `paper/main_v2.tex` successfully with XeLaTeX + biber after the new citation.
- Evaluated the remaining obvious compute gap:
  - `n=29` LC + near-miss
  - estimated Modal cost is on the order of `$1k` before credits at the established `1024`-worker setting
- Dispatched then stopped a trial Modal app for `n=29` LC/NM:
  - app id: `ap-D5xjmI0DeGpzgxaxo1wQUG`
  - final state: `stopped`
  - no result artifact collected; `results/analysis_n29_modal_lc_nm.json` still does not exist
- Project disposition:
  - exhaustive unimodality result through `n=29` stands as the final computational frontier
  - no further large-compute work is planned under current budget
  - treat the project as shelved unless a cheap paper-only cleanup task arises

## Current state (2026-03-06)

### Session notes (2026-03-06)
- Repaired Modal shard fanout by adding `dispatch_missing` entrypoints to:
  - `search_modal_exhaustive.py`
  - `search_modal_exhaustive_n29.py`
  - `analyze_modal_lc_nm.py`
  - `analyze_modal_lc_nm_n28.py`
- Collected final `n=28` Modal artifacts:
  - `results/analysis_n28_modal_unimodality.json`
  - `results/analysis_n28_modal_lc_nm.json`
- Final `n=28` exhaustive unimodality result:
  - `2,023,443,032` trees
  - `0` unimodality failures
  - tree count matches OEIS A000055
- Final `n=28` exhaustive LC + near-miss result:
  - `2,023,443,032` trees
  - `0` non-unimodal trees
  - `19` log-concavity failures, all at `k = 14`
  - worst LC ratio `1.5027777777777778`
  - best near-miss ratio `0.8565665724120973` at `k = 13`
- Exhaustive unimodality frontier is now:
  - `8,691,747,673` trees through `n = 29`
  - `0` unimodality failures
- Collected final `n=29` exhaustive unimodality artifact:
  - `results/analysis_n29_modal_unimodality.json`
  - `5,469,566,585` trees
  - `0` counterexamples
  - tree count matches OEIS A000055
- The last `n=29` shard showed pathological binary-residue skew; completion was diagnosed and independently verified via an odd-mod rescue split before the main dict closed.

## Current state (2026-03-05)

### Session notes (2026-03-05)
- Cleaned up stale docs after the committed `n=27` LC + near-miss artifact:
  - `README.md`, `STATUS.md`, `paper/main_v2.tex`
  - `results/analysis_n27_modal_lc_nm.json` now reflected in top-level summaries
- Added a dict-backed collector for Modal shard jobs:
  - `scripts/collect_modal_results.py`
  - supports `status` snapshots and `collect` merges for both unimodality and LC/NM runs
- Hardened Modal wrappers for `n>=28`:
  - raised worker timeout from `2h` to `12h`
  - fixed stale wrapper defaults in `search_modal_exhaustive_n29.py` and `analyze_modal_lc_nm_n28.py`
  - added `launch_partitions` to those wrappers for parity with the base scripts
- Started new large-`n` runs:
  - `n=28` exhaustive unimodality
  - `n=28` exhaustive LC + near-miss
  - `n=29` exhaustive unimodality
- Live snapshot at logging time:
  - `n=28` unimodality: `396/1024` partitions, `614,576,414` trees, `0` counterexamples so far
  - `n=28` LC/NM: `75/1024` partitions, `97,539,110` trees, `0` non-unimodal, `2` LC failures so far, best `nm=0.8565665724120973`
  - `n=29` unimodality: `36/1024` partitions, `111,044,351` trees, `0` counterexamples so far
- Session log:
  - `notes/modal_jobs_2026-03-05.md`

## Current state (2026-03-04)

### Session notes (2026-03-04)
- Folded in the committed `n=27` LC + near-miss Modal artifact:
  - `results/analysis_n27_modal_lc_nm.json`
  - 751,065,460 trees, 0 unimodality failures, 0 log-concavity failures
  - best near-miss ratio `0.8571425274916726` at `k=13`
- Implemented and pushed deterministic runtime packages:
  - `pi_n/` (exact-rational `Pi(n)` with transcript certificates + replay verifiers)
  - `orchestrator_v13/` (authoritative gating, partition selection, queueing, obligations, replay checker)
- Added fixture-based conformance harness for orchestrator v13:
  - `orchestrator_v13/tests/fixtures/T0..T7`
  - byte-stable golden artifacts under `orchestrator_v13/tests/goldens/`
  - replay and golden tests passing (`test_fixtures.py`, `test_golden_bytes.py`)
- Added direct data extractor for orchestrator input from Modal lambda-frontier outputs:
  - `scripts/build_orchestrator_input_from_lambda_frontiers.py`
  - uses script-native semantics (`X = Lambda-D`, `R = R_shift`, `sum_all = Σ err_s`)
- Produced and committed reproducible machine artifacts:
  - mixed-source normalization run:
    - `results/orchestrator_v13_input_from_modal_n24_n25.json`
    - `results/orchestrator_v13_run_n24_n25/*.json`
    - `results/pi_n_cert_from_modal_n24_n25.json`
  - direct lambda-source run:
    - `results/orchestrator_v13_input_direct_lambda_n22_n25.json`
    - `results/orchestrator_v13_run_direct_lambda_n22_n25/*.json`
    - `results/pi_n_cert_from_direct_lambda_n22_n25.json`
  - provenance notes:
    - `notes/orchestrator_v13_n24_n25_run_2026-03-03.md`
    - `notes/orchestrator_v13_direct_lambda_run_2026-03-03.md`
- Direct lambda-source artifact status:
  - orchestrator v13 selected `PI0`, `Phi=0`, global closure `CLOSED`
  - `Pi(n)` certificate status `PASS` with replay verification successful

## Current state (2026-02-16)

### Session notes (2026-02-16, night)
- **n=27 exhaustive search COMPLETED** via Modal cloud compute
  - 751,065,460 trees, all unimodal, 0 counterexamples
  - Tree count matches OEIS A000055 exactly
  - Modal: 1024 workers, 78 minutes wall time, app ID: ap-T9RkZ9fGOtXyvZgPEuYBkZ
  - Results saved to `results/analysis_n27_modal.json`
- Set up Modal account (workspace: brettrey), applied for academic credits ($25K program)
- Created `search_modal.py` (persistent Dict, streaming progress, --detach mode)
- Updated `paper/main_v2.tex`: Modal app ID in reproducibility appendix, companion biology paper mention in Discussion
- Killed orphaned multiprocessing workers from earlier local search attempt
- **Lessons learned**: never kill running processes without asking; Modal --detach can leave stale local clients; Python multiprocessing can orphan workers; check logs before speculating

### Session notes (2026-02-16, afternoon)
- Reviewed Gemini 3's unsolicited patent application ("Hub Exclusion Scheduling")
- Assessment: not patentable (prior art: crown reduction is decades old; math theorems aren't patentable; internal inconsistencies; thin evidence)
- No patentable applications of the project's math results identified
- Deleted all 8 patent-related files (PDF, tex, draft, scripts, figures)

### Session notes (2026-02-16, morning)
- Rewrote `plot_roots_n26.py` using house style (EB Garamond via `text.usetex`, house palette, frameless legend, zoomed clip at |z| < 1.4)
- Integrated root plot as Figure 1 in Section 5 of `main_v2.tex`
- Drafted email to David Galvin (dgalvin1@nd.edu) for feedback on the paper
- Created biology paper folder: `papers/Tree_Independence_Polynomials_and_Biological_Network_Motifs/`
- n=27 exhaustive search launched locally (8 geng workers, est. 30-40 hours) -- superseded by Modal

## Previous state (2026-02-15)

### Recent Progress (Proof Push)

- **Analytic Subdivision Lemma:**
    - Analytically proved that $i_k(T/e) \le i_k(T)$ for all $k$.
    - Verified this bound on 123,867 trees ($n=18$).
    - This bound, combined with ECMS, proves that subdivision preserves unimodality.
- **Edge Contraction Mode Stability (ECMS):**
    - Verified the conjecture $|\text{mode}(I(T)) - \text{mode}(I(T/e))| \le 1$ on all trees up to $n=18$.
    - This conjecture implies that minimal counterexamples have no degree-2 vertices.

### Paper v2: "A subdivision-contraction identity and structural reductions"

The paper has been completely rewritten to foreground proved theorems rather than just computational verification. Target venue: Experimental Mathematics.

**Proved results in the paper:**
1. **Subdivision-contraction identity** (Theorem 3.1): I(T_e) = I(T) + x I(T/e)
2. **Conditional subdivision lemma** (Theorem 3.3): ECMS implies subdivision preserves unimodality
3. **Private Neighbor Bound** (Lemma 4.1): P >= n - 2k + 1
4. **Hub Exclusion** (Lemma 4.3): d_leaf(v) >= 2 and 1-Private => v not in S, leaves in S
5. **Transfer Lemma** (Lemma 4.4): 1-Private transfers to residual tree
6. **Spider mean bound** (Proposition 4.6): S(2^k,1) has n/3 - mu -> 1/6
7. **Edge bound** (Theorem 4.7): P(u) + P(v) < 2/3 for every tree edge
8. **Leaf-attachment asymptotics** (Theorem 6.2): nm(s) = 1 - C/s + O(1/s^2), C in [4,8)

**Open conjectures in the paper:**
- **ECMS** (Conjecture 3.2): |mode(I(T)) - mode(I(T/e))| <= 1. Verified 24.7M edges, 0 violations.
- **Conjecture A** (Conjecture 4.5): d_leaf <= 1 => mode <= floor(n/3)+1. Verified 528K trees through n=23.

### Literature search (2026-02-15)

All major claims of originality vetted against published literature:
- Subdivision-contraction identity: **Novel** (not found in surveys by Levit-Mandrescu or edge elimination polynomial framework)
- PNP reduction (Transfer Lemma + d_leaf <= 1): **Novel**
- ECMS: **Novel**
- Edge bound P(u)+P(v) < 2/3: **Novel** as stated (related to hard-core model literature)
- Multi-arm stars + near-miss ratio: **Novel**
- Leaf-attachment asymptotics: **Novel**
- n=26 exhaustive verification: Kadrawi & Levit (2023) checked LC at n=26 (implicitly confirming unimodality); our contribution is explicit unimodality check + near-miss metrics. Paper now acknowledges this.

### Exhaustive verification through n = 29

No unimodality violations among all 8,691,747,673 non-isomorphic trees on n <= 29 vertices. Tree counts match OEIS A000055.

| n | Trees | Time / notes |
|---|------:|-------------:|
| 1--15 | 13,188 | <1s |
| 16 | 19,320 | 1s |
| 17 | 48,629 | 3s |
| 18 | 123,867 | 9s |
| 19 | 317,955 | 23s |
| 20 | 823,065 | 68s |
| 21 | 2,144,505 | 55s |
| 22 | 5,623,756 | 1m 44s |
| 23 | 14,828,074 | 4m 41s |
| 24 | 39,299,897 | 12m 5s |
| 25 | 104,636,890 | 38m 33s |
| 26 | 279,793,450 | 4h 51m |
| 27 | 751,065,460 | 78m (Modal, 1024 workers) |
| 28 | 2,023,443,032 | Modal dict-backed (1024 workers) |
| 29 | 5,469,566,585 | Modal dict-backed + odd-mod tail rescue |
| **Total** | **8,691,747,673** | |

n=26 details:
- Exactly 2 log-concavity failures (both at k = 13), matching Kadrawi & Levit (2023).
- Best near-miss ratio nm = 0.845.

n=27 details:
- Computed on Modal cloud (app ID: ap-T9RkZ9fGOtXyvZgPEuYBkZ)
- LC + near-miss follow-up completed (`results/analysis_n27_modal_lc_nm.json`)
- 0 log-concavity failures
- Best near-miss ratio `nm = 0.8571425274916726` (first tail rise candidate at `k = 13`)

n=28 details:
- Unimodality artifact collected (`results/analysis_n28_modal_unimodality.json`)
- LC + near-miss artifact collected (`results/analysis_n28_modal_lc_nm.json`)
- 0 unimodality failures
- 19 log-concavity failures, all at `k = 14`
- Worst LC ratio `1.5027777777777778`
- Best near-miss ratio `nm = 0.8565665724120973` (first tail rise candidate at `k = 13`)

n=29 details:
- Unimodality artifact collected (`results/analysis_n29_modal_unimodality.json`)
- 0 unimodality failures
- Tree count matches OEIS exactly: `5,469,566,585`
- LC + near-miss metrics have not yet been run at `n = 29`

### Computational certificates

| Check | Count | n range | Failures |
|-------|-------|---------|----------|
| Unimodality (exhaustive) | 8,691,747,673 trees | <= 29 | 0 |
| LC + near-miss (exhaustive) | 3,222,181,088 trees | <= 28 | 21 LC failures |
| I(T_e) = I(T) + x I(T/e) | 66,697 edges | <= 14 | 0 |
| ECMS | 24,710,099 edges | <= 20 | 0 |
| A(x) unimodal and LC | 9,071,864 edges | <= 19 | 0 |
| Combined tail | 9,071,864 edges | <= 19 | 0 |
| xR_uR_v ascending before mode | 24,710,099 edges | <= 20 | 0 |
| |delta(T,v)| <= 1 | 26,056,121 pairs | <= 20 | 0 |
| Conjecture A | 931,596 trees | <= 23 | 0 |
| mu < n/3 (d_leaf <= 1) | 931,596 trees | <= 23 | 0 |
| Case B bound | 8,710,881 trees | <= 22 | 0 |

### Targeted search on structured families (n up to 500)

145,362 tested trees across five families, no unimodality violations.

| Family | Trees | LC failures | Best nm |
|--------|------:|----------:|---------|
| Galvin SST T_{m,t,1} | 571 | 108 | 0.936 |
| Generalized SST T_{m,t,d} | 680 | 268 | 0.981 |
| Caterpillars | 5,196 | 0 | -- |
| Spiders and brooms | 133,915 | 0 | 0.992 |
| Random Ramos-Sun-style | 5,000 | 2 | 0.804 |
| **Total** | **145,362** | **378** | |

Multi-arm stars surpass brooms as the true extremal family. Champion at n >= 200: M(s; 5,5,4,2).

## Artifacts

- `paper/main_v2.tex` -- current manuscript (compiles cleanly)
- `paper/figures/roots_n26_lc_failures.pdf` -- root plot (Figure 1)
- `plot_roots_n26.py` -- generates the root plot
- `email_galvin.md` -- draft email to David Galvin for feedback
- `paper/main.tex` -- previous version (9pp, computational verification + broom asymptotics)
- `search_modal.py` -- Modal cloud search script (n=27)
- `results/analysis_n26.json` -- n=26 exhaustive LC and near-miss summary
- `results/analysis_n27_modal.json` -- n=27 exhaustive unimodality check (Modal)
- `results/analysis_n27_modal_lc_nm.json` -- n=27 exhaustive LC and near-miss summary (Modal)
- `results/analysis_n28_modal_unimodality.json` -- n=28 exhaustive unimodality check (Modal)
- `results/analysis_n28_modal_lc_nm.json` -- n=28 exhaustive LC and near-miss summary (Modal)
- `results/analysis_n29_modal_unimodality.json` -- n=29 exhaustive unimodality check (Modal)
- `results/targeted_n500.json` -- targeted search summary + top near-misses
- `results/targeted_families.json` -- per-family summary for the targeted search
- `notes/shutdown_2026-03-18.md` -- end-of-session handoff note with current DOI, release tag, and submission posture
- `subdivision_correct.py` -- definitive subdivision identity analysis
- `verify_subdivision_formula.py` -- subdivision-contraction identity verifier
- `notes/one_private_status.md` -- definitive PNP framework
- `notes/conjecture_A_analysis.md` -- Conjecture A reduction details

## Immediate next actions

1. Submit the manuscript to `European Journal of Combinatorics` or `Discrete Mathematics`.
2. Treat any Kadrawi or other reader feedback as parallel input; do not block submission waiting on it.
3. If another paper-only DOI snapshot is needed, use Zenodo's direct `New version` flow instead of GitHub-release retries.
4. Leave the project shelved unless there is a genuinely new idea on ECMS or Conjecture A.

## Dead ends (do NOT revisit)

See `notes/scc_false_n28_2026-03-01.md` and the other obsolete notes in `notes/` for dead ends not to revisit. Key ones:
- Core avg < 1/3: FAILS at n=7+
- Matching cover / SDR / Hall's condition: all FAIL
- Product mode domination: FALSE
- delta(T,v) >= 0: FALSE
- A <= (1+x)I for gap cases: too loose
- Leaf-Mode Inequality: FALSE at n=62
