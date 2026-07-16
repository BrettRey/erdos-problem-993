# Decisions Log
<!-- SUMMARY: Append-only decisions governing proof scope, verification standards, release posture, and bounded AI-assisted work on Erdős #993 · status: active decision log · updated: 2026-07-15 -->

Append-only record of project decisions. Agents: add an entry whenever a non-trivial decision is made during a session (structural changes, venue choices, theoretical commitments, scope changes, reviewer feedback acted on). Keep entries short.

Format: `## YYYY-MM-DD` then bullet points with **bold topic** and brief rationale.

---

## 2026-07-16

- **Flat-mode contamination found in the raw valley margin; metric corrected**
  Component diagnosis of the top depth-3 scaling ray showed its scoring index
  b coincides with the global mode (no prior descent; root-included component
  nowhere dominant), so the raw margin `V` was measuring mode flatness there,
  not a dip. At a smooth mode the deficit factorizes as (local ratio slope,
  ~C/n) x (fractional lattice offset of the r=1 crossing, quasi-random in
  (0,1)), which fully explains the apparent deficit "acceleration" at
  n=19,777 and the non-monotone n(1-V) constants in the dumbbell rays:
  lattice wobble, no new mechanism. Plateau-approach can push V arbitrarily
  close to 1 without ever crossing, so V near 1 is not evidence of a nearby
  counterexample. Corrected ranking metric: descent-thresholded rebound
  `R_gap(theta)` (rebound suffix-max/i_b over dips with prior descent at
  least theta, theta = 0.001 and 0.01, with rise distance c-b reported);
  witness detection (V > 1 exactly) is unaffected. Killed the mix
  re-optimization and ray-extension runs that were optimizing the artifact;
  re-scored the depth-3 grid and all 2026-07-15 champions under R_gap
  (`scratch_rgap_rescore_20260715.py`).
- **Corrected barrier is stronger than the raw one** Under R_gap(0.001) the
  best genuine-descent rebounds are n(1-R) ~ 10.9-13.5 (depth-3 grid best
  10.87 at n=1378; two-phase shoulders 13.3-13.5 at n~320), the constant
  grows with n (20.5 at n=2537, 16.8 at n=6563) and with theta (23-42 at
  theta=0.01), and every achieving rebound in ~3,500 rescored configs plus
  all champions has rise distance c-b=1: locally slowed decay, never a
  cross-gap recovery. No witnesses. The disproof-relevant gap to V=1 is
  roughly 11/n at best and widening with size, so no family in the searched
  classes is on a crossing trajectory.
- **Aperiodic composition attack closed** Trees grown by iterating five
  pairs of rooted maps under Fibonacci, Thue-Morse, periodic, and random
  control words (103 constructions, exact map/adjacency cross-check):
  certified root angle 0 at every dominance ratio <= 1.5 for every
  control class, best 0.095 at ratio 2 (below T_{6,6,1}). Aperiodic
  ordering does not produce exotic root geometry when each composition
  map individually smooths. No witnesses; c-b = 1.
- **Dominance-angle frontier measured; root-space attacks closed**
  Pareto evolution with per-threshold champions and subtree crossover
  (85,114 certified evals): phi*(tau) = 0, 0, 0.054, 0.103, 0.181,
  0.324, 0.449 at tau = 1.05, 1.1, 1.2, 1.35, 1.5, 1.75, 2.0. The two
  strictest slots never left zero; the curve extrapolates to zero near
  tau ~ 1.1, i.e., an empirical exclusion zone: no certified non-real
  root pair within ~10% of dominant modulus was ever produced by any
  search. R <= 0.957 and c-b = 1 throughout. The single-champion
  tradeoff of the 611k-eval run is therefore a universal frontier at
  this search power, and the spectral route to a counterexample requires
  entering a region no tree has ever exhibited. Frontier trees saved to
  `results/pareto_root_frontier_20260716.json`.
- **Bump-interference attack (exact pathology blocks) closed** 540,492 exact
  joins pairing the 21 exhaustive LC-failure trees (n=26/28) with each other
  at every vertex pair and with smooth partners sized to translate the k=14
  bump into the composite mid-band
  (`scratch_bump_interference_20260716.py`): zero witnesses, zero rebounds
  with c-b > 1. The bump does survive embedding as the campaign's best
  constant (n(1-R)=5.3 at n~104, failure tree x broom), but it is a fixed
  28-vertex object that dilutes under scaling; not a counterexample seed.
- **Root-herding attack: float64 pipeline invalidated, certified frontier
  established** A claw-free control (P_101 must be real-rooted) exposed
  numpy-roots noise: the float pipeline reported phi_dom=0.685 for P_101
  and 0.949 for an evolved champion whose certified value is 0.045.
  Certified Arb root isolation (python-flint 0.9.0, scratch venv;
  `scratch_root_herding_20260716.py` + certified recheck) gives the true
  near-dominant-root angular frontier: phi_dom <= 0.256, maximized by
  T_{6,6,1} (the LC-fragile class), failure trees 0.06-0.11, smooth
  families 0. Root geometry is the right search coordinate (angle tracks
  the known convexity pathology exactly), but any future root-herding
  evolution must score with certified roots in the loop (~1-2 s/eval at
  degree ~70). Same lesson as the repo's standing rule: do not trust
  floating-point combinatorics.
- **Certified root-herding evolution run and closed: the dominance-angle
  tradeoff is the obstruction** 611,897 certified evals
  (`scratch_certified_root_evolution_20260716.py`, Arb in the scoring
  loop): phi_dom(|z| <= 2 z_min) rose 0.256 -> 0.444 in the first ~10
  seconds, then zero improvements for 49+ minutes; no witness; R <= 0.90
  and rise distance c-b = 1 on every champion. The 0.26 "ceiling" was
  family bias, but the champion's own certified spectrum shows the real
  invariant: root angle falls monotonically as modulus dominance rises
  (ratio 2.0 -> phi 0.44; 2.5 -> 0.26; 2.7 -> 0.13; 2.76 -> 0.03), heading
  to phi -> 0 at amplitude parity. A counterexample needs a near-dominant
  pair at substantial angle, i.e. the opposite corner of the measured
  frontier; oscillations from subdominant pairs are exponentially damped
  in k. This is the spectral form of the same barrier the coefficient-space
  campaign measured. Champion saved to
  `results/certified_root_evolution_20260716.json`.
- **Subtractive carving channel tested and closed** The last
  composition mechanism outside the D17/D19/D20/D23 smoothing obstructions
  was the edge-join correction `-x^2 I(T1-N[u]) I(T2-N[v])`, which subtracts
  rather than mixes and so could in principle carve a dip. A 1,406-join scan
  (`scratch_join_carving_20260716.py`) over 11 component types at
  structurally distinct join vertices (root/gadget center/mid-leg/leaf,
  connector lengths 1-2) produced no witness, no rebound with rise distance
  c-b > 1, and at most a ~20% improvement in the deficit constant
  (best n(1-R)=8.8 at n=338). Join-vertex choice shifts R by ~0.2%: the
  correction is too smooth to carve, consistent with the manuscript's
  edge bound P(u)+P(v) < 2/3. With additive mixing, depth stacking, free
  mutation, and subtractive carving all measured, no composition mechanism
  reachable from spider/star/path blocks shows a crossing trajectory; a
  disproof would need a qualitatively different source of coefficient
  interference than tree composition provides.

- **Opened a user-directed valley-first disproof campaign** Brett explicitly
  asked for a disproof attempt, overriding the parked posture. Honoring the
  2026-07-11 "future disproof posture" rule, the campaign used direct
  valley-producing architectures scored by the exact witness margin
  `V(T) = max_b min(prefix-max, suffix-max)/i_b` (V > 1 iff counterexample),
  not log-concavity defect. New tools: `scripts/valley_search.py`
  (spider-bouquet grammar sweeps + hill climb, self-tested against the
  generic DP), `scripts/product_valley_search.py` (forest/product lane with
  path-join tree realization), `scripts/valley_scaling_probe.py`
  (Kronecker-packed exact arithmetic for n in the thousands), plus dated
  scratch probes for free-form mutation, deficit-constant mapping, and
  depth-3 phase stacking.
- **No counterexample; ~150k exact evaluations** Zero valley witnesses across
  bouquet sweeps (107,253 specs, n <= 320), 6,441 forest products with
  dumbbell realizations, scaling probes to n = 7,299, a 263k-mutation
  free-form optimizer, and an independent Codex (GPT-5.6 sol) search at
  n <= 3,000 whose top candidates were replayed exactly and matched to
  10 decimals.
- **Empirical barrier law identified** Every architecture (perturbed
  T_{3,M,N}, hub-bouquets, two-type hybrids, subdivided stars, dumbbells,
  depth-3/4 block trees, engineered connectors) obeys
  `1 - V(n) ~ C/n` with best-achieved C converging to about 4, matching the
  lower end of the manuscript's proved leaf-attachment constant C in [4,8).
  The optimized-mix constant is flat in the phase-gap parameter (t = 6..90),
  so widening phase contrast buys nothing. The descent/rebound decomposition
  shows the binding tradeoff: genuine descents (depth 1.014) pair with weak
  rebounds (~1 - 13/n); near-unit rebounds (~1 - 4/n) sit on plateau-scale
  descents (1 + 1e-5). T_{3,M,N}-type shoulders additionally slide outside
  the legal rebound window (rise index passes ceil((2*alpha-1)/3)) as n grows,
  so that mechanism is asymptotically dead regardless of depth.
- **Exhaustive-data motif audit came back empty** The 19 exact n=28
  LC-failure trees are all small hub-with-spider-star variants (the same
  late-shoulder mechanism) and the n=28 near-miss champions are brooms
  (binomial C(25,k) tails); no third mechanism is hiding in the exhaustive
  artifacts.

- **Use Dobriban (2026) only as proof-workflow precedent** The author's reported
  GPT-5.6 Pro discovery of a BH counterexample illustrates a useful separation
  among AI exploration, a narrow analytic theorem, exact replayable
  certification, finite simulation, and human audit. It supplies no
  combinatorial lemma, certificate, or evidence for Erdős #993. Do not reopen
  the search, alter the submitted manuscript, or imply progress on #993 from
  it. The one-shot and cross-version capability story is an uncontrolled
  anecdote and does not displace the project's bounded-target discipline.
- **Require clean-clone verification of the whole release bundle** A local
  clean-clone check replayed Dobriban's standalone `python-flint`/Arb
  certificate, but the advertised full verifier failed because the manifest
  and repository tree disagree about the README, manuscript, and simulation
  artifact paths. Certificate replay is not an independent audit of the
  analytic reduction. For any future #993 release, test the proof reduction,
  exact certificate, simulations, and complete package separately from a clean
  checkout before making a reproduction claim.

## 2026-07-11

- **Reduce Erdős #993 to prefix GSB and treat `r <= 3` as the proved base** The
  variance/moment lane (D6/D13) reduces tree unimodality to the prefix
  inequality `Var(e_r) <= 2 mu_r + 2 eta_r` for `r <= L-2`, proves it
  unconditionally for `r=1,2` (D6, discriminant argument) and as an audited
  symbolic theorem for `r=3` (D12, `alpha>=7 ==> 5 i_3 i_5 <= 4 i_4^2 + i_3 i_4`).
  A minimal failure now has `r>=4`. This is a reduction plus a base, not a
  solution; the resolution gate is unmet.
- **Record the rank-four limit of the rank-three radius-two tuple** Two order-13
  trees with identical recorded radius-two moments differ at `i_6` (170 vs 168),
  so the rank-three certificate cannot be extended from those recorded
  statistics alone. This does not establish that every `r=4` proof requires
  radius-three data.
- **Block the D14 down--up spectral route at the separated energy estimate** The
  separated aggregate energy inequality is false at an exact audited `n=210`
  witness, though that tree is still unimodal. Reopen D14 only with a corrected
  non-separated energy inequality, not with broader spherical scans.
- **Treat the whole solve-today round as falsification and localization, not
  resolution** The numbered round through D27 is predominantly obstruction and
  localization work, mostly through exact counterexamples to proposed
  invariants (D1 at `n=29`, unedged `Var<=2 mu` at
  `T_(6,6,1)`, unrestricted GSB at `T_(14,8,1)`, D7 flows at small order). No
  manuscript change, release, or Erdős Problems announcement follows.
- **Close the bounded augmented-fork cycle as refuted** Exact small-tree and
  random scans initially left `AF_lambda` viable, but a structured same-branch
  search produced an exact last-prefix obstruction on a center joined to five
  copies of `G(60,18)`. Its gap is negative at `lambda=0` and decreases with
  `lambda`, so every `lambda in [0,1]` fails. This blocks only that fixed-parameter
  local allocation, not GSB. Aggregate `q_far<=0` remains separate and parked.
- **Do not promote full marked flow as a reduction** Its total target-capacity
  cut is exactly GSB. Saturation through order 11 is finite diagnostic evidence
  only.
- **Park direct #993 search after the D25 closeout** The remaining aggregate
  `q_far<=0` proposition is theorem-strength and has no live local mechanism.
  Do not open D28. The next high-value activity is publication triage of the
  universal Poisson-binomial reserve theorem; if proof work later resumes, use
  the bounded fixed-arm `A+xR` perturbation target.
- **Separate disproof engineering from proof-route falsification** A
  counterexample to an auxiliary lemma is not evidence against #993. Any future
  negative program must begin with a concept-level architecture for an exact
  descent followed by a later ascent, explain what escapes D17/D19/D20/D23,
  and only then authorize computation.

## 2026-07-10

- **Treat the universal finite Poisson-binomial reserve bound as a proved private theorem, not as a solution of Erdős #993** Exact symbolic replay, independent mathematical audits, and the completed conditional Lean core support the theorem itself. It closes the finite signed nonterminal bridge, but it does not control the additive `xR` perturbation in the hub-bouquet argument and therefore does not close Case B, PNP, issue #5, or the original problem.
- **Keep the current manuscript unchanged pending a deliberate publication decision** The result may justify a standalone note or a later manuscript revision, but migration requires a Hillion--Johnson citation and manuscript-grade presentation of the exact certificates. No release, Erdős Problems announcement, or solution claim should be made from the current repository state.
- **Use fixed-arm `A+xR` perturbation as the next proof target** The next bounded research lane is separate pre-descent and post-descent increment control for fixed arm length. Growing brooms and end-to-end Lean formalization are later targets; the latter is tracked separately in GitHub issue #7, while issue #5 remains open.

## 2026-07-04

- **Switch issue #5 signed-reserve work to audit-hardening mode** External audits found real errors in the recent proof route (missing X-side boundary term, fair-binomial constant error, strict-transfer overstatement, float plateau risk). The durable stance is now: patch audit findings before extending the route, keep theorem/algebra/computation/conjecture/disproof separated, regression-test every concrete error, and do not close issue #5 or claim signed reserve / hub-bouquet reserve / Erdos #993 from these notes.
- **Use `notes/literature/audit_hardening_claims_ledger_2026-07-04.md` as the claim-status surface** The ledger records what is theorem-level, algebraic, computational, superseded, or still blocked. It is a reliability surface, not a proof note; no manuscript import from this lane without a later audit pass.

## 2026-07-01

- **Guard Lean LC/STP2 at `k >= 1`** The old Nat-subtraction boundary created artificial `k = 0` obligations. `Formal/STP2Closure.lean` now treats LC/STP2 as guarded non-boundary conditions and records a regression for `leafI * leafI = [1,2,1]`.
- **Quarantine abstract STP2 closure behind tree-DP realizability** The contiguous-support toy pair `I=[1,4,1]`, `E=[1,1,1]` satisfies the abstract shape package, including `noGaps`, but its self-convolution violates guarded STP2. The broad closure shells are now deliberately vacuous via `treeDPPair`; the next mathematical target is a genuine tree-DP realizability invariant.
- **Expose three public work tracks** Public coordination is now through issues for tree-DP realizability, exact fixed-r certificate emission, and forest/product valley search. The README states that the repo does not contain a proof of Erdos #993 and lists these as July 2026 targets.

## 2026-06-12

- **Close the abstract fixed-`r` Lean bridge lane** Fable packets were useful only after local Lean replay. The checked interface now runs from one-instance bridge lemmas through `Route2FamilyCertificate`; do not keep asking Fable to re-audit the same gap. The next Lean work is either script emission into these records or the larger spider-polynomial modeling block.
- **Treat Li's symmetric-function method as a research lead, not a current bridge dependency** Li's 2026 KL-family result is relevant background, but 2-Schur positivity is essentially a log-concavity certificate. The reusable possibility is partial 2-Schur positivity plus tail monotonicity, not an immediately available projectible invariant for this Lean packet.

## 2026-06-05

- **Reopen only the fixed-`r` formalization lane after LEAP** The LEAP paper increases confidence in bounded Lean-agent workflows, not in broad #993 proof search. Added a checked mean-shift bridge in `gpt_attack/axiom_fixed_r_certificate/problem.lean`: finite-support variance lemmas plus a derivative-identity theorem replacing the opaque Lipschitz certificate with explicit smaller obligations.

## 2026-05-27

- **Aim AxiomProver-style runs at fixed-`r` certificates first** Recent AxiomProver papers suggest the highest-yield pattern is narrow, formalizable algebraic subgoals with exact certificate replay. For this project, the best match is the fixed-`r` Route-2 certificate bridge rather than a direct attack on Erdős #993, ECMS, or generic log-concavity. Added `gpt_attack/axiom_fixed_r_certificate/` as the handoff packet.
- **Mothball main project again after packet creation** The manuscript is already under E-JC review and no immediate high-return proof/computation path is open. Keep the project parked unless E-JC reviews arrive, AxiomProver or another strong Lean agent is available for `gpt_attack/axiom_fixed_r_certificate/`, or a deliberately bounded session targets only the finite-support Gibbs mean-shift lemma.

---

## 2026-05-24

- **Treat DeepMind AlphaProof Nexus paper as methodology, not #993 input** Reviewed `/Users/brettreynolds/Downloads/2605.22763v1.md`. The solved Erdős list does not include #993, and the highlighted graph-theory result is not about independence-polynomial modes. The durable takeaway is workflow-level: if #993 proof-search resumes, aim agents at narrow Lean/certificate subgoals with frozen statements, not at broad informal proof discovery or manuscript revision.

## 2026-05-21

- **Ship fixed-`r` Route-2 work as a proof workspace, not a submission target** The session produced a substantial fixed-`r`/Route-2 attack package, including scripts, notes, and certificate snapshots, and pushed it as commit `6d2c29c`. The mathematical posture remains conservative: `notes/fixed_r_proof_note_clean_2026-05-21.md` states a fixed-`r` computable-threshold criterion for structured lanes `S(2^a,r)`, not an unconditional theorem for all fixed `r` or for all `d_leaf<=1` trees.
- **Do not submit the fixed-`r` material now** The current manuscript remains under review at E-JC. The new fixed-`r` material should be kept as a private technical reserve unless the fragile proof points are audited, `r=2,3` are handled, and the certificate machinery is made manuscript-grade or extended to a broader tree class.
- **Exclude large/noisy artifacts from the shipped commit** The shipped commit intentionally left untracked `erdos_993_gpt_attack_2026-04-24.zip`, `gpt_attack/main_v2.pdf`, and the 43 MB `results/route2_spider_lane_j0_1_a200_r80_float.json`; these are not needed for the replayable fixed-`r` certificate snapshot.

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
