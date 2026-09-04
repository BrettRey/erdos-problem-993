# Decisions Log
<!-- SUMMARY: Append-only decisions governing proof scope, verification standards, release posture, and bounded AI-assisted work on Erdős #993 · status: active decision log · updated: 2026-08-08 -->

Append-only record of project decisions. Agents: add an entry whenever a non-trivial decision is made during a session (structural changes, venue choices, theoretical commitments, scope changes, reviewer feedback acted on). Keep entries short.

Format: `## YYYY-MM-DD` then bullet points with **bold topic** and brief rationale.

---

## 2026-08-08

- **Retire EADS as a proof route; preserve its counterexample as an auxiliary
  certificate only.** Exact independent replay verifies an order-60 tree for
  which every vertex split has modal-interval distance two. Its full
  independence polynomial is unimodal and log-concave, so this refutes EADS
  but does not refute Erdős #993. Certificate and verifier:
  `results/eads_counterexample_n60_20260808.json` and
  `verify_eads_counterexample_20260808.py`.
- **Discard the rooted-product 0.9999999769 score as a first-descent artifact.**
  Exact coefficient replay shows monotone decay after that index and no later
  rebound. Future searches must score a rise strictly after a genuine prior
  descent; plateau proximity at the first descent is not a valley signal.
- **Pause the campaign unresolved and decaffeinated.** At Brett's morning
  stopping point, three remaining exploratory searches were interrupted and
  all `caffeinate` processes stopped. No negative search is promoted to a
  theorem. Full handoff: `notes/session_handoff_2026-08-08.md`.
- **Record the fixed-arm threshold as an effective refinement, not a new
  family theorem.** For
  each fixed multi-arm vector, finitely many early differences are dominated
  explicitly by the `q_0=1` binomial contribution; after the fixed correction
  support ends, real-rootedness of the path-product term supplies a unimodal
  tail. This proves eventual unimodality with a computable threshold and,
  after exact finite replay, all-`s` unimodality for four hard vectors.
  However, every multi-arm star is a spider and the stronger all-spider
  log-concavity theorem is already known. The genuinely open extension is
  multi-hub structure. Record:
  `notes/fixed_arm_effective_unimodality_2026-08-08.md`.
- **Reject pairwise synchronization as the two-hub bridge.** Exact double-star
  witnesses show that both corrected synchronization/partial synchronization
  for the adjacent-hub summands and synchronization of the two terms in
  `I(T_e)=I(T)+xI(T/e)` can fail while the resulting polynomial remains
  log-concave. Any `C_2` proof must retain aggregate Turán-gap compensation.
  Certificate: `results/two_hub_ratio_order_obstructions_20260808.json`.
- **Treat the aggregate factorization as a hand reproof, not a new base
  theorem.** Galvin--Hilyard Theorem 1.9(1) already includes all double stars.
  After factoring off the smaller leaf side, the polynomial is a log-concave
  product plus one `x`; only the index-2 Turán inequality becomes harder, and
  its gap is a concave quadratic with nonnegative endpoint values. The short
  argument replaces a computational inequality check in the known base case.
  Record:
  `notes/double_star_log_concavity_2026-08-08.md`.
- **Close arbitrary hub-connector subdivision for double stars.** The
  Fibonacci recurrence `F_t=F_{t-1}+xF_{t-2}` and partial-synchronization
  closure reduce the whole unit-arm connector family to two base relations.
  Product-binomial-basis expansions prove both relations coefficientwise for
  all leaf counts, so every double star with an arbitrarily subdivided hub
  connector is log-concave. Ordinary synchronization still fails at
  `p=q=2,t=3`; partial synchronization is essential. Exact finite replay is
  an audit of the general identities, not the proof. Subdivided pendant arms
  remain open. Record:
  `notes/connector_partial_sync_route_2026-08-08.md`.
- **Promote the bounded pendant-core transfer to a computer-assisted
  theorem, not to full `C_2`.** For arbitrary pendant-arm multisets, the
  connector recurrence reduces all connector lengths to log-concavity of two
  base polynomials and two partial-synchronization relations. Exact integer
  enumeration verifies those four hypotheses for all 163,523 unordered arm
  pairs of total pendant weight at most 24. Hence every such core is
  log-concave for an arbitrary number of internal connector vertices. The
  pendant bound is computational and must not be extrapolated; the universal
  base invariant and Erdős #993 remain open. Record and certificate:
  `notes/c2_bounded_pendant_core_2026-08-08.md` and
  `results/c2_bounded_pendant_core_20260808.json`.
- **Publish the double-broom result as a repository note, not a paper.** Brett
  declined a journal article or standalone preprint. The curated public
  package is the proof note, exact replay code, two certificates, focused
  tests, and a README entry, with substantive generative-AI assistance and
  lack of peer review stated explicitly. It was pushed to public `master` as
  commit `51885d3`. The submitted manuscript remains untouched; no Zenodo
  deposit or journal submission is planned, and the exploratory overnight
  artifacts remain local and uncommitted.

## 2026-07-16

- **Target Electronic Communications in Probability for the standalone Poisson--binomial theorem** Brett approved ECP as the first venue after a Roughdraft review of the full `venue-selection` record. The journal-reader contract is an intrinsic variance-normalized curvature theorem for finite Bernoulli sums, supported by a current ECP neighborhood in independent-sum inequalities, Bernoulli structure, and real-rooted probability-generating polynomials. Draft in the ECP class and pivot early to CPC if the complete proof does not fit legibly within 12 pages or gains a genuine combinatorial lead theorem. No public preprint is required. Proof/certificate audit, cautious novelty language, the page-fit test, AI disclosure, and Brett's prediction-ledger forecast remain pre-submission gates. Record: `submission/venue-decision-2026-07-16.md`.
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
- **(Parallel literature-audit session) Uniform positive-axis sector
  refuted at the source; logged here at shutdown** Jerrum--Patel
  (arXiv:2510.01466v2, Lemma 17) supplies the asymptotic counterexample
  class (subdivided complete binary trees, max degree 3, zeros
  accumulating on the positive real axis); a 41-case certified stress
  corpus (`scripts/stress_literature_root_families.py`,
  `results/literature_root_stress_20260716.json`) confirmed the
  mechanism is not near-dominant at accessible scale (min positive-axis
  angle 1.103 at modulus ratio 20.41 through height 8). Conjecture C is
  replaced by target C' (quantified saddle-sector profile) in
  notes/real_collar_conjecture_2026-07-16.md; the audit also flagged an
  apparent index shift in the lemma's appeal to Buys, worth
  clarification before any citation-bearing use.
- **Bridge packet RESOLVED same-day: REFUTED(T3), zero-location
  hypotheses provably insufficient** Codex (GPT-5.6) returned a
  refutation within its 60-minute budget:
  `P_s = ((1+x)^6+20x)(1+x^2)^{2s}` (base = the repo's split-graph
  control polynomial) satisfies dominance, U_24-avoidance, the 0.9
  sector, and the cusp envelope, yet has a mid-window valley at
  b = 2s ~ alpha/2 for every s. Every lemma is hand-checkable (Rouche
  collar, sine-rule sector, binomial valley algebra) and was
  independently replayed in exact arithmetic by the parent session
  (valleys at s = 3, 4, 10, 50; six exact root boxes; s=3 sequence
  digit-identical). Filed under
  `gpt_attack/bridge_window_unimodality/outcome_2026-07-16/`.
  Consequences: (i) the bridge to #993 cannot be proved from zero
  locations plus positivity alone; the missing input is the tree-DP
  structure that excludes parity-lacunary factors like (1+x^2)^m;
  (ii) issue #1's realizability invariant is now provably necessary,
  not just desirable; (iii) the reformulated bridge target is to prove
  trees admit no cross-gap rebound (c-b > 1), the empirical law the
  whole campaign observed without exception. This is the packet
  discipline working as designed: refutation-as-success redirected the
  program in one dispatch.
- **Bridge-lemma packet frozen and dispatched (mode-1 transfer from the
  AI-Erdos claim wave)** Methods audit of the Chojecki/GPT Erdos claims
  identified the transferable pattern: cross-domain reframing of a
  tightly-scoped statement (the #1196 Markov-chain solve), not one-shot
  problem solving. Applied here as `gpt_attack/bridge_window_unimodality/`:
  frozen target (effective no-valley theorem in the window
  k < 2*alpha/3 from zero-geometry hypotheses, graded sub-targets
  T1 pilot spiders / T2 bounded degree / T3 refutation mode / T4
  degree-uniform), binding verification rules (theorem-vs-evidence
  labels, no float roots, effective constants, refutation counts as
  success), certified data bundle, and the full literature txt set from
  the 2026-07-16 intake. Dispatched to Codex GPT-5.6 with a strict
  deliverable format and self-grading bar. Replay discipline unchanged:
  any returned proof gets hand/exact verification before any claim.
- **Spectral conjectures distilled, literature-placed, and same-day
  kill-tested; constant collar refuted, cusp envelope and sector survive**
  Following the campaign, the spectral content was stated precisely
  (`notes/real_collar_conjecture_2026-07-16.md`) and placed against the
  literature: Csikvari (CPC 2013) proves the smallest-modulus root of
  I(G,x) is unique and real for every graph; Prakash--Sharma
  (FSTTCS 2025, arXiv:2510.09197) quantify the all-roots gap as
  exponentially small in general (paths: 1 + Omega(1/n^2), so a uniform
  all-roots collar is false for trees); Hibi--Kara--Vien
  (arXiv:2604.18824) is complementary gamma-positivity work. The
  first-draft conjecture (uniform real collar for non-real zeros) was
  killed within the hour by its own kill-test
  (`scratch_collar_stress_20260716.py`, 28,998 certified trees): random
  trees reach ratio 1.054 and adversarial minimization reaches 1.01839,
  via real-root collisions whose emergent pairs have angle ~2e-5.
  Surviving conjectures: B' (cusp envelope: angle of non-real zeros
  vanishes as dominance ratio -> 1) and C (uniform positive-axis
  sector: |arg z| >= 0.90 measured over all certified trees). A
  non-unimodal control (split graph K_20 v E_6, valley at k=2 with
  non-real zeros at ratio 30+) shows any spectral bridge to
  unimodality must be asymptotic (window coefficients are governed by
  positive-axis saddle analysis, not the dominant zero), with n <= 29
  exhaustive covering small cases.
- **Cusp envelope and sector both survive dedicated adversarial
  kill-tests; DP invariant made concrete** B' attack (maximize certified
  angle in the modulus band <= 1.2): 25,163 mutations, zero improvement
  over 0.05404. C attack (minimize |arg z|): 17,526 mutations, zero
  improvement over 0.90284. The refuted collar moved within seconds
  under identical pressure, so adversarial rigidity is functioning as a
  truth signal. The issue #1 invariant is now stated concretely in
  notes/real_collar_conjecture_2026-07-16.md section 4: occupation-ratio
  reachability for the extend/merge semigroup (R -> x/(1+R);
  R1R2/x), zeros = parameters where -1 is reachable; verified
  consistent with the normality-region characterization in
  arXiv:2111.06451 (abstract level; full read is the open task, for the
  negative-axis cusp geometry and the exact positive-boundary formula).
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

## 2026-07-17

- **Public outreach page for the Poisson--binomial theorem deployed to brettreynolds.ca.** Brett supplied `first-descent.html` (interactive explainer/game for the variance-scaled Turán inequality). Claude audited every mathematical claim against `paper/poisson_binomial/main.tex` before deployment: theorem, corollary, tail bound, the (n+1)/(3(n-1)) family, the c★ bracket, and the proof-map stations all matched; one directional error in the BMM maximal-mass gloss ("high-variance hill cannot have a tall peak") was corrected to the lower-bound reading (p ≥ (1+12V)^{-1/2}). Deployed at https://brettreynolds.ca/first-descent.html with the current manuscript build hosted as https://brettreynolds.ca/variance-scaled-turan.pdf (PDF is gitignored here; site copy must be refreshed when the manuscript changes). Website publications page now lists the manuscript with no venue claim; ECP submission gates unchanged.
- **Independent AI audit of the Poisson--binomial proof recorded as passed (2026-07-17).** Exact replay, fresh symbolic re-derivation of the scalar section, exact end-to-end tests of Prop 3.1 and Theorem 1.1, HJ source-quotation check, and a clean Lean build found no defects; record and reusable audit scripts committed under `notes/` and `scripts/` (see `notes/poisson_binomial_independent_audit_2026-07-17.md`). This does not by itself discharge the "human proof/certificate audit" submission gate; Brett decides whether to count it. MathSciNet/expert novelty checks, durable certificate DOI, and the prediction-ledger forecast remain open.

## 2026-07-18

- **Valley Hunt outreach page deployed for the tree conjecture itself.** Companion explainer/game at https://brettreynolds.ca/valley-hunt.html: exact BigInt independence-sequence engine, near-miss-ratio scoring per the paper's nm metric, challenges calibrated by fresh local computation (the multi-arm champion configs reproduce Table values 0.9437/0.9575/0.9792 exactly; the two Kadrawi n=26 LC-failure polynomials reproduce results/analysis_n26.json coefficient-for-coefficient). All page claims audited against paper/main_v2.tex. Cross-linked with the First Descent page; publications.html tree entry now carries the companion link. Calibration script: website session scratchpad; engine tests replayed in node before deploy.

## 2026-07-22

- **Open-problem handoff packet built for other models.** `gpt_attack/erdos993_handoff.py` (single self-contained file: brief in docstring, exact std-lib engine, embedded grounded data, self-test) plus a JSON-free zip and source dir under `gpt_attack/erdos993_open_handoff_2026-07-22/`. Framing follows the proven packet pattern: graded sub-targets, refutation-as-success, self-grade bar, binding verification rules (exact integer arithmetic, no float roots, theorem-vs-evidence labels, certificate-or-it-doesn't-count). Data is source-grounded: reproduces both Kadrawi n=26 LC-failure polynomials and all five multi-arm champion ratios (n=75..1000). JSON was dropped because some model uploaders (Kimi, z.ai) reject `.json`.
- **Star-of-stars (GLM-5.2) evaluated and dropped — not original, not important.** GLM proposed "star of stars" trees (central vertex joined to k star-centers) as a new unimodal class. Replayed exactly: the tree is always log-concave over 6500+ configs incl. n>=30 (so on the easy side of the log-concave/unimodal fence), GLM's pass-1 proof was false (its |mode(A)-mode(B)|<=1 claim fails, max gap 6), and pass-2 correctly reduced to a cross-term inequality C_j>=0 that is empirically robust (no negatives in 3082 configs, global min 2) but left unproven (the residual is u_j+w_j<=2, which log-concavity gives only as u_j*w_j<=1, not the sum bound). Literature check: the class is substantially covered by the 2025-26 spider / chromatic-symmetric-function program (Li 2026 arXiv:2603.03025; Li-Yang-Zhang 2-Schur-positivity Corollary 2.4), whose general criterion already targets exactly these families; the all-leaf-count-1 case is the spider S(2^n), a labeled figure there. Verdict: a correct-but-subsumed special case; do not pursue.

## 2026-07-23

- **Reconsider a long-horizon #993 proof run only if a new model arrives with a meaningful allowance reset.** The CDC-style experiment is not untried: on 2026-07-11, Codex session `019f50b2-d9d6-7381-a975-737224dd5d0b` ran GPT-5.6 Sol at Ultra effort with a root agent and three dynamically recycled adversarial agents until the allowance stopped it. It produced the D1--D27 registry, prefix-GSB reduction, audited rank-three proof, and exact obstruction packets, but no full proof or counterexample. With the present weekly allowance nearly exhausted, do not relaunch now. If a new model ships with a reset, a further run may be worth testing, but it must ingest the July 11 and July 16 frontier, blacklist the recorded closed lanes, use a persistent goal and checkpoints, and seek a genuinely new continuation rather than replaying the original prompt from scratch. A model release without a reset is not the trigger.

## 2026-07-31

- **Cite Lin & Li (arXiv:2607.27199) as the AI-assisted-combinatorics context in the introduction.** The 2026-03-11 entry records adding the Nagda--Raghavan--Thakurta Ramsey/AlphaEvolve paper to `paper/references.bib` and citing it in `paper/main_v2.tex` for exactly this purpose. Neither file contained it (no match on `nagda`, `ramsey`, `evolve`), so the slot was empty; either it was dropped in a later revision or the log entry overstated. Filled with `linli2026` after the `ramos2025` sentence in the intro, flagged as outside the tree problem. Lin & Li settled whether the exponent $1/2$ in $\sigma(A)^{1/2} \le \delta(A) \le \sigma(A)^2$ can be improved (it cannot) from a construction developed with Hyra, and released a Lean~4 / mathlib formalization of the sharp statement. Bib entry added directly to `paper/references.bib`, not `references-local.bib`: `main_v2.tex:12` records that local entries were merged into `references.bib` for submission, and `/push-bib` would route a math entry into the linguistics central bibliography. Build replayed clean (30 pages, no undefined citations).
  - **This does not reach E-JC.** `main_v2.tex` is the manuscript submitted 2026-04-25 (OJS 15526); the change lands only if there is a revision round.

- **Lin & Li's formal-verification claim is more careful than "zero `sorry`" implies, and the distinction is the one to copy.** Per the repo's own axiom footprint (`scripts/CheckAxioms.lean`): `growthExponent_lt_two` and `supremum_not_attained` depend only on `propext`, `Classical.choice`, `Quot.sound`, so the analytic upper bound is fully kernel-checked; the supremum, convergence, and quantitative-witness theorems additionally depend on three `native_decide` certificates about a 12-element base-39 block. `native_decide` moves the Lean compiler and native execution into the trusted base. Relevant to `Formal/` if the exhaustive $n \le 29$ verification is ever formalized, where `native_decide` is the obvious shortcut and carries the same asterisk. Standing audit rule (replay locally, check for `sorry` **and** unexpected axioms) is what catches this and needs no change.

- **Hyra does not meet the 2026-07-23 relaunch trigger, and isn't usable anyway.** The agent framework is unreleased; only result artifacts are public (`Tencent-Hunyuan/hyra-results`). The base model Hy3 (295B A21B MoE) is open-weights under Apache 2.0 and available through third-party APIs, but that is a model, not the agent that produced the construction. No allowance reset for this project, so no relaunch. Recorded as evidence bearing on the underlying bet, not as a trigger.

- **Noted for possible outreach: Hyra's published results include an Erdős problem.** `hyra-results` lists an `erdos_min_overlap` track (C₅ improved from 0.380868 to 0.380859 against the SimpleTES baseline, arXiv:2604.19341) and 100 record-breaking packings credited on Erich Friedman's Packing Center as "Found by Haowei Lin". Lin is the natural contact if #993 is ever pitched as a target for that pipeline. Not acted on.

## 2026-08-11

- **Ship gate cleared via override, not clean pass, after `passes.py record` proved unable to log the run correctly.** `/ship`'s pass gate blocked on `build-integrity` and `house-style` (never run). Both procedures were actually run against `paper/main_v2.tex`, the manuscript pinned via `passes.py pin`: (1) build-integrity — includes resolve (`manuscript_fingerprint.py`), no hardcoded TeX Live path, `references.bib` confirmed a documented deliberate regular file rather than a central-bib symlink (this is a standalone math project, not routed through the linguistics central bib), full 4-pass `xelatex`+`biber` build clean with zero undefined refs/citations and zero overfull boxes; (2) house-style — `check-style.py` run directly plus an independent `/check-style` subagent pass, both concluding all 55 findings are false positives (underscore-subscript rule firing on math-mode subscripts, `\path{}` filenames, and Modal/shell command lines with underscored filenames, none of it prose) plus one mis-flagged formula (line 739, a genuine vertex-deletion identity, not an AI false-range tic). Both `passes.py record` calls still logged the run against `paper/main.tex` and immediately came back stale: `cmd_record()` (passes.py:734) calls `fingerprint(project)` directly, without the pin-aware resolution `cmd_status`/`cmd_gate` use (passes.py:465), so `record` ignores a pin. Gate cleared with `passes.py override build-integrity/house-style --gate ship --reason ...` documenting the above rather than re-running against a false pin bug. The `passes.py record` pin bug should be fixed in `Project-Management/tools/passes.py`, not worked around per-project.

- **External partial-sync obstruction packet accepted into the repo after exact replay and independent audit.** Brett supplied a note, verifier script, and results JSON from an external session mapping where the HWZZ partial-synchronization induction fails on the Bautista-Ramos pattern-graph families. Replay of the shipped script was byte-identical. The shipped artifacts under-covered the note's claims (`T12` only to `m = 14`, sampled `T13`/`T17`), so the sweep was extended to the stated ranges before integration (158 instances, 134/134 closed-form checks against arXiv:2603.14204 Cor. 4.8/4.12/4.13/4.14, read from `notes/literature/arxiv_2603_14204.txt`); an independent audit (fresh violation/Turán/break code, brute-force subset enumeration on four small trees, including an independently reconstructed `T_{5,5,2}`) passed. Two findings were corrected before publishing the note: F4's Levit–Mandrescu escape threshold is `m = 21` under the manuscript's `ceil((2 alpha - 1)/3)` convention, not 22 (boundary touched exactly at `m = 20`), and F5's contiguous-top-block claim has a real exception (three-break family breaks at reflected index {2} alone for `3 <= m <= 7`, consistent with Cor. 4.14's reflected-index-1 deficit staying negative below `m = 8`). F2 was sharpened: both larger m-direction families fail synchronization already at `m = 2`. Files: `notes/partial_sync_obstruction_2026-08-11.md`, `verify_partial_sync_obstruction_20260811.py` (output path adjusted to `results/`), `results/partial_sync_obstruction_20260811.json`. Downloads originals left in place.
- **Three revision-plan notes stay local, per Brett (2026-08-11 ship).** `notes/main_v2_referee_revision_plan_2026-07-16.md`, `notes/poisson_binomial_external_review_revision_plan_2026-07-16.md`, `notes/poisson_binomial_review_board_revision_plan_2026-07-16.md` are internal deliberation on live/in-progress submissions; excluded from the public push. Revisit after the E-JC decision.
- **Absorption-margin kill-test recorded; family-growth disproof route closed (2026-08-11, later).** The probe (Sonnet subagent, replay-verified in the main session: byte-identical JSON regeneration, exact re-derivation of the closest-approach record, independent TG builder reproducing arXiv:2511.00334's exact break sets) found the absorption factor rho never below ~3.37 across 9,961 deficit indices in four families, with closest approaches at the smallest parameters and orders-of-magnitude growth thereafter. Decision: stop treating parameter growth in the known pattern families as a counterexample route; treat bounded absorption as a lemma target for the banded-sync + head-reserve proof lane. Evidence is finite-corpus, not a theorem. Literature sweep same day: nothing since the 2026-07-16 intake bears on tree unimodality; 2511.00334 confirmed head-anchored by two independent reads (main session + agy Claude Sonnet 4.6).
- **C2 sync facts recast as single-index laws plus a head strip (2026-08-11 evening).** Decision: pursue the C2 base facts through the proved assembly lemma (C-LC cited from Li-Li-Yang-Zhang arXiv:2501.04245 + V,U,W single-index forms) rather than direct two-index partial-synchronization arguments. Reason: V and W are exceptionless and adversarially rigid over 8,166 exact instances, U fails only at reflected depth <= 3, and the lemma (proved, notes/c2_single_index_laws_2026-08-11.md) shows these plus the cited spider theorem imply both sync relations off the head strip. Bottleneck-family head window proved by polynomial-in-h Sturm certificates. Termwise bilinear positivity checked and refuted (aggregate cancellation is essential), so proof attempts should target aggregate dominance of the two dirty cross terms, not generator-pairwise sync.
- **Universal C2 kill-test recorded: target survives to pendant weight 6,002 (2026-08-11 night).** Sonnet-built harness (validated against the known (1,1)x(9,9) margin-6 closed form before running), run twice — subagent and main-session replay — with identical verdicts: zero violations of B-LC, C ~_p B, xC ~_p B over 16.5k pairs (structured, random, hill-climbed); global bottleneck is exactly the proved h+2 top-corner margin of the (1,1)x(2h+1,2h+1) family through h = 1460, and nothing hill-climbing found within two orders of magnitude of it. Decision: treat the universal C2 base facts as the standing conjecture of record and direct proof effort at (in order) LAW-V/LAW-W, the depth<=3 U-strip, and a body argument for B's log-concavity (its margin minimum is interior, k=1269 at m=3000, so head certificates alone cannot close it). Banded/sparse tiers are exact where they look but not full-grid at every parameter; recorded as evidence, not proof.
- **Contraction band law adopted as the campaign's universal working conjecture (2026-08-11 night).** The single-index forms U, V, W generalize from the C2 arm-pair setting to arbitrary tree edge contraction (identical A/D four-state shape at the contracted edge), and after the first, stronger version (V and W clean everywhere) was killed within the hour by non-LC composites, the surviving statement is: all three forms nonnegative at reflected depth >= 3, every tree, every edge. Evidence: ~478k contraction pairs exact (exhaustive n<=16; random n<=18; KL non-LC composites), no violation at depth >= 3; verifier runs in 13 s. Decision: this, not the C2-specific laws, is the statement to freeze for proof dispatch; the C2 clean laws remain the special case needed for the two-hub theorem. Next adversarial step before any dispatch: composites with multi-break pattern-family components, which could deepen the strip.
- **Both universal contraction conjectures refuted; freeze nothing universal, keep the C2 laws (2026-08-11, late night).** The stated next adversarial step (multi-break pattern components) ran before any dispatch and killed the absolute band law: U-failures to reflected depth 23, V to 5, W to 4 on composites over the 1/2/3-break families, with failure depth tracking the components' own sync-band depth. This supersedes the same-evening decision to adopt the band law as the working conjecture. Standing state: the C2-restricted clean laws (path pieces) remain the proof target of record for the two-hub theorem; no universal contraction statement is frozen; the next session should measure how contraction-form failure depth scales with the component parameter (candidate normalizations: component band depth, Levit-Mandrescu zone) before conjecturing again. All four sweeps regenerate in ~25 s via verify_contraction_band_law_20260811.py.

## 2026-08-12

- **Outward survey run (light tier) plus Karlin read plus LAW-V packet built.** (1) The prefix-certification survey (hyperresearch light run, ship gate passed; report in the gitignored research vault) found: the prefix-plus-tail reduction is PUBLISHED as Levit-Kadrawi arXiv:2603.17114 Lemma 2.19 / Corollary 2.20 (a paper already in notes/literature under the unicyclic-counterexamples title); the Chan-Pak atlas has never been index-truncated (evidenced negative; its engine is the matroid exchange property, so truncation for trees stacks two open problems); Lorentzian cone/truncation variants restrict domains, not coefficient indices; Bendjeddou-Hardiman propose the unimodality-direct weakening as future work; Bjorner 1994 is the older partial-unimodality lineage. Decision: cite LK Lemma 2.19/Cor 2.20 as the frame for all prefix-LC work; do not pursue atlas truncation ahead of the Karlin/TP and head-window lanes (a cheap local-conditions kill-test is the reopening trigger). (2) HWZZ MIA 2017 fetched (ele-math open access) into notes/literature; exact closure statements now quotable. (3) TP2 toolkit note written (notes/karlin_lr_order_toolkit_2026-08-12.md): LAW-V is consecutive-TP2 of [C; xB]; common-LC-convolution closure proved via Cauchy-Binet and verified on thousands of exact instances; consequence recorded that dirty-term failure geometry reduces to single-hub pairs, making the single-hub case load-bearing. One memory-only item (Shaked-Shanthikumar mixture statement) explicitly marked UNVERIFIED. (4) Frozen proof packet built at gpt_attack/law_v_packet_2026-08-12/ (graded targets G0-G6, refutation-as-success, sourced toolkit, 21-check grounded self-test passing). Not dispatched; awaiting Brett's choice of vehicle per the 2026-07-23 allowance rule.
- **Atlas truncation closed as a lane; reopening trigger fired negative (2026-08-12, follow-up).** The cheap kill-test named in the morning's survey ran: the Chan-Pak atlas's entire direct verification burden sits at sink vertices carrying the bordered compatible-extensions matrix A(empty,1), whose required one-positive-eigenvalue property is, per the paper's own Lemma 4.6, equivalent at the bottom level to transitivity of the conflict relation — the matroid exchange property. For graph independence complexes the conflict relation is adjacency; exact computation (harness validated on two matroid controls) shows every tree on 3..9 vertices fails (93/93, with (Hyp) <=> cluster conflict graph exactly across the sweep), paths failing WORSE (3 positive eigenvalues) than claws or spiders (2). The failing vertex lies in the atlas of every certifying index k, so index truncation removes nothing. Verifier: verify_atlas_bottom_level_20260812.py (also records the sympy count_roots endpoint-inclusivity gotcha: strip zero roots before counting positives on rank-deficient matrices). Note: notes/atlas_truncation_killtest_2026-08-12.md. Proof effort stays on the ordered queue: LAW-V/LAW-W packet, U head strip, B body log-concavity.
- **Three comments posted to erdosproblems.com/993 (2026-08-12).** Comment 1: literature block (Levit-Mandrescu 2006 tail theorem, LLYZ spiders, the Galvin -> Bautista-Ramos multiple-breaks chain, Li's unimodality of the KL families). Comment 2: progress update to the 12 Mar comment (preprint DOI, the C2 four-base-facts reduction, weight-24 exhaustive + weight-6,002 adversarial evidence, h+2 bottleneck; AI disclosure included). Comment 3: follow-up literature (Levit-Kadrawi Lemma 2.19 / Cor 2.20 prefix-plus-tail reduction; Bendjeddou-Hardiman with their unimodality-direct future-work quote; Bjorner 1994). Drafts with source-verification notes: outreach/erdosproblems_993_comments_2026-08-12.md. Deliberately NOT posted: any campaign-internal results since Comment 2 (single-index laws, contraction kills, atlas kill-test); next self-update waits for a landed theorem or the E-JC decision.
- **Beurling-majorant lead opened and killed same day; lane closed (2026-08-12, evening).** A new lead (notes/beurling_majorant_difference_test_lead_2026-08-12.md, prompted by the toolbox of Cohen's 2025 Annals FUP paper, filed to central literature/) proposed testing window coefficient differences against one-sided band-limited majorants instead of local-limit estimates. Its own toy kill-test (scratch_beurling_band_killtest_20260812.py; results/beurling_band_killtest_20260812.json; exact coefficients, exact rational saddle bisection, Arb-certified roots at max ball radius 1.9e-26, P_101 real-rooted control passing) killed the pure lane structurally: Selberg interval slack ~2pi/w in count units means isolating one coefficient needs w > 2pi, while saddle-circle test kernels are torus-capped at w <= pi ~-- so a width-1 valley is invisible on EVERY tree, zero geometry aside; the real-rooted control saturates the bound at 2-coefficient resolution. Hybrid residual (interval pinches at measured Delta_k 2.2-4.7 on the pareto champions and JP phases, n <= 189, plus an undesigned local step) inherits the Jerrum--Patel degradation (w -> 0) and is not pursued. Proof queue unchanged: LAW-V/LAW-W packet, U head strip, B body log-concavity. Reopening bar recorded in the note: beat the 2pi/w slack or escape the torus cap.

## 2026-08-14

- **Rota flat-count disproof (arXiv:2608.07342) read and transfer-probed; recipe anatomy adopted as the campaign's frame.** Divoux-Larson-Lowen-Wang refute Rota's unimodality conjecture via severe LC failure (theta graphs, unbounded violation ratio) + direct-sum replication + Whittle q-lift (the in-category geometric tilt; LC is exactly tilt-invariant unimodality). Exact probe (scripts/probe_matroid_recipe_transfer_20260814.py) showed trees already have steps 1-2 — the Bautista-Ramos TG family carries unbounded violation ratios (2.4e8 in-grid) and m replicated breaks — but every break sits in the Levit-Kadrawi-immunized tail, structurally so within the family (depth 2m-1 vs threshold depth ~(t+1)m), and trees lack any tilt: their mode-raising operations are all binomial smoothers (the D17/D20/D23 content). Decision: treat break DEPTH, not severity or count, as the binding resource for any #993 counterexample; note the cross-category rhyme that the authors' method also cannot touch points-lines-planes (the bottom of the rank sequence — the analogue of the tree prefix). Paper filed to notes/literature/arxiv_2608_07342.{pdf,txt}; analysis in notes/matroid_nonunimodality_2026-08-14.md. No manuscript change; erdosproblems comment drafted as a follow-up option only.
- **Census comment posted to erdosproblems.com/993 (2026-08-14, Brett).** One comment, census-update style per the 12 Aug conventions: Rota-disproof motivation (arXiv:2608.07342), the complete n <= 32 LC census with failure counts 2, 0, 19, 7, 121, 159, 922 and the single-geometry finding (every failure a lone break at k = alpha-1), unimodality record left attributed to the thread, the k - q window arithmetic pinning any distance-3 break to depth >= 3, alpha in {17,18,19}, n in [33,38] (framed as an open bounded target), AI-disclosure line included. Draft with per-claim verification appendix: outreach/erdosproblems_993_comment_2026-08-14.md (gitignored, local). Decision owner: Brett (posted personally). Supporting data public at commit 64b8b9f.
- **Full LC census of all trees n <= 32 run same day (Brett's direction, after calling out the stale "day of C work" estimate); depth-2 window closed unconditionally.** New C engine (scripts/lc_census.c: exact uint64 full-polynomial DP fused to gentreeg parent arrays, 1.38M trees/s/core; validated by byte-matching all 106 n=10 polys against the Python DP, the 2 stored n=26 failures, 0 at n=27, and the 19 stored n=28 failures by exact polynomial multiset; every rung's total equals an independent gentreeg -u count). Result over 173,378,186,130 trees (n=27..32, 7.5 h laptop wall): zero non-unimodal sequences (independent re-verification of unimodality through n=32); 1,228 LC failures; EVERY failure is a lone depth-1 break at k = alpha-1 with alpha in {ceil(n/2)+1..+3}; failure counts by order 26..32: 2, 0, 19, 7, 121, 159, 922 (the n>=29 values believed new). Decision consequences: dist >= 4 is now unconditional for n <= 32 by complete enumeration; combined with the window arithmetic, no tree of ANY size can reach dist 3 via a depth-2 break; the dist-3 frontier is depth >= 3, alpha in {17,18,19}, n in [33,38], best attacked via a matching-constrained enumeration (alpha <= 19 means matching number >= n-19). Artifacts: scripts/lc_census.c, scripts/run_lc_census_20260814.py, scripts/run_census_ladder_20260814.sh, results/lc_census_20260814/ (shards + per-order summaries). The n>=29 failure counts are candidate erdosproblems.com/OEIS material, Brett's call.
- **Break-depth lane opened and first campaign run (Brett's direction); dist >= 4 wall recorded; dist = 3 reduced to a finite window.** New objective, never previously optimized: minimize dist = k_break - ceil((2*alpha-1)/3). Results (all exact arithmetic, scripts/search_break_depth_20260814.py + three results JSONs): every known break has dist = 4 (21 exhaustive n<=28 failure trees recomputed; 982-tree structured sweep over TG/KL-pattern/mixed-ks grammars; complete 1- and 2-mutation balls, 10.6k distinct polys; 58,476,064-eval six-island evolutionary run, zero unimodality alarms). Exact window arithmetic: a depth-d break has dist = 3 iff alpha in {3d+8, 3d+9, 3d+10}, so with n <= 28 exhausted the first possible dist-3 tree is a depth-2 break at alpha in {15,16}, n in [29,32] — a bounded Modal-shaped enumeration (n=29 unimodality sweep precedent: 5.47e9 trees, 1024 workers, 2026-03) with a cheap alpha-prefilter. Smallest alpha carrying any depth-2 break in the whole corpus: 19 (window needs <= 16). Decision: record "every tree LC break has dist >= 4" as an adversarially-tested working observation (evidence, not theorem); the named next steps are the n in [29,32] alpha-filtered enumeration and, if wanted, a diversity-forced second evolutionary campaign (first campaign's deepen islands all collapsed into the n=35 alpha=19 basin). Note: notes/break_depth_lane_2026-08-14.md.

## 2026-08-19

- **Depth-3 window sized exactly by counting, not enumerating; the dist-3 enumeration ladder provably ends at rung 3.** The named next lane decision (consequence 4, notes/break_depth_lane_2026-08-14.md) was "sizing the matching-constrained enumeration". Done without generating a tree: trees are bipartite, so alpha = n - mu (Konig), and an Otter-dissymmetry DP with bivariate weight x^n y^mu and the classical root-exposable state counts unlabeled trees by matching number exactly (scripts/count_trees_by_matching_20260819.py; validated against all 20 grounded totals incl. the six census rungs, plus exhaustive brute-force (n,mu) histograms n=4..14 with per-tree Konig assertion). Results: rung 3 (alpha in {17,18,19}, n in [33,38]) = 531,192,598,181 trees inside 74.2T generated (7.2e-3; filter worthless at n=33 at 30% pass, decisive at n=38 at 3e-4); rung 4 (alpha in {20,21,22}, n in [33,44]) = 129.65T trees inside 34.06 QUADRILLION — 244x/459x growth per rung, so rung 3 is the last enumeration-reachable point of the dist-3 question, and beyond it only proof. If rung 3 comes back empty, dist >= 4 becomes unconditional for EVERY tree with alpha <= 19 (quotable, erdosproblems-comment-shaped). Filter engine built and validated (scripts/lc_census_alpha.c: n <= 38, O(n) greedy-matching prefilter, per-tree greedy-vs-DP alpha self-check aborting on mismatch, PASS_ALPHA shard reconciliation against the exact counts, __uint128_t Turan products since coefficients reach C(38,19) ~ 2^34; full-enumeration per-alpha counts match the DP exactly at n=10, 20, 26). Measured pipeline rates 1.96M trees/s (n=33) to 4.2M trees/s (n=38, producer-bound) price the rung-3 run at ~5,300 pipeline-hours as-is: not a laptop job. Decision: run NOT launched — recommendation recorded as fuse-the-filter-into-gentreeg first (heavy orders are producer-bound; nauty plugin hook, expected ~8x, bringing the job to ~700-800 core-hours = 3-4 laptop-days or ~1 Modal-hour at 1024 workers), with Modal-as-is the fallback; Modal spend is Brett's call. Note: notes/depth3_window_sizing_2026-08-19.md; results/depth3_window_sizing_20260819.json. Decision owner: Brett (launch gate); session was the 4pm Fable-credit burn directive, lane chosen autonomously as the named next step.
- **Sizing cross-validated against OEIS after the fact; window counts turn out to be A339829 territory.** The trees-by-independence-number triangle is OEIS A339829 (b-file through n = 50); the perfect-matching diagonal is A005751. All 990 b-file terms with n <= 44 match the Otter DP exactly — full independent external validation by a different method. Recorded candidly as the inward-search lesson, OEIS edition: the window counts could have been read off A339829 directly; the DP still supplies independent verification, the mu-state machinery the filter engine reuses, and the rung analysis. No OEIS contribution warranted below n = 51.
- **Session shipped as 0238bf8 (Brett's "ship"), payload constraints passed in the /ship launch args per the 2026-08-14 near-miss lesson.** Exactly the seven session files staged; the three 2026-07-16 revision-plan notes and notes/source-hooks/ deliberately left untracked. Verified post-push: HEAD == origin/master, exclusions still untracked.

## 2026-08-21

- **Replace the four queued August Erdős Problems comments with one concise theorem update.** The moderator's 2 July rejection of an earlier comment as incomprehensible showed that the 12 and 14 August literature, conditional-reduction, and census comments were too dense for the forum even where mathematically relevant. Brett deleted all four while they were still awaiting approval. He then submitted one replacement comment stating only the new structural result—every tree with at most two degree-at-least-three vertices has a log-concave independence polynomial, so any counterexample to #993 must have at least three such vertices—together with the Li--Li--Yang--Zhang spider context and the commit-pinned formal result note at `95f86d96dd89e5ddfff16b65f500fa9c85cb661d`. The replacement is awaiting moderator approval. Future forum updates should lead with one unconditional advance and omit campaign notation, proof-search history, and census detail unless directly requested.

## 2026-08-26

- **Schweitzer arXiv:2608.23262 adopted as the citation for the Lorentzian-route boundary; two-sentence paragraph added to the main_v2 intro.** The arXiv-watch harvest (triage: notes/arxiv-watch-triage-2026-08-26.md) surfaced Schweitzer's note proving Mason-type ULC for the common independent sets of any two matroids and killing it at three partition matroids. Since graph independence systems are intersections of one partition matroid per edge, this turns "the Lorentzian lane doesn't reach trees" from examples into a theorem boundary, and structurally explains why Bendjeddou-Hardiman needed a special subclass. Pre-use checks all run before citing (source hook, notes/source-hooks/schweitzer-2026-matroid-intersection-ulc.md): derivation chain replayed by hand against AMCD+26 v2; the three-matroid witness, its split-graph reading, and the per-edge matroid reduction replayed in exact arithmetic (scripts/verify_schweitzer_witness_20260826.py, all pass); still v1, no replacement, as of 14:35 UTC. The added paragraph also states the derived remark that the order-26 LC failures certify trees whose independence systems are not two-matroid intersections (ULC implies LC for downward-closed families). Still unrefereed: re-run the forward-check before any submission carrying the citation. Companion arXiv:2601.02547 skimmed same day: matroid-side (ultrametric distance matrices, not tree independence), filed to literature/, not cited.
- **arxiv_watch.py hardened after the first wholesale failure (06:15, all 21 queries 429/503, day silently unscreened).** Retry-with-backoff per request (15 s / 60 s, Retry-After honored to a 120 s cap, retryable 429/500/502/503) and, the part the morning actually needed, loud wholesale failure: exit 1, stderr warning, macOS notification under --notify, instead of "nothing new". Partial failure warns coverage is partial. Mocked tests for every retry path plus a clean live --replay run. Portfolio-level file (Project-Management/tools/), logged here because the failure and repair belong to this project's watch.
- **Matroid-intersection-number lane opened; subcubic LC conjecture stated and supported to n=33; root-sector route killed (2026-08-26 afternoon, Fable-credit burn directive).** Sharpening the morning's Schweitzer citation gave a lemma: the least number of matroids on V whose common independent sets are IS(G) equals the least number of disjoint-clique subgraphs covering E(G), so chi'(G) for triangle-free G and exactly Delta(T) for trees (lower bound holds for arbitrary matroids, not just partition matroids; constructive half verified on all 200 trees n=2..10). Consequence: Schweitzer's two-matroid ULC theorem covers exactly the paths, and subcubic trees are the three-matroid rung where his counterexample lives. A scan of all 1,230 stored LC failures pins three tight frontiers: every failure has Delta >= 4 (81 at exactly 4), >= 4 branch vertices, and >= 8 leaves. New computation: complete subcubic census at n=33 (2,076,217,086 trees, zero LC failures, total matching OEIS A000672, offset pre-validated at n=10/15/20/24), behind three gates (stored n=26 failures reproduced as a positive control; subcubic n=26 and n=28 complete and clean); n=34..36 left running, restartable. Route killed the same session: the pi/3 root-sector lemma cannot prove the subcubic conjecture, because 269 of 20,512 subcubic trees swept at n=14..20 have Arb-certified deviation > pi/3 while every one is log-concave, and complete binary trees reach 1.815 rad at n=127 still LC (scripts/sector_angle_probe_20260826.py; python-flint venv created at venv/, per the mandatory root-finding rule). Reusable by-product: binary caterpillars are real-rooted to n=90 and binary spiders sit near 0.31 rad, so path-like subcubic trees are inside the sector and bushy ones far outside, both LC — any future route must be insensitive to that split. Note: notes/matroid_intersection_number_lane_2026-08-26.md. NOT done, and explicitly not claimable: the prior-art search for the representation lemma (probably the classical equivalence covering number eq(G), marked UNVERIFIED from memory) never ran, so no novelty may be claimed; and the three-hub thread produced no artifact. Manuscript unchanged by this lane; the sentence added this morning keeps the safe crude per-edge count.
- **Two Fable subagent forks died mid-flight on a usage limit, writing nothing (2026-08-26).** An adversarial critic (representation-lemma attack, prior-art search, subcubic-graph LC search, manuscript referee-proofing) and a constructive proof-target developer (three-hub obstruction map, subcubic mechanisms, leaf-bound topological-type census) both terminated on "reached your Fable 5 limit" before writing their output files. Their interim one-line summaries are visible in the transcript and are NOT evidence: in particular the claim that the C2 single-index laws replicate at three hubs has no artifact and must be re-derived from scratch if wanted. Recorded so a later session does not mistake a dead agent's teaser for a result. Operational lesson, consistent with the 2026-07-17 18-agent loss: a wide fan-out launched near a limit loses everything, including work already done, because output files are written at the end.
- **First subcubic log-concavity failure found, at n=34; the session's own subcubic conjecture refuted hours after being stated, and the matroid frame confirmed (2026-08-26, late afternoon).** The subcubic census commissioned by the matroid lane returned exactly one LC failure among 4,790,669,518 trees at n=34 (total matching OEIS A000672; n=33's 2,076,217,086 trees were clean). Independently replayed in exact arithmetic from the parent array (tree rebuilt, polynomial recomputed by a separate DP, coefficients matching the census engine exactly): Delta = 3, 7 branch vertices, 9 leaves, diameter 10, alpha = 18, lone break at k = 17 where i_17^2 = 5184 < i_16 i_18 = 5325 (deficit 141), and STILL UNIMODAL, so not a #993 counterexample. Consequences: (1) "every tree with Delta <= 3 is LC", stated in notes/matroid_intersection_number_lane_2026-08-26.md §2 that afternoon, is false; (2) the Delta >= 4 frontier over the 1,230 previously known failures is dead, while the two hub-counting frontiers (>= 4 branch vertices, >= 8 leaves) survive this witness; (3) the matroid frame is CONFIRMED rather than damaged, because the tree LC-failure threshold and the matroid-representation threshold now coincide exactly at three matroids, so the count does predictive work instead of merely redescribing known facts; (4) the break has dist 5, the first observed at that depth against dist 4 for all 1,230 earlier failures. Recorded caution: an exhaustive verification one rung short of the answer (n=33) reads exactly like a theorem, and the sector probe's earlier finding that subcubic trees violate the pi/3 sector freely was evidence against the conjecture hours before the census killed it, misread at the time as a fact about the method rather than about the class.
- **Matroid-representation lower bound sent to Aristotle (2026-08-26).** Packet formalization/matroid_representation_lower_bound_aristotle_input_20260826.md, project 546ab2b8-189c-4a4b-8316-3bbe6bfac967, task d8843604-e4ab-4ac2-ad3c-3f7b706446cd. Target G0: if matroids M_1..M_k on V(G) have IS(G) as their common independent sets, then k >= deg(v) at every vertex whose neighbourhood is independent; G1 the triangle-free/tree corollary k >= Delta; G2 (adversarial) whether ground-set enlargement can beat the bound. Scope chosen against the actual Mathlib surface, checked in the repo checkout rather than from memory: the matroid library (Basic/Circuit/Loop/Closure/Minor) and SimpleGraph.IsIndepSet exist, but there is no parallel-class API, no partition matroids, and no edge colouring, chromatic index, Koenig or Vizing anywhere in Combinatorics. So the packet deliberately excludes the upper bound (needs partition-matroid infrastructure, and is the half already verified computationally on all 200 trees n <= 10) and the chi' route to parts (b)/(c) (would mean formalizing Koenig, expensive and classical). The lower bound is the half with no computational check behind it. Unaffected by the n=34 refutation, which concerns log-concavity, not representation. Replay locally and audit for sorry/axioms before recording any verdict.
- **Independent attack on the representation lemma returned one usable result (2026-08-26).** Packet gpt_attack/matroid_representation_lemma_2026-08-26.md run through ocx against the Z.AI/Zhipu GLM-series preview initially released as "Ox Alpha" (provenance confirmed by Z.AI Co. per Bloomberg, shared-memory entry b6490c13, 2026-08-26; the route label opencode/x-preview-f-free and the name "Ox Alpha" are kept verbatim in the packet and route manifest as historical labels). Verdict T1 HOLDS-AS-STATED, with a clean reduction on the ground-set worry: auxiliary elements plus restriction cannot beat the bound, because IS(G) = intersection of (M_i restricted to V) and each restriction is a matroid on V. Verified concretely here on partition matroids with auxiliary elements before acceptance, and folded into the Aristotle packet's G2 along with the still-open contraction case. The no-fabricated-citations rule in the packet held: it returned concept names with confidence labels and 12 search strings, zero invented references, which is the intended way to use a model weak at exact citations. Its T3 independently identified the subcubic-triangle-free question as the fork in the road and proposed searching to n ~ 16; the census answered it at n = 34 while the model was writing. Prior-art search for the lemma STILL NOT RUN; no novelty claimable.
- **Second-round external return: one theorem kept with a corrected proof, one claim refuted (2026-08-26).** The same GLM-series preview route, unprompted, closed the ground-set/contraction gap left open in its first pass: mi*(G) = mi(G) even allowing auxiliary elements and contraction. The conclusion is correct and adopted, but its proof is not sound (it argues contraction can add dependences, then asserts the restriction reproduces the same dependence pattern); the correct one-line proof, already in the Aristotle packet's G2, is that a minor of a matroid is a matroid, so ((M_i)/F)|_V is a matroid on V. Three structure facts it offered were verified by brute-force eqcov and kept: mi of a disjoint union is the max; mi = 1 iff the graph is a disjoint union of cliques; mi <= chi' always, with equality when triangle-free and a large gap otherwise (mi(K_4) = 1 vs chi' = 3). One claim refuted: that the triangle-free identity covers all odd cycles giving mi(C_2m+1) = Delta+1 — false at C_3, which is itself a triangle, so mi(C_3) = 1; the claim holds from C_5 on. Its sweep code does not parse and its design imposes min-degree 2, excluding trees, so it could not have found the n=34 witness. Standing lesson reinforced: this route generates strong candidates and unsound derivations at the same confidence, so conclusions are worth checking and proofs are worth rewriting.
- **Matroid-representation lower bound formally proved and locally replayed (2026-08-26).** Aristotle returned project 546ab2b8-189c-4a4b-8316-3bbe6bfac967 complete on all three graded targets. Verified here before recording: escape-hatch grep empty (no sorry/admit/axiom/native_decide/implemented_by), statements read and compared declaration-by-declaration against the packet, lake build clean at 8,027 jobs, and #print axioms reporting only propext, Classical.choice, Quot.sound across all 13 declarations. So the lower-bound half of the representation lemma is now machine-checked: deg v <= k at any vertex with independent neighbourhood, hence k >= Delta for triangle-free graphs and for trees. G2 came back STRONGER than requested: instead of just the restriction case it proves the bound for representations by arbitrary minors ((P i).contract (C i)).restrict (R i) via Matroid.comap along an embedding, which independently settles the contraction question the external route had reached by an unsound argument. Aristotle also improved the statement by dropping the ground-set hypothesis, having observed that the representation condition on singletons already forces every vertex into every (M i).E. Remaining unproved: the upper bound (deliberately out of scope, needs partition-matroid infrastructure absent from Mathlib) and the chi'=Delta step (would need Koenig). Prior-art search STILL NOT RUN, so novelty remains unclaimable.
- **Subcubic census n=35 clean; the n=34 failure was isolated at that checkpoint (2026-08-26).** 11,077,270,335 subcubic trees at n=35, zero LC failures, total matching OEIS A000672. With n=33 clean (2.08e9), n=34 carrying exactly one failure, and n=35 clean, subcubic LC failures are extremely sparse rather than absent, sharpening the open question from the refuted boolean conjecture to a density one. At this checkpoint n=36 was still running; the later n=36 entry below supersedes the provisional density discussion.
- **Prior-art search finally run on the representation lemma; two of three components are classical (2026-08-26).** The invariant I called eqcov is Alon's equivalence number eq(G) (Combinatorica 1986, DOI 10.1007/BF02579381), and the matroid count is the standard mu(I), with 'every independence system is a finite intersection of matroids' textbook. Rename accordingly: mi -> mu, eqcov -> eq. No source found for the identity mu(IS(G)) = eq(G) itself, which is the only part that could be a contribution; status is 'searched, not found', NOT 'novel'. The suggestive finding: mu for matching systems (arXiv:math/0212235) and eq of line graphs (Esperet-Gimbel-King, arXiv:1006.3692) are, under the lemma, the same quantity, since matchings of G are the independent sets of L(G) — and a full-text fetch of math/0212235 shows it cites neither Alon nor equivalence covers, so the two literatures appear not to know they overlap. Their numbers are mutually consistent (mu(15)=5 against the EGK upper bound ~5.9 at chi=15), weak evidence the identity holds and is unnoticed. Five search strings recorded in the lane note; zbMATH, MathSciNet and a forward-citation sweep of Alon 1986 NOT run, so the novelty question stays open.
- **The n=34 subcubic witness motivated a three-fold symmetric slice (2026-08-26).** Vertex 1 is a centre with three pairwise-isomorphic 11-vertex arms (automorphism group of order 48, arm diameter 8, two branch vertices per arm). That structure justified a cheap targeted sweep, not a claim that symmetry is enriched or causal. The sweep and its later unequal-arm correction are recorded below.
- **External route retracted its own contraction proof under challenge; independently matched the defect already recorded (2026-08-26).** It supplied the counterexample to its own reasoning (a triangle circuit {e,u,v} leaves {u,v} independent in M but parallel in M/e, so contraction manufactures dependences) and replaced it with the uniform argument that a 2-set is dependent only by a loop or parallelism, parallelism being transitive among nonloops. That defect was already recorded here hours earlier from an independent read, which is the result that justifies the checking discipline. Status unchanged either way: the question was closed by machine proof (le_card_of_represents_minor), and the Lean development already takes the loop-or-parallel route. Kept as the campaign's cleanest instance of a right conclusion carried on a wrong proof.
- **Symmetric subcubic slice swept: 11,892 observed LC failures through n=67, all irrelevant to #993 in the tested range (2026-08-26 evening).** The exhaustive three-fold symmetric sweep (unordered rooted binary arms; $I(T)=(F+G)^3+xF^3$; exact integers; positive control re-finds the n=34 witness) finds 11,892 failures and zero unimodality alarms. Every observed break is at $k=\alpha-1$ and inside the Levit--Mandrescu decreasing zone; depths range from 5 to 12. This is finite negative evidence and a reason to deprioritize the slice for shallow-break searches, not a proof that the infinite family is structurally incapable of producing a #993 counterexample. The two-fold symmetric slice has zero failures through the same arm-size bound. Within the tested three-fold slice, 0 of 135,451 gap-0 arms fail; the named candidate extremizer $\mathrm{unit}^j(P_2)$ remains clean through $n=358$, but neither the extremality claim beyond $m=20$ nor the gap condition is a theorem. The completed n=36 census below shows that the gap condition fails immediately for unequal arms. Notes and exact command: notes/symmetric_subcubic_family_2026-08-26.md; `python3 scripts/symmetric_subcubic_sweep_20260826.py 22`.
- **A code bug briefly reported three #993 counterexamples (2026-08-26).** The first 2-fold construction multiplied the subtracted term by x and produced three NON-UNIMODAL alarms at n=24. They were artifacts: the polynomials ended '...,174,0,1', an internal zero, impossible for an independence polynomial since deleting a vertex from an independent set leaves one. The sweep scripts now assert i_0=1, i_1=n, nonnegativity and no internal zeros on EVERY computed polynomial and treat a violation as a code error. Standing lesson: a search whose success condition is 'find an anomaly' will report its own bugs as discoveries, so the cheap structural invariants of the computed object must be asserted per instance, not spot-checked.
- **Symmetry-enrichment lead opened and KILLED by its own control (2026-08-26 evening).** All 21 known LC failures at n=26 and n=28 have nontrivial automorphism groups (orders 48 to 82,944), none asymmetric, which looked like it would generalize the symmetric approach past subcubic. The control refutes it: 94.74%, 95.60%, 96.55% and 97.38% of ALL trees at n=12, 14, 16, 18 have a nontrivial automorphism, rising with n. At a ~97% base rate, 21/21 symmetric has probability ~0.53 and carries no information. 'LC failures are symmetric' is a fact about trees, not about failures. What survives: the n=34 witness's specific three-isomorphic-arm structure still made the productive family sweep possible, and whether failures have unusually LARGE automorphism groups is a separate untested question. Method note: the first control used GraphMatcher.isomorphisms_iter(), which enumerates the group (a star has (n-1)! automorphisms) and hung until killed; deciding whether any nontrivial automorphism exists needs AHU canonical forms of sibling subtrees instead, which runs the whole n<=18 census in seconds.
- **Subcubic census completed at n=36: 15 extremely sparse failures, all independently replayed (2026-08-26).** The 32 shards total 25,664,800,714 trees, exactly A000672, with 15 LC failures and zero unimodality alarms; the sequence for n=33..36 is 0,1,0,15. A separate tree DP rebuilds all 15 parent arrays, matches every stored polynomial, finds no duplicate graphs, and confirms $\alpha=19$, one break at $k=18$, and depth 5. All are unicentral with a degree-3 centre; five have two rooted-arm shapes and ten have three, so the unequal three-arm class is a justified search target but not an established mechanism. Arm sizes range from 2 to 25, invalidating the proposed size-11--12 restriction. More decisively, six witnesses have arm-gap vector $(1,1,1)$ and nine have $(0,0,0)$, so the symmetric gap-0 condition does not survive unequal arms. The witness-seeded default perturbation tests 969 arm multisets/921 distinct trees, reproduces the exact 15-tree n=36 set, and finds 126 LC failures through n=76 with depths 5..10 and no unimodality alarms. The n=36 rate is $5.84\times10^{-10}$ (one per 1.71 billion), 2.8 times the n=34 rate but nowhere near enough data for a density trend. Replay: `python3 scripts/verify_subcubic_n36_failures_20260826.py`; seeded follow-up: `python3 scripts/witness_perturbation_20260826.py`.
- **Log-concavity convolution and the $\pi/3$ root-sector criterion formally proved and locally replayed (2026-08-26).** Aristotle project 23421e64-d858-4065-be32-71ffbc8f64be, task db1608bb-ed2f-4bbc-b353-aebec33f15fd, proves G0 convolution closure with exact no-internal-zeros hypotheses and a polynomial corollary; the linear/quadratic G1 statements including the sharp iff; the full optional G2 sector theorem; and an explicit counterexample showing why no-internal-zeros is necessary. Local `lake build` completed 8,034 jobs against pinned Mathlib v4.28.0; the escape-hatch scan is empty and all target axiom audits report only `propext`, `Classical.choice`, and `Quot.sound`. Imported at `formalization/lc_convolution_sector_aristotle/`.

## 2026-08-28

- **Adopt a validated local proof-frontier graph, using TheoremGraph as a candidate source rather than an authority.** The current depth-3 matching-bag lane is encoded in `proof_graph/erdos993_frontier.json`; its validator rejects incomplete bridge records, missing grounded evidence, and dependency cycles. Search packets must begin from a ready local claim, expand at most one dependency hop, and keep source statements separate from project adaptations.
- **Split the blocked correction into blocked-profile log-concavity and a mixed cross-reserve.** This is the exact decomposition of the remaining depth-3 inequality after the extendable Pascal reserve. An exhaustive probe of all 13,179 relevant non-isomorphic trees through `n <= 15` finds zero failures of either combined target and zero full depth-3 failures. This licenses the split as the next attack surface, not either inequality as a theorem; the negative mixed term in 223 trees also rules out proving its nonnegativity in isolation.
- **Prioritize ordered cover deletion / first violated comparison as the next bridge kill test.** It retains the ordered-forest code and directly targets the blocked mass. Minimal-unsatisfiable tree CSPs, dependency-graph/Bonferroni bounds, and nonpure simplicial faces stay in the queue as speculative representations. A bridge earns proof-development time only after its dictionary is complete, its source theorem survives adaptation, and a bounded exact test beats the current baseline. No manuscript edit or large enumeration was authorized by this decision.
- **Keep Vatter-style residual-tail descent as a capped, high-upside side bet; do not spend the banked reset on it now.** For a leaf-peeling residual `H_j`, the candidate `i_k(H_j) <= i_{k-1}(H_j)` for every `k >= d(T)` would rule out every post-descent rebound by the least-vertex coefficient decomposition. It survives every leaf-peeling order for all 13,188 trees through `n <= 15` (46,829,967 exact residual checks), 1,002 orders on each known `n=26` log-concavity witness, and 711 adversarial family cases through 180 vertices; the worst residual ratio through `n <= 15` is `357/415`. This is bounded evidence, not a theorem. The route bypasses rather than proves the tie-fugacity mean inequality, and its natural induction intersects the existing leaf-step descent, mode-alignment, and forest-product bottlenecks. Stop the exploratory block here; resume only for one narrowly scoped attack on the tight residual boundary case, preferably using ordinary tokens near reset. A banked-reset sprint requires an actual local lemma or a reduction to a standard inequality.

## 2026-09-04

- **Run the canonical first-violated-comparison partition as the next bounded matching-bag experiment.** The weekly arXiv/Astra audit strengthened the rationale without supplying a theorem. The exact `n<=15` corpus supports `C<0 => b_3^2 e_2 e_4 >= b_2 b_4 e_3^2` in all 223 adverse-cross cases; the normalized inequality fails 49 times outside that regime, so retain the implication rather than a universal strengthening. Partition blocked assignments disjointly by their first violated comparison/path, require exact reconstruction of `b_2,b_3,b_4`, test down-closure before importing Guo--Kang ballot/O-sequence machinery, and seek Pahari-style summandwise ratio comparisons before addition. Astra #126's direct CND/sign-kernel route is closed by the five-vertex-star kill test. Baldi--Kummer tropicalization is authorized only as a downstream discovery step after a fixed semialgebraic obstruction-profile model exists. This is an experiment and proof-strategy decision, not a manuscript or all-tree claim; full record: `notes/arxiv_astra_transfer_2026-09-04.md`.
- **Complete only the proportional first gate, then stop.** The canonical first-obstruction implementation exactly partitions all 4,788,857 blocked defect-2/3/4 assignments through `n<=15` and reconstructs the profiles on all 13,179 relevant trees. The 1,400,970 multiply certified assignments show that disjointification is nontrivial. Pair certificates range from path length 2 to 6, so any later structural pass should test a shortest-path-first or minimal-obstruction order. Runtime was 73.89 seconds wall time. No down-closure, class-ratio, proof-search, or tropicalization work was started; resumption is a separate decision. Exact certificate: `results/first_violation_partition_probe_20260904.json`.
- **After Brett resumed the experiment, keep the grouped summand route and close the direct individual-class down-set route.** With unary-first and shortest-path-first pair order, individual first-certificate classes fail one-bag-filling closure 694,084 times, first at `n=6`; cumulative prefixes have zero failures for the formal reason that filling preserves existing violations and may only add earlier ones. All 46,557 nonempty class profiles through `n<=15` are log-concave, and all 41,427 profiles spanning defects 2 and 4 satisfy `c_3^2 e_2 e_4 >= c_2 c_4 e_3^2`. Cross-class pair positivity is too strong (9,991 failures overall and one on an adverse tree), but natural unary/path-length grouping resolves the bounded adverse corpus: 220 of the 223 `C<0` trees are unary-only and three pair-only; all 224 group diagonals and the sole cross-group pair term are nonnegative. Prioritize explicit formulas or injections for those natural groups if this lane resumes. Do not invoke Guo--Kang for individual classes, and do not tropicalize without a fixed semialgebraic group-profile model. This is exact bounded evidence, not a theorem or manuscript claim. Certificate: `results/obstruction_class_structure_probe_20260904.json`.
- **Replace opaque unary certificate counts by an exact induced-forest formula, then stop at the unsupported inequality.** Let `U` be the vertices occurring in no maximum independent set and `J_d=[x^(alpha-d)]I(T-U)`. Then the unary and pair groups are exactly `u_d=s_d-J_d` and `p_d=J_d-e_d`. Ordering `U=(v_1,...,v_r)`, the `i`-th unary class is `[x^(alpha-d-1)]I(T-N[v_i]-{v_j:j<i})`. These identities are elementary set partitions; the exact replay checks 41,685 unary certificates and both coarse identities on all 13,179 trees with zero mismatches. The new frozen target is the relative depth-3 coefficient inequality for these specially arising residual forests. General forest log-concavity is unavailable, so do not record it as proved or spend more compute merely extending the census. Pair-only adverse cases remain a separate fixed-path calculation. No manuscript claim changed.
- **Close the conditional relative-ratio and grouped-summand route after the proportional `n<=17` kill test.** The optimized 46.98-second extension checks 81,128 relevant trees and finds the first exact tree with both `C<0` and `R<0`: graph6 ``PpG`A?@O??g??@O???I??O??``, `n=17`, `alpha=11`, `e_2..e_4=(584,1408,2241)`, `b_2..b_4=(1,8,44)`, `C=-5409`, and `R=-3468800`. The original target survives comfortably (`extendable reserve + C = 668311`; full depth-3 margin 668331). The witness's two unary classes each pass the individual ratio test, but their cross term defeats the grouped ratio, so neither classwise inequalities nor the natural grouping prove the aggregate statement. Retain the exact filtration and induced-forest formulas, but use them only for a direct cross-reserve attack. Also retain the numerical kill test `(e_2,e_3,e_4)=(10,20,20)`, `(b_2,b_3,b_4)=(1,10,100)`, which satisfies the proved extendable reserve and blocked LC while giving `C=-620`, `R=-20000`; tree-specific structure is indispensable. Exact certificate: `results/blocked_profile_depth3_probe_n17_20260904.json`. No manuscript claim changed.
- **Adopt the conditional blocked-density bound as the live direct cross-reserve bridge after the `n<=18` kill test.** The universal Pascal reserve shows it is sufficient, whenever `C<0`, to prove `b_2/e_2+b_4/e_4 <= (5 alpha+17)/(27(alpha-3))`; this discards the favorable `2e_3b_3` term entirely. The optimized scan checks 204,995 relevant trees and 1,539 adverse cases through `n<=18` in 159.83 seconds, with zero failures of this stronger bound, the actual cross-reserve, blocked-profile depth-three LC, or the full depth-three inequality. The worst observed adverse `-C` consumes `9/35` of the guaranteed Pascal reserve. This candidate retains the obstruction-path charging problem in a simpler weighted-density form and is the next proof target; it is not a theorem. No auxiliary failure found in this session is a refutation of Erdős #993, and no manuscript claim changed. Certificate: `results/blocked_profile_depth3_probe_n18_20260904.json`.
- **Stop untargeted random sampling in the actual finite-reduction rung.** The `n<=18` census does not cover the live `alpha=17..19`, `deficiency<=5`, `n=33..38` cells. A conditioned random port-labeled bag-tree probe accepted 50 trees in each of the 12 admissible cells (600 trees from 1,560 trials in 240.37 seconds), but none had `C<0`; consequently the run did not exercise the conditional blocked-density candidate. Zero blocked-LC, cross-reserve, full depth-three, and unimodality failures are only biased random evidence. Do not scale this sampler further. A future computational test must target negative `C` inside the rung, and any auxiliary failure remains distinct from a refutation of Erdős #993. Certificate: `results/cross_reserve_target_scope_probe_20260904.json`.
- **Retire blocked-profile log-concavity and replace the old assembly route by a joint blocked-plus-cross charge.** Eight pendant matched-pair lifts take the exact `n=17`, `alpha=11`, deficiency-5 negative-cross witness into the live `(n,alpha,deficiency)=(33,19,5)` cell. A 56.66-second targeted search evaluates 500 nonisomorphic lifts, including 57 with `C<0`. The conditional cross-endpoint bound survives all 57, with worst load 0.0250234 of the Pascal reserve, but nine lifts refute `b_3^2>=b_2b_4`. The first has `b_2..b_4=(1,16,4216)` and blocked margin `-3960`; its full margin is still positive at 462,192,427,400, and independent graph6 replay is exact. The new sufficient target, needed only when `D=(b_3^2-b_2b_4)+C<0`, is `e_2b_4+e_4b_2+b_2b_4 <= ((5alpha+17)/(27(alpha-3)))e_2e_4`; it survives all 57 adverse combined-correction lifts with worst load 0.0250394. It fails outside the adverse regime and is not a universal inequality. Update the proof frontier accordingly. This auxiliary refutation is not a refutation of Erdős #993; no non-unimodal tree or negative full depth-three margin was found. Certificate: `results/cross_reserve_witness_lifts_20260904.json`.
