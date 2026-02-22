# Proof Closure Checklist (2026-02-21)

## Scope
This file is the single closure tracker for the current symbolic endgame:
- Route 1 (`mu_P(lambda_m(T)) >= m-2`)
- STRONG C2 (`lc_surplus + mismatch >= 0`)
- Final bridge to Conjecture A closure language

It intentionally avoids adding new scanners that duplicate existing verifiers.

## Existing Resources (Do Not Rebuild)
- Triage wrapper (orchestration only):
  - `scripts/closure_triage.py`
- Route 1 direct verifier + extrema:
  - `verify_muP_sum_of_means_2026_02_19.py`
  - `results/verify_muP_sum_of_means_n23_2026_02_19.json`
- STRONG C2 split verifier + mismatch witness list:
  - `verify_strong_c2_rise_identity_2026_02_19.py`
  - `results/verify_strong_c2_rise_identity_2026_02_19.json`
- Route-B pairwise decomposition minima:
  - `verify_strong_c2_route_b_pair_bounds_2026_02_19.py`
  - `results/verify_strong_c2_route_b_pair_bounds_2026_02_19_n23.json`
- Unified STRONG C2 certificate:
  - `prove_strong_c2_unified.py`
  - `results/prove_strong_c2_unified_n23.json`
- Tight witness packs:
  - `results/route1_transfer_failures_all4_2026_02_19.json`
  - `results/analyze_strong_c2_qdrop_witnesses_2026_02_19.json`

## Remaining Symbolic Gaps (Dependency Order)
1. Route 1 universal symbolic bound:
   - Prove `exact_slack_B >= exact_excess_D` for all sizes
   - Equivalent target: `sum_c mu_{T_c}(lambda_m(T)) >= m-2`
   - Reference note: `notes/path1_direct_muP_proof_2026-02-19.md`

2. STRONG C2 universal symbolic bound:
   - Prove `lc_surplus + mismatch >= 0`
   - Hard regime reduction to close: `(-dq/db) <= p1/b1` when `dq < 0`
   - Route-B alternative: structural proof of `R = cross + mismatch >= 0`
   - Reference notes:
     - `notes/strong_c2_algebraic_proof_2026-02-19.md`
     - `notes/strong_c2_route_b_proof_2026-02-19.md`
     - `notes/strong_c2_route_c_proof_2026-02-19.md`

3. Route-C fallback lemma (optional but useful):
   - Prove `mode(P) >= m-1` (or `p_{m-1} >= q_{m-1}`) structurally
   - This is a fallback closure path if Route-B remains stubborn

## Triage Protocol for New Lemmas (from 5.2 or manual work)
1. Classify lemma as one of:
   - `route1_exact_slack`, `strong_c2_hard_ratio`, `routeB_R_nonneg`, `routeC_modeP`
2. Test on existing full frontier first (no new custom scanner):
   - Route 1: `verify_muP_sum_of_means_2026_02_19.py`
   - STRONG C2 split/hard regime: `verify_strong_c2_rise_identity_2026_02_19.py`
   - Route-B decomposition: `verify_strong_c2_route_b_pair_bounds_2026_02_19.py`
3. If failure appears, immediately map to known tight witnesses:
   - Route 1: `results/route1_transfer_failures_all4_2026_02_19.json`
   - STRONG C2 Q-drop: `results/analyze_strong_c2_qdrop_witnesses_2026_02_19.json`
4. Promote lemma only after both:
   - symbolic derivation in notes
   - zero-failure computational check on the canonical frontier

## Exit Criteria
Closure is ready for manuscript integration when all are true:
1. A symbolic proof note exists for Route 1 universal bound.
2. A symbolic proof note exists for STRONG C2 universal bound.
3. The final notes explicitly show how these imply the Conjecture A closure step used in `paper/main_v2.tex`.
4. Existing verification scripts still report zero failures on the same frontier.

## Guardrails (Known Dead Ends)
- Do not rely on `mode(P') >= m-1` (false).
- Do not rely on TP2 multiplicative closure for 3+ factors.
- Do not rely on exact transfer `D <= 1/(1+lambda)` (fails on 4 known witnesses).
- Do not spin up new brute-force infrastructure unless a candidate lemma cannot be encoded in existing verifiers.

## Triage Command
Use this command shape for every incoming lemma candidate:

```bash
python3 scripts/closure_triage.py \
  --lemma-class <route1_exact_slack|strong_c2_hard_ratio|routeB_R_nonneg|routeC_modeP> \
  --lemma-name <short_id> \
  --note "<candidate summary>"
```

For fast artifact-only precheck (no full rerun), add:

```bash
--skip-full
```

## Triage Log

### 2026-02-22 00:06:02Z - route1_exact_slack - baseline_snapshot_2026_02_22
- Tight check: PASS (lower_bound_from_chain=0.186241122304, implies=True)
- Full verifier: SKIPPED
- Note: artifact baseline before first external candidate

### 2026-02-22 00:06:04Z - strong_c2_hard_ratio - baseline_snapshot_2026_02_22
- Tight check: PASS (hard_ratio_fail=0, rise_fail=0, combined_neg=0, min_hard_ratio_gap=0.901605807263)
- Full verifier: SKIPPED
- Note: artifact baseline before first external candidate

### 2026-02-22 00:06:08Z - routeB_R_nonneg - baseline_snapshot_2026_02_22
- Tight check: PASS (R_neg=0, R_min=4, cross_minus_abs_mismatch_neg=0, cross_minus_abs_mismatch_min=2)
- Full verifier: SKIPPED
- Note: artifact baseline before first external candidate

### 2026-02-22 00:06:11Z - routeC_modeP - baseline_snapshot_2026_02_22
- Tight check: PASS (p_ge_q=77141/77141, min_gap=1, qdrop_worst_p1_over_b1=0.9283765824998251)
- Full verifier: SKIPPED
- Note: artifact baseline before first external candidate

### 2026-02-22 00:09:39Z - route1_exact_slack - candidate_52_framework_intake_2026_02_22
- Tight check: PASS (lower_bound_from_chain=0.186241122304, implies=True)
- Full verifier: SKIPPED
- Note: 5.2 framework response; awaiting concrete derived inequality with filled packet

### 2026-02-22 00:09:39Z - strong_c2_hard_ratio - candidate_52_framework_intake_2026_02_22
- Tight check: PASS (hard_ratio_fail=0, rise_fail=0, combined_neg=0, min_hard_ratio_gap=0.901605807263)
- Full verifier: SKIPPED
- Note: 5.2 framework response; identifies hard-regime slope as core gap

### 2026-02-22 01:06:36Z - route1_exact_slack - candidate_52_primary_route_2026_02_22
- Tight check: PASS (lower_bound_from_chain=0.186241122304, implies=True)
- Full verifier: PASS (muP_fail=0, sum_identity_fail=0, chain_identity_fail=0, transfer_cap_violation=0)
- Full artifact: `results/triage/triage_route1_exact_slack_candidate_52_primary_route_2026_02_22.json`
- Note: full triage on first concrete 5.2 route

### 2026-02-22 01:07:45Z - strong_c2_hard_ratio - candidate_52_primary_route_2026_02_22
- Tight check: PASS (hard_ratio_fail=0, rise_fail=0, combined_neg=0, min_hard_ratio_gap=0.901605807263)
- Full verifier: PASS (combined_neg=0, rise_fail=0, hard_ratio_fail=0, identity_fail=0)
- Full artifact: `results/triage/triage_strong_c2_hard_ratio_candidate_52_primary_route_2026_02_22.json`
- Note: full triage on first concrete 5.2 route

### 2026-02-22 01:30:18Z - formal_lean_bridge_hardratio_route1 - local_formalization_progress
- Lean build: PASS (`lake build`)
- Formalized in `Formal/Basic.lean`:
  - bridge algebra identities (`mismatch`, `lcSurplus`, `combined` decompositions)
  - coefficient bridge to fugacity ratio (`lambda = (b1+b0+p0)/(b2+b1+p1)`)
  - Route-1 closure interface (`route1_iff_E`: `m-2<=muP` iff `exact_excess_D<=exact_slack_B` under chain identity)
  - hard-ratio interface (`rise_nonneg_iff_hard_ratio` and combined+sign -> hard-ratio theorem)
  - STRONG C2 stitch theorem (`strong_c2_of_lc_and_abs`)
- Note: This is symbolic scaffolding; remaining open content is proving the missing input lemmas (E/sign lemmas), not algebraic rewrites.

### 2026-02-22 01:32:59Z - hard_regime_sign_lemma_check - local_symbolic_triage
- Checked claimed sign lemma `(dq<0 and db>0) => (b1-b2<=0)` against certified hard-regime witnesses in `results/verify_strong_c2_rise_identity_2026_02_19.json`.
- Result: FAIL (2/2 hard-regime witnesses have `b1-b2>0`).
- Witness 1: n=20, m=7, `b1-b2=796`, `db=658`, `dq=-14`.
- Witness 2: n=23, m=8, `b1-b2=2983`, `db=1793`, `dq=-48`.
- Action: hard-regime target reset to minimal lemma `E_hard := rise-neg>=0` (equiv. `p1*db+b1*dq>=0`), and new focused prompt prepared: `notes/prompt_for_52pro_round5_2026-02-22.md`.
- Lean status: PASS (`lake build`) after adding `route1_and_hard_ratio_of_minimal_E` in `Formal/Basic.lean`.

### 2026-02-22 01:36:46Z - hard_regime_counterexample_formalized - local_formalization_progress
- Lean build: PASS (`lake build`).
- Added explicit formal counterexamples in `Formal/Basic.lean`:
  - `hard_regime_not_imply_b1_sub_b2_nonpos_example1`
  - `hard_regime_not_imply_b1_sub_b2_nonpos_example2`
- Purpose: lock in that `dq<0 and db>0` does not imply `b1-b2<=0`, preventing regression to invalid route.
- Closure chain retained via minimal hard lemma `E_hard := p1*db + b1*dq >= 0`.

### 2026-02-22 02:08:24Z - 52_round5_audit_and_round6_pivot - local_symbolic_triage
- Audited 5.2 round-5 output:
  - `E_hard*` candidate is structurally useful (`p1>=q1` and `db>=-2dq` => `rise-neg>=0`).
  - Route-1 replacement `lambda P'(lambda) >= (m-2)P(lambda)` is circular (equivalent to target), rejected.
- Verified hard-witness consistency from `results/verify_strong_c2_rise_identity_2026_02_19.json`:
  - hard_count=2, violations of `p1>=q1`: 0, violations of `db>=-2dq`: 0.
- Lean update: added theorem `rise_nonneg_of_star_conditions` in `Formal/Basic.lean`; build PASS.
- Next external prompt prepared with strict non-circularity constraint: `notes/prompt_for_52pro_round6_2026-02-22.md`.

### 2026-02-22 02:11:23Z - 52_analogical_output_audit - local_symbolic_triage
- Audited “out-of-domain analogies” response: selected winner (convex flux/entropy) rejected as circular.
- Reason: core assumption (existence of convex potential with secant `-dq/db` and derivative `p1/b1`) encodes target inequality directly.
- Action: prepared stricter follow-up prompt forbidding latent-potential circularity and requiring coefficient-level lemmas only:
  - `notes/prompt_for_52pro_round7_2026-02-22.md`

### 2026-02-22 02:34:39Z - constructor_theory_schema_and_prompt - local_formalization_progress
- Lean update (build PASS): added `ConstructorSchema` section in `Formal/Basic.lean` with explicit bad-task predicates:
  - `HardBad p1 b1 db dq := dq<0 ∧ db>0 ∧ p1*db+b1*dq<0`
  - `RouteBad exact_slack_B exact_excess_D := exact_slack_B < exact_excess_D`
  and impossibility lemmas from nonnegative monotones.
- Prepared constructor-theory focused external prompt:
  - `notes/prompt_for_52pro_round8_constructor_2026-02-22.md`
- Intent: force primitive-task + invariant-preservation lemma output, avoiding metaphor/circular convexity assumptions.

### 2026-02-22 02:44:54Z - 52_round8_constructor_output_audit - local_symbolic_triage
- Audited latest 5.2 output.
- Kept: hard-side decomposition to local lemmas
  - `H1_local: p0>=q0`
  - `H2_local: p1-p0 >= 3*(q0-q1)`
  and derivation `H1_local+H2_local => E_hard`.
- Rejected: Route-1 replacement `R1*` (false).
- Counterexample found in canonical regime:
  - `n=8`, `g6='G?`@F_'`, `m=3`, `lambda=21/23`
  - `p0=5,p1=8,pm=4,b1=10,b2=5`
  - `R1*_LHS=328965/12167 < R1*_RHS=529/2`.
  - Route-1 target still true on same witness: `mu_P-(m-2)=0.769580... > 0`.
- Action: prepared stricter follow-up prompt with explicit R1* counterexample and hard-side salvage:
  - `notes/prompt_for_52pro_round9_2026-02-22.md`.

### 2026-02-22 03:35:57Z - factorization_route_finalpass_setup - local_symbolic_triage
- Audited latest factorization-themed output: still missing explicit `J_attach(lambda)` and `J_boundary` formulas; no no-go theorem for candidate family.
- Status: does not meet binary success criterion for final pass.
- Action taken:
  - Prepared strict final-pass prompt with determinant normal form and explicit success/fail conditions:
    `notes/prompt_for_52pro_round11_factorization_finalpass_2026-02-22.md`
  - Added Lean validation hooks for matrix-factorization route in `Formal/Basic.lean`:
    `rise_nonneg_of_matrix_factorization`, `hard_ratio_of_det_factorization`.
  - Lean build PASS.

### 2026-02-22 11:16:29Z - matrix_shape_no_go_sweep - local_symbolic_triage
- Ran decisive shape-family sweep for final-pass 2x2 factorization on canonical data (`n<=14`).
- Candidate family tested (36 total):
  - `J_boundary` rows from `{(b1,q1),(b0,q0),(db,dq)}`
  - optional per-row transform `q -> q-b`
- No-go result summary:
  - 6 candidates are vacuous (always singular `J_boundary`, det=0).
  - 28 non-singular candidates fail gadget-only criterion: for same `lambda`, implied `J_attach = M_target * J_boundary^{-1}` differs across trees.
  - Only 2 candidates survive same-`lambda` consistency:
    1) `('hi',q-b ; 'diff',q)` gives `J_boundary = M_target`, `J_attach = I`.
    2) row-swapped variant gives constant swap matrix.
  - Both survivors are degenerate/circular: determinant sign obligation on `J_boundary` is equivalent to (or sign-flip equivalent to) target `rise-neg` sign.
- Conclusion: this plausible finite candidate family yields no non-circular factorization route.

### 2026-02-22 11:20:16Z - hard_local_chain_strengthened - local_formalization_progress
- Lean update (build PASS): added direct theorem
  - `hard_ratio_of_h1_h2_hard`
  giving `H1_local + H2_local + hard-regime signs => (-dq/db)<=p1/b1`.
- Empirical profile on certified `mismatch<0` witness set (`results/verify_strong_c2_rise_identity_2026_02_19.json`):
  - total mismatch-negative: 129
  - hard subset (`dq<0,db>0`): 2
  - H1_local failures: 0 (all), 0 (hard)
  - H2_local failures: 0 (all), 0 (hard)
  - min margins:
    - `min(p0-q0)=28` (all), `2842` (hard)
    - `min((p1-p0)-3*(q0-q1))=30` (all), `630` (hard)
- Route status: hard side now has a concrete local-lemma target; route-side (`E_route1`) still needs a non-circular replacement lemma.

### 2026-02-22 11:29:04Z - round12_stale_instance_audit - local_symbolic_triage
- Audited output from instance that previously ran factorization rounds.
- Hard side: accepted algebraic chain (matches local-lemma path already formalized):
  - `H1_local: p0>=q0`, `H2_local: p1-p0>=3*(q0-q1)` + hard signs => `E_hard`.
- Route-1 side: rejected.
  - Uses untrusted identification `exact_excess_D = p0*b2 - p1*b1` (not in packet).
  - Claimed `R1_new` fails immediately on known canonical witness (`n=8`, `g6='G?`@F_'`):
    `p1*b1-p0*b2 = 55 < 70 = (m-2)(b1^2-b0*b2)`.
- Additional evidence for structural obstruction under strict local symbol set:
  - Exhaustive canonical scan (`n<=14`) found many collisions where
    `(m,lambda,p0,p1,pm,q0,q1,qm,b0,b1,b2)` are identical but `mu_P(lambda)` differs.
  - Hence local tuple does not uniquely determine Route-1 witness value; strict local-only replacement is likely underdetermined.

### 2026-02-22 11:39:08Z - route1_tail_scaffold_added - local_formalization_progress
- Lean update (build PASS): added Route-1 tail-certificate scaffolding in `Formal/Basic.lean`:
  - `mu_ge_of_tail_certificate`
  - `route1_of_tail_certificate`
- Purpose: immediate slot-in point for 5.2 `R1_tail` output.
- Interface shape:
  - provide balance identity `lam*dP - (m-2)*Pval = headGain - tailDef`
  - prove certificate `tailDef <= headGain`
  - conclude `E_route1` via trusted chain identity.

### 2026-02-22 12:02:46Z - round13_candidate_audit - local_symbolic_triage
- Audited proposed `R1_tail: TailDef <= 2*pm*lambda^m` on canonical `d_leaf<=1` through `n<=20` (77,140 checked cases).
- Result: `R1_tail` fails (5 explicit counterexamples; all at `n=20`, `m=7`).
- Primitive obligations audit on same set:
  - Obligation (2): no failures observed.
  - Obligation (3): no failures observed.
  - Obligation (4): fails in all checked cases (77,140/77,140), including small trees; thus invalid as a proof step.
- Additional structural issue: proof text assumes `p_k=0` for `k>m`; this is generally false in canonical data.
- Conclusion: reject this `R1_tail` candidate and its obligation set.

### 2026-02-22 12:07:08Z - round13_reaudit_and_tail2_scan - local_symbolic_triage
- Re-ran a fresh canonical scan (`d_leaf<=1`, canonical `deg(s)=2` bridge, `n<=20`) with direct recomputation of `P=dp_B[u][0]`, `lambda=i_{m-1}(T)/i_m(T)`, and
  - `TailDef := sum_{k=0}^{m-3} (m-2-k) p_k lambda^k`.
- Reaudit result for prior candidate:
  - `R1_tail: TailDef <= 2*p_m*lambda^m` fails in `5/77,141` checked cases.
  - Explicit failing witnesses (all `n=20`, `m=7`):  
    `S???????????_?O?C??o?@_?@_??oFig?`,  
    `S?????????O?O?G?A??O?@??B??@_F|O?`,  
    `S?????????O?O?G?A??O?@??B?C@_B|O?`,  
    `S?????????O?O?G?A??O?@??B?E@_@|O?`,  
    `S?????????O?O?G?A??O?@??B??D_FxO?`.
  - Obligation `(4)` from that proof attempt (`p_{m-1} lambda^{m-1} <= p_m lambda^m`) fails in `77,141/77,141` checked cases.
- New stronger candidate scan:
  - `R1_tail2: TailDef <= p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m`.
  - Observed `0/77,141` failures through `n<=20` on the same canonical domain.
- Next action prepared:
  - Prompt file `notes/prompt_for_52pro_round14_route1_tail2_2026-02-22.md` created for 5.2 with strict non-circular constraints and forced `SUCCESS/BLOCKED` output.
- Lean interface updated (build PASS):
  - Added `mu_ge_of_tail2_certificate` and `route1_of_tail2_certificate` in `Formal/Basic.lean`.
  - `lake build` completed successfully after update.

### 2026-02-22 13:55:51Z - round14_output_audit - route1_tail2_triage
- Audited the round-14 proof output:
  - The closure chain `R1_tail2 => mu_P(lambda)>=m-2 => E_route1` is algebraically valid **if** `R1_tail2` is granted.
  - The proposed primitive proof of `R1_tail2` is invalid: it relies on one-step obligation
    `p_{k+2} lambda^{k+2} >= p_{k+1} lambda^{k+1} + p_k lambda^k` for all `0<=k<=m-2`.
- Fresh exhaustive canonical scan (`d_leaf<=1`, canonical `deg(s)=2`, `n<=20`, 77,141 checked cases):
  - Obligation above fails in `77,140/77,141` cases.
  - Earliest witness: `n=5`, `g6=DQo`, `m=2`, `lambda=5/6`,
    `p0=1, p1=2, p2=0`, so `lhs=0 < 2.666...=rhs`.
- Conclusion:
  - Keep `R1_tail2` as an empirically plausible candidate.
  - Reject the round-14 obligation set/proof method; request a repaired proof with new primitive obligations or a no-go.

### 2026-02-22 14:28:37Z - round15_output_audit - drift3_obstruction_rejected
- Audited round-15 response:
  - Correctly returns `BLOCKED`.
  - Introduces minimal obstruction lemma candidate `Drift3`:
    ` (m-2-k) p_k λ^k <= 2 p_{k+3} λ^{k+3} - p_{k+2} λ^{k+2} - p_{k+1} λ^{k+1}` for `0<=k<=m-3`.
- Fresh exhaustive canonical check (`d_leaf<=1`, canonical `deg(s)=2`, `n<=20`):
  - `Drift3` fails in `77,130/77,141` checked cases.
  - Failing instances (quantified `k`) fail in `77,130/295,675` evaluated `(tree,k)` pairs.
  - Earliest witness: n=8, g6 string = G?`@F_, m=3, lambda≈0.913043, k=0,
    LHS `=1`, RHS `≈-5.145`.
- Decision:
  - Reject `Drift3` as a valid obstruction lemma.
  - Continue with `R1_tail2` branch, but require any next obligation set to be DP-derived and witness-validated before claiming progress.
- Next action prepared:
  - `notes/prompt_for_52pro_round16_route1_tail2_dpstep_2026-02-22.md`
  - Enforces exact one-child DP-step derivation plus mandatory witness-pack falsification table (DQo, G?`@F_, and a known n=20 tail-failure witness for earlier candidates).

### 2026-02-22 15:36:26Z - round16_output_audit - identity_ok_obstruction_tautology
- Audited round-16 output:
  - One-step DP product identity is algebraically correct (independently spot-checked by random symbolic tests).
  - Closure chain from `R1_tail2` to `E_route1` remains valid.
  - Returned `BLOCKED`, which is appropriate.
- Critical issue:
  - Proposed “minimal obstruction lemma” is tautological after rearranging the telescoped identity:
    it is equivalent to `R1_tail2` rather than an independent local step lemma.
- Next action prepared:
  - `notes/prompt_for_52pro_round17_route1_tail2_nontaut_2026-02-22.md`
  - Adds an explicit non-tautology gate: local lemma `L` must exclude final-target symbols and must survive witness-pack checks before `SUCCESS`.

### 2026-02-22 16:45:41Z - round17_output_audit - local_bound_insufficient_and_M_rejected
- Audited round-17 response:
  - Correct `BLOCKED` status.
  - Acceptable: step-local identity for decrement `D(A,F)` and resulting lower bound `L` (degree-window bound).
  - Main issue: bound is too weak to close `R1_tail2` (as reported on W3).
- Additional audit of proposed “minimal additional local lemma” `M`:
  - `M`: weighted tail2 bound on partial products  
    `sum_{i>=t} a_i lambda^i <= a_t lambda^t + a_{t+1} lambda^{t+1}`.
  - Canonical scan (`d_leaf<=1`, canonical `deg(s)=2`, `n<=20`) over root-product partial steps:
    - case failures: `11,920 / 77,141` trees,
    - partial-step failures: `12,295 / 106,706` tested partial products.
  - Earliest failing witness appears at `n=8`, `g6=G?`@F_`, step `3`, cutoff `t=0`:
    LHS tail `≈ 7.9868`, RHS `≈ 4.6522`.
- Decision:
  - Reject lemma `M`.
  - Continue branch only with stronger non-tautological local state, or terminate local-tail route via a formal no-go statement.

### 2026-02-22 17:34:30Z - round18_output_audit - c_window_nogo_accepted_markov_insufficient
- Audited round-18 response:
  - Returns `BLOCKED` with an explicit no-go class (`C_window`) based on degree-window-only decrement bounds.
  - Internal numerics are consistent with direct recomputation on W3:
    - `max B_window ≈ 3.1263772916`,
    - required threshold for `R1_tail2` implication `≈ 4.3472138532`.
  - So `C_window` no-go is accepted.
- Tested proposed “next minimal lemma” (Markov-prefix bound using only `mu_A(lambda)`):
  - Per-step lower bound `prefix >= max(0, 1 - mu_A/(t+1))` is valid but too weak.
  - On W3, even best factor order gives aggregate lower bound only `≈ 2.9195 < 4.3472`.
- Decision:
  - Reject first-moment (Markov-prefix) strengthening as sufficient for closing `R1_tail2`.
  - Next step must use stronger local information than degree-window + first moment.
- Next action prepared:
  - `notes/prompt_for_52pro_round19_route1_localclass_nogo_2026-02-22.md`
  - Forces either:
    - a formal no-go for class `C_{deg,mu}` (degree + first moment only), or
    - a genuinely stronger local lemma outside rejected classes with witness-pack validation.

### 2026-02-22 18:16:44Z - round19_output_audit - c_deg_mu_nogo_provisionally_accepted_data_corrections
- Audited round-19 output:
  - Core claim (no-go for class `C_{deg,mu}`) is directionally correct and consistent with prior computed gap on W3:
    - best first-moment aggregate bound `≈ 2.9195`,
    - required threshold `≈ 4.3472`.
  - Proof idea (Markov-prefix extremality under `(deg,mu)`-only access) is acceptable as the working no-go mechanism.
- Required corrections:
  - Witness table had data errors:
    - W2 (`g6=G?`@F_`) true mode ratio is `i_{m-1}/i_m = 21/23 ≈ 0.913043`, not `7/8`.
    - W3 row used corrupted `g6` string with `*` characters; correct witness is
      `S???????????_?O?C??o?@_?@_??oFig?` (tree, `m=7`, `i_6=8088`, `i_7=8142`, `lambda=8088/8142=1348/1357`).
- Decision:
  - Treat `C_{deg,mu}` no-go as established for current workflow.
  - Proceed to next strictly stronger class (must include shape data beyond first moment), with hard witness-gated numerics using the corrected pack.
- Additional local verification:
  - On corrected W3 factors, class-optimal bound in `C_{deg,mu,mu2,mu3}` exceeds threshold:
    - `B_max ≈ 4.3923839874`,
    - threshold `≈ 4.3472138532`,
    - gap `≈ +0.04517`.
  - This suggests second-moment no-go is not the end-state; a concrete `mu3` local construction may be viable.
- Next action prepared:
  - `notes/prompt_for_52pro_round21_route1_mu3_dual_2026-02-22.md`
  - Forces explicit LP-dual one-step lemma in class `C_{deg,mu,mu2,mu3}` (or hard block), with corrected witness decoding required.

### 2026-02-22 21:14:24Z - round21_output_audit - mu3_witness_success_not_global
- Audited round-21 response:
  - Local LP-dual lemma in class `C_{deg,mu,mu2,mu3}` is structurally valid and non-tautological.
  - Witness-pack success claims are numerically consistent on corrected W3-old witness.
- Critical gap:
  - Round-21 claimed `SUCCESS`, but only witness-level validation was provided; no universal class guarantee.
- Full canonical audit performed (canonical `d_leaf<=1`, `n<=20`, all 77,141 checked trees):
  - Computed class-optimal `B_max(C_{deg,mu,mu2,mu3})` vs threshold per tree (all child-factor orders, canonical root decomposition).
  - Found a concrete global failure:
    - witness g6: `S???????C?G?G?C?@??G??_?@??@?F~_?`
    - `m=7`, `lambda≈0.9469606675`, child count `d=9`
    - threshold `≈ 4.3286125513`
    - `B_max ≈ 4.3046782725`
    - gap `≈ -0.0239342788` (insufficient).
- Decision:
  - Reject universal `SUCCESS` for class `C_{deg,mu,mu2,mu3}`.
  - Promote this witness as a formal no-go certificate for the mu3 class.
- Next action prepared:
  - `notes/prompt_for_52pro_round22_route1_mu3_nogo_2026-02-22.md`
  - Forces explicit mu3-class no-go statement (or strictly stronger class beyond mu3 with full witness checks).
- Next action prepared:
  - `notes/prompt_for_52pro_round20_route1_localclass_mu2_2026-02-22.md`
  - Resolves class `C_{deg,mu,mu2}` (degree + first + second factorial moment only) with mandatory exact g6 decode table before any bound claims.
