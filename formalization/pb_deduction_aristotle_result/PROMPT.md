# Aristotle request: Sections 2–3 (Theorem 1.1 from Proposition 3.1) of the Poisson–binomial paper

Work only in this Lean 4 project: Lean `v4.28.0`, Mathlib `v4.28.0`. Read
`PROOF_CONTEXT.md` first. It has a roadmap, Sections 2–3 of the paper verbatim,
and the Hillion–Johnson appendix.

## Goal

Fill the `sorry`s in `PBDeduction/Statement.lean` and make `lake build` pass.

`PBScalar/` is a **verified dependency**: Proposition 3.1 is
`PBScalar.scalar_inequality`, already kernel-checked. **Do not modify any file
in `PBScalar/`**. Import from it freely.

## Graded targets

- **G1 `deduction`:** the paper's Sections 2–3, conditional on `DeductionHyp`
  (positivity, support, total mass, strict log-concavity, Hillion–Johnson cubic
  inequalities) and on the maximal-mass bound `hmax`. This is the main target.
  For partial credit, state each proved step as its own theorem:
  - the recurrence (2.4);
  - the reciprocal bound (2.5);
  - Lemma 2.1;
  - the mass bounds;
  - `M = maxMass`;
  - `pairVar ≥ M² · A δ K`.
- **G2 `max_mass_bound`:** `(M⁻² − 1)/12 ≤ pairVar` for any mass function on
  `{0,…,n}`.
- **G3 `pb_*`:** the Poisson–binomial basics:
  - `pb_pos`, `pb_zero_out`, `pb_sum_one`, `pb_pairVar`;
  - `pb_strict_lc` (Newton's inequalities, strict);
  - `pb_first_descent_exists`.
- **G4 `pb_hj_left`, `pb_hj_right`:** Hillion and Johnson's Theorem A.2 and
  Corollary A.3. This is the hardest target, and a stretch goal.
- **G0 `theorem_1_1`:** assemble G1–G4.

Expected difficulty: G1, G2 and G3 should be within reach; G4 is a genuine
formalization project. **Do G1 first.**

## Refutation is a success mode

Every statement was tested in exact rational arithmetic on random instances
before submission. That covers 300 random Poisson–binomial laws for the
`DeductionHyp` properties and Theorem 1.1, and 2000 random mass functions for
G2, with no violations. The statements are believed true.

Report a refutation, with the exact counterexample, if you find:
- a statement false as written (for example an off-by-one in an index range,
  or an endpoint convention);
- that `DeductionHyp` plus `hmax` is too weak for G1.

A correct refutation is worth more than a proof.

## Requirements

- **No escape hatches:** no `sorry`, `admit`, new `axiom`, `implemented_by`, or
  `native_decide`. Kernel-checked computation (`decide`, `decide +kernel`,
  `norm_num`, `ring`, `nlinarith`, `positivity`) is fine.
- **Don't weaken statements:**
  - Don't weaken the theorem statements or change the definitions `deficit`,
    `pairVar`, `maxMass`, `IsFirstDescent`, `DeductionHyp`, `pgf`, `pbPmf`,
    `pbVar`.
  - If a definition needs an elaboration-only change, record it in a README
    traceability table, with the mathematics unchanged.
  - You may add files, definitions and lemmas freely.
- **Build:** run `lake build` before returning, and give `#print axioms` for
  every theorem in `PBDeduction/Statement.lean`.

## Deliverable and grading

- A compiling project.
- A README with a traceability table: each Lean declaration, the paper
  equation or step it formalizes, and its status.
- A self-assessment with one grade:
  - **COMPLETE:** G0 proved.
  - **PARTIAL:** list exactly which of G1–G4 and which sub-steps are proved.
  - **REFUTED:** a statement is false, with the counterexample.

An honest PARTIAL with G1 proved beats a grandiose FAILED in disguise. A result
that quietly strengthens a hypothesis or specializes a statement does not
count.
