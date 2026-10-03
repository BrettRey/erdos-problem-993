# Summary of changes for run eccc6851-19fc-4712-9d4d-613cbb6d15ff
**Grade: COMPLETE.** `theorem_1_1` (G0) is proved, along with everything it depends on (G1–G4). `lake build` passes. There are no `sorry`, `admit`, new `axiom`, `implemented_by` or `native_decide` in `PBDeduction/`, and `PBScalar/` is unchanged. I found no refutation: every statement is proved exactly as written, and `DeductionHyp` plus `hmax` is enough for G1.

**Axioms.** Running `lake env lean PBDeduction/Axioms.lean` shows that each of `deduction`, `max_mass_bound`, `pb_pos`, `pb_zero_out`, `pb_sum_one`, `pb_pairVar`, `pb_strict_lc`, `pb_first_descent_exists`, `pb_hj_left`, `pb_hj_right` and `theorem_1_1` depends only on `[propext, Classical.choice, Quot.sound]`.

**What was proved, by target**

- **G1 `deduction`** (in `PBDeduction/Basic.lean`, `PBDeduction/Propagation.lean` and `PBDeduction/Deduction.lean`). Each partial-credit step is its own theorem:
  - (2.4) `deficit_recurrence`
  - (2.5) `reciprocal_bound`
  - Lemma 2.1 `propagation`
  - mass bounds `ratio_D_ge`, `right_mass_bound` (2.9), `left_mass_bound` (2.10)
  - `M = maxMass`: `maxMass_eq`
  - (3.2) `pairVar_ge_A`
  - the final scalar step `final_algebra`, which closes the argument using `PBScalar.scalar_inequality`.
- **G2 `max_mass_bound`** (`PBDeduction/MaxMass.lean`): a purely discrete proof. It compares each mass with a clipped antiderivative \(\Phi(y)=r^2y-y^3/3\), \(r=1/(2M)\), then telescopes the sum.
- **G3** (`PBDeduction/PBBasics.lean`, `PBDeduction/PBLogConcave.lean`): each result is proved by induction on the step that adds one Bernoulli summand.
  - Positivity, support, total mass, and variance via the first and second moments.
  - Strict Newton inequalities, using a Turán expansion together with a "cross" inequality.
  - Existence of the first descent, via \(f_{n-1}=f_n\sum(1-p_i)/p_i\) and \(V<\sum(1-p_i)/p_i\).
- **G4** (`PBDeduction/PBHillionJohnson.lean`): Hillion–Johnson (78), by induction using the cubic Bernstein expansion in the new parameter. The two middle coefficients are shown nonnegative by exact identities (\(c\,X = 2b\,C_1(k)+e\,C_1(k-1)\), and a similar one for \(c^2d\,Y\)), with the boundary cases handled through the support. (79) then follows from the reflection `pbPmf n p (n-k) = pbPmf n (1-p) k`.
- **G0**: `DeductionHyp` is assembled from G3 and G4, `hmax` comes from G2, and the result follows from `deduction`.

**One change outside the proofs.** The eight definitions (`deficit`, `pairVar`, `maxMass`, `IsFirstDescent`, `DeductionHyp`, `pgf`, `pbPmf`, `pbVar`) were moved, unchanged, into `PBDeduction/Defs.lean` so the proof files can import them; I checked mechanically that the text matches the original exactly. No theorem statement in `PBDeduction/Statement.lean` was changed; only the `sorry`s were replaced.

**Deliverables.** `README.md` contains the self-assessment, the change table, the file layout, and a traceability table mapping each Lean declaration to its paper equation or step and its status.