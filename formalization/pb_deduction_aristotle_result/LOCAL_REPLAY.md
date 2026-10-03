# Local replay of the Aristotle return for Sections 2–3 and Theorem 1.1
<!-- SUMMARY: local replay of Aristotle's Lean proof of Theorem 1.1 of the PB paper (Sections 2–3 plus cited inputs) · status: replay passed; Theorem 1.1 kernel-checked, standard axioms only · updated: 2026-10-03 -->

**Return.**
- Aristotle (Harmonic) project `e7a03302-49a3-4402-b8a8-9f88ed293b55`, run `eccc6851-19fc-4712-9d4d-613cbb6d15ff`.
- Downloaded 2026-10-03 with `scripts/aristotle_cli.py result`. The archive, kept outside git, has sha256 `65a75564f5be8325d445d4cf5d31deb41d3313b155d884989931d6677cedb667`.
- The packet sent was `formalization/pb_deduction_aristotle/`, prompt copy `formalization/pb_deduction_aristotle_input_20261003.md`.

**Self-assessment (Aristotle's claim):** COMPLETE (G0–G4).

## What was replayed here (Claude Code, Opus 5.5, 2026-10-03)

1. **Build.** `lake build` on Lean 4.28.0 with Mathlib v4.28.0, using the same
   Mathlib packages as the Proposition 3.1 replay: "Build completed
   successfully (8049 jobs)".
2. **Axioms.** `lake env lean PBDeduction/Axioms.lean` gives
   `[propext, Classical.choice, Quot.sound]` for each of the following, with
   neither `Lean.ofReduceBool` nor `Lean.trustCompiler` appearing:
   - the abstract deduction and its inputs: `deduction`, `max_mass_bound`;
   - the Poisson–binomial basics: `pb_pos`, `pb_zero_out`, `pb_sum_one`,
     `pb_pairVar`, `pb_strict_lc`, `pb_first_descent_exists`;
   - the Hillion–Johnson inequalities: `pb_hj_left`, `pb_hj_right`;
   - the main result: `theorem_1_1`.
3. **Escape hatches.** I grepped `PBDeduction/` and `PBDeduction.lean` for
   `sorry`, `admit`, `axiom`, `implemented_by`, `native_decide`,
   `ofReduceBool`, `trustCompiler`, `unsafe`, `@[extern`, `opaque`,
   `skipKernelTC`, `macro_rules`, `elab` and `set_option`. There were no
   matches.
4. **Statements.**
   - The seven definitions (`deficit`, `pairVar`, `maxMass`, `IsFirstDescent`,
     `pgf`, `pbPmf`, `pbVar`) and the structure `DeductionHyp` in
     `PBDeduction/Defs.lean` are identical, after whitespace normalization, to
     the packet's `PBDeduction/Statement.lean`.
   - All 11 theorem signatures in the returned `Statement.lean` are identical
     to the packet's up to `:=`.
   - The returned `Statement.lean` defines nothing.
5. **Dependency unchanged.** The returned `PBScalar/` and `PBScalar.lean` are
   byte-identical to the replayed Proposition 3.1 project
   (`formalization/pb_scalar_inequality_aristotle_result/`).
6. **Inputs unchanged.** `PROMPT.md`, `PROOF_CONTEXT.md`, `lakefile.toml`,
   `lake-manifest.json`, `lean-toolchain` and
   `data/PBReserve_Core_reference.lean` are byte-identical to the packet's.
7. **Statements tested before submission.** Exact rational tests on 300
   random Poisson–binomial laws covered the `DeductionHyp` properties and
   Theorem 1.1, and 2000 random mass functions covered `max_mass_bound`. No
   violations.

## Fidelity of `theorem_1_1` to the paper's Theorem 1.1

| Paper | Lean |
|---|---|
| $W=\sum_{i=1}^nB_i$, $B_i\sim\mathrm{Bernoulli}(p_i)$ independent, standing assumption $0<p_i<1$ | `p : Fin n → ℝ`, `hp : ∀ i, 0 < p i ∧ p i < 1` |
| $f_k=\mathbb P(W=k)$, extended by $f_{-1}=f_{n+1}=0$ | `pbPmf n p k` = coefficient of $X^k$ in $\prod_i(1-p_i+p_iX)$ for $k\geq0$, and $0$ for $k<0$ (and $k>n$ by degree) |
| $V=\sum p_i(1-p_i)\geq1$ | `pbVar n p`, `hV : 1 ≤ pbVar n p` |
| $D=\min\{1\leq k\leq n: f_k<f_{k-1}\}$ | `IsFirstDescent (pbPmf n p) n D`: $1\leq D\leq n$, $f_D<f_{D-1}$, and $f_{k-1}\leq f_k$ for $1\leq k<D$. Existence: `pb_first_descent_exists` |
| $V\delta_D\geq1/4$, $\delta_D=1-f_{D-1}f_{D+1}/f_D^2$ | `1 / 4 ≤ pbVar n p * deficit (pbPmf n p) D` |

## Verdict

**Replay passed.** Theorem 1.1 is kernel-checked in Lean 4 with Mathlib,
starting from the definition of a Poisson–binomial law through its
probability-generating polynomial, and depends only on the three standard
axioms. The formal proof includes:
- the inputs the paper cites: the Hillion–Johnson cubic inequalities, the
  strict log-concavity of the pmf (the consequence of Newton's inequalities
  that the proof uses, not the binomially normalized Newton inequalities),
  and the maximal-mass bound of Bobkov,
  Marsiglietti and Melbourne (proved discretely, by a different route from the
  paper's §3);
- Proposition 3.1, via `PBScalar`.

**Not formalized:**
- the reduction for $p_i\in\{0,1\}$ (§2.1, an elementary translation);
- the extension to differences $W-W'$;
- Corollary 1.3, Propositions 1.2, 1.4 and 1.5, Table 1, Example 1.6.
