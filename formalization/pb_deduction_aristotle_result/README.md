This project was edited by [Aristotle](https://aristotle.harmonic.fun).

To cite Aristotle:
- Tag @Aristotle-Harmonic on GitHub PRs/issues
- Add as co-author to commits:
```
Co-authored-by: Aristotle (Harmonic) <aristotle-harmonic@harmonic.fun>
```

# Sections 2–3 of "Variance and local log-concavity of Poisson–binomial laws"

## Self-assessment: **COMPLETE**

`PBDeduction.theorem_1_1` (G0, Theorem 1.1) is proved, and so are all of its
ingredients G1–G4. There are no `sorry`, `admit`, new `axiom`,
`implemented_by` or `native_decide` anywhere in `PBDeduction/`. `PBScalar/`
was not modified.

`lake build` succeeds. `lake env lean PBDeduction/Axioms.lean` prints:

```
'PBDeduction.deduction' depends on axioms: [propext, Classical.choice, Quot.sound]
'PBDeduction.max_mass_bound' depends on axioms: [propext, Classical.choice, Quot.sound]
'PBDeduction.pb_pos' depends on axioms: [propext, Classical.choice, Quot.sound]
'PBDeduction.pb_zero_out' depends on axioms: [propext, Classical.choice, Quot.sound]
'PBDeduction.pb_sum_one' depends on axioms: [propext, Classical.choice, Quot.sound]
'PBDeduction.pb_pairVar' depends on axioms: [propext, Classical.choice, Quot.sound]
'PBDeduction.pb_strict_lc' depends on axioms: [propext, Classical.choice, Quot.sound]
'PBDeduction.pb_first_descent_exists' depends on axioms: [propext, Classical.choice, Quot.sound]
'PBDeduction.pb_hj_left' depends on axioms: [propext, Classical.choice, Quot.sound]
'PBDeduction.pb_hj_right' depends on axioms: [propext, Classical.choice, Quot.sound]
'PBDeduction.theorem_1_1' depends on axioms: [propext, Classical.choice, Quot.sound]
```

No refutation was found. Every statement is proved exactly as written.

## Changes to the original files (elaboration-only)

| Change | Reason | Mathematics |
|---|---|---|
| The definitions `deficit`, `pairVar`, `maxMass`, `IsFirstDescent`, `DeductionHyp`, `pgf`, `pbPmf`, `pbVar` moved from `PBDeduction/Statement.lean` to `PBDeduction/Defs.lean` | This lets the proof files import them without import cycles (the same arrangement as `PBScalar/Defs.lean`) | Unchanged: the text is copied character for character, which was checked mechanically against the original file |
| In `PBDeduction/Statement.lean`, each `:= by sorry` is replaced by a proof term | — | No theorem statement was changed |
| New files under `PBDeduction/` and `PBDeduction/Axioms.lean` | Proofs and the axiom report | — |

## File layout

| File | Content |
|---|---|
| `PBDeduction/Defs.lean` | The original definitions, verbatim |
| `PBDeduction/Basic.lean` | Basic facts about deficits, recurrence (2.4), reciprocal bound (2.5) |
| `PBDeduction/Propagation.lean` | Lemma 2.1 and the mass bounds (2.8)–(2.10) |
| `PBDeduction/Deduction.lean` | `M = maxMass`, pairwise variance (3.2), final step, `deduction_main` |
| `PBDeduction/MaxMass.lean` | G2: the maximal-mass bound (discrete proof) |
| `PBDeduction/PBBasics.lean` | Bernoulli recursion, support, total mass, mean, second moment, variance |
| `PBDeduction/PBLogConcave.lean` | Positivity, strict Newton inequalities, existence of the first descent |
| `PBDeduction/PBHillionJohnson.lean` | G4: Hillion–Johnson (78) and (79) |
| `PBDeduction/Statement.lean` | The graded statements, now proved |

## Traceability table

All entries are **proved**, with no `sorry`.

### G1: the deduction (Sections 2–3, conditional on `DeductionHyp` and `hmax`)

| Lean declaration | Paper step | Status |
|---|---|---|
| `DeductionHyp.deficit_pos`, `deficit_le_one`, `deficit_zero`, `deficit_top` | `0 < δ_k ≤ 1` on the support; `δ_0 = δ_n = 1` | proved |
| `deficit_recurrence` | (2.4): `δ_{k-1}(1-δ_k) ≤ δ_k`, `δ_{k+1}(1-δ_k) ≤ δ_k` for `1 ≤ k ≤ n-1` | proved |
| `reciprocal_bound` | (2.5): `\|1/δ_{j+1} - 1/δ_j\| ≤ 1` for `0 ≤ j ≤ n-1` | proved |
| `propagation` | Lemma 2.1: for `(r+1)δ < 1`, `D±r ∈ {1,…,n-1}` and `δ_{D±r} ≤ δ/(1-rδ) < 1` | proved |
| `exists_K` | Definition of `K`: `(K+1)δ < 1 ≤ (K+2)δ` and `K ≥ 3` | proved |
| `window_bounds` | `D - K ≥ 1` and `D + K ≤ n-1` | proved |
| `ratio_succ` | (2.6): `q_{k+1} = q_k(1-δ_k)` | proved |
| `ratio_D_ge` | `q_D ≥ a` (using `q_{D-1} ≥ 1`) | proved |
| `ratio_right_ge` | `q_{D+j} ≥ a(1-jδ)` for `0 ≤ j ≤ K-1` | proved |
| `ratio_left_le` | `q_{D-j}(1-(j+1)δ)/(1-δ) ≤ q_D`, hence `1/q_{D-j} > (1-(j+1)δ)/(1-δ)` | proved |
| `right_mass_bound` | (2.9): `f(c+r) ≥ M R_r` for `0 ≤ r ≤ K` | proved |
| `left_mass_bound` | (2.10): `f(c-r) ≥ M L_r` for `0 ≤ r ≤ K` | proved |
| `maxMass_eq` | `M = f_c = max_k f_k` (unimodality) | proved |
| `pairVar_ge_window`, `pairVar_ge_A` | (3.2): `pairVar f n ≥ M² · A δ K` | proved |
| `final_algebra` | Proof of Theorem 1.1 from Proposition 3.1 (end of Section 3) | proved |
| `deduction_main` / **`deduction`** | **G1** | **proved** |

### G2

| Lean declaration | Paper step | Status |
|---|---|---|
| `pairVar_eq_central` | Pairwise variance equals the central second moment (when `∑ f = 1`) | proved |
| `PhiR`, `PhiC`, `T_ge`, `PhiC_le`, `PhiC_ge`, `sum_Icc_telescope` | A discrete version of the uniform-smoothing comparison in Section 3 | proved |
| `max_mass_bound_proof` / **`max_mass_bound`** | (3.3): `V ≥ (M⁻² - 1)/12` (Bobkov–Marsiglietti–Melbourne, Cor. 3.2) | **proved** |

The G2 proof is purely discrete. Let `μ` be the mean, `r = 1/(2M)` and
`Φ(y) = r²y - y³/3`, and write `Φc = Φ ∘ clip_{[-r,r]}`. For every real `y`,
the quantity `T(y) = Φc(y+½) - Φc(y-½)` is nonnegative and at least
`r² - y² - 1/12`. Since `0 ≤ f_k ≤ M`, this gives
`f_k((k-μ)² + 1/12 - r²) ≥ -M·T(k-μ)`. Summing over `k` and telescoping gives
`V + 1/12 - r² ≥ -4Mr³/3`, which is exactly the bound.

### G3

| Lean declaration | Statement | Status |
|---|---|---|
| `pbPmf_succ` | `f^{[n+1]}_k = (1-p_0) f^{[n]}_k + p_0 f^{[n]}_{k-1}` | proved |
| `pbPmf_nonneg_pos` / **`pb_pos`** | G3a | **proved** |
| `pbPmf_zero_out` / **`pb_zero_out`** | G3b | **proved** |
| `pbPmf_sum_one` / **`pb_sum_one`** | G3c | **proved** |
| `pbPmf_mean`, `pbPmf_second`, `pbPmf_pairVar` / **`pb_pairVar`** | G3d | **proved** |
| `cross_ineq`, `conv_turan`, `pbPmf_lc` / **`pb_strict_lc`** | G3e, Newton (strict), proved by induction on `n` | **proved** |
| `pbPmf_penultimate`, `pbPmf_first_descent_exists` / **`pb_first_descent_exists`** | G3f, via `f_{n-1} = f_n ∑(1-p_i)/p_i` and `V < ∑(1-p_i)/p_i` | **proved** |

### G4 (Hillion–Johnson, Appendix A)

| Lean declaration | Statement | Status |
|---|---|---|
| `hjC1`, `hjC1_conv` | `C1(k)` (78) and its cubic Bernstein expansion in the new parameter | proved |
| `hj_middle_nonneg` | The coefficients of `p(1-p)²` (`D1(k)`, (84)) and `p²(1-p)` (Prop. A.4) are nonnegative | proved |
| `pbPmf_hjC1` | Theorem A.2: `C1(k) ≥ 0` for all `k` | proved |
| `pbPmf_reflect`, `pbPmf_hjC1_dual` | Duality `k ↦ n-k`, `p ↦ 1-p`; Corollary A.3 (79) | proved |
| **`pb_hj_left`**, **`pb_hj_right`** | G4a, G4b | **proved** |

The G4 certificates were checked with exact `ring` identities. Write
`a,…,e = g_{k-3},…,g_{k+1}` and `D_{k-1} = c² - bd`. Then:

- `c·X = 2b·C1(k) + e·C1(k-1)`, where `X = ade - 3bce + 2bd²`.
- `c²d·Y = b D_{k-1} C1(k) + e D_{k-1} C1(k-1) + D_{k-1}²(cd - be) + acd·C1(k) + 2c²e·C1(k-1)`,
  where `Y = ace + ad² - 2b²e - bcd + c³`.

When `c = 0` or `d = 0`, the support structure is used instead.

### G0

| Lean declaration | Statement | Status |
|---|---|---|
| **`theorem_1_1`** | Theorem 1.1: `V·δ_D ≥ 1/4`. It combines G1 with G2 (for `hmax`), G3 and G4 (for `DeductionHyp`) | **proved** |
