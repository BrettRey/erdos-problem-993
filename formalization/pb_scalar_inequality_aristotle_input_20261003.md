# Aristotle request: the scalar inequality (Proposition 3.1) of the Poisson–binomial paper

Work only in this minimal Lean 4 project. The toolchain is Lean `v4.28.0` with
Mathlib `v4.28.0`. Read `PROOF_CONTEXT.md` first. It gives the definitions, a
roadmap, and Section 4 of the paper verbatim.

## Goal

Fill the `sorry`s in `PBScalar/Statement.lean` and make `lake build` pass. The
target is a single explicit real inequality: for `0 < δ < 1/4`, with
`(K+1)δ < 1 ≤ (K+2)δ`,

```text
(3 + δ) / (4 δ^2) ≤ A δ K,   where   A δ K = (1/2) ∑_{i,j=-K}^{K} w_i w_j (i - j)^2
```

and the weights `w` are built from `R_r` and `L_r` as defined in the file.

This step is the computer-assisted part of the paper. A Lean proof would
replace "two independent programs and a third-party reconstruction agree" with
a kernel-checked statement. That is the point of the request.

## Graded targets

- **G1 (compact range):** `scalar_inequality_compact`, for `1/17 ≤ δ < 1/4`.
  This consists of 13 cells in `H = 1/δ - 1`:
  - `[3,4]`, where the weights stay unsymmetrized;
  - `[m, m+1]` for `m = 4, …, 15`, which use the symmetrized `S_m T_m ≥ Q`.

  **Partial credit is per cell.** If you finish some cells and not others,
  state each finished cell as its own theorem and list which cells remain.
- **G2 (large range):** `scalar_inequality_large`, for `0 < δ ≤ 1/17`. This
  uses the triangular cells `J ≥ 5`, the closed forms `S̃_J` and `T̃_J`, and
  the explicit quartic `N_J`. The `J ≥ 6` case is a polynomial identity in
  `u` (and `t`) with positive coefficients.
- **G0:** `scalar_inequality`. Assemble G1 and G2.
- **Shared lemmas worth proving once:**
  - `L_r = ∏_{s=1}^r (1 - s/H)`;
  - `R_r ≥ L_r` for `r ≤ K`;
  - monotonicity of `A` in nonnegative weights;
  - the pairwise identity `A_sym = S·T`;
  - the Weierstrass product inequality.

## Refutation is a success mode

The target is believed true. Two independent exact programs and an independent
reconstruction found every certificate coefficient positive. If you find any of
the following false, report it with the exact counterexample (the cell, the
value of `H` or `t`, the failing identity or coefficient) and stop that branch:
- a printed identity;
- a Bernstein coefficient's sign;
- the coverage of the cells;
- the `K = m` versus `K = m - 1` boundary handling;
- `R_r ≥ L_r`.

A correct refutation is worth more than a proof.

## Requirements

- **No escape hatches:** no `sorry`, `admit`, new `axiom`, `implemented_by`,
  or `native_decide`. Kernel-checked computation (`decide`, `norm_num`,
  `ring`, `nlinarith`, `positivity`, `polyrith`-produced certificates) is fine.
- **Statements:** do not weaken the three theorem statements. The definitions
  may be adjusted only for elaboration (coercions, `Finset` idioms) with the
  mathematics unchanged. Record any such change in a README traceability
  table.
- **Freedom of route:** you may restructure the proof freely. For example:
  - encode the 13 cell polynomials from the data file as Lean `Polynomial`s;
  - prove `4H^{2m}(S_m T_m - Q) = P_m(H)` by `ring`;
  - prove `P_m(m+t) = ∑ β_i binom(d,i) t^i (1-t)^{d-i}` coefficientwise;
  - or use any other exact method (interval subdivision with polynomial
    bounds, sum-of-squares certificates) that avoids `native_decide`.
- **Build:** run `lake build` before returning, and give `#print axioms` for
  the main theorems.

## Deliverable and grading

- A compiling project.
- A README with a traceability table: each Lean declaration, the paper
  equation or step it formalizes, and its status.
- A self-assessment with one grade:
  - **COMPLETE:** G0 proved.
  - **PARTIAL:** say exactly which of G1's cells, which of G2's sub-steps, and
    which shared lemmas are proved.
  - **REFUTED:** a statement is false, with the counterexample.

An honest PARTIAL with one real proved cell beats a grandiose FAILED in
disguise. A COMPLETE that quietly specializes the statement does not count:
for example, proving the inequality at sampled `δ` only, or assuming the
Bernstein coefficients without deriving the polynomial identities.
