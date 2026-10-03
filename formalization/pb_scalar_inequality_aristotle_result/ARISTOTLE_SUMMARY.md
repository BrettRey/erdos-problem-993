# Summary of changes for run 16660e17-944d-4715-a182-78bdea524f43
**Grade: COMPLETE.** G0 (`scalar_inequality`), G1 (`scalar_inequality_compact`) and G2 (`scalar_inequality_large`) are proved with their statements unchanged. `lake build` passes, and there is no `sorry`, `admit`, new `axiom`, `implemented_by` or `native_decide`.

`#print axioms` (run with `lake env lean PBScalar/Axioms.lean`) gives only `[propext, Classical.choice, Quot.sound]` for all three theorems. The finite rational checks use `decide +kernel`, so the kernel itself checks them and `Lean.ofReduceBool` does not appear.

**No refutation.** Every printed identity I checked holds exactly, all 275 certificate coefficients are positive, the cells cover the whole range, the `K = m` versus `K = m − 1` boundary works out, and `R_r ≥ L_r` holds.

**How each part was proved:**
- **G2 (`H ≥ 16`)** follows the paper step by step:
  - Lean picks the triangular cell `J ≥ 5` itself and proves `J ≤ K`.
  - `b_r ≥ λ_r` comes from the Weierstrass inequality.
  - The closed forms `S̃_J` and `T̃_J`, and the expansion of the quartic `N_J`, are proved.
  - The `J = 5` Bernstein coefficients and the `J ≥ 6` identity `β_i = μ_i π_i(u)/2880` are each proved as exact identities by `ring`.
- **G1 (`3 < H ≤ 16`)** covers the same 13 cells but uses a different exact method from the paper's Bernstein route:
  - For fixed `K`, all the weights are nonnegative and shrink as `δ` grows. Because `A` is monotone in its weights, `A δ K` decreases in `δ`, and so does the target.
  - So on each sub-interval, one exact rational check at the endpoints is enough, evaluated in `ℚ` by the kernel.
  - Cell `[3,4]` (unsymmetrized weights) needs 20 sub-intervals; each cell `m = 4…15` needs one check.
  - The hypotheses `(K+1)δ < 1 ≤ (K+2)δ` put `H` in `(m, m+1]` when `K = m`, which settles the boundary question without a special case.
- **Shared lemmas, all proved:**
  - `L_eq_prod`: `L_r = ∏(1 − s/H)`
  - `R_ge_L`: `R_r ≥ L_r`
  - `Aw_mono`: monotonicity of `A` in nonnegative weights
  - `Aw_pairwise`: the pairwise identity, with `Aw_wsym` giving `A_sym = S·T`
  - `weierstrass_prod`: the Weierstrass product inequality

**Paper's certificate, checked separately (not needed for G0):**
- The printed identity `A − Q = P(H)/(4H⁵(H+1)³)` on `[3,4]` holds exactly.
- For all 13 cells in the data file, Lean proves the cleared-numerator identity, the Bernstein-expansion identity and positivity of every coefficient. From these it derives each cell by the paper's own route (`cell3_paper`, `cellm_paper`).
- One limit: for `m ≥ 4`, `cellm_paper` states `Q ≤ S_m T_m` and is not connected to `A` in Lean. G1 does not rely on it.

**Change to the given files:** the definitions `a`, `R`, `L`, `w`, `A` were moved word-for-word into `PBScalar/Defs.lean`, so the proof files can import them while the theorems stay in `PBScalar/Statement.lean`. This is recorded in the README.

**Files:**
- `README.md`: self-assessment, traceability table and change record
- `PBScalar/Statement.lean`: the three graded theorems
- `PBScalar/Defs.lean`: the definitions
- `PBScalar/Weights.lean`, `PBScalar/Quadratic.lean`: weight lemmas and the quadratic form
- `PBScalar/Compact.lean`: G1
- `PBScalar/Large.lean`: G2
- `PBScalar/PaperIdentities.lean`, `PBScalar/PaperCells/*.lean`: the certificate check, with the `PaperCells` files generated from the JSON by `scripts/regen_paper_cells.sh`
- `PBScalar/Axioms.lean`: the axiom printout