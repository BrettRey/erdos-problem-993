import Mathlib

/-!
# Definitions for the cross-check of the paper's cells `[m, m+1]`, `m = 4, …, 15`

Not used by the main proof (which uses a different exact method on the compact range).
For each cell, generated from
`data/universal_pb_finite_bernstein_full_certificate_2026-07-16.json`:

* `cellm_numerator`: `4 H^{2m} (S_m T_m - Q(H))` equals the certificate's numerator
  polynomial at `H = m + t` (exact identity, `field_simp; ring`);
* `cellm_bernstein`: that polynomial equals the certificate's Bernstein expansion
  `∑ β_i binom(d,i) t^i (1-t)^{d-i}` (exact identity, `ring`);
* `cellm_coeffs_pos`: every `β_i` is positive (`norm_num`).
-/

namespace PBScalar.Certificate

/-- `b_r = ∏_{s=1}^r (1 - s/H)`. -/
noncomputable def bH (H : ℝ) (r : ℕ) : ℝ :=
  ∏ s ∈ Finset.range r, (1 - ((s : ℝ) + 1) / H)

/-- `S_m = 1 + 2 ∑_{r=1}^m b_r`. -/
noncomputable def Sm (H : ℝ) (m : ℕ) : ℝ := 1 + 2 * ∑ r ∈ Finset.range m, bH H (r + 1)

/-- `T_m = 2 ∑_{r=1}^m r² b_r`. -/
noncomputable def Tm (H : ℝ) (m : ℕ) : ℝ :=
  2 * ∑ r ∈ Finset.range m, ((r : ℝ) + 1) ^ 2 * bH H (r + 1)

/-- `Q(H) = (3H+4)(H+1)/4`. -/
noncomputable def Qf (H : ℝ) : ℝ := (3 * H + 4) * (H + 1) / 4

end PBScalar.Certificate
