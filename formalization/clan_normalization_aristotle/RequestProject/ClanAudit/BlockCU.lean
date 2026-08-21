import RequestProject.ClanAudit.DoubleSpider

/-!
# Central unimodality of the blocks on the derived range of scalars

`RequestProject/ClanAudit/Blocks.lean` proved

* `Fblock_not_centrally_unimodal` — the claimed central unimodality of the four-map
  block is *false* for general `c, d > 0`;
* `Fblock_decr` — but the decreasing half is true as soon as `c, d ≥ 1`.

This file completes the two remaining halves (symmetry and nonnegativity, which hold for
all `c, d ≥ 0`) and packages the result as `IsCU`, and does the same for the one-hub
block `A(r,c;z)`.  These are exactly the two facts consumed by the global partition:
every block of the partition is one of `Ablock r c`, `Fblock r s c d` (with derived
scalars `c, d` that are powers of two, hence `≥ 1`), or a product of "good" factors.
-/

namespace ClanAudit

open Finset LaurentPolynomial

theorem cbin_neg (n : ℕ) (m : ℤ) : cbin n (-m) = cbin n m := by
  rw [← coeffL_zz_pow, ← coeffL_zz_pow, (isCU_zz_pow n).symm]

theorem cbin_nonneg (n : ℕ) (m : ℤ) : 0 ≤ cbin n m := by
  rw [← coeffL_zz_pow]
  exact (isCU_zz_pow n).nonneg m

/-! ### The one-hub block -/

/-- **The one-hub block is centrally unimodal on the derived range `c ≥ 1`.** -/
theorem Ablock_isCU {r : ℕ} (hr : 1 ≤ r) {c : ℚ} (hc : 1 ≤ c) : IsCU (Ablock r c) := by
  have hsplit : Ablock r c = (c - 1) • zz ^ r + (zz ^ r + Pw r) := by
    rw [Ablock, sub_smul, one_smul]
    abel
  rw [hsplit]
  exact ((isCU_zz_pow r).smul (by linarith)).add (isCU_zz_pow_add_Pw r hr)

/-! ### The four-map block -/

theorem Fblock_comm (r s : ℕ) (c d : ℚ) : Fblock r s c d = Fblock s r d c := by
  rw [Fblock, Fblock, mul_comm, Nat.add_comm]

private theorem Fblock_symm_aux {r s : ℕ} (h : s ≤ r) (c d : ℚ) (m : ℤ) :
    coeffL (Fblock r s c d) (-m) = coeffL (Fblock r s c d) m := by
  rw [coeffL_Fblock h, coeffL_Fblock h, cbin_neg,
    show -m - (s : ℤ) = -(m + s) by ring, show -m + (s : ℤ) = -(m - s) by ring,
    show -m - (r : ℤ) = -(m + r) by ring, show -m + (r : ℤ) = -(m - r) by ring,
    cbin_neg, cbin_neg, cbin_neg, cbin_neg]
  have h1 : (-m = ((r - s : ℕ) : ℤ)) ↔ (m = -((r - s : ℕ) : ℤ)) := by omega
  have h2 : (-m = -((r - s : ℕ) : ℤ)) ↔ (m = ((r - s : ℕ) : ℤ)) := by omega
  simp only [h1, h2]
  ring

private theorem Fblock_nonneg_aux {r s : ℕ} (h : s ≤ r) {c d : ℚ} (hc : 0 ≤ c) (hd : 0 ≤ d)
    (m : ℤ) : 0 ≤ coeffL (Fblock r s c d) m := by
  rw [coeffL_Fblock h]
  have e1 : (0:ℚ) ≤ c * d * cbin (r + s) m :=
    mul_nonneg (mul_nonneg hc hd) (cbin_nonneg _ _)
  have e2 : (0:ℚ) ≤ c * (cbin r (m - s) + cbin r (m + s)) :=
    mul_nonneg hc (by have := cbin_nonneg r (m - s); have := cbin_nonneg r (m + s); linarith)
  have e3 : (0:ℚ) ≤ d * (cbin s (m - r) + cbin s (m + r)) :=
    mul_nonneg hd (by have := cbin_nonneg s (m - r); have := cbin_nonneg s (m + r); linarith)
  have e4 : (0:ℚ) ≤ (if m = ((r - s : ℕ) : ℤ) then 1 else 0)
      + (if m = -((r - s : ℕ) : ℤ) then 1 else 0) := by
    positivity
  linarith

/-- **The four-map block is centrally unimodal on the derived range `c, d ≥ 1`.**
Together with `Fblock_not_centrally_unimodal` this delimits exactly the range on which
the published claim is true. -/
theorem Fblock_isCU {r s : ℕ} (hr : 1 ≤ r) (hs : 1 ≤ s) {c d : ℚ} (hc : 1 ≤ c) (hd : 1 ≤ d) :
    IsCU (Fblock r s c d) := by
  rcases le_total s r with h | h
  · exact ⟨Fblock_symm_aux h c d, Fblock_nonneg_aux h (by linarith) (by linarith),
      Fblock_decr hr hs hc hd⟩
  · rw [Fblock_comm]
    exact ⟨Fblock_symm_aux h d c, Fblock_nonneg_aux h (by linarith) (by linarith),
      Fblock_decr hs hr hd hc⟩

end ClanAudit
