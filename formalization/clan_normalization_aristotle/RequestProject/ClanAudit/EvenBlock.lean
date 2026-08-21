import RequestProject.EvenConnector
import RequestProject.ClanAudit.Blocks

/-!
# The even-connector block in the notation of the clan development

The even-connector Laurent inequality is proved in `RequestProject/EvenConnector.lean`
(with its supporting files `RequestProject/Binomial.lean` and
`RequestProject/Coeff.lean`) as `EvenConnector.Gpoly_centrallyUnimodal`.  That
development uses its own names `XL`, `HL`, `lcoeff` and `CentrallyUnimodal`; this file
identifies them with `zz`, `Pw`, `coeffL` and `IsCU` of the clan development, so that the
even-connector block can be used exactly like `Fblock`.

* `Gblock r s c d` — `c d (z+z⁻¹)^(r+s+1) + d (z+z⁻¹)^s (z^(r+1)+z^(-(r+1)))
  + c (z+z⁻¹)^r (z^(s+1)+z^(-(s+1))) + z^(r+s+1) + z^(-(r+s+1))`;
* `Gblock_eq` — it *is* `EvenConnector.Gpoly r s c d`;
* `Gblock_isCU` — for `r, s ≥ 1` and `c, d ≥ 1` it is centrally unimodal.
-/

namespace ClanAudit

open LaurentPolynomial

/-- The two notions of central unimodality agree. -/
theorem isCU_of_centrallyUnimodal {f : LaurentPolynomial ℚ}
    (h : EvenConnector.CentrallyUnimodal f) : IsCU f where
  symm := h.symm
  nonneg := h.nonneg
  decr := h.decr

theorem zz_eq_XL : zz = EvenConnector.XL := by rw [zz_eq, EvenConnector.XL]

theorem Pw_eq_HL (t : ℕ) : Pw t = EvenConnector.HL (t : ℤ) := by rw [Pw, EvenConnector.HL]

/-- The even-connector block, in the notation of the clan development. -/
noncomputable def Gblock (r s : ℕ) (c d : ℚ) : LaurentPolynomial ℚ :=
  (c * d) • zz ^ (r + s + 1) + d • (zz ^ s * Pw (r + 1)) + c • (zz ^ r * Pw (s + 1))
    + Pw (r + s + 1)

theorem Gblock_eq (r s : ℕ) (c d : ℚ) : Gblock r s c d = EvenConnector.Gpoly r s c d := by
  rw [Gblock, EvenConnector.Gpoly, zz_eq_XL, Pw_eq_HL, Pw_eq_HL, Pw_eq_HL]
  push_cast
  simp only [smul_eq_C_mul, map_mul]
  ring

/-- **The even-connector block is centrally unimodal** for `r, s ≥ 1` and `c, d ≥ 1`. -/
theorem Gblock_isCU {r s : ℕ} {c d : ℚ} (hr : 1 ≤ r) (hs : 1 ≤ s) (hc : 1 ≤ c) (hd : 1 ≤ d) :
    IsCU (Gblock r s c d) := by
  rw [Gblock_eq]
  exact isCU_of_centrallyUnimodal (EvenConnector.Gpoly_centrallyUnimodal r s c d hr hs hc hd)

end ClanAudit
