import PBScalar.Quadratic

/-!
# The large range `0 < δ ≤ 1/17`, i.e. `H ≥ 16` (G2)

Following Section 4.2 of the paper: triangular cells `J(J+1)/2 ≤ H ≤ (J+1)(J+2)/2`
with `J ≥ 5`, the lower bounds `λ_r = 1 - r(r+1)/(2H)`, the closed forms `S̃_J`, `T̃_J`,
the explicit quartic `N_J`, and its degree-four Bernstein expansions
(`J = 5` numerically, `J ≥ 6` as a polynomial identity in `u = J - 6` and `t`).
-/

namespace PBScalar

open Finset

/-- `λ_r = 1 - r(r+1)/(2H)` (paper, eq. `bonferroni`), indexed so that `lam H r = λ_{r+1}`. -/
noncomputable def lam (H : ℝ) (r : ℕ) : ℝ := 1 - ((r : ℝ) + 1) * ((r : ℝ) + 1 + 1) / (2 * H)

/-- `S̃_J = 1 + 2 ∑_{r=1}^J λ_r`. -/
noncomputable def St (J : ℕ) (H : ℝ) : ℝ := 1 + 2 * ∑ r ∈ Finset.range J, lam H r

/-- `T̃_J = 2 ∑_{r=1}^J r² λ_r`. -/
noncomputable def Tt (J : ℕ) (H : ℝ) : ℝ :=
  2 * ∑ r ∈ Finset.range J, ((r : ℝ) + 1) ^ 2 * lam H r

/-- `σ_3 + σ_4` as a function of a real `J`. -/
noncomputable def sig34 (J : ℝ) : ℝ :=
  J ^ 2 * (J + 1) ^ 2 / 4 + J * (J + 1) * (2 * J + 1) * (3 * J ^ 2 + 3 * J - 1) / 30

/-- Paper, eq. `S0`: `S̃_J = 2J + 1 - J(J+1)(J+2)/(3H)`. -/
theorem St_closed (J : ℕ) (H : ℝ) (hH : H ≠ 0) :
    St J H = 2 * J + 1 - (J : ℝ) * (J + 1) * (J + 2) / (3 * H) := by
  have : ∑ r ∈ Finset.range J, lam H r = J - (J : ℝ) * (J + 1) * (J + 2) / (6 * H) := by
    induction J with
    | zero => simp
    | succ J ih =>
      rw [Finset.sum_range_succ, ih, lam]
      push_cast
      field_simp
      ring
  rw [St, this]
  field_simp
  ring

/-- Paper, eq. `T0`: `T̃_J = J(J+1)(2J+1)/3 - (σ_3 + σ_4)/H`. -/
theorem Tt_closed (J : ℕ) (H : ℝ) (hH : H ≠ 0) :
    Tt J H = (J : ℝ) * (J + 1) * (2 * J + 1) / 3 - sig34 J / H := by
  have : ∑ r ∈ Finset.range J, ((r : ℝ) + 1) ^ 2 * lam H r =
      (J : ℝ) * (J + 1) * (2 * J + 1) / 6 - sig34 J / (2 * H) := by
    induction J with
    | zero => simp [sig34]
    | succ J ih =>
      rw [Finset.sum_range_succ, ih, lam, sig34, sig34]
      push_cast
      field_simp
      ring
  rw [Tt, this]
  field_simp
  ring

/-- The explicit quartic `N_J(H)` displayed in the paper. -/
noncomputable def NJ (J H : ℝ) : ℝ :=
  -(3 / 4) * H ^ 4 - (7 / 4) * H ^ 3 + (J * (J + 1) * (2 * J + 1) ^ 2 / 3 - 1) * H ^ 2
    - ((2 * J + 1) * sig34 J + J ^ 2 * (J + 1) ^ 2 * (J + 2) * (2 * J + 1) / 9) * H
    + J * (J + 1) * (J + 2) / 3 * sig34 J

/-- Paper: `N_J(H) = H² (S̃_J T̃_J - Q(H))` expands to the displayed quartic. -/
theorem NJ_eq (J H : ℝ) (hH : H ≠ 0) :
    H ^ 2 * ((2 * J + 1 - J * (J + 1) * (J + 2) / (3 * H)) *
        (J * (J + 1) * (2 * J + 1) / 3 - sig34 J / H) - (3 * H + 4) * (H + 1) / 4) =
      NJ J H := by
  unfold NJ
  field_simp
  ring

/-- `J = 5`: Bernstein expansion of `N_5(16 + 5t)` with coefficients
`2360, 7500, 25055/2, 254205/16, 31115/2`. -/
theorem NJ_five_bernstein (t : ℝ) :
    NJ 5 (16 + 5 * t) =
      2360 * (1 - t) ^ 4 + 4 * 7500 * t * (1 - t) ^ 3 + 6 * (25055 / 2) * t ^ 2 * (1 - t) ^ 2
        + 4 * (254205 / 16) * t ^ 3 * (1 - t) + (31115 / 2) * t ^ 4 := by
  unfold NJ sig34
  ring

/-- `π_0, …, π_4` from the paper. -/
def pi0 (u : ℝ) : ℝ := 121 * u ^ 4 + 2474 * u ^ 3 + 17431 * u ^ 2 + 46014 * u + 25400
def pi1 (u : ℝ) : ℝ :=
  121 * u ^ 5 + 3442 * u ^ 4 + 37967 * u ^ 3 + 199405 * u ^ 2 + 480475 * u + 384510
def pi2 (u : ℝ) : ℝ :=
  121 * u ^ 6 + 4410 * u ^ 5 + 65863 * u ^ 4 + 512764 * u ^ 3 + 2172926 * u ^ 2
    + 4668316 * u + 3829440
def pi3 (u : ℝ) : ℝ :=
  121 * u ^ 5 + 3684 * u ^ 4 + 43735 * u ^ 3 + 249961 * u ^ 2 + 671859 * u + 644400
def pi4 (u : ℝ) : ℝ := 121 * u ^ 4 + 2958 * u ^ 3 + 25579 * u ^ 2 + 88782 * u + 91440

/-- `J ≥ 6`: with `J = u + 6` and `H = J(J+1)/2 + (J+1)t`, the degree-four Bernstein
coefficients of `N_J(H)` in `t` are `β_i = μ_i π_i(u)/2880` (polynomial identity in `u`, `t`). -/
theorem NJ_bernstein (u t : ℝ) :
    let J := u + 6
    NJ J (J * (J + 1) / 2 + (J + 1) * t) =
      (J ^ 2 * (J + 1) ^ 2 * pi0 u / 2880) * (1 - t) ^ 4
      + 4 * (J * (J + 1) ^ 2 * pi1 u / 2880) * t * (1 - t) ^ 3
      + 6 * ((J + 1) ^ 2 * pi2 u / 2880) * t ^ 2 * (1 - t) ^ 2
      + 4 * ((J + 1) ^ 2 * (J + 2) * pi3 u / 2880) * t ^ 3 * (1 - t)
      + ((J + 1) ^ 2 * (J + 2) ^ 2 * pi4 u / 2880) * t ^ 4 := by
  intro J
  simp only [J, NJ, sig34, pi0, pi1, pi2, pi3, pi4]
  ring

/-- `N_J(H) ≥ 0` on the triangular cell (`J ≥ 5`, `H ≥ 16`). -/
theorem NJ_nonneg (J : ℕ) (hJ : 5 ≤ J) (H : ℝ) (h16 : 16 ≤ H)
    (hlo : (J : ℝ) * (J + 1) / 2 ≤ H) (hhi : H ≤ ((J : ℝ) + 1) * (J + 2) / 2) :
    0 ≤ NJ J H := by
  rcases Nat.eq_or_lt_of_le hJ with h5 | h6
  · subst h5
    push_cast at hhi
    set t := (H - 16) / 5 with ht
    have hH : H = 16 + 5 * t := by rw [ht]; ring
    have ht0 : 0 ≤ t := by rw [ht]; linarith
    have ht1 : 0 ≤ 1 - t := by rw [ht]; linarith
    push_cast
    rw [hH, NJ_five_bernstein]
    generalize 1 - t = s at ht1 ⊢
    positivity
  · have hJ6 : (6 : ℝ) ≤ J := by exact_mod_cast h6
    set u : ℝ := (J : ℝ) - 6 with hu
    have hJu : (J : ℝ) = u + 6 := by rw [hu]; ring
    have hu0 : 0 ≤ u := by rw [hu]; linarith
    have hJ1 : (0 : ℝ) < J + 1 := by positivity
    set t := (H - (J : ℝ) * (J + 1) / 2) / (J + 1) with ht
    have hH : H = (J : ℝ) * (J + 1) / 2 + (J + 1) * t := by
      rw [ht]; field_simp; ring
    have ht0 : 0 ≤ t := by rw [ht]; apply div_nonneg _ hJ1.le; linarith
    have ht1 : 0 ≤ 1 - t := by
      rw [ht, sub_nonneg, div_le_one hJ1]; nlinarith
    rw [hH, hJu, NJ_bernstein u t]
    simp only [pi0, pi1, pi2, pi3, pi4]
    generalize 1 - t = s at ht1 ⊢
    positivity

/-- The scalar reduction: `Q(H) ≤ S̃_J T̃_J` on the triangular cell. -/
theorem Q_le_StTt (J : ℕ) (hJ : 5 ≤ J) (H : ℝ) (h16 : 16 ≤ H)
    (hlo : (J : ℝ) * (J + 1) / 2 ≤ H) (hhi : H ≤ ((J : ℝ) + 1) * (J + 2) / 2) :
    (3 * H + 4) * (H + 1) / 4 ≤ St J H * Tt J H := by
  have hH : H ≠ 0 := by positivity
  rw [St_closed J H hH, Tt_closed J H hH]
  have h := NJ_eq J H hH
  have hN := NJ_nonneg J hJ H h16 hlo hhi
  rw [← h] at hN
  have hH2 : 0 < H ^ 2 := by positivity
  have := (mul_nonneg_iff_of_pos_left hH2).mp hN
  linarith

/-- `S ≥ S̃_J` and `T ≥ T̃_J` for `J ≤ K`, `J ≤ H`. -/
theorem ST_ge_tilde (δ : ℝ) (K J : ℕ) (hδ0 : 0 < δ) (hK : ((K : ℝ) + 1) * δ < 1)
    (hJK : J ≤ K) (hJH : (J : ℝ) ≤ 1 / δ - 1) :
    St J (1 / δ - 1) ≤ Ssym δ K ∧ Tt J (1 / δ - 1) ≤ Tsym δ K := by
  have hδ1 : δ < 1 := by
    have : (0 : ℝ) ≤ K := K.cast_nonneg
    nlinarith
  have hLnn : ∀ r ∈ Finset.range K, 0 ≤ L δ (r + 1) := by
    intro r hr
    have hr' : r + 1 ≤ K := Finset.mem_range.mp hr
    have : ((r + 1 : ℕ) : ℝ) ≤ K := by exact_mod_cast hr'
    exact (L_nonneg_le_R δ hδ0 (r + 1) (by nlinarith)).1
  have hlam : ∀ r ∈ Finset.range J, lam (1 / δ - 1) r ≤ L δ (r + 1) := by
    intro r hr
    have hr' : r + 1 ≤ J := Finset.mem_range.mp hr
    have : ((r + 1 : ℕ) : ℝ) ≤ J := by exact_mod_cast hr'
    have := L_ge_lambda δ hδ0 hδ1 (r + 1) (by linarith)
    simpa [lam] using this
  have hsub : Finset.range J ⊆ Finset.range K := Finset.range_subset_range.mpr hJK
  constructor
  · unfold St Ssym
    have h1 : ∑ r ∈ Finset.range J, lam (1 / δ - 1) r ≤ ∑ r ∈ Finset.range J, L δ (r + 1) :=
      Finset.sum_le_sum hlam
    have h2 : ∑ r ∈ Finset.range J, L δ (r + 1) ≤ ∑ r ∈ Finset.range K, L δ (r + 1) :=
      Finset.sum_le_sum_of_subset_of_nonneg hsub (fun r hr _ => hLnn r hr)
    linarith
  · unfold Tt Tsym
    have h1 : ∑ r ∈ Finset.range J, ((r : ℝ) + 1) ^ 2 * lam (1 / δ - 1) r ≤
        ∑ r ∈ Finset.range J, ((r : ℝ) + 1) ^ 2 * L δ (r + 1) :=
      Finset.sum_le_sum (fun r hr => mul_le_mul_of_nonneg_left (hlam r hr) (sq_nonneg _))
    have h2 : ∑ r ∈ Finset.range J, ((r : ℝ) + 1) ^ 2 * L δ (r + 1) ≤
        ∑ r ∈ Finset.range K, ((r : ℝ) + 1) ^ 2 * L δ (r + 1) :=
      Finset.sum_le_sum_of_subset_of_nonneg hsub
        (fun r hr _ => mul_nonneg (sq_nonneg _) (hLnn r hr))
    linarith

/-- `S̃_J ≥ 0` and `T̃_J ≥ 0` when `J(J+1)/2 ≤ H`. -/
theorem StTt_nonneg (J : ℕ) (H : ℝ) (hH : 0 < H) (hlo : (J : ℝ) * (J + 1) / 2 ≤ H) :
    0 ≤ St J H ∧ 0 ≤ Tt J H := by
  have hlam : ∀ r ∈ Finset.range J, 0 ≤ lam H r := by
    intro r hr
    have hr' : r + 1 ≤ J := Finset.mem_range.mp hr
    have hrJ : (r : ℝ) + 1 ≤ J := by exact_mod_cast hr'
    unfold lam
    rw [sub_nonneg, div_le_one (by positivity)]
    have : ((r : ℝ) + 1) * ((r : ℝ) + 1 + 1) ≤ (J : ℝ) * (J + 1) :=
      mul_le_mul hrJ (by linarith) (by positivity) (by positivity)
    linarith
  constructor
  · unfold St
    have := Finset.sum_nonneg hlam
    linarith
  · unfold Tt
    have := Finset.sum_nonneg (fun (r : ℕ) hr => mul_nonneg (sq_nonneg ((r : ℝ) + 1)) (hlam r hr))
    linarith

/-- Existence of the triangular cell index `J ≥ 5` for `H ≥ 16`. -/
theorem exists_triangular_cell (H : ℝ) (h16 : 16 ≤ H) :
    ∃ J : ℕ, 5 ≤ J ∧ (J : ℝ) * (J + 1) / 2 ≤ H ∧ H ≤ ((J : ℝ) + 1) * (J + 2) / 2 := by
  classical
  have hex : ∃ n : ℕ, H < ((n : ℝ) + 1) * ((n : ℝ) + 2) / 2 := by
    obtain ⟨n, hn⟩ := exists_nat_gt H
    refine ⟨n, ?_⟩
    have : (0 : ℝ) ≤ n := n.cast_nonneg
    nlinarith
  set J := Nat.find hex with hJ
  have hspec : H < ((J : ℝ) + 1) * ((J : ℝ) + 2) / 2 := Nat.find_spec hex
  have hJ5 : 5 ≤ J := by
    by_contra hcon
    have : (J : ℝ) ≤ 4 := by exact_mod_cast (by omega : J ≤ 4)
    have : (0 : ℝ) ≤ J := J.cast_nonneg
    nlinarith
  refine ⟨J, hJ5, ?_, hspec.le⟩
  have hmin := Nat.find_min hex (show J - 1 < J by omega)
  push_neg at hmin
  have hc : ((J - 1 : ℕ) : ℝ) = (J : ℝ) - 1 := by
    rw [Nat.cast_sub (by omega)]; simp
  rw [hc] at hmin
  linarith [hmin, show ((J : ℝ) - 1 + 1) * ((J : ℝ) - 1 + 2) / 2 = (J : ℝ) * (J + 1) / 2 by ring]

/-- The target in terms of `H = 1/δ - 1`: `(3+δ)/(4δ²) = (3H+4)(H+1)/4`. -/
theorem target_eq_Q (δ : ℝ) (hδ0 : 0 < δ) :
    (3 + δ) / (4 * δ ^ 2) = (3 * (1 / δ - 1) + 4) * ((1 / δ - 1) + 1) / 4 := by
  field_simp
  ring

/-- **G2** in its final form (same statement as `scalar_inequality_large`). -/
theorem large_main (δ : ℝ) (K : ℕ)
    (hδ0 : 0 < δ) (hδ17 : δ ≤ 1 / 17)
    (hK : ((K : ℝ) + 1) * δ < 1) (hK' : 1 ≤ ((K : ℝ) + 2) * δ) :
    (3 + δ) / (4 * δ ^ 2) ≤ A δ K := by
  set H := 1 / δ - 1 with hHdef
  have h16 : 16 ≤ H := by
    rw [hHdef, le_sub_iff_add_le, le_div_iff₀ hδ0]; linarith
  have hHK : H ≤ (K : ℝ) + 1 := by
    rw [hHdef, sub_le_iff_le_add, div_le_iff₀ hδ0]; linarith
  obtain ⟨J, hJ5, hlo, hhi⟩ := exists_triangular_cell H h16
  have hJ5r : (5 : ℝ) ≤ J := by exact_mod_cast hJ5
  have hJH : (J : ℝ) ≤ H := by nlinarith
  have hJK : J ≤ K := by
    have : (J : ℝ) ≤ K := by nlinarith
    exact_mod_cast this
  obtain ⟨hS, hT⟩ := ST_ge_tilde δ K J hδ0 hK hJK hJH
  obtain ⟨hS0, hT0⟩ := StTt_nonneg J H (by linarith) hlo
  have hQ := Q_le_StTt J hJ5 H h16 hlo hhi
  rw [target_eq_Q δ hδ0]
  calc (3 * H + 4) * (H + 1) / 4 ≤ St J H * Tt J H := hQ
    _ ≤ Ssym δ K * Tsym δ K := mul_le_mul hS hT hT0 (hS0.trans hS)
    _ ≤ A δ K := A_ge_ST δ K hδ0 hK

end PBScalar
