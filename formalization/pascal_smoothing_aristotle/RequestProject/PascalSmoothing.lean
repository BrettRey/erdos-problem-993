import Mathlib

/-!
# A universal Pascal-smoothing lemma for erasure shadows

Formalization of `pascal_smoothing_shadow_lemma_2026-08-20.md`.

Let `Δ` be a downward-closed family of subsets of a set of `M` labelled coordinates
(`Fin M`), with `a j = faceCount Δ j` the number of its `j`-faces.  Put

  `E d = ∑ j, a j * (M - j).choose d`.

The main result (Pascal-smoothing lemma) is, for `1 ≤ d ≤ M - 1` and `m = M - d`,

  `E d ^ 2 / (E (d-1) * E (d+1)) ≥ (8/9) * ((d+1)*(m+1)) / (d*m)`,

stated with cleared denominators over `ℕ`.  Its consequences (log-concavity through
defect depth eight, log-concavity for `M ≤ 33`, and the quantitative depth-three
reserve) are derived at the end of the file.
-/

namespace PascalSmoothing

open Finset
open scoped FinsetFamily

variable {M : ℕ}

/-- `faceCount Δ j` is the number of `j`-element members ("`j`-faces") of the family `Δ`. -/
def faceCount (Δ : Finset (Finset (Fin M))) (j : ℕ) : ℕ :=
  (Δ.filter (fun s => s.card = j)).card

/-- The erasure shadow profile `E d = ∑ j a_j * C(M - j, d)` of a family on `M` coordinates. -/
def erasureShadow (Δ : Finset (Finset (Fin M))) (d : ℕ) : ℕ :=
  ∑ j ∈ range (M + 1), faceCount Δ j * (M - j).choose d

/-- The face numbers normalized by the Boolean layers: `b j = a j / C(M, j)`. -/
noncomputable def layerRatio (Δ : Finset (Finset (Fin M))) (j : ℕ) : ℚ :=
  (faceCount Δ j : ℚ) / (M.choose j : ℚ)

/-- The Pascal transform `S m = ∑ j C(m, j) * b j` of a sequence `b`. -/
noncomputable def pascalTransform (b : ℕ → ℚ) (m : ℕ) : ℚ :=
  ∑ j ∈ range (m + 1), (m.choose j : ℚ) * b j

/-- The shifted Pascal transform `U k r = ∑ j C(k, j) * b (j + r)`. -/
noncomputable def shiftedTransform (b : ℕ → ℚ) (k r : ℕ) : ℚ :=
  ∑ j ∈ range (k + 1), (k.choose j : ℚ) * b (j + r)

/-! ### The abstract Pascal-smoothing inequality for decreasing sequences -/

lemma shiftedTransform_zero (b : ℕ → ℚ) (k : ℕ) :
    shiftedTransform b k 0 = pascalTransform b k := by
  simp [shiftedTransform, pascalTransform]

lemma shiftedTransform_succ (b : ℕ → ℚ) (k r : ℕ) :
    shiftedTransform b (k + 1) r = shiftedTransform b k r + shiftedTransform b k (r + 1) := by
  unfold shiftedTransform
  rw [Finset.sum_range_succ' (fun j => ((k + 1).choose j : ℚ) * b (j + r)) (k + 1)]
  rw [Finset.sum_range_succ' (fun j => (k.choose j : ℚ) * b (j + r)) k]
  have hcongr : ∀ i ∈ range (k + 1), ((k + 1).choose (i + 1) : ℚ) * b (i + 1 + r)
      = (k.choose i : ℚ) * b (i + (r + 1)) + (k.choose (i + 1) : ℚ) * b (i + 1 + r) := by
    intro i _
    rw [Nat.choose_succ_succ' k i]
    push_cast
    have h : i + 1 + r = i + (r + 1) := by omega
    rw [h]; ring
  rw [Finset.sum_congr rfl hcongr, Finset.sum_add_distrib]
  have h2 : ∑ i ∈ range (k + 1), (k.choose (i + 1) : ℚ) * b (i + 1 + r)
      = ∑ i ∈ range k, (k.choose (i + 1) : ℚ) * b (i + 1 + r) := by
    rw [Finset.sum_range_succ]
    simp
  rw [h2]
  simp only [Nat.choose_zero_right, Nat.cast_one, Nat.zero_add, one_mul]
  ring

lemma shiftedTransform_nonneg {b : ℕ → ℚ} (hb : ∀ j, 0 ≤ b j) (k r : ℕ) :
    0 ≤ shiftedTransform b k r :=
  Finset.sum_nonneg fun _ _ => mul_nonneg (by positivity) (hb _)

lemma shiftedTransform_antitone {b : ℕ → ℚ} (hb : ∀ j, b (j + 1) ≤ b j) (k r : ℕ) :
    shiftedTransform b k (r + 1) ≤ shiftedTransform b k r := by
  refine Finset.sum_le_sum fun j _ => ?_
  have hstep : b (j + (r + 1)) ≤ b (j + r) := by
    simpa [Nat.add_assoc] using hb (j + r)
  exact mul_le_mul_of_nonneg_left hstep (by positivity)

/-- The core inequality: three adjacent Pascal transforms of a nonnegative decreasing
sequence satisfy `9 S_m ^ 2 ≥ 8 S_{m-1} S_{m+1}`. -/
lemma pascalTransform_smoothing {b : ℕ → ℚ} (hb0 : ∀ j, 0 ≤ b j)
    (hb : ∀ j, b (j + 1) ≤ b j) (k : ℕ) :
    8 * (pascalTransform b k * pascalTransform b (k + 2)) ≤ 9 * pascalTransform b (k + 1) ^ 2 := by
  set A := shiftedTransform b k 0 with hA
  set B := shiftedTransform b k 1 with hB
  set C := shiftedTransform b k 2 with hC
  have h0 : pascalTransform b k = A := (shiftedTransform_zero b k).symm
  have h1 : pascalTransform b (k + 1) = A + B := by
    rw [← shiftedTransform_zero, shiftedTransform_succ]
  have h2 : pascalTransform b (k + 2) = A + 2 * B + C := by
    rw [← shiftedTransform_zero, show k + 2 = (k + 1) + 1 from rfl, shiftedTransform_succ,
      shiftedTransform_succ, shiftedTransform_succ]
    ring
  have hA0 : 0 ≤ A := shiftedTransform_nonneg hb0 k 0
  have hCB : C ≤ B := shiftedTransform_antitone hb k 1
  rw [h0, h1, h2]
  nlinarith [sq_nonneg (A - 3 * B), mul_nonneg hA0 (sub_nonneg.2 hCB)]

/-! ### Local LYM for the normalized face numbers -/

/-- Local LYM: `(j+1) * a_{j+1} ≤ (M - j) * a_j` for a downward-closed family. -/
lemma faceCount_local_lym (Δ : Finset (Finset (Fin M)))
    (hdc : ∀ s ∈ Δ, ∀ t ⊆ s, t ∈ Δ) (j : ℕ) :
    (j + 1) * faceCount Δ (j + 1) ≤ (M - j) * faceCount Δ j := by
  rcases le_or_gt (j + 1) M with hj | hj
  · set 𝒜 := Δ.filter (fun s => s.card = j + 1) with h𝒜def
    have hsized : ((𝒜 : Finset (Finset (Fin M))) : Set (Finset (Fin M))).Sized (j + 1) := by
      intro s hs
      simp only [h𝒜def, Finset.coe_filter, Set.mem_setOf_eq] at hs
      exact hs.2
    have hlym := Finset.card_mul_le_card_shadow_mul hsized
    have hsub : ∂ 𝒜 ⊆ Δ.filter (fun s => s.card = j) := by
      intro t ht
      rw [Finset.mem_shadow_iff] at ht
      obtain ⟨s, hs, a, ha, rfl⟩ := ht
      simp only [h𝒜def, Finset.mem_filter] at hs ⊢
      refine ⟨hdc s hs.1 _ (Finset.erase_subset _ _), ?_⟩
      rw [Finset.card_erase_of_mem ha, hs.2]
      rfl
    have hcard : (∂ 𝒜).card ≤ faceCount Δ j := Finset.card_le_card hsub
    have hM : Fintype.card (Fin M) - (j + 1) + 1 = M - j := by
      simp only [Fintype.card_fin]; omega
    rw [hM] at hlym
    calc (j + 1) * faceCount Δ (j + 1) = 𝒜.card * (j + 1) := by rw [faceCount]; ring
      _ ≤ (∂ 𝒜).card * (M - j) := hlym
      _ ≤ faceCount Δ j * (M - j) := Nat.mul_le_mul_right _ hcard
      _ = (M - j) * faceCount Δ j := by ring
  · have hzero : faceCount Δ (j + 1) = 0 := by
      rw [faceCount, Finset.card_eq_zero, Finset.filter_eq_empty_iff]
      intro s _ hs
      have hcard := Finset.card_le_univ s
      simp only [Fintype.card_fin, hs] at hcard
      omega
    simp [hzero]

lemma faceCount_mul_choose_le (Δ : Finset (Finset (Fin M)))
    (hdc : ∀ s ∈ Δ, ∀ t ⊆ s, t ∈ Δ) (j : ℕ) :
    faceCount Δ (j + 1) * M.choose j ≤ faceCount Δ j * M.choose (j + 1) := by
  have h := faceCount_local_lym Δ hdc j
  have key : (j + 1) * (faceCount Δ (j + 1) * M.choose j)
      ≤ (j + 1) * (faceCount Δ j * M.choose (j + 1)) := by
    calc (j + 1) * (faceCount Δ (j + 1) * M.choose j)
        = ((j + 1) * faceCount Δ (j + 1)) * M.choose j := by ring
      _ ≤ ((M - j) * faceCount Δ j) * M.choose j := Nat.mul_le_mul_right _ h
      _ = faceCount Δ j * (M.choose j * (M - j)) := by ring
      _ = faceCount Δ j * (M.choose (j + 1) * (j + 1)) := by rw [Nat.choose_succ_right_eq]
      _ = (j + 1) * (faceCount Δ j * M.choose (j + 1)) := by ring
  exact Nat.le_of_mul_le_mul_left key (Nat.succ_pos j)

lemma layerRatio_nonneg (Δ : Finset (Finset (Fin M))) (j : ℕ) : 0 ≤ layerRatio Δ j := by
  unfold layerRatio; positivity

lemma layerRatio_antitone (Δ : Finset (Finset (Fin M)))
    (hdc : ∀ s ∈ Δ, ∀ t ⊆ s, t ∈ Δ) (j : ℕ) :
    layerRatio Δ (j + 1) ≤ layerRatio Δ j := by
  have h := faceCount_mul_choose_le Δ hdc j
  unfold layerRatio
  rcases Nat.eq_zero_or_pos (M.choose (j + 1)) with h1 | h1
  · rw [h1]
    simp only [Nat.cast_zero, div_zero]
    positivity
  · have hj : j + 1 ≤ M := by
      by_contra hc
      rw [Nat.choose_eq_zero_of_lt (by omega)] at h1
      omega
    have h2 : 0 < M.choose j := Nat.choose_pos (by omega)
    rw [div_le_div_iff₀ (by exact_mod_cast h1) (by exact_mod_cast h2)]
    exact_mod_cast h

/-! ### The absorption identity `E_d = C(M,d) * S_{M-d}` -/

lemma choose_mul_choose_sub (M j d : ℕ) :
    M.choose j * (M - j).choose d = M.choose d * (M - d).choose j := by
  have h1 := Nat.choose_mul (n := M) (k := j + d) (s := j) (by omega)
  have h2 := Nat.choose_mul (n := M) (k := j + d) (s := d) (by omega)
  have h3 : (j + d).choose j = (j + d).choose d := by
    rw [← Nat.choose_symm (by omega : d ≤ j + d)]
    congr 1
    omega
  simp only [Nat.add_sub_cancel_left, Nat.add_sub_cancel] at h1 h2
  rw [h3] at h1
  rw [← h1, h2]

lemma erasureShadow_eq (Δ : Finset (Finset (Fin M))) {d : ℕ} (hd : d ≤ M) :
    (erasureShadow Δ d : ℚ) = (M.choose d : ℚ) * pascalTransform (layerRatio Δ) (M - d) := by
  unfold erasureShadow layerRatio
  push_cast
  rw [pascalTransform, Finset.mul_sum]
  rw [← Finset.sum_subset (s₁ := range (M - d + 1)) (s₂ := range (M + 1))
    (by intro x hx; simp only [Finset.mem_range] at *; omega)
    (by
      intro x hx hx2
      simp only [Finset.mem_range] at hx hx2
      have hz : (M - x).choose d = 0 := Nat.choose_eq_zero_of_lt (by omega)
      rw [hz]; simp)]
  refine Finset.sum_congr rfl fun j hj => ?_
  simp only [Finset.mem_range] at hj
  have hjM : j ≤ M := by omega
  have hpos : (0 : ℚ) < (M.choose j : ℚ) := by exact_mod_cast Nat.choose_pos hjM
  have habs : ((M.choose j : ℚ)) * ((M - j).choose d : ℚ)
      = (M.choose d : ℚ) * ((M - d).choose j : ℚ) := by
    exact_mod_cast congrArg (fun n : ℕ => (n : ℚ)) (choose_mul_choose_sub M j d)
  field_simp
  nlinarith [habs]

/-! ### The Pascal-smoothing lemma -/

lemma choose_mul_identity (e k : ℕ) :
    (e + 2) * (k + 2) * ((e + k + 2).choose e * (e + k + 2).choose (e + 2))
      = (e + 1) * (k + 1) * (e + k + 2).choose (e + 1) ^ 2 := by
  have he : (e + k + 2) - e = k + 2 := by omega
  have he2 : (e + k + 2) - (e + 1) = k + 1 := by omega
  have h1 := Nat.choose_succ_right_eq (e + k + 2) e
  have h2 := Nat.choose_succ_right_eq (e + k + 2) (e + 1)
  rw [he] at h1
  rw [he2] at h2
  zify at h1 h2 ⊢
  linear_combination ((((e + k + 2).choose (e + 1) : ℕ) : ℤ) * ((e : ℤ) + 1)) * h2
    - ((((e + k + 2).choose (e + 2) : ℕ) : ℤ) * ((e : ℤ) + 2)) * h1

/-- The Pascal-smoothing lemma, in the parametrized form `M = e + k + 2`, `d = e + 1`,
`m = M - d = k + 1`. -/
theorem pascal_smoothing_core (e k : ℕ) (Δ : Finset (Finset (Fin (e + k + 2))))
    (hdc : ∀ s ∈ Δ, ∀ t ⊆ s, t ∈ Δ) :
    8 * ((e + 2) * (k + 2)) * (erasureShadow Δ e * erasureShadow Δ (e + 2))
      ≤ 9 * ((e + 1) * (k + 1)) * erasureShadow Δ (e + 1) ^ 2 := by
  set T := pascalTransform (layerRatio Δ) with hT
  have hb0 : ∀ j, 0 ≤ layerRatio Δ j := layerRatio_nonneg Δ
  have hb : ∀ j, layerRatio Δ (j + 1) ≤ layerRatio Δ j := layerRatio_antitone Δ hdc
  have hE0 : (erasureShadow Δ e : ℚ) = ((e + k + 2).choose e : ℚ) * T (k + 2) := by
    rw [erasureShadow_eq Δ (by omega : e ≤ e + k + 2)]
    congr 2
    omega
  have hE1 : (erasureShadow Δ (e + 1) : ℚ) = ((e + k + 2).choose (e + 1) : ℚ) * T (k + 1) := by
    rw [erasureShadow_eq Δ (by omega : e + 1 ≤ e + k + 2)]
    congr 2
    omega
  have hE2 : (erasureShadow Δ (e + 2) : ℚ) = ((e + k + 2).choose (e + 2) : ℚ) * T k := by
    rw [erasureShadow_eq Δ (by omega : e + 2 ≤ e + k + 2)]
    congr 2
    omega
  have hcore : 8 * (T k * T (k + 2)) ≤ 9 * T (k + 1) ^ 2 :=
    pascalTransform_smoothing hb0 hb k
  have hid : ((e + 2) * (k + 2) * ((e + k + 2).choose e * (e + k + 2).choose (e + 2)) : ℚ)
      = ((e + 1) * (k + 1) * ((e + k + 2).choose (e + 1) : ℚ) ^ 2) := by
    exact_mod_cast congrArg (fun n : ℕ => (n : ℚ)) (choose_mul_identity e k)
  rw [← Nat.cast_le (α := ℚ)]
  push_cast
  rw [hE0, hE1, hE2]
  have hP : (0 : ℚ) ≤ (e + 1) * (k + 1) * ((e + k + 2).choose (e + 1) : ℚ) ^ 2 := by positivity
  calc 8 * ((e + 2) * (k + 2))
        * (((e + k + 2).choose e : ℚ) * T (k + 2) * (((e + k + 2).choose (e + 2) : ℚ) * T k))
      = ((e + 2) * (k + 2) * (((e + k + 2).choose e : ℚ) * ((e + k + 2).choose (e + 2) : ℚ)))
          * (8 * (T k * T (k + 2))) := by ring
    _ = ((e + 1) * (k + 1) * ((e + k + 2).choose (e + 1) : ℚ) ^ 2) * (8 * (T k * T (k + 2))) := by
        rw [hid]
    _ ≤ ((e + 1) * (k + 1) * ((e + k + 2).choose (e + 1) : ℚ) ^ 2) * (9 * T (k + 1) ^ 2) :=
        mul_le_mul_of_nonneg_left hcore hP
    _ = 9 * ((e + 1) * (k + 1)) * (((e + k + 2).choose (e + 1) : ℚ) * T (k + 1)) ^ 2 := by ring

/-- **Pascal-smoothing lemma.**  For any nonempty downward-closed family `Δ` on `M` labelled
coordinates and any `1 ≤ d ≤ M - 1`, writing `m = M - d`,
`E d ^ 2 / (E (d-1) * E (d+1)) ≥ (8/9) * (d+1)(m+1) / (dm)`, stated with cleared denominators.
(The nonemptiness hypothesis `hΔ`, present in the source statement, is not needed for this
inequality; it is kept for faithfulness.) -/
theorem pascal_smoothing_shadow_lemma (Δ : Finset (Finset (Fin M))) (hΔ : Δ.Nonempty)
    (hdc : ∀ s ∈ Δ, ∀ t ⊆ s, t ∈ Δ) {d : ℕ} (hd1 : 1 ≤ d) (hd2 : d ≤ M - 1) :
    8 * ((d + 1) * (M - d + 1)) * (erasureShadow Δ (d - 1) * erasureShadow Δ (d + 1))
      ≤ 9 * (d * (M - d)) * erasureShadow Δ d ^ 2 := by
  obtain ⟨e, rfl⟩ : ∃ e, d = e + 1 := ⟨d - 1, by omega⟩
  obtain ⟨k, rfl⟩ : ∃ k, M = e + k + 2 := ⟨M - e - 2, by omega⟩
  have h1 : e + 1 + 1 = e + 2 := by omega
  have h2 : e + k + 2 - (e + 1) + 1 = k + 2 := by omega
  have h3 : e + 1 - 1 = e := by omega
  have h4 : e + k + 2 - (e + 1) = k + 1 := by omega
  rw [h1, h2, h3, h4]
  exact pascal_smoothing_core e k Δ hdc

/-- Log-concavity of the erasure shadow profile under the numerical condition
`d (M - d) ≤ 8 (M + 1)`. -/
theorem erasureShadow_log_concave_of (Δ : Finset (Finset (Fin M))) (hΔ : Δ.Nonempty)
    (hdc : ∀ s ∈ Δ, ∀ t ⊆ s, t ∈ Δ) {d : ℕ} (hd1 : 1 ≤ d) (hd2 : d ≤ M - 1)
    (hnum : d * (M - d) ≤ 8 * (M + 1)) :
    erasureShadow Δ (d - 1) * erasureShadow Δ (d + 1) ≤ erasureShadow Δ d ^ 2 := by
  have hmain := pascal_smoothing_shadow_lemma Δ hΔ hdc hd1 hd2
  have hM : d + 1 ≤ M := by omega
  obtain ⟨m, hm⟩ : ∃ m, M - d = m := ⟨M - d, rfl⟩
  rw [hm] at hmain hnum
  have hm1 : 1 ≤ m := by omega
  have hsum : d + m = M := by omega
  have hcoef : 9 * (d * m) ≤ 8 * ((d + 1) * (m + 1)) := by nlinarith [hnum, hsum]
  have hpos : 0 < 9 * (d * m) := by positivity
  refine Nat.le_of_mul_le_mul_left ?_ hpos
  calc 9 * (d * m) * (erasureShadow Δ (d - 1) * erasureShadow Δ (d + 1))
      ≤ 8 * ((d + 1) * (m + 1)) * (erasureShadow Δ (d - 1) * erasureShadow Δ (d + 1)) :=
        Nat.mul_le_mul_right _ hcoef
    _ ≤ 9 * (d * m) * erasureShadow Δ d ^ 2 := hmain

/-- Every upper-tail log-concavity inequality through defect depth eight holds. -/
theorem erasureShadow_log_concave_depth_le_eight (Δ : Finset (Finset (Fin M))) (hΔ : Δ.Nonempty)
    (hdc : ∀ s ∈ Δ, ∀ t ⊆ s, t ∈ Δ) {d : ℕ} (hd1 : 1 ≤ d) (hd2 : d ≤ M - 1) (hd8 : d ≤ 8) :
    erasureShadow Δ (d - 1) * erasureShadow Δ (d + 1) ≤ erasureShadow Δ d ^ 2 := by
  refine erasureShadow_log_concave_of Δ hΔ hdc hd1 hd2 ?_
  calc d * (M - d) ≤ 8 * (M - d) := Nat.mul_le_mul_right _ hd8
    _ ≤ 8 * (M + 1) := Nat.mul_le_mul_left _ (by omega)

/-- The whole coefficient sequence is log-concave whenever `M ≤ 33`. -/
theorem erasureShadow_log_concave_of_M_le_33 (Δ : Finset (Finset (Fin M))) (hΔ : Δ.Nonempty)
    (hdc : ∀ s ∈ Δ, ∀ t ⊆ s, t ∈ Δ) (hM : M ≤ 33) {d : ℕ} (hd1 : 1 ≤ d) (hd2 : d ≤ M - 1) :
    erasureShadow Δ (d - 1) * erasureShadow Δ (d + 1) ≤ erasureShadow Δ d ^ 2 := by
  refine erasureShadow_log_concave_of Δ hΔ hdc hd1 hd2 ?_
  have key : ∀ x N : ℕ, 4 * x ≤ N * N → N ≤ 33 → x ≤ 8 * (N + 1) := by
    intro x N h hN
    interval_cases N <;> omega
  refine key _ M ?_ hM
  obtain ⟨m, hm⟩ : ∃ m, M - d = m := ⟨M - d, rfl⟩
  have hsum : d + m = M := by omega
  rw [hm, ← hsum]
  zify
  nlinarith [sq_nonneg ((d : ℤ) - m)]

/-- The depth-three quantitative reserve:
`27 (M - 3) E_3 ^ 2 ≥ 32 (M - 2) E_2 E_4`. -/
theorem erasureShadow_depth_three (Δ : Finset (Finset (Fin M))) (hΔ : Δ.Nonempty)
    (hdc : ∀ s ∈ Δ, ∀ t ⊆ s, t ∈ Δ) (hM : 4 ≤ M) :
    32 * (M - 2) * (erasureShadow Δ 2 * erasureShadow Δ 4)
      ≤ 27 * (M - 3) * erasureShadow Δ 3 ^ 2 := by
  have hmain := pascal_smoothing_shadow_lemma Δ hΔ hdc (d := 3) (by omega) (by omega)
  have h1 : M - 3 + 1 = M - 2 := by omega
  rw [h1] at hmain
  calc 32 * (M - 2) * (erasureShadow Δ 2 * erasureShadow Δ 4)
      = 8 * ((3 + 1) * (M - 2)) * (erasureShadow Δ (3 - 1) * erasureShadow Δ (3 + 1)) := by
        simp only [show (3 : ℕ) - 1 = 2 from rfl, show (3 : ℕ) + 1 = 4 from rfl]
        ring
    _ ≤ 9 * (3 * (M - 3)) * erasureShadow Δ 3 ^ 2 := hmain
    _ = 27 * (M - 3) * erasureShadow Δ 3 ^ 2 := by ring

/-- The depth-three log-concavity, strict whenever `E_2 E_4 > 0`. -/
theorem erasureShadow_depth_three_strict (Δ : Finset (Finset (Fin M))) (hΔ : Δ.Nonempty)
    (hdc : ∀ s ∈ Δ, ∀ t ⊆ s, t ∈ Δ) (hM : 4 ≤ M)
    (hpos : 0 < erasureShadow Δ 2 * erasureShadow Δ 4) :
    erasureShadow Δ 2 * erasureShadow Δ 4 < erasureShadow Δ 3 ^ 2 := by
  have hmain := erasureShadow_depth_three Δ hΔ hdc hM
  obtain ⟨q, rfl⟩ : ∃ q, M = q + 4 := ⟨M - 4, by omega⟩
  have h1 : q + 4 - 2 = q + 2 := by omega
  have h2 : q + 4 - 3 = q + 1 := by omega
  rw [h1, h2] at hmain
  by_contra hcon
  push_neg at hcon
  nlinarith [hmain, hpos, hcon]

/-! ### The generating-function description of the profile -/

section GenFun

open Polynomial

/-- The polynomial `F(t) = (1+t)^M A_Δ(t/(1+t)) = ∑ j a_j t^j (1+t)^{M-j}`, where
`A_Δ(u) = ∑ j a_j u^j` is the face-generating polynomial of `Δ`. -/
noncomputable def shadowGenFun (Δ : Finset (Finset (Fin M))) : ℚ[X] :=
  ∑ j ∈ range (M + 1), C (faceCount Δ j : ℚ) * ((1 + X) ^ (M - j) * X ^ j)

/-- `E_d` is the coefficient of `t^{M-d}` in `F(t)`. -/
theorem erasureShadow_eq_coeff (Δ : Finset (Finset (Fin M))) {d : ℕ} (hd : d ≤ M) :
    (shadowGenFun Δ).coeff (M - d) = (erasureShadow Δ d : ℚ) := by
  unfold shadowGenFun erasureShadow
  rw [Polynomial.finset_sum_coeff]
  push_cast
  refine Finset.sum_congr rfl fun j hj => ?_
  simp only [Finset.mem_range] at hj
  rw [← mul_assoc, Polynomial.coeff_mul_X_pow']
  rcases le_or_gt j (M - d) with h | h
  · rw [if_pos h, Polynomial.coeff_C_mul, Polynomial.coeff_one_add_X_pow]
    have hsub : M - d - j = (M - j) - d := by omega
    rw [hsub, Nat.choose_symm (by omega : d ≤ M - j)]
  · rw [if_neg (by omega)]
    have hzero : (M - j).choose d = 0 := Nat.choose_eq_zero_of_lt (by omega)
    rw [hzero]
    simp

end GenFun

end PascalSmoothing
