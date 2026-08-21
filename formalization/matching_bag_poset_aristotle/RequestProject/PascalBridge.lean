import RequestProject.Antichain
import RequestProject.PascalSmoothing

/-!
# Bridge: the poset erasure profile is the erasure shadow of a downward-closed family

This file connects Section 3 of `matching_bag_poset_reduction_2026-08-20.md` with the
universal Pascal-smoothing lemma of `RequestProject/PascalSmoothing.lean`, and thereby
proves equations (7) and (7a) of the note.

Fix a finite poset `P` with `r = |P|`, a number `c` of constant coordinates, and put
`M = r + c`.  The erasure profile of the note is

  `E d = [t^{M-d}] ((1+t)^c I(B(P);t))`,

which by formula (10) equals `∑_j a_j C(M-j, d)` with `a_j = antichainCount P j`.  We
record this natural-number sequence as `MatchingBag.erasureProfile P c d`.

The Pascal development is stated for downward-closed families of subsets of `Fin M`.
Since `P` is an arbitrary finite type, we *transport*: choosing an embedding
`P ↪ Fin M` (which exists because `r ≤ M`), the antichains of `P` are carried to a
downward-closed family `pascalFamily P c` of subsets of `Fin M` with the same face
numbers.  The transport is proved, not assumed: `faceCount_pascalFamily` and
`erasureShadow_pascalFamily` establish that `E` is *precisely* the erasure shadow of
`pascalFamily P c`.

Main results:

* `antichainsOf_downward_closed` — the antichains form a downward-closed family;
* `erasureProfile_eq_coeff` — `E d` is the coefficient of `t^{M-d}` of `(1+t)^c I(B(P);t)`;
* `erasureProfile_eq_sum` — `E d = ∑_j a_j C(M-j,d)` (formula (10), in `ℕ`);
* `erasureShadow_pascalFamily` — `E` is the erasure shadow of a downward-closed family;
* `erasureProfile_pascal_smoothing` — the denominator-cleared Pascal-smoothing inequality;
* `erasureProfile_log_concave_depth_le_eight`, `erasureProfile_log_concave_of_M_le_33` —
  log-concavity through defect depth eight, and at every interior defect for `M ≤ 33`;
* `erasureProfile_depth_three` (7a), `erasureProfile_depth_three_strict` (7).
-/

open scoped BigOperators
open Finset Polynomial

namespace MatchingBag

variable {P : Type*} [Fintype P] [DecidableEq P] [PartialOrder P]
  [DecidableRel ((· ≤ ·) : P → P → Prop)]

/-! ### The erasure profile as a natural-number sequence -/

/-- The erasure profile `E_d = ∑_j a_j C(M - j, d)` with `M = |P| + c`, as a natural number.
`erasureProfile_eq_coeff` identifies it with `[t^{M-d}]((1+t)^c I(B(P);t))`. -/
def erasureProfile (P : Type*) [Fintype P] [DecidableEq P] [PartialOrder P]
    [DecidableRel ((· ≤ ·) : P → P → Prop)] (c d : ℕ) : ℕ :=
  ∑ j ∈ range (Fintype.card P + 1), antichainCount P j * (Fintype.card P + c - j).choose d

/-- **Formula (10)** over `ℕ`: `E_d = ∑_j a_j C(M-j, d)`. -/
theorem erasureProfile_eq_sum (c d : ℕ) :
    erasureProfile P c d
      = ∑ j ∈ range (Fintype.card P + 1),
          antichainCount P j * (Fintype.card P + c - j).choose d := rfl

/-- The erasure profile is the coefficient of `t^{M-d}` in `(1+t)^c I(B(P);t)`. -/
theorem erasureProfile_eq_coeff (c : ℕ) {d : ℕ} (hd : d ≤ Fintype.card P + c) :
    (((1 : Polynomial ℤ) + X) ^ c * indepPoly P).coeff (Fintype.card P + c - d)
      = (erasureProfile P c d : ℤ) := by
  rw [coeff_erasure c d hd, erasureProfile]
  push_cast
  rfl

/-! ### The antichains form a downward-closed family -/

/-- Any subset of an antichain is an antichain: the antichains of `P` are downward closed. -/
theorem antichainsOf_downward_closed :
    ∀ A ∈ antichainsOf P, ∀ B ⊆ A, B ∈ antichainsOf P := by
  intro A hA B hBA
  rw [antichainsOf, Finset.mem_filter] at hA ⊢
  exact ⟨Finset.mem_univ _, fun i hi j hj hij => hA.2 i (hBA hi) j (hBA hj) hij⟩

lemma antichainCount_eq_zero_of_lt {j : ℕ} (hj : Fintype.card P < j) :
    antichainCount P j = 0 := by
  rw [antichainCount, Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  intro A _ hcard
  have h := Finset.card_le_univ A
  rw [hcard] at h
  omega

/-! ### Transport to a family of subsets of `Fin M` -/

/-- An embedding of the poset `P` into the `M = |P| + c` labelled coordinates. -/
noncomputable def posetEmb (P : Type*) [Fintype P] (c : ℕ) :
    P ↪ Fin (Fintype.card P + c) :=
  (Fintype.equivFin P).toEmbedding.trans (Fin.castLEEmb (by omega))

/-- The family of antichains of `P`, transported into `Finset (Fin M)` along `posetEmb`. -/
noncomputable def pascalFamily (P : Type*) [Fintype P] [DecidableEq P] [PartialOrder P]
    [DecidableRel ((· ≤ ·) : P → P → Prop)] (c : ℕ) :
    Finset (Finset (Fin (Fintype.card P + c))) :=
  (antichainsOf P).image (fun A => A.map (posetEmb P c))

lemma mem_pascalFamily {c : ℕ} {s : Finset (Fin (Fintype.card P + c))} :
    s ∈ pascalFamily P c ↔ ∃ A ∈ antichainsOf P, A.map (posetEmb P c) = s := by
  simp [pascalFamily]

/-- The transported family is downward closed. -/
theorem pascalFamily_downward_closed (c : ℕ) :
    ∀ s ∈ pascalFamily P c, ∀ t ⊆ s, t ∈ pascalFamily P c := by
  intro s hs t hts
  obtain ⟨A, hA, rfl⟩ := mem_pascalFamily.1 hs
  obtain ⟨B, hBA, rfl⟩ := Finset.subset_map_iff.1 hts
  exact mem_pascalFamily.2 ⟨B, antichainsOf_downward_closed A hA B hBA, rfl⟩

/-- The transported family is nonempty (the empty antichain belongs to it). -/
theorem pascalFamily_nonempty (c : ℕ) : (pascalFamily P c).Nonempty :=
  ⟨(∅ : Finset P).map (posetEmb P c),
    mem_pascalFamily.2 ⟨∅, by simp [antichainsOf], rfl⟩⟩

/-- Transport preserves face numbers: the `j`-faces of `pascalFamily P c` are exactly the
antichains of `P` with `j` elements. -/
theorem faceCount_pascalFamily (c j : ℕ) :
    PascalSmoothing.faceCount (pascalFamily P c) j = antichainCount P j := by
  classical
  have hset : (pascalFamily P c).filter (fun s => s.card = j)
      = ((antichainsOf P).filter (fun A => A.card = j)).image
          (fun A => A.map (posetEmb P c)) := by
    ext s
    simp only [Finset.mem_filter, Finset.mem_image, pascalFamily]
    constructor
    · rintro ⟨⟨A, hA, rfl⟩, hcard⟩
      exact ⟨A, ⟨hA, by simpa using hcard⟩, rfl⟩
    · rintro ⟨A, ⟨hA, hcard⟩, rfl⟩
      exact ⟨⟨A, hA, rfl⟩, by simpa using hcard⟩
  rw [PascalSmoothing.faceCount, hset,
    Finset.card_image_of_injective _ (Finset.map_injective (posetEmb P c)), antichainCount]

/-- **The bridge**: the erasure profile of the poset is precisely the erasure shadow of the
downward-closed family `pascalFamily P c` on `M = |P| + c` coordinates. -/
theorem erasureShadow_pascalFamily (c d : ℕ) :
    PascalSmoothing.erasureShadow (pascalFamily P c) d = erasureProfile P c d := by
  rw [PascalSmoothing.erasureShadow, erasureProfile]
  rw [Finset.sum_congr rfl (fun j _ => by rw [faceCount_pascalFamily])]
  refine (Finset.sum_subset (by
    intro x hx
    simp only [Finset.mem_range] at *
    omega) ?_).symm
  intro x _ hx
  simp only [Finset.mem_range, not_lt] at hx
  rw [antichainCount_eq_zero_of_lt (by omega), Nat.zero_mul]

/-! ### Equations (7) and (7a) for the poset erasure profile -/

/-- **The Pascal-smoothing inequality for the poset erasure profile**, with denominators
cleared: for `1 ≤ d ≤ M - 1` and `m = M - d`,
`E_d^2 / (E_{d-1} E_{d+1}) ≥ (8/9)·(d+1)(m+1)/(dm)`. -/
theorem erasureProfile_pascal_smoothing (c : ℕ) {d : ℕ} (hd1 : 1 ≤ d)
    (hd2 : d ≤ Fintype.card P + c - 1) :
    8 * ((d + 1) * (Fintype.card P + c - d + 1))
        * (erasureProfile P c (d - 1) * erasureProfile P c (d + 1))
      ≤ 9 * (d * (Fintype.card P + c - d)) * erasureProfile P c d ^ 2 := by
  have h := PascalSmoothing.pascal_smoothing_shadow_lemma (pascalFamily P c)
    (pascalFamily_nonempty c) (pascalFamily_downward_closed c) hd1 hd2
  simpa only [erasureShadow_pascalFamily] using h

/-- Log-concavity of the erasure profile under the numerical condition
`d (M - d) ≤ 8 (M + 1)`. -/
theorem erasureProfile_log_concave_of (c : ℕ) {d : ℕ} (hd1 : 1 ≤ d)
    (hd2 : d ≤ Fintype.card P + c - 1)
    (hnum : d * (Fintype.card P + c - d) ≤ 8 * (Fintype.card P + c + 1)) :
    erasureProfile P c (d - 1) * erasureProfile P c (d + 1) ≤ erasureProfile P c d ^ 2 := by
  have h := PascalSmoothing.erasureShadow_log_concave_of (pascalFamily P c)
    (pascalFamily_nonempty c) (pascalFamily_downward_closed c) hd1 hd2 hnum
  simpa only [erasureShadow_pascalFamily] using h

/-- **Log-concavity through defect depth eight**, for arbitrary `M = |P| + c`. -/
theorem erasureProfile_log_concave_depth_le_eight (c : ℕ) {d : ℕ} (hd1 : 1 ≤ d)
    (hd2 : d ≤ Fintype.card P + c - 1) (hd8 : d ≤ 8) :
    erasureProfile P c (d - 1) * erasureProfile P c (d + 1) ≤ erasureProfile P c d ^ 2 := by
  have h := PascalSmoothing.erasureShadow_log_concave_depth_le_eight (pascalFamily P c)
    (pascalFamily_nonempty c) (pascalFamily_downward_closed c) hd1 hd2 hd8
  simpa only [erasureShadow_pascalFamily] using h

/-- **Log-concavity at every interior defect when `M ≤ 33`.** -/
theorem erasureProfile_log_concave_of_M_le_33 (c : ℕ) (hM : Fintype.card P + c ≤ 33)
    {d : ℕ} (hd1 : 1 ≤ d) (hd2 : d ≤ Fintype.card P + c - 1) :
    erasureProfile P c (d - 1) * erasureProfile P c (d + 1) ≤ erasureProfile P c d ^ 2 := by
  have h := PascalSmoothing.erasureShadow_log_concave_of_M_le_33 (pascalFamily P c)
    (pascalFamily_nonempty c) (pascalFamily_downward_closed c) hM hd1 hd2
  simpa only [erasureShadow_pascalFamily] using h

/-- **Equation (7a)**, the quantitative depth-three reserve:
`E_3^2 ≥ (32(M-2) / 27(M-3)) E_2 E_4`, with denominators cleared. -/
theorem erasureProfile_depth_three (c : ℕ) (hM : 4 ≤ Fintype.card P + c) :
    32 * (Fintype.card P + c - 2) * (erasureProfile P c 2 * erasureProfile P c 4)
      ≤ 27 * (Fintype.card P + c - 3) * erasureProfile P c 3 ^ 2 := by
  have h := PascalSmoothing.erasureShadow_depth_three (pascalFamily P c)
    (pascalFamily_nonempty c) (pascalFamily_downward_closed c) hM
  simpa only [erasureShadow_pascalFamily] using h

/-- **Equation (7)**, in strict form: `E_3^2 > E_2 E_4` whenever `E_2 E_4 > 0` and `M ≥ 4`. -/
theorem erasureProfile_depth_three_strict (c : ℕ) (hM : 4 ≤ Fintype.card P + c)
    (hpos : 0 < erasureProfile P c 2 * erasureProfile P c 4) :
    erasureProfile P c 2 * erasureProfile P c 4 < erasureProfile P c 3 ^ 2 := by
  have h := PascalSmoothing.erasureShadow_depth_three_strict (pascalFamily P c)
    (pascalFamily_nonempty c) (pascalFamily_downward_closed c) hM
    (by simpa only [erasureShadow_pascalFamily] using hpos)
  simpa only [erasureShadow_pascalFamily] using h

/-- **Equation (7)**: `E_3^2 ≥ E_2 E_4`, for every finite poset `P` and every `c`, with
`M = |P| + c ≥ 4`. -/
theorem erasureProfile_depth_three_log_concave (c : ℕ) (hM : 4 ≤ Fintype.card P + c) :
    erasureProfile P c 2 * erasureProfile P c 4 ≤ erasureProfile P c 3 ^ 2 := by
  rcases Nat.eq_zero_or_pos (erasureProfile P c 2 * erasureProfile P c 4) with h | h
  · rw [h]; exact Nat.zero_le _
  · exact le_of_lt (erasureProfile_depth_three_strict c hM h)

/-- Equation (7) in the coefficient form in which the note states it:
`([t^{M-3}] (1+t)^c I(B(P);t))^2 ≥ ([t^{M-2}]…) · ([t^{M-4}]…)`. -/
theorem coeff_erasure_depth_three_log_concave (c : ℕ) (hM : 4 ≤ Fintype.card P + c) :
    (((1 : Polynomial ℤ) + X) ^ c * indepPoly P).coeff (Fintype.card P + c - 2)
        * (((1 : Polynomial ℤ) + X) ^ c * indepPoly P).coeff (Fintype.card P + c - 4)
      ≤ ((((1 : Polynomial ℤ) + X) ^ c * indepPoly P).coeff (Fintype.card P + c - 3)) ^ 2 := by
  rw [erasureProfile_eq_coeff c (by omega), erasureProfile_eq_coeff c (by omega),
    erasureProfile_eq_coeff c (by omega)]
  exact_mod_cast erasureProfile_depth_three_log_concave (P := P) c hM

end MatchingBag
