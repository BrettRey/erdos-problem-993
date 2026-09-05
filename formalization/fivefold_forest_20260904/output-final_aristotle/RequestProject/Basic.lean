import RequestProject.Defs

/-!
# Basic properties of independent sets, `alpha`, and the layer counts.
-/

open Finset

namespace FivefoldForest

open scoped Classical

variable {V : Type*} {G : SimpleGraph V} {S A B M : Finset V}

@[simp] lemma mem_indFam : A ∈ indFam G S ↔ A ⊆ S ∧ G.IsIndepSet (A : Set V) := by
  simp [indFam, Finset.mem_filter, Finset.mem_powerset]

lemma indep_of_subset (h : G.IsIndepSet (A : Set V)) (hBA : B ⊆ A) :
    G.IsIndepSet (B : Set V) :=
  Set.Pairwise.mono (by exact_mod_cast hBA) h

lemma indFam_subset (hA : A ∈ indFam G S) (hBA : B ⊆ A) : B ∈ indFam G S := by
  rw [mem_indFam] at hA ⊢
  exact ⟨hBA.trans hA.1, indep_of_subset hA.2 hBA⟩

lemma empty_mem_indFam : (∅ : Finset V) ∈ indFam G S := by
  simp

lemma indFam_nonempty : (indFam G S).Nonempty := ⟨∅, empty_mem_indFam⟩

lemma card_le_alpha (hA : A ∈ indFam G S) : A.card ≤ alpha G S :=
  Finset.le_sup hA

lemma exists_max_indFam (G : SimpleGraph V) (S : Finset V) :
    ∃ M ∈ indFam G S, M.card = alpha G S := by
  obtain ⟨M, hM, hcard⟩ := Finset.exists_mem_eq_sup (indFam G S) indFam_nonempty Finset.card
  exact ⟨M, hM, hcard.symm⟩

lemma alpha_le_card (G : SimpleGraph V) (S : Finset V) : alpha G S ≤ S.card := by
  obtain ⟨M, hM, hcard⟩ := exists_max_indFam G S
  rw [← hcard]
  exact Finset.card_le_card (mem_indFam.1 hM).1

lemma indFam_mono (h : S ⊆ B) : indFam G S ⊆ indFam G B := by
  intro A hA
  rw [mem_indFam] at hA ⊢
  exact ⟨hA.1.trans h, hA.2⟩

lemma alpha_mono (h : S ⊆ B) : alpha G S ≤ alpha G B :=
  Finset.sup_mono (indFam_mono h)

/-- Every extendable set is independent and inside `S`. -/
lemma Ext.mem_indFam (h : Ext G S A) : A ∈ indFam G S := by
  obtain ⟨M, hM, _, hAM⟩ := h
  exact indFam_subset hM hAM

lemma Ext.mono (h : Ext G S A) (hBA : B ⊆ A) : Ext G S B := by
  obtain ⟨M, hM, hc, hAM⟩ := h
  exact ⟨M, hM, hc, hBA.trans hAM⟩

lemma ext_of_max (hM : M ∈ indFam G S) (hc : M.card = alpha G S) : Ext G S M :=
  ⟨M, hM, hc, Subset.rfl⟩

lemma ecnt_add_bcnt (G : SimpleGraph V) (S : Finset V) (k : ℕ) :
    ecnt G S k + bcnt G S k = icnt G S k := by
  classical
  have h := Finset.card_filter_add_card_filter_not
    (s := (indFam G S).filter fun A => A.card = k) (p := fun A => Ext G S A)
  simpa [ecnt, bcnt, icnt, Finset.filter_filter, and_comm] using h

/-- Independent sets of size larger than `alpha` do not exist. -/
lemma icnt_eq_zero_of_gt (h : alpha G S < k) : icnt G S k = 0 := by
  rw [icnt, Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  rintro A hA rfl
  exact absurd (card_le_alpha hA) (by omega)

lemma ecnt_eq_zero_of_gt (h : alpha G S < k) : ecnt G S k = 0 := by
  have := ecnt_add_bcnt G S k
  have h0 := icnt_eq_zero_of_gt (G := G) (S := S) (k := k) h
  omega

lemma bcnt_eq_zero_of_gt (h : alpha G S < k) : bcnt G S k = 0 := by
  have := ecnt_add_bcnt G S k
  have h0 := icnt_eq_zero_of_gt (G := G) (S := S) (k := k) h
  omega

/-- Every independent set of size `alpha` is extendable, so at defect `0` nothing is blocked. -/
lemma bD_zero (G : SimpleGraph V) (S : Finset V) : bD G S 0 = 0 := by
  rw [bD, if_pos (Nat.zero_le _), bcnt, Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  intro A hA
  rintro ⟨hc, hne⟩
  exact hne (ext_of_max hA (by simpa using hc))

lemma eD_zero_pos (G : SimpleGraph V) (S : Finset V) : 1 ≤ eD G S 0 := by
  obtain ⟨M, hM, hc⟩ := exists_max_indFam G S
  rw [eD, if_pos (Nat.zero_le _)]
  refine Finset.card_pos.2 ⟨M, ?_⟩
  simp only [Finset.mem_filter]
  exact ⟨hM, by simpa using hc, ext_of_max hM hc⟩

lemma eD_eq_zero_of_gt (h : alpha G S < d) : eD G S d = 0 := by
  rw [eD, if_neg (by omega)]

lemma bD_eq_zero_of_gt (h : alpha G S < d) : bD G S d = 0 := by
  rw [bD, if_neg (by omega)]

lemma margin_eq_zero_of_gt (h : alpha G S < d) : margin G S d = 0 := by
  simp [margin, eD_eq_zero_of_gt h, bD_eq_zero_of_gt h]

/-- Every size below `alpha` is realised by an independent set. -/
lemma icnt_pos (h : k ≤ alpha G S) : 1 ≤ icnt G S k := by
  obtain ⟨M, hM, hcard⟩ := exists_max_indFam G S
  obtain ⟨A, hAM, hA⟩ := Finset.exists_subset_card_eq (show k ≤ M.card by omega)
  exact Finset.card_pos.2 ⟨A, by simp only [Finset.mem_filter]; exact ⟨indFam_subset hM hAM, hA⟩⟩

/-- `alpha S ≥ 1` for nonempty `S`. -/
lemma one_le_alpha (hS : S.Nonempty) : 1 ≤ alpha G S := by
  obtain ⟨v, hv⟩ := hS
  have : ({v} : Finset V) ∈ indFam G S := by
    rw [mem_indFam]
    refine ⟨Finset.singleton_subset_iff.2 hv, ?_⟩
    intro x hx y hy hxy
    simp only [Finset.coe_singleton, Set.mem_singleton_iff] at hx hy
    exact absurd (hx.trans hy.symm) hxy
  simpa using card_le_alpha this

/-- Every layer at defect `d ≤ alpha` contains at least one extendable set. -/
lemma eD_pos (h : d ≤ alpha G S) : 1 ≤ eD G S d := by
  obtain ⟨M, hM, hcard⟩ := exists_max_indFam G S
  obtain ⟨A, hAM, hA⟩ := Finset.exists_subset_card_eq (show alpha G S - d ≤ M.card by omega)
  rw [eD, if_pos h]
  refine Finset.card_pos.2 ⟨A, ?_⟩
  simp only [Finset.mem_filter]
  exact ⟨indFam_subset hM hAM, hA, ⟨M, hM, hcard, hAM⟩⟩

/-- Weak incidence bound between consecutive layers of extendable sets. -/
lemma ecnt_incidence (G : SimpleGraph V) (S : Finset V) (k : ℕ) :
    (k + 1) * ecnt G S (k + 1) ≤ S.card * ecnt G S k := by
  classical
  set E1 := (indFam G S).filter (fun A => A.card = k + 1 ∧ Ext G S A) with hE1
  set E0 := (indFam G S).filter (fun A => A.card = k ∧ Ext G S A) with hE0
  have hstep : ∀ A ∈ E1, ∀ x ∈ A, A.erase x ∈ E0 := by
    intro A hA x hx
    rw [hE1, Finset.mem_filter] at hA
    rw [hE0, Finset.mem_filter]
    exact ⟨indFam_subset hA.1 (Finset.erase_subset _ _),
      by rw [Finset.card_erase_of_mem hx, hA.2.1]; rfl,
      hA.2.2.mono (Finset.erase_subset _ _)⟩
  have hmap : ∀ q ∈ E1.sigma (fun A => A),
      (⟨q.1.erase q.2, q.2⟩ : (_ : Finset V) × V) ∈ E0.sigma (fun B => S \ B) := by
    rintro ⟨A, x⟩ hq
    rw [Finset.mem_sigma] at hq
    have hxA : x ∈ A := hq.2
    have hAS : A ⊆ S := by
      have := hq.1
      rw [hE1, Finset.mem_filter, mem_indFam] at this
      exact this.1.1
    rw [Finset.mem_sigma]
    exact ⟨hstep A hq.1 x hxA, by
      simp only [Finset.mem_sdiff, Finset.mem_erase]
      exact ⟨hAS hxA, by simp⟩⟩
  have hinj : Set.InjOn (fun q : (_ : Finset V) × V => (⟨q.1.erase q.2, q.2⟩ : (_ : Finset V) × V))
      (E1.sigma (fun A => A)) := by
    rintro ⟨A, x⟩ hA ⟨B, y⟩ hB hEq
    simp only [Finset.coe_sigma, Set.mem_sigma_iff, Finset.mem_coe] at hA hB
    simp only [Sigma.mk.injEq] at hEq
    obtain ⟨h1, h2⟩ := hEq
    subst h2
    have : A = insert x (A.erase x) := (Finset.insert_erase hA.2).symm
    have hB' : B = insert x (B.erase x) := (Finset.insert_erase hB.2).symm
    rw [this, hB', h1]
  have hle := Finset.card_le_card_of_injOn _ hmap hinj
  rw [Finset.card_sigma, Finset.card_sigma] at hle
  have hL : ∑ A ∈ E1, A.card = (k + 1) * E1.card := by
    rw [Finset.sum_congr rfl (fun A hA => ?_), Finset.sum_const, smul_eq_mul, mul_comm]
    rw [hE1, Finset.mem_filter] at hA
    exact hA.2.1
  have hR : ∑ B ∈ E0, (S \ B).card ≤ S.card * E0.card := by
    calc ∑ B ∈ E0, (S \ B).card ≤ ∑ _B ∈ E0, S.card :=
          Finset.sum_le_sum fun B _ => Finset.card_le_card Finset.sdiff_subset
      _ = S.card * E0.card := by rw [Finset.sum_const, smul_eq_mul, mul_comm]
  rw [hL] at hle
  calc (k + 1) * ecnt G S (k + 1) = (k + 1) * E1.card := by rw [ecnt, hE1]
    _ ≤ ∑ B ∈ E0, (S \ B).card := hle
    _ ≤ S.card * E0.card := hR
    _ = S.card * ecnt G S k := by rw [ecnt, hE0]

/-- The forest hypothesis enters only through this bound. -/
lemma e0_le_two_e1 (hS : S.Nonempty) (hcard : S.card ≤ 2 * alpha G S) :
    eD G S 0 ≤ 2 * eD G S 1 := by
  have ha : 1 ≤ alpha G S := one_le_alpha hS
  have hinc := ecnt_incidence G S (alpha G S - 1)
  have h1 : alpha G S - 1 + 1 = alpha G S := by omega
  rw [h1] at hinc
  have he0 : eD G S 0 = ecnt G S (alpha G S) := by rw [eD, if_pos (Nat.zero_le _), Nat.sub_zero]
  have he1 : eD G S 1 = ecnt G S (alpha G S - 1) := by rw [eD, if_pos (by omega)]
  rw [he0, he1]
  have h2 : alpha G S * ecnt G S (alpha G S) ≤ (2 * alpha G S) * ecnt G S (alpha G S - 1) :=
    le_trans hinc (Nat.mul_le_mul_right _ hcard)
  have h3 : alpha G S * ecnt G S (alpha G S) ≤ alpha G S * (2 * ecnt G S (alpha G S - 1)) := by
    calc alpha G S * ecnt G S (alpha G S) ≤ (2 * alpha G S) * ecnt G S (alpha G S - 1) := h2
      _ = alpha G S * (2 * ecnt G S (alpha G S - 1)) := by ring
  exact Nat.le_of_mul_le_mul_left h3 (by omega)

lemma e1_le_four_e2 (h2 : 2 ≤ alpha G S) (hcard : S.card ≤ 2 * alpha G S) :
    eD G S 1 ≤ 4 * eD G S 2 := by
  have hinc := ecnt_incidence G S (alpha G S - 2)
  have h1 : alpha G S - 2 + 1 = alpha G S - 1 := by omega
  rw [h1] at hinc
  have he1 : eD G S 1 = ecnt G S (alpha G S - 1) := by rw [eD, if_pos (by omega)]
  have he2 : eD G S 2 = ecnt G S (alpha G S - 2) := by rw [eD, if_pos (by omega)]
  rw [he1, he2]
  have hstep : (alpha G S - 1) * ecnt G S (alpha G S - 1)
      ≤ (alpha G S - 1) * (4 * ecnt G S (alpha G S - 2)) := by
    calc (alpha G S - 1) * ecnt G S (alpha G S - 1)
        ≤ S.card * ecnt G S (alpha G S - 2) := hinc
      _ ≤ (2 * alpha G S) * ecnt G S (alpha G S - 2) := Nat.mul_le_mul_right _ hcard
      _ ≤ (alpha G S - 1) * (4 * ecnt G S (alpha G S - 2)) := by
          have : 2 * alpha G S ≤ (alpha G S - 1) * 4 := by omega
          calc (2 * alpha G S) * ecnt G S (alpha G S - 2)
              ≤ ((alpha G S - 1) * 4) * ecnt G S (alpha G S - 2) :=
                Nat.mul_le_mul_right _ this
            _ = (alpha G S - 1) * (4 * ecnt G S (alpha G S - 2)) := by ring
  exact Nat.le_of_mul_le_mul_left hstep (by omega)

lemma e1_le_three_e2 (h3 : 3 ≤ alpha G S) (hcard : S.card ≤ 2 * alpha G S) :
    eD G S 1 ≤ 3 * eD G S 2 := by
  have hinc := ecnt_incidence G S (alpha G S - 2)
  have h1 : alpha G S - 2 + 1 = alpha G S - 1 := by omega
  rw [h1] at hinc
  have he1 : eD G S 1 = ecnt G S (alpha G S - 1) := by rw [eD, if_pos (by omega)]
  have he2 : eD G S 2 = ecnt G S (alpha G S - 2) := by rw [eD, if_pos (by omega)]
  rw [he1, he2]
  have hstep : (alpha G S - 1) * ecnt G S (alpha G S - 1)
      ≤ (alpha G S - 1) * (3 * ecnt G S (alpha G S - 2)) := by
    calc (alpha G S - 1) * ecnt G S (alpha G S - 1)
        ≤ S.card * ecnt G S (alpha G S - 2) := hinc
      _ ≤ (2 * alpha G S) * ecnt G S (alpha G S - 2) := Nat.mul_le_mul_right _ hcard
      _ ≤ (alpha G S - 1) * (3 * ecnt G S (alpha G S - 2)) := by
          have : 2 * alpha G S ≤ (alpha G S - 1) * 3 := by omega
          calc (2 * alpha G S) * ecnt G S (alpha G S - 2)
              ≤ ((alpha G S - 1) * 3) * ecnt G S (alpha G S - 2) :=
                Nat.mul_le_mul_right _ this
            _ = (alpha G S - 1) * (3 * ecnt G S (alpha G S - 2)) := by ring
  exact Nat.le_of_mul_le_mul_left hstep (by omega)

/-- The deepest layer consists of the empty set alone. -/
lemma eD_alpha_eq_one (G : SimpleGraph V) (S : Finset V) : eD G S (alpha G S) = 1 := by
  rw [eD, if_pos le_rfl, Nat.sub_self, ecnt]
  rw [Finset.card_eq_one]
  refine ⟨∅, ?_⟩
  ext A
  simp only [Finset.mem_filter, Finset.mem_singleton, Finset.card_eq_zero]
  constructor
  · rintro ⟨-, h, -⟩; exact h
  · rintro rfl
    exact ⟨empty_mem_indFam, rfl, (exists_max_indFam G S).elim
      fun M hM => ⟨M, hM.1, hM.2, Finset.empty_subset _⟩⟩

lemma bD_alpha_eq_zero (G : SimpleGraph V) (S : Finset V) : bD G S (alpha G S) = 0 := by
  rw [bD, if_pos le_rfl, Nat.sub_self, bcnt, Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  intro A hA
  rintro ⟨hc, hne⟩
  rw [Finset.card_eq_zero] at hc
  subst hc
  exact hne ((exists_max_indFam G S).elim fun M hM => ⟨M, hM.1, hM.2, Finset.empty_subset _⟩)

/-- Beyond the deepest layer the margin is trivially nonnegative. -/
lemma margin_nonneg_of_alpha_le (h : alpha G S ≤ d) : 0 ≤ margin G S d := by
  rcases eq_or_lt_of_le h with rfl | h'
  · simp [margin, bD_alpha_eq_zero, eD_alpha_eq_one]
  · simp [margin_eq_zero_of_gt h']

/-- If `alpha S = 1` then the unique layer at defect one is `{∅}`. -/
lemma eD_one_of_alpha_one (h : alpha G S = 1) : eD G S 1 = 1 := by
  have := eD_alpha_eq_one G S
  rwa [h] at this

/-- A convenient combination of the two incidence bounds, valid for all nonempty forests. -/
lemma e1_le_e0_add_six_e2 (hS : S.Nonempty) (hcard : S.card ≤ 2 * alpha G S) :
    eD G S 1 ≤ eD G S 0 + 6 * eD G S 2 := by
  rcases lt_or_ge (alpha G S) 2 with h | h
  · have h1 : alpha G S = 1 := le_antisymm (by omega) (one_le_alpha hS)
    have := eD_one_of_alpha_one (G := G) (S := S) h1
    have h0 := eD_pos (G := G) (S := S) (d := 0) (Nat.zero_le _)
    omega
  · have := e1_le_four_e2 (G := G) (S := S) h hcard
    omega


/-- Removing a vertex drops the independence number by at most one. -/
lemma alpha_erase_succ_ge (G : SimpleGraph V) (S : Finset V) (v : V) :
    alpha G S ≤ alpha G (S.erase v) + 1 := by
  obtain ⟨M, hM, hcard⟩ := exists_max_indFam G S
  have hsub : M.erase v ∈ indFam G (S.erase v) := by
    rw [mem_indFam] at hM ⊢
    refine ⟨?_, indep_of_subset hM.2 (Finset.erase_subset _ _)⟩
    intro x hx
    rw [Finset.mem_erase] at hx ⊢
    exact ⟨hx.1, hM.1 hx.2⟩
  have := card_le_alpha hsub
  have h2 := Finset.pred_card_le_card_erase (s := M) (a := v)
  omega

lemma alpha_erase_le (G : SimpleGraph V) (S : Finset V) (v : V) :
    alpha G (S.erase v) ≤ alpha G S :=
  alpha_mono (Finset.erase_subset _ _)

lemma icnt_mono {S T : Finset V} (h : S ⊆ T) (k : ℕ) : icnt G S k ≤ icnt G T k := by
  refine Finset.card_le_card ?_
  intro A hA
  rw [Finset.mem_filter] at hA ⊢
  exact ⟨indFam_mono h hA.1, hA.2⟩

/-- Eligible sets: nonempty, with no blocked sets at defect one or two. -/
def Eligible (G : SimpleGraph V) (S : Finset V) : Prop :=
  S.Nonempty ∧ bD G S 1 = 0 ∧ bD G S 2 = 0

/-- The fivefold bounds at defects three and four. -/
def Good (G : SimpleGraph V) (S : Finset V) : Prop :=
  0 ≤ margin G S 3 ∧ 0 ≤ margin G S 4

/-- The exceptional profile, that of the four-leaf star `K₁,₄`. -/
def StarProf (G : SimpleGraph V) (S : Finset V) : Prop :=
  alpha G S = 4 ∧ eD G S 0 = 1 ∧ eD G S 1 = 4 ∧ eD G S 2 = 6 ∧ eD G S 3 = 4 ∧
    eD G S 4 = 1 ∧ bD G S 3 = 1 ∧ bD G S 4 = 0

lemma StarProf.margin3 (h : StarProf G S) : margin G S 3 = -1 := by
  simp [margin, h.2.2.2.2.1, h.2.2.2.2.2.2.1]

lemma StarProf.margin4 (h : StarProf G S) : margin G S 4 = 1 := by
  simp [margin, h.2.2.2.2.2.1, h.2.2.2.2.2.2.2]

end FivefoldForest
