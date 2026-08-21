import RequestProject.ClanAudit.GlobalPartition

/-!
# The explicit collision of the proof notes, resolved

`notes/two_hub_clan_cancellation_attack_2026-08-20.md` records a concrete collision that
blocks a naive combination of the published one-hub pairing with the new two-hub pairing:
for five unit arms at `v` and three unit arms at `u`,

* the **double-hub source** puts multiplicity one on both hubs and on all eight leaves;
* the **single-hub source** puts multiplicity one on `v` and its five leaves, zero on `u`,
  and `(2,1,1)` on the three `u`-leaves;

and the two "sources" have the same image under the naive transform.

This file formalizes that pair and proves that the block construction of
`GlobalPartition.lean` *separates* it: the "single-hub source" is not a source at all —
its `u`-side has an inactive hub, hence is not an `ActiveSide` — it is exactly the image
`transf` of the `u`-side of the double-hub source, so both maps lie in **one and the same**
four-element block.  No partner is used twice, and no map is missing.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph
open scoped Classical

/-! ### A spider with `k` unit arms -/

/-- `k` pendant paths of length one. -/
def unitArms (k : ℕ) : Fin k → ℕ := fun _ => 1

theorem unitArms_apply (k : ℕ) (i : Fin k) : unitArms k i = 1 := rfl

/-- The all-ones clan map of a spider with `k` unit arms: the "bad" source state. -/
def onesSide (k : ℕ) : SpiderV k (unitArms k) → ℕ := fun _ => 1

/-- The image of `onesSide` under the canonical transformation: the hub is switched off
and the first leaf is cloned. -/
def twoOnesSide (k : ℕ) : SpiderV k (unitArms k) → ℕ
  | none => 0
  | some p => if p.1.val = 0 then 2 else 1

theorem pref_onesSide (k : ℕ) (i : Fin k) : pref (onesSide k) i = 1 := by
  refine pref_eq _ i 1 (le_refl 1) (fun j hj => ?_) (Or.inl rfl)
  rw [armVal, dif_pos (show j < unitArms k i from by rw [unitArms_apply]; omega)]
  rfl

theorem admRun_onesSide (k : ℕ) : AdmRun (onesSide k) := by
  intro i j hj
  have hpos := hj j le_rfl
  by_cases h : j < unitArms k i
  · rw [armVal, dif_pos h]
    rfl
  · rw [armVal_of_ge _ _ (by omega)] at hpos
    omega

theorem pNum_onesSide (k : ℕ) : pNum (onesSide k) = k := by
  rw [pNum, Finset.filter_true_of_mem (fun i _ => by rw [pref_onesSide]), Finset.card_univ,
    Fintype.card_fin]

theorem activeSide_onesSide {k : ℕ} (hk : 2 ≤ k) : ActiveSide (onesSide k) :=
  ⟨rfl, admRun_onesSide k, by rw [pNum_onesSide]; exact hk⟩

/-! ### The canonical choices on `onesSide` -/

theorem armKey_onesSide (k : ℕ) (i : Fin k) : armKey (onesSide k) i = i.val + k := by
  rw [armKey, pref_onesSide, one_mul]

theorem mem_keySet_onesSide (k : ℕ) (i : Fin k) : i.val + k ∈ keySet (onesSide k) :=
  ⟨i, by rw [pref_onesSide], by rw [armKey_onesSide]⟩

theorem key0_onesSide {k : ℕ} (hk : 0 < k) : key0 (onesSide k) = k := by
  have hmem : k ∈ keySet (onesSide k) := by
    have := mem_keySet_onesSide k ⟨0, hk⟩
    simpa using this
  refine le_antisymm (Nat.sInf_le hmem) ?_
  obtain ⟨i, -, hi⟩ := Nat.sInf_mem (s := keySet (onesSide k)) ⟨k, hmem⟩
  rw [armKey_onesSide] at hi
  show k ≤ sInf (keySet (onesSide k))
  omega

theorem key1_onesSide {k : ℕ} (hk : 1 < k) : key1 (onesSide k) = k + 1 := by
  rw [key1, key0_onesSide (by omega)]
  have hmem : k + 1 ∈ keySet (onesSide k) \ {k} := by
    refine ⟨?_, ?_⟩
    · have := mem_keySet_onesSide k ⟨1, hk⟩
      simpa [Nat.add_comm] using this
    · simp
  refine le_antisymm (Nat.sInf_le hmem) ?_
  obtain ⟨⟨i, -, hi⟩, hne⟩ := Nat.sInf_mem (s := keySet (onesSide k) \ {k}) ⟨k + 1, hmem⟩
  rw [armKey_onesSide] at hi
  simp only [Set.mem_singleton_iff] at hne
  show k + 1 ≤ sInf (keySet (onesSide k) \ {k})
  omega

theorem idx0_onesSide {k : ℕ} (hk : 0 < k) : idx0 (onesSide k) = 0 := by
  rw [idx0, key0_onesSide hk, Nat.mod_self]

theorem plen_onesSide {k : ℕ} (hk : 0 < k) : plen (onesSide k) = 1 := by
  rw [plen, key0_onesSide hk, Nat.div_self hk]

theorem idx1_onesSide {k : ℕ} (hk : 1 < k) : idx1 (onesSide k) = 1 := by
  rw [idx1, key1_onesSide hk, Nat.add_mod_left, Nat.mod_eq_of_lt hk]

/-- **The canonical transformation of the all-ones state.**  With unit arms the shortest
odd prefix has length one, so the published rule replaces the first leaf by a cloned `K₂`
and switches the hub off. -/
theorem transf_onesSide {k : ℕ} (hk : 1 < k) : transf (onesSide k) = twoOnesSide k := by
  rw [transf, idx0_onesSide (by omega), idx1_onesSide hk, plen_onesSide (by omega)]
  funext x
  match x with
  | none => rfl
  | some ⟨i, j⟩ =>
      have hj : j.val = 0 := by
        have h1 : (j : ℕ) < 1 := j.isLt
        omega
      show (if i.val = 0 ∧ j.val < 1 then altVal j.val
            else if i.val = 1 ∧ j.val < 1 - 1 then altVal j.val
            else onesSide k (some ⟨i, j⟩)) = twoOnesSide k (some ⟨i, j⟩)
      show _ = (if i.val = 0 then 2 else 1)
      by_cases h0 : i.val = 0
      · rw [if_pos ⟨h0, by omega⟩, if_pos h0, hj, altVal_zero]
      · rw [if_neg (fun hc => h0 hc.1), if_neg (fun hc => by omega), if_neg h0]
        rfl

/-! ### The collision of the notes -/

/-- The double-hub source of the notes: three unit arms at `u`, five at `v`, everything of
multiplicity one. -/
noncomputable def collisionSource : DV 3 (unitArms 3) 5 (unitArms 5) → ℕ :=
  Sum.elim (onesSide 3) (onesSide 5)

/-- The "single-hub source" of the notes: `u` switched off with leaf multiplicities
`(2,1,1)`, and the all-ones state at `v`. -/
noncomputable def collisionPartner : DV 3 (unitArms 3) 5 (unitArms 5) → ℕ :=
  Sum.elim (twoOnesSide 3) (onesSide 5)

theorem collisionPartner_eq :
    collisionPartner = Sum.elim (transf (onesSide 3)) (onesSide 5) := by
  rw [collisionPartner, transf_onesSide (by omega)]

/-- The two maps of the collision are distinct. -/
theorem collisionSource_ne_collisionPartner : collisionSource ≠ collisionPartner := by
  intro h
  have := congrFun h (Sum.inl none)
  exact absurd this (by decide)

/-- Both maps have the same total order `10`, as the notes state. -/
theorem sum_collisionSource : ∑ x, collisionSource x = 10 := by decide

theorem sum_collisionPartner : ∑ x, collisionPartner x = 10 := by decide

/-- **The collision is resolved by the block construction.**  The "single-hub source" is
the canonical image of the `u`-side of the double-hub source, so the two maps have the
same representative and lie in one and the same block. -/
theorem collision_same_block : rep collisionPartner = rep collisionSource := by
  rw [collisionPartner_eq, collisionSource]
  exact collision_resolved_left (activeSide_onesSide (by omega)) (onesSide 5)

/-- The block in question is a genuine four-element block, and it contains both maps of
the collision. -/
theorem collision_block :
    collisionSource ∈ block 3 (unitArms 3) 5 (unitArms 5) 10 collisionSource
      ∧ collisionPartner ∈ block 3 (unitArms 3) 5 (unitArms 5) 10 collisionSource
      ∧ (block 3 (unitArms 3) 5 (unitArms 5) 10 collisionSource).card = 4 := by
  have hact3 : ActiveSide (onesSide 3) := activeSide_onesSide (by omega)
  have hact5 : ActiveSide (onesSide 5) := activeSide_onesSide (by omega)
  have hfix : rep collisionSource = collisionSource := by
    show Sum.elim (repSide (onesSide 3)) (repSide (onesSide 5)) = collisionSource
    rw [repSide_of_active hact3, repSide_of_active hact5]
    rfl
  have hmem : collisionSource ∈ mapsOfOrder 3 (unitArms 3) 5 (unitArms 5) 10 := by
    rw [mem_mapsOfOrder]
    exact sum_collisionSource
  have hmemp : collisionPartner ∈ mapsOfOrder 3 (unitArms 3) 5 (unitArms 5) 10 := by
    rw [mem_mapsOfOrder]
    exact sum_collisionPartner
  refine ⟨Finset.mem_filter.mpr ⟨hmem, hfix⟩, ?_, ?_⟩
  · refine Finset.mem_filter.mpr ⟨hmemp, ?_⟩
    rw [collision_same_block, hfix]
  · have hl : ActiveSide (fun x => collisionSource (Sum.inl x)) := hact3
    have hr : ActiveSide (fun y => collisionSource (Sum.inr y)) := hact5
    rw [card_block hfix hmem, card_sideBlock_active hl, card_sideBlock_active hr]

end ClanAudit
