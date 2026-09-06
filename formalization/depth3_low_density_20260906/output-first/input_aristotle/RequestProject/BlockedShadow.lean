import RequestProject.ExtendableCount
import RequestProject.BlockedCertificate

/-!
# The blocked shadow bound `(α - 4) b₂ ≤ 6 b₃`

Double counting the incidences between blocked independent `(α-2)`-sets and their blocked
one-vertex deletions.

* `MatchingBag.TreeMatching.card_extVerts_le_six`: an independent `(α-3)`-set leaves three
  matching bags empty, each of size at most two, hence has at most six one-vertex
  extensions.
* `MatchingBag.TreeMatching.card_blocked_deletions_ge`: a blocked independent `(α-2)`-set
  admits at least `α-4` blocked one-vertex deletions, namely those outside a chosen
  unary-or-pair certificate.
* `MatchingBag.TreeMatching.blocked_shadow`: the resulting inequality.
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable {D : TreeMatching V}

lemma blocked_mono {T U : Finset V} (hTU : T ⊆ U)
    (hT : ¬ ∃ I ∈ D.maxIndepSets, T ⊆ I) : ¬ ∃ I ∈ D.maxIndepSets, U ⊆ I := by
  rintro ⟨I, hI, hUI⟩
  exact hT ⟨I, hI, hTU.trans hUI⟩

/-! ### At most six one-vertex extensions -/

variable (D)

/-- The one-vertex extensions of an independent set. -/
noncomputable def extVerts (S : Finset V) : Finset V :=
  Finset.univ.filter (fun v => v ∉ S ∧ ∀ u ∈ insert v S, ∀ w ∈ insert v S, ¬ D.G.Adj u w)

variable {D}

lemma mem_extVerts {S : Finset V} {v : V} :
    v ∈ D.extVerts S ↔ v ∉ S ∧ ∀ u ∈ insert v S, ∀ w ∈ insert v S, ¬ D.G.Adj u w := by
  simp [extVerts]

/-- An independent set of size `α - 3` leaves exactly three bags empty, so it has at most
six one-vertex extensions. -/
theorem card_extVerts_le_six {S : Finset V} (hind : ∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v)
    (hcard : S.card + 3 = Fintype.card D.Bag) :
    (D.extVerts S).card ≤ 6 := by
  classical
  have himg : (D.extVerts S).image D.bagOf ⊆ (Finset.univ \ S.image D.bagOf) := by
    intro b hb
    obtain ⟨v, hv, rfl⟩ := Finset.mem_image.1 hb
    rw [mem_extVerts] at hv
    rw [Finset.mem_sdiff]
    refine ⟨Finset.mem_univ _, ?_⟩
    intro hc
    obtain ⟨u, hu, hbu⟩ := Finset.mem_image.1 hc
    have huv : u = v := bagOf_injOn_indep hv.2
      (by exact_mod_cast Finset.mem_insert_of_mem hu)
      (by exact_mod_cast Finset.mem_insert_self v S) hbu
    exact hv.1 (huv ▸ hu)
  have hcompl : (Finset.univ \ S.image D.bagOf).card = 3 := by
    rw [Finset.card_sdiff_of_subset (Finset.subset_univ _), Finset.card_univ,
      card_image_bagOf hind]
    omega
  have hfib : ∀ b ∈ (D.extVerts S).image D.bagOf,
      ((D.extVerts S).filter (fun v => D.bagOf v = b)).card ≤ 2 := by
    intro b _
    refine le_trans (Finset.card_le_card ?_) (card_bagSet_le b)
    intro v hv
    rw [Finset.mem_filter] at hv
    exact bagOf_eq_iff.1 hv.2
  have h := Finset.card_le_mul_card_image (D.extVerts S) 2 hfib
  have h2 : ((D.extVerts S).image D.bagOf).card ≤ 3 := by
    rw [← hcompl]; exact Finset.card_le_card himg
  omega

/-! ### At least `α - 4` blocked deletions -/

/-- A blocked independent `(α-2)`-set has at least `α - 4` blocked one-vertex deletions. -/
theorem card_blocked_deletions_ge {S : Finset V} {a : ℕ} (ha : 4 ≤ a)
    (hS : S ∈ D.blkSets (a - 2)) :
    a - 4 ≤ (S.filter (fun v => S.erase v ∈ D.blkSets (a - 3))).card := by
  classical
  rw [mem_blkSets] at hS
  obtain ⟨hcard, hind, hblk⟩ := hS
  obtain ⟨T, hTS, hT2, hTblk⟩ := exists_small_blocked_subset hind hblk
  have hsub : S \ T ⊆ S.filter (fun v => S.erase v ∈ D.blkSets (a - 3)) := by
    intro v hv
    rw [Finset.mem_sdiff] at hv
    rw [Finset.mem_filter]
    refine ⟨hv.1, ?_⟩
    rw [mem_blkSets]
    refine ⟨?_, ?_, ?_⟩
    · rw [Finset.card_erase_of_mem hv.1, hcard]
      omega
    · intro u hu w hw
      exact hind u (Finset.mem_of_mem_erase hu) w (Finset.mem_of_mem_erase hw)
    · refine blocked_mono ?_ hTblk
      intro z hz
      exact Finset.mem_erase.2 ⟨fun hzv => hv.2 (hzv ▸ hz), hTS hz⟩
  have hcards : (S \ T).card = S.card - T.card := Finset.card_sdiff_of_subset hTS
  have := Finset.card_le_card hsub
  omega

/-! ### The double count -/

/-- **Blocked shadow bound.**  For a forest with `α ≥ 4`, `(α - 4) b₂ ≤ 6 b₃`. -/
theorem blocked_shadow (D : TreeMatching V) (ha : 4 ≤ Fintype.card D.Bag) :
    (Fintype.card D.Bag - 4) * (D.blkSets (Fintype.card D.Bag - 2)).card
      ≤ 6 * (D.blkSets (Fintype.card D.Bag - 3)).card := by
  classical
  set a := Fintype.card D.Bag with hadef
  set A := D.blkSets (a - 2) with hA
  set B := D.blkSets (a - 3) with hB
  set Inc : Finset (V × Finset V) :=
    A.biUnion (fun S => (S.filter (fun v => S.erase v ∈ B)).image (fun v => (v, S))) with hInc
  have hmemInc1 : ∀ p ∈ Inc, p.1 ∈ p.2 := by
    intro p hp
    rw [hInc, Finset.mem_biUnion] at hp
    obtain ⟨S, -, hp2⟩ := hp
    obtain ⟨v, hv, rfl⟩ := Finset.mem_image.1 hp2
    rw [Finset.mem_filter] at hv
    exact hv.1
  have hmemInc2 : ∀ p ∈ Inc, p.2 ∈ A := by
    intro p hp
    rw [hInc, Finset.mem_biUnion] at hp
    obtain ⟨S, hSA, hp2⟩ := hp
    obtain ⟨v, -, rfl⟩ := Finset.mem_image.1 hp2
    exact hSA
  have hmemInc3 : ∀ p ∈ Inc, p.2.erase p.1 ∈ B := by
    intro p hp
    rw [hInc, Finset.mem_biUnion] at hp
    obtain ⟨S, -, hp2⟩ := hp
    obtain ⟨v, hv, rfl⟩ := Finset.mem_image.1 hp2
    rw [Finset.mem_filter] at hv
    exact hv.2
  -- count the incidences by the blocked `(α-2)`-sets
  have hdisj : ∀ S ∈ A, ∀ S' ∈ A, S ≠ S' →
      Disjoint ((S.filter (fun v => S.erase v ∈ B)).image (fun v => (v, S)))
        ((S'.filter (fun v => S'.erase v ∈ B)).image (fun v => (v, S'))) := by
    intro S _ S' _ hne
    rw [Finset.disjoint_left]
    rintro p hp hp'
    obtain ⟨v, -, rfl⟩ := Finset.mem_image.1 hp
    obtain ⟨v', -, hv'⟩ := Finset.mem_image.1 hp'
    exact hne (congrArg Prod.snd hv').symm
  have hcount : Inc.card = ∑ S ∈ A, (S.filter (fun v => S.erase v ∈ B)).card := by
    rw [hInc, Finset.card_biUnion hdisj]
    refine Finset.sum_congr rfl fun S _ => ?_
    exact Finset.card_image_of_injOn fun x _ y _ hxy => congrArg Prod.fst hxy
  have hlow : (a - 4) * A.card ≤ Inc.card := by
    rw [hcount]
    calc (a - 4) * A.card = ∑ _S ∈ A, (a - 4) := by
          rw [Finset.sum_const, smul_eq_mul, mul_comm]
      _ ≤ ∑ S ∈ A, (S.filter (fun v => S.erase v ∈ B)).card :=
          Finset.sum_le_sum fun S hS => card_blocked_deletions_ge ha (by rwa [hA] at hS)
  -- count the incidences by the blocked `(α-3)`-sets
  have hfiber : ∀ S' ∈ B, (Inc.filter (fun p => p.2.erase p.1 = S')).card ≤ 6 := by
    intro S' hS'
    have hS'mem := hS'
    rw [hB, mem_blkSets] at hS'mem
    obtain ⟨hc, hi, -⟩ := hS'mem
    refine le_trans (Finset.card_le_card_of_injOn Prod.fst ?_ ?_)
      (card_extVerts_le_six hi (by rw [hc]; omega))
    · intro p hp
      simp only [Finset.mem_coe, Finset.mem_filter] at hp
      simp only [Finset.mem_coe]
      have hvS : p.1 ∈ p.2 := hmemInc1 p hp.1
      have hSA := hmemInc2 p hp.1
      rw [hA, mem_blkSets] at hSA
      have hSeq : p.2 = insert p.1 S' := by rw [← hp.2, Finset.insert_erase hvS]
      rw [mem_extVerts]
      refine ⟨?_, ?_⟩
      · rw [← hp.2]; exact Finset.notMem_erase p.1 p.2
      · rw [← hSeq]; exact hSA.2.1
    · intro p hp q hq hpq
      simp only [Finset.mem_coe, Finset.mem_filter] at hp hq
      have e1 : p.2 = insert p.1 S' := by
        rw [← hp.2, Finset.insert_erase (hmemInc1 p hp.1)]
      have e2 : q.2 = insert q.1 S' := by
        rw [← hq.2, Finset.insert_erase (hmemInc1 q hq.1)]
      have hsnd : p.2 = q.2 := by rw [e1, e2, hpq]
      exact Prod.ext hpq hsnd
  have hhigh : Inc.card ≤ 6 * B.card := by
    rw [Finset.card_eq_sum_card_fiberwise (t := B) (fun p hp => hmemInc3 p hp)]
    refine le_trans (Finset.sum_le_sum hfiber) ?_
    rw [Finset.sum_const, smul_eq_mul, mul_comm]
  omega

end TreeMatching

end MatchingBag
