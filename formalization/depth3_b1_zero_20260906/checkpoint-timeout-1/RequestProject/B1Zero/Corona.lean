import RequestProject.B1Zero.Reach
import RequestProject.RootedForestCertificate

/-!
# Coronas and their codimension profiles

A `CoronaData` records a set `base` of vertices together with an injective assignment of a
private pendant neighbour `pend b` to each of them: the pendant is adjacent to `b` and to no
other vertex of the corona.  The *corona* over `X ⊆ base` is `X ∪ pend '' X`.

For such a corona we count independent subsets by codimension below `|X|` — avoiding a set
`A` of vertices — and package the first five counts into the `TopFive` vector `vec`.  The
two structural identities proved here are

* `vec_union`: the profile of a disjoint corona with no edges across is the truncated
  product of the profiles;
* `vec_peel` / `vec_peel_root`: peeling a base vertex `v` off a corona.

These are exactly the two clauses of `RootedForestCertificate.profile`.
-/

open Finset SimpleGraph

namespace MatchingBag

namespace B1Zero

open DepthThree.RootedForestCertificate

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V] {G : SimpleGraph V}

/-- A corona structure: every `base` vertex carries a private pendant neighbour. -/
structure CoronaData (G : SimpleGraph V) where
  /-- The base vertices. -/
  base : Finset V
  /-- The pendant attached to a base vertex. -/
  pend : V → V
  /-- Pendants are not base vertices. -/
  pend_notMem : ∀ b ∈ base, pend b ∉ base
  /-- Distinct base vertices have distinct pendants. -/
  pend_inj : ∀ b ∈ base, ∀ b' ∈ base, pend b = pend b' → b = b'
  /-- A pendant is adjacent to its base vertex. -/
  adj_pend : ∀ b ∈ base, G.Adj b (pend b)
  /-- A pendant has no other neighbour in the corona. -/
  pend_nbr : ∀ b ∈ base, ∀ w ∈ base ∪ base.image pend, G.Adj (pend b) w → w = b

namespace CoronaData

variable (C : CoronaData G)

/-- The corona over a set of base vertices. -/
noncomputable def cor (X : Finset V) : Finset V := X ∪ X.image C.pend

variable {C}

lemma mem_cor {X : Finset V} {x : V} : x ∈ C.cor X ↔ x ∈ X ∨ ∃ b ∈ X, C.pend b = x := by
  simp [CoronaData.cor]

lemma subset_cor {X : Finset V} : X ⊆ C.cor X := fun _ hx => mem_cor.2 (Or.inl hx)

lemma pend_mem_cor {X : Finset V} {b : V} (hb : b ∈ X) : C.pend b ∈ C.cor X :=
  mem_cor.2 (Or.inr ⟨b, hb, rfl⟩)

lemma cor_mono {X Y : Finset V} (h : X ⊆ Y) : C.cor X ⊆ C.cor Y := by
  intro x hx
  rcases mem_cor.1 hx with hx | ⟨b, hb, rfl⟩
  · exact mem_cor.2 (Or.inl (h hx))
  · exact pend_mem_cor (h hb)

lemma cor_union (X Y : Finset V) : C.cor (X ∪ Y) = C.cor X ∪ C.cor Y := by
  simp only [CoronaData.cor, Finset.image_union]
  ext x
  simp only [Finset.mem_union]
  tauto

@[simp] lemma cor_empty : C.cor ∅ = ∅ := by
  rw [CoronaData.cor]
  simp

lemma cor_insert (v : V) (X : Finset V) :
    C.cor (insert v X) = insert v (insert (C.pend v) (C.cor X)) := by
  simp only [CoronaData.cor, Finset.image_insert]
  ext x
  simp only [Finset.mem_union, Finset.mem_insert]
  tauto

/-- The base vertex under a corona vertex. -/
noncomputable def pinv (C : CoronaData G) (x : V) : V :=
  if h : ∃ b, b ∈ C.base ∧ C.pend b = x then h.choose else x

lemma pinv_pend {b : V} (hb : b ∈ C.base) : C.pinv (C.pend b) = b := by
  have h : ∃ b', b' ∈ C.base ∧ C.pend b' = C.pend b := ⟨b, hb, rfl⟩
  rw [pinv, dif_pos h]
  exact C.pend_inj _ h.choose_spec.1 _ hb h.choose_spec.2

lemma pinv_base {x : V} (hx : x ∈ C.base) : C.pinv x = x := by
  rw [pinv, dif_neg]
  rintro ⟨b, hb, rfl⟩
  exact C.pend_notMem b hb hx

lemma pinv_mem {X : Finset V} (hX : X ⊆ C.base) {x : V} (hx : x ∈ C.cor X) : C.pinv x ∈ X := by
  rcases mem_cor.1 hx with hx | ⟨b, hb, rfl⟩
  · rwa [pinv_base (hX hx)]
  · rwa [pinv_pend (hX hb)]

lemma eq_or_eq_pend {X : Finset V} (hX : X ⊆ C.base) {x : V} (hx : x ∈ C.cor X) :
    x = C.pinv x ∨ x = C.pend (C.pinv x) := by
  rcases mem_cor.1 hx with hx | ⟨b, hb, rfl⟩
  · exact Or.inl (pinv_base (hX hx)).symm
  · exact Or.inr (by rw [pinv_pend (hX hb)])

/-- Independent subsets of a corona have at most one vertex per base vertex. -/
theorem card_le_of_indep {X : Finset V} (hX : X ⊆ C.base) {S : Finset V}
    (hS : S ⊆ C.cor X) (hind : IndepOn G S) : S.card ≤ X.card := by
  refine Finset.card_le_card_of_injOn C.pinv (fun x hx => pinv_mem hX (hS hx)) ?_
  intro x hx y hy hxy
  rw [Finset.mem_coe] at hx hy
  by_contra hne
  have hx' := eq_or_eq_pend hX (hS hx)
  have hy' := eq_or_eq_pend hX (hS hy)
  have hpb : C.pinv x ∈ C.base := hX (pinv_mem hX (hS hx))
  rcases hx' with hx' | hx' <;> rcases hy' with hy' | hy'
  · exact hne (by rw [hx', hy', hxy])
  · refine hind x hx y hy ?_
    rw [hx', hy', hxy]
    exact C.adj_pend _ (hxy ▸ hpb)
  · refine hind y hy x hx ?_
    rw [hy', hx', ← hxy]
    exact C.adj_pend _ hpb
  · exact hne (by rw [hx', hy', hxy])

/-- No edges between the coronas over disjoint, mutually non-adjacent base sets. -/
theorem cor_no_cross {X Y : Finset V} (hX : X ⊆ C.base) (hY : Y ⊆ C.base)
    (hdisj : Disjoint X Y) (hcross : ∀ u ∈ X, ∀ w ∈ Y, ¬ G.Adj u w) :
    ∀ u ∈ C.cor X, ∀ w ∈ C.cor Y, ¬ G.Adj u w := by
  have hbase : ∀ Z : Finset V, Z ⊆ C.base → C.cor Z ⊆ C.cor C.base := fun Z hZ => cor_mono hZ
  intro u hu w hw hadj
  rcases mem_cor.1 hu with hu' | ⟨b, hb, rfl⟩
  · rcases mem_cor.1 hw with hw' | ⟨b', hb', rfl⟩
    · exact hcross u hu' w hw' hadj
    · have := C.pend_nbr b' (hY hb') u (hbase X hX hu) hadj.symm
      exact Finset.disjoint_left.1 hdisj hu' (this ▸ hb')
  · rcases mem_cor.1 hw with hw' | ⟨b', hb', rfl⟩
    · have := C.pend_nbr b (hX hb) w (hbase Y hY hw) hadj
      exact Finset.disjoint_left.1 hdisj (this ▸ hb) hw'
    · have := C.pend_nbr b (hX hb) (C.pend b') (hbase Y hY (pend_mem_cor hb')) hadj
      exact C.pend_notMem b' (hY hb') (this ▸ hX hb)

/-- A base vertex outside `X` is outside the corona over `X`. -/
lemma notMem_cor_of_base {X : Finset V} (hX : X ⊆ C.base) {v : V} (hv : v ∈ C.base)
    (hvX : v ∉ X) : v ∉ C.cor X := by
  intro h
  rcases mem_cor.1 h with h | ⟨b, hb, rfl⟩
  · exact hvX h
  · exact C.pend_notMem b (hX hb) hv

/-- The pendant of a base vertex outside `X` is outside the corona over `X`. -/
lemma pend_notMem_cor {X : Finset V} (hX : X ⊆ C.base) {v : V} (hv : v ∈ C.base)
    (hvX : v ∉ X) : C.pend v ∉ C.cor X := by
  intro h
  rcases mem_cor.1 h with h | ⟨b, hb, hbv⟩
  · exact C.pend_notMem v hv (hX h)
  · exact hvX (C.pend_inj b (hX hb) v hv hbv ▸ hb)

/-- A pendant is adjacent to nothing in a corona not containing its base vertex. -/
lemma not_adj_pend_of_notMem {X : Finset V} (hX : X ⊆ C.base) {v : V} (hv : v ∈ C.base)
    (hvX : v ∉ X) : ∀ w ∈ C.cor X, ¬ G.Adj (C.pend v) w := by
  intro w hw hadj
  have hwv := C.pend_nbr v hv w (cor_mono hX hw) hadj
  exact notMem_cor_of_base hX hv hvX (hwv ▸ hw)

end CoronaData

end B1Zero

end MatchingBag
