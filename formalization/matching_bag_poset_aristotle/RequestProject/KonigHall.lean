import Mathlib
import RequestProject.MatchingCover

/-!
# König's theorem for bipartite graphs, from Hall's marriage theorem

Section 1 of `matching_bag_poset_reduction_2026-08-20.md` uses the fact that a minimum
vertex cover of a tree with maximum matching `M` has exactly `|M|` vertices.  This is
König's theorem; it is proved here for arbitrary bipartite graphs (acyclicity is not
needed) from Mathlib's Hall marriage theorem, and *not* postulated.

* `MatchingBag.IsMatchingSet G S`: `S` is a set of pairwise disjoint edges of `G`, recorded
  as ordered pairs.
* `MatchingBag.IsVertexCover G C`: `C` meets every edge of `G`.
* `MatchingBag.card_matching_le_cover`: every cover is at least as large as every matching.
* `MatchingBag.konig`: for a bipartite `G` with a maximum matching `M`, some vertex cover has
  exactly `|M|` elements.

The proof of `konig` is the classical one: a minimum cover `C` splits by colour as
`C = C_L ⊔ C_R`; Hall's condition holds for the neighbourhoods of `C_L` inside the
complement of `C`, since otherwise replacing a subset `A ⊆ C_L` by that neighbourhood would
give a smaller cover; the two matchings produced by Hall use disjoint sets of vertices, so
their union is a matching with `|C|` edges.
-/

open Finset

namespace MatchingBag

variable {V : Type*}

/-- `S` is a matching of `G`: a set of edges of `G`, recorded as ordered pairs, which are
pairwise disjoint. -/
def IsMatchingSet (G : SimpleGraph V) (S : Finset (V × V)) : Prop :=
  (∀ e ∈ S, G.Adj e.1 e.2) ∧
    ∀ e ∈ S, ∀ e' ∈ S, e ≠ e' → e.1 ≠ e'.1 ∧ e.1 ≠ e'.2 ∧ e.2 ≠ e'.1 ∧ e.2 ≠ e'.2

/-- `C` is a vertex cover of `G`. -/
def IsVertexCover (G : SimpleGraph V) (C : Finset V) : Prop :=
  ∀ ⦃u v⦄, G.Adj u v → u ∈ C ∨ v ∈ C

variable [DecidableEq V]

/-- Every vertex cover has at least as many vertices as any matching has edges. -/
theorem card_matching_le_cover {G : SimpleGraph V} {S : Finset (V × V)} {C : Finset V}
    (hS : IsMatchingSet G S) (hC : IsVertexCover G C) : S.card ≤ C.card := by
  classical
  set f : V × V → V := fun e => if e.1 ∈ C then e.1 else e.2 with hf
  have hfC : ∀ e ∈ S, f e ∈ C := by
    intro e he
    by_cases h : e.1 ∈ C
    · simp [hf, h]
    · rcases hC (hS.1 e he) with h1 | h2
      · exact absurd h1 h
      · simpa [hf, h] using h2
  have hinj : Set.InjOn f (S : Set (V × V)) := by
    intro e he e' he' hee'
    by_contra hne
    obtain ⟨h1, h2, h3, h4⟩ := hS.2 e he e' he' hne
    by_cases ha : e.1 ∈ C <;> by_cases hb : e'.1 ∈ C <;>
      simp only [hf, ha, hb, if_true, if_false] at hee'
    · exact h1 hee'
    · exact h2 hee'
    · exact h3 hee'
    · exact h4 hee'
  calc S.card = (S.image f).card := (Finset.card_image_of_injOn hinj).symm
    _ ≤ C.card := Finset.card_le_card (fun v hv => by
        obtain ⟨e, he, rfl⟩ := Finset.mem_image.1 hv
        exact hfC e he)

variable [Fintype V]

/-- Hall's condition for the `b`-coloured part of a minimum cover: every subset `A` of it has
at least `|A|` neighbours outside the cover.  Consequently there is an injective choice of
such a neighbour for each vertex of that part. -/
theorem hall_side {G : SimpleGraph V} {col : V → Bool}
    (hcol : ∀ ⦃u v⦄, G.Adj u v → col u ≠ col v) {C : Finset V}
    (hC : IsVertexCover G C)
    (hmin : ∀ C' : Finset V, IsVertexCover G C' → C.card ≤ C'.card) (b : Bool) :
    ∃ f : {v // v ∈ C.filter (fun v => col v = b)} → V, Function.Injective f ∧
      ∀ x, G.Adj x.1 (f x) ∧ f x ∉ C := by
  classical
  set r : {v // v ∈ C.filter (fun v => col v = b)} → V → Prop :=
    fun x y => G.Adj x.1 y ∧ y ∉ C with hr
  have hall : ∀ A : Finset {v // v ∈ C.filter (fun v => col v = b)},
      A.card ≤ (Finset.univ.filter (fun y => ∃ a ∈ A, r a y)).card := by
    intro A
    by_contra hlt
    push_neg at hlt
    set N := Finset.univ.filter (fun y => ∃ a ∈ A, r a y) with hN
    set CA := A.image (fun a => a.1) with hCA
    have hCAsub : CA ⊆ C := by
      intro w hw
      obtain ⟨a, _, rfl⟩ := Finset.mem_image.1 hw
      exact (Finset.mem_filter.1 a.2).1
    have hCAcol : ∀ w ∈ CA, col w = b := by
      intro w hw
      obtain ⟨a, _, rfl⟩ := Finset.mem_image.1 hw
      exact (Finset.mem_filter.1 a.2).2
    have hCAcard : CA.card = A.card :=
      Finset.card_image_of_injective _ (fun x y h => Subtype.ext h)
    set C' := (C \ CA) ∪ N with hC'
    have hcov : IsVertexCover G C' := by
      have key : ∀ w w' : V, w ∈ C → G.Adj w w' → w ∈ C' ∨ w' ∈ C' := by
        intro w w' hwC hadj
        by_cases hw : w ∈ CA
        · obtain ⟨a, haA, hav⟩ := Finset.mem_image.1 hw
          have hcolw : col w = b := hCAcol w hw
          have hcolw' : col w' ≠ b := by
            have h := hcol hadj
            rw [hcolw] at h
            exact fun hb => h hb.symm
          by_cases hw'C : w' ∈ C
          · right
            refine Finset.mem_union_left _ (Finset.mem_sdiff.2 ⟨hw'C, ?_⟩)
            intro hmem
            exact hcolw' (hCAcol w' hmem)
          · right
            refine Finset.mem_union_right _
              (Finset.mem_filter.2 ⟨Finset.mem_univ _, a, haA, ?_, hw'C⟩)
            rw [hav]; exact hadj
        · exact Or.inl (Finset.mem_union_left _ (Finset.mem_sdiff.2 ⟨hwC, hw⟩))
      intro u v huv
      rcases hC huv with hu | hv
      · exact key u v hu huv
      · exact (key v u hv huv.symm).symm
    have hcard : C'.card < C.card := by
      have h1 : C'.card ≤ (C \ CA).card + N.card := Finset.card_union_le _ _
      have h2 : (C \ CA).card = C.card - CA.card := Finset.card_sdiff_of_subset hCAsub
      have h3 : CA.card ≤ C.card := Finset.card_le_card hCAsub
      omega
    exact absurd (hmin C' hcov) (by omega)
  exact (Fintype.all_card_le_filter_rel_iff_exists_injective r).1 hall

/-- **König's theorem** (the direction that needs Hall): in a bipartite graph, a maximum
matching `M` admits a vertex cover with exactly `|M|` vertices. -/
theorem konig {G : SimpleGraph V} {col : V → Bool}
    (hcol : ∀ ⦃u v⦄, G.Adj u v → col u ≠ col v) {M : Finset (V × V)}
    (hM : IsMatchingSet G M) (hmax : ∀ S : Finset (V × V), IsMatchingSet G S → S.card ≤ M.card) :
    ∃ C : Finset V, IsVertexCover G C ∧ C.card = M.card := by
  classical
  obtain ⟨C, hCmem, hCmin⟩ := Finset.exists_min_image
    ((Finset.univ : Finset (Finset V)).filter (fun C => IsVertexCover G C)) Finset.card
    ⟨Finset.univ, Finset.mem_filter.2 ⟨Finset.mem_univ _, fun u v _ => Or.inl (Finset.mem_univ u)⟩⟩
  have hC : IsVertexCover G C := (Finset.mem_filter.1 hCmem).2
  have hmin : ∀ C' : Finset V, IsVertexCover G C' → C.card ≤ C'.card := fun C' hC' =>
    hCmin C' (Finset.mem_filter.2 ⟨Finset.mem_univ _, hC'⟩)
  obtain ⟨f, hfinj, hf⟩ := hall_side hcol hC hmin true
  obtain ⟨g, hginj, hg⟩ := hall_side hcol hC hmin false
  have hcolL : ∀ a : {v // v ∈ C.filter (fun v => col v = true)}, col a.1 = true :=
    fun a => (Finset.mem_filter.1 a.2).2
  have hmemL : ∀ a : {v // v ∈ C.filter (fun v => col v = true)}, a.1 ∈ C :=
    fun a => (Finset.mem_filter.1 a.2).1
  have hcolR : ∀ b : {v // v ∈ C.filter (fun v => col v = false)}, col b.1 = false :=
    fun b => (Finset.mem_filter.1 b.2).2
  have hmemR : ∀ b : {v // v ∈ C.filter (fun v => col v = false)}, b.1 ∈ C :=
    fun b => (Finset.mem_filter.1 b.2).1
  have hcolf : ∀ a : {v // v ∈ C.filter (fun v => col v = true)}, col (f a) = false := by
    intro a
    have h := hcol (hf a).1
    rw [hcolL a] at h
    cases hfa : col (f a)
    · rfl
    · exact absurd hfa.symm h
  have hcolg : ∀ b : {v // v ∈ C.filter (fun v => col v = false)}, col (g b) = true := by
    intro b
    have h := hcol (hg b).1
    rw [hcolR b] at h
    cases hgb : col (g b)
    · exact absurd hgb.symm h
    · rfl
  have hTF : ∀ x y : V, col x = true → col y = false → x ≠ y := by
    intro x y hx hy h
    rw [h, hy] at hx
    exact Bool.noConfusion hx
  set S : Finset (V × V) :=
    ((C.filter (fun v => col v = true)).attach.image (fun a => (a.1, f a))) ∪
      ((C.filter (fun v => col v = false)).attach.image (fun b => (g b, b.1))) with hSdef
  have hmemS : ∀ e ∈ S, (∃ a, e = (a.1, f a)) ∨ (∃ b, e = (g b, b.1)) := by
    intro e he
    rcases Finset.mem_union.1 he with h | h
    · obtain ⟨a, -, rfl⟩ := Finset.mem_image.1 h
      exact Or.inl ⟨a, rfl⟩
    · obtain ⟨b, -, rfl⟩ := Finset.mem_image.1 h
      exact Or.inr ⟨b, rfl⟩
  have hSmatching : IsMatchingSet G S := by
    refine ⟨?_, ?_⟩
    · intro e he
      rcases hmemS e he with ⟨a, rfl⟩ | ⟨b, rfl⟩
      · exact (hf a).1
      · exact (hg b).1.symm
    · intro e he e' he' hne
      rcases hmemS e he with ⟨a, rfl⟩ | ⟨b, rfl⟩ <;>
        rcases hmemS e' he' with ⟨a', rfl⟩ | ⟨b', rfl⟩
      · have haa : a ≠ a' := by rintro rfl; exact hne rfl
        exact ⟨fun h => haa (Subtype.ext h), hTF _ _ (hcolL a) (hcolf a'),
          fun h => hTF _ _ (hcolL a') (hcolf a) h.symm, fun h => haa (hfinj h)⟩
      · refine ⟨fun h => (hg b').2 ?_, hTF _ _ (hcolL a) (hcolR b'),
          fun h => hTF _ _ (hcolg b') (hcolf a) h.symm, fun h => (hf a).2 ?_⟩
        · have h' : (a : V) = g b' := h
          rw [← h']; exact hmemL a
        · have h' : f a = (b' : V) := h
          rw [h']; exact hmemR b'
      · refine ⟨fun h => (hg b).2 ?_, fun h => hTF _ _ (hcolg b) (hcolf a') h,
          fun h => hTF _ _ (hcolL a') (hcolR b) h.symm, fun h => (hf a').2 ?_⟩
        · have h' : g b = (a' : V) := h
          rw [h']; exact hmemL a'
        · have h' : (b : V) = f a' := h
          rw [← h']; exact hmemR b
      · have hbb : b ≠ b' := by rintro rfl; exact hne rfl
        exact ⟨fun h => hbb (hginj h), fun h => hTF _ _ (hcolg b) (hcolR b') h,
          fun h => hTF _ _ (hcolg b') (hcolR b) h.symm, fun h => hbb (Subtype.ext h)⟩
  have hScard : S.card = C.card := by
    have hdisj : Disjoint ((C.filter (fun v => col v = true)).attach.image (fun a => (a.1, f a)))
        ((C.filter (fun v => col v = false)).attach.image (fun b => (g b, b.1))) := by
      rw [Finset.disjoint_left]
      rintro e he he'
      obtain ⟨a, -, rfl⟩ := Finset.mem_image.1 he
      obtain ⟨b, -, hb⟩ := Finset.mem_image.1 he'
      have h1 : g b = a.1 := congrArg Prod.fst hb
      exact (hg b).2 (by rw [h1]; exact hmemL a)
    rw [hSdef, Finset.card_union_of_disjoint hdisj,
      Finset.card_image_of_injective _ (fun x y h => Subtype.ext (congrArg Prod.fst h)),
      Finset.card_image_of_injective _ (fun x y h => Subtype.ext (congrArg Prod.snd h)),
      Finset.card_attach, Finset.card_attach]
    have hsplit := Finset.card_filter_add_card_filter_not (s := C) (p := fun v => col v = true)
    have heq : C.filter (fun v => ¬ (col v = true)) = C.filter (fun v => col v = false) := by
      apply Finset.filter_congr
      intro v _
      cases col v <;> simp
    rw [heq] at hsplit
    omega
  refine ⟨C, hC, le_antisymm ?_ (card_matching_le_cover hM hC)⟩
  rw [← hScard]
  exact hmax S hSmatching

end MatchingBag
