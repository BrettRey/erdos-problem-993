import RequestProject.Cases
import RequestProject.Paths

/-!
# What the forest hypothesis provides

The induction uses acyclicity of `G` only through the two properties collected in
`ForestLike`:

* every finite vertex set `T` satisfies `|T| ≤ 2 α(T)` (a forest is bipartite);
* every nonempty finite vertex set is a star, or splits into two non-adjacent
  nonempty parts, or carries a pendant configuration with a base of size at least two.
-/

open Finset SimpleGraph

namespace FivefoldForest

open scoped Classical

variable {V : Type*} {G : SimpleGraph V}

/-- The consequences of acyclicity used by the induction. -/
structure ForestLike (G : SimpleGraph V) : Prop where
  card_le : ∀ T : Finset V, T.card ≤ 2 * alpha G T
  decomp : ∀ S : Finset V, S.Nonempty →
    (∃ c, IsStarAt G S c) ∨
    (∃ S₁ S₂, Split G S S₁ S₂ ∧ S₁.Nonempty ∧ S₂.Nonempty) ∨
    (∃ p w L, Pendant G S p w L ∧ 2 ≤ (pBase S p L).card)

/-- If all neighbours of `c` are leaves, then `S` is a star or splits. -/
lemma star_or_split {S : Finset V} {c : V} (hc : c ∈ S)
    (hleaf : ∀ x ∈ S, G.Adj c x → ∀ y ∈ S, G.Adj x y → y = c) :
    (∃ c, IsStarAt G S c) ∨ (∃ S₁ S₂, Split G S S₁ S₂ ∧ S₁.Nonempty ∧ S₂.Nonempty) := by
  classical
  have hKS : insert c (S.filter fun x => G.Adj c x) ⊆ S := by
    intro x hx
    rcases Finset.mem_insert.1 hx with rfl | hx
    · exact hc
    · exact (Finset.mem_filter.1 hx).1
  by_cases hKeq : insert c (S.filter fun x => G.Adj c x) = S
  · left
    refine ⟨c, hc, fun x hx hxc => ?_⟩
    have hcx : G.Adj c x := by
      rw [← hKeq] at hx
      rcases Finset.mem_insert.1 hx with rfl | hx
      · exact absurd rfl hxc
      · exact (Finset.mem_filter.1 hx).2
    exact ⟨hcx, fun y hy hadj => hleaf x hx hcx y hy hadj⟩
  · right
    obtain ⟨z, hzS, hzK⟩ := Finset.exists_of_ssubset (lt_of_le_of_ne hKS hKeq)
    refine ⟨insert c (S.filter fun x => G.Adj c x), S \ insert c (S.filter fun x => G.Adj c x),
      ⟨(Finset.union_sdiff_of_subset hKS).symm, Finset.disjoint_sdiff, ?_⟩,
      ⟨c, Finset.mem_insert_self _ _⟩, ⟨z, Finset.mem_sdiff.2 ⟨hzS, hzK⟩⟩⟩
    intro x hx y hy hadj
    obtain ⟨hyS, hyK⟩ := Finset.mem_sdiff.1 hy
    rcases Finset.mem_insert.1 hx with rfl | hx
    · exact hyK (Finset.mem_insert_of_mem (Finset.mem_filter.2 ⟨hyS, hadj⟩))
    · obtain ⟨hxS, hcx⟩ := Finset.mem_filter.1 hx
      exact hyK ((hleaf x hxS hcx y hyS hadj) ▸ Finset.mem_insert_self _ _)

/-- Every acyclic graph is `ForestLike`. -/
theorem forestLike_of_isAcyclic (hac : G.IsAcyclic) : ForestLike G := by
  classical
  refine ⟨fun T => card_le_two_alpha hac T.card T le_rfl, ?_⟩
  intro S hS
  have wrap : ∀ c ∈ S, (∀ x ∈ S, G.Adj c x → ∀ y ∈ S, G.Adj x y → y = c) →
      (∃ c, IsStarAt G S c) ∨
      (∃ S₁ S₂, Split G S S₁ S₂ ∧ S₁.Nonempty ∧ S₂.Nonempty) ∨
      (∃ p w L, Pendant G S p w L ∧ 2 ≤ (pBase S p L).card) :=
    fun c hc hl => (star_or_split hc hl).imp id Or.inl
  obtain ⟨u, v, p, hp, hlen⟩ := tPathLen_max (G := G) hS
  cases p with
  | nil =>
    refine wrap u (hp.2 u (by simp)) ?_
    intro x hx hadj
    exfalso
    have hpath : IsTPath G S (Walk.cons hadj Walk.nil) := by
      refine ⟨(Walk.cons_isPath_iff _ _).2 ⟨Walk.IsPath.nil, by simp [hadj.ne]⟩, ?_⟩
      intro z hz
      rw [Walk.support_cons, Walk.support_nil, List.mem_cons, List.mem_singleton] at hz
      rcases hz with rfl | rfl
      · exact hp.2 z (by simp)
      · exact hx
    have hlt := hpath.length_lt
    rw [Walk.length_cons, Walk.length_nil] at hlt
    exact not_tPathLen_of_gt (m := 1) (by rw [← hlen]; simp) (by omega)
      ⟨_, _, _, hpath, by simp⟩
  | cons h₁ q =>
    cases q with
    | nil =>
      -- a longest path with a single edge `u - v`
      refine wrap u (hp.2 u (by simp)) ?_
      intro x hx hadj
      have hxb : x = v := max_nbr_eq hac h₁ hp hlen hx hadj
      subst hxb
      intro y hy hadj'
      have hp' : IsTPath G S (Walk.cons h₁.symm (Walk.nil : G.Walk u u)) := by
        refine ⟨(Walk.cons_isPath_iff _ _).2 ⟨Walk.IsPath.nil, by simp [h₁.ne']⟩, ?_⟩
        intro z hz
        rw [Walk.support_cons, Walk.support_nil, List.mem_cons, List.mem_singleton] at hz
        rcases hz with rfl | rfl
        · exact hx
        · exact hp.2 z (by simp)
      exact max_nbr_eq hac h₁.symm hp'
        (by simp only [Walk.length_cons, Walk.length_nil] at hlen ⊢; omega) hy hadj'
    | cons h₂ t =>
      cases t with
      | nil =>
        rename_i b
        -- a longest path `u - b - v`: the middle vertex is a star centre
        obtain ⟨hpath1, hu⟩ := (Walk.cons_isPath_iff _ _).1 hp.1
        refine wrap b (hp.2 b (by simp)) ?_
        intro x hx hadj
        by_cases hxc : x = v
        · subst hxc
          intro y hy hadj'
          have hp' : IsTPath G S
              (Walk.cons h₂.symm (Walk.cons h₁.symm (Walk.nil : G.Walk u u))) := by
            refine ⟨(Walk.cons_isPath_iff _ _).2
              ⟨(Walk.cons_isPath_iff _ _).2 ⟨Walk.IsPath.nil, by simp [h₁.ne']⟩, ?_⟩, ?_⟩
            · simp only [Walk.support_cons, Walk.support_nil, List.mem_cons,
                List.not_mem_nil, or_false]
              push_neg
              refine ⟨h₂.ne', fun hh => hu ?_⟩
              rw [← hh]
              simp
            · intro z hz
              simp only [Walk.support_cons, Walk.support_nil, List.mem_cons,
                List.not_mem_nil, or_false] at hz
              rcases hz with rfl | rfl | rfl
              · exact hx
              · exact hp.2 z (by simp)
              · exact hp.2 z (by simp)
          exact max_nbr_eq hac h₂.symm hp'
            (by simp only [Walk.length_cons, Walk.length_nil] at hlen ⊢; omega) hy hadj'
        · exact max_leaf hac h₁ h₂ hp hlen hx hadj hxc
      | cons h₃ t' =>
        rename_i b c d
        -- a longest path with at least three edges: a pendant configuration
        obtain ⟨hpath1, hu⟩ := (Walk.cons_isPath_iff _ _).1 hp.1
        obtain ⟨hpath2, hb⟩ := (Walk.cons_isPath_iff _ _).1 hpath1
        obtain ⟨hpath3, hc⟩ := (Walk.cons_isPath_iff _ _).1 hpath2
        right; right
        refine ⟨b, c, S.filter (fun x => G.Adj b x ∧ x ≠ c),
          ⟨hp.2 b (by simp), hp.2 c (by simp), h₂, ?_, Finset.filter_subset _ _, ?_, ?_, ?_, ?_,
            ?_⟩, ?_⟩
        · intro hmem
          exact (Finset.mem_filter.1 hmem).2.2 rfl
        · intro hmem
          exact G.irrefl (Finset.mem_filter.1 hmem).2.1
        · refine ⟨u, Finset.mem_filter.2 ⟨hp.2 u (by simp), h₁.symm, fun hh => hu ?_⟩⟩
          rw [hh]
          simp
        · exact fun x hx => (Finset.mem_filter.1 hx).2.1
        · intro x hx
          exact max_leaf hac h₁ h₂ hp hlen (Finset.mem_filter.1 hx).1
            (Finset.mem_filter.1 hx).2.1 (Finset.mem_filter.1 hx).2.2
        · intro y hy hadj
          by_cases hyc : y = c
          · exact Or.inl hyc
          · exact Or.inr (Finset.mem_filter.2 ⟨hy, hadj, hyc⟩)
        · -- the base contains the third and fourth vertices of the path
          have hcB : c ∈ pBase S b (S.filter (fun x => G.Adj b x ∧ x ≠ c)) := by
            simp only [pBase, Finset.mem_sdiff, Finset.mem_insert, Finset.mem_filter]
            refine ⟨hp.2 c (by simp), ?_⟩
            push_neg
            exact ⟨fun hh => G.irrefl (hh ▸ h₂), fun _ _ => rfl⟩
          have hdB : d ∈ pBase S b (S.filter (fun x => G.Adj b x ∧ x ≠ c)) := by
            simp only [pBase, Finset.mem_sdiff, Finset.mem_insert, Finset.mem_filter]
            refine ⟨hp.2 d (by simp), ?_⟩
            push_neg
            refine ⟨fun hh => hb (by rw [← hh]; simp), fun _ hbd => ?_⟩
            by_cases hdc : d = c
            · exact hdc
            · exfalso
              have hcb := max_leaf hac h₁ h₂ hp hlen (hp.2 d (by simp)) hbd hdc
                c (hp.2 c (by simp)) h₃.symm
              exact G.ne_of_adj h₂ hcb.symm
          refine Finset.one_lt_card.2 ⟨c, hcB, d, hdB, ?_⟩
          intro hh
          subst hh
          exact hc (Walk.start_mem_support t')

end FivefoldForest
