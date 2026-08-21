import RequestProject.ClanAudit.IsoTransport
import RequestProject.ClanAudit.SpiderLogConcave
import RequestProject.ClanAudit.ConnectorLogConcave

/-!
# The `C₂` family of models and its log-concavity

A finite tree with at most two vertices of degree at least three is a path, a spider, or a
two-hub connector tree.  All three shapes are realized by the explicit models of this
development:

* a single vertex is `spider 0 Fin.elim0`;
* a path on `N + 1` vertices is `spider 1 (fun _ => N)` (`pathGraphIsoSpider`);
* a spider is `spider n len`;
* a two-hub connector tree is `connGraph t m a n b`.

`IsC2Model G` says that `G` is isomorphic to one of the two models `spider n len` or
`connGraph t m a n b`.  The main results are that every `IsC2Model` graph has a
log-concave independence-count sequence (`indepCount_logConcave_of_isC2Model`) and a
log-concave independence polynomial (`indepPoly_logConcave_of_isC2Model`).
-/

namespace ClanAudit

open Finset SimpleGraph

open scoped Classical

/-! ### Paths are spiders -/

/-- A path on `N + 1` vertices is the one-armed spider with an arm of length `N`. -/
def pathGraphIsoSpider (N : ℕ) :
    pathGraph (N + 1) ≃g spider 1 (fun _ => N) where
  toFun i := if h : i.val = 0 then none else some ⟨0, ⟨i.val - 1, show i.val - 1 < N by omega⟩⟩
  invFun x :=
    match x with
    | none => ⟨0, Nat.succ_pos N⟩
    | some p => ⟨p.2.val + 1, Nat.succ_lt_succ p.2.isLt⟩
  left_inv := by
    intro i
    by_cases h : i.val = 0
    · simp only [h, dif_pos]
      exact Fin.ext (by simpa using h.symm)
    · simp only [dif_neg h]
      exact Fin.ext (by simp; omega)
  right_inv := by
    rintro (_ | ⟨i, j⟩)
    · rfl
    · have hi : i = 0 := Subsingleton.elim _ _
      subst hi
      dsimp only
      rw [dif_neg (by simp : ¬ ((⟨j.val + 1, Nat.succ_lt_succ j.isLt⟩ : Fin (N + 1)).val = 0))]
      simp
  map_rel_iff' := by
    intro i j
    simp only [Equiv.coe_fn_mk]
    rcases Nat.eq_zero_or_pos i.val with hi | hi <;>
      rcases Nat.eq_zero_or_pos j.val with hj | hj
    · rw [dif_pos hi, dif_pos hj, pathGraph_adj]
      exact ⟨fun h => absurd h (by simp [spider]), fun h => by omega⟩
    · rw [dif_pos hi, dif_neg hj.ne', pathGraph_adj]
      show (j.val - 1 = 0) ↔ _
      omega
    · rw [dif_pos hj, dif_neg hi.ne', pathGraph_adj]
      show (i.val - 1 = 0) ↔ _
      omega
    · rw [dif_neg hi.ne', dif_neg hj.ne', pathGraph_adj]
      show ((0 : Fin 1) = 0 ∧ _) ↔ _
      simp only [true_and]
      omega

/-- Log-concavity of the independence polynomial of a path, for every number of
vertices. -/
theorem pathGraph_indepPoly_logConcave (N j : ℕ) :
    (indepPoly (pathGraph N)).coeff j * (indepPoly (pathGraph N)).coeff (j + 2)
      ≤ (indepPoly (pathGraph N)).coeff (j + 1) * (indepPoly (pathGraph N)).coeff (j + 1) := by
  cases N with
  | zero =>
      have h2 : indepCount (pathGraph 0) (j + 2) = 0 := by
        have hemp : indepSets (pathGraph 0) (j + 2) = ∅ := by
          ext S
          simp only [mem_indepSets, Finset.notMem_empty, iff_false, not_and]
          intro _
          have hS : S = ∅ := Finset.eq_empty_of_forall_notMem (fun x _ => x.elim0)
          simp [hS]
        simp [indepCount, hemp]
      simp [coeff_indepPoly, h2]
  | succ M =>
      exact indepPoly_logConcave_of_iso (pathGraphIsoSpider M).symm
        (spider_indepCount_logConcave 1 (fun _ => M)) j

/-! ### The `C₂` family -/

/-- `G` is a `C₂` model: it is isomorphic to a spider (which includes single vertices and
paths) or to a two-hub connector tree. -/
def IsC2Model {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V) : Prop :=
  (∃ (n : ℕ) (len : Fin n → ℕ), Nonempty (G ≃g spider n len)) ∨
    (∃ (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ), Nonempty (G ≃g connGraph t m a n b))

theorem isC2Model_spider {n : ℕ} (len : Fin n → ℕ) : IsC2Model (spider n len) :=
  Or.inl ⟨n, len, ⟨Iso.refl⟩⟩

theorem isC2Model_connGraph (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    IsC2Model (connGraph t m a n b) :=
  Or.inr ⟨t, m, a, n, b, ⟨Iso.refl⟩⟩

theorem isC2Model_pathGraph (N : ℕ) : IsC2Model (pathGraph (N + 1)) :=
  Or.inl ⟨1, fun _ => N, ⟨pathGraphIsoSpider N⟩⟩

/-- Every `C₂` model has a log-concave independence-count sequence. -/
theorem indepCount_logConcave_of_isC2Model {V : Type*} [Fintype V] [DecidableEq V]
    {G : SimpleGraph V} (hG : IsC2Model G) (j : ℕ) :
    indepCount G j * indepCount G (j + 2) ≤ indepCount G (j + 1) * indepCount G (j + 1) := by
  rcases hG with ⟨n, len, ⟨e⟩⟩ | ⟨t, m, a, n, b, ⟨e⟩⟩
  · exact indepCount_logConcave_of_iso e.symm (spider_indepCount_logConcave n len) j
  · exact indepCount_logConcave_of_iso e.symm
      (connGraph_indepCount_logConcave (t := t) (m := m) (a := a) (n := n) (b := b)) j

/-- Every `C₂` model has a log-concave independence polynomial. -/
theorem indepPoly_logConcave_of_isC2Model {V : Type*} [Fintype V] [DecidableEq V]
    {G : SimpleGraph V} (hG : IsC2Model G) (j : ℕ) :
    (indepPoly G).coeff j * (indepPoly G).coeff (j + 2)
      ≤ (indepPoly G).coeff (j + 1) * (indepPoly G).coeff (j + 1) := by
  simpa [coeff_indepPoly] using indepCount_logConcave_of_isC2Model hG j

end ClanAudit
