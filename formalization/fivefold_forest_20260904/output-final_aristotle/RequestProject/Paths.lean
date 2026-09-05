import RequestProject.Basic

/-!
# Longest paths inside a vertex set

The forest hypothesis enters the argument through longest paths of the induced subgraph.
This file develops the two facts that are needed:

* the initial vertex of a longest path has at most one neighbour inside the set;
* a neighbour of the *second* vertex of a longest path, other than the third vertex,
  has no neighbour inside the set other than that second vertex.

Both are consequences of acyclicity: a shortcut would close a cycle.
-/

open Finset SimpleGraph

namespace FivefoldForest

open scoped Classical

variable {V : Type*} {G : SimpleGraph V} {T : Finset V}

/-- A path all of whose vertices lie in `T`. -/
def IsTPath (G : SimpleGraph V) (T : Finset V) {u v : V} (p : G.Walk u v) : Prop :=
  p.IsPath ∧ ∀ x ∈ p.support, x ∈ T

/-- There is a path of length `n` all of whose vertices lie in `T`. -/
def TPathLen (G : SimpleGraph V) (T : Finset V) (n : ℕ) : Prop :=
  ∃ (u v : V) (p : G.Walk u v), IsTPath G T p ∧ p.length = n

lemma IsTPath.length_lt {u v : V} {p : G.Walk u v} (hp : IsTPath G T p) : p.length < T.card := by
  have hnd : p.support.Nodup := hp.1.support_nodup
  have hsub : p.support.toFinset ⊆ T := fun x hx => hp.2 x (List.mem_toFinset.1 hx)
  have h1 : p.support.toFinset.card = p.support.length := List.toFinset_card_of_nodup hnd
  have h2 := Finset.card_le_card hsub
  rw [h1, Walk.length_support] at h2
  omega

/-- The maximal length of a path inside `T`. -/
noncomputable def maxTPathLen (G : SimpleGraph V) (T : Finset V) : ℕ :=
  Nat.findGreatest (TPathLen G T) T.card

lemma tPathLen_max (hT : T.Nonempty) : TPathLen G T (maxTPathLen G T) := by
  obtain ⟨a, ha⟩ := hT
  refine Nat.findGreatest_spec (Nat.zero_le _) ⟨a, a, Walk.nil, ⟨Walk.IsPath.nil, ?_⟩, rfl⟩
  intro x hx
  rw [Walk.support_nil, List.mem_singleton] at hx
  exact hx ▸ ha

lemma not_tPathLen_of_gt {m : ℕ} (h : maxTPathLen G T < m) (hm : m ≤ T.card) :
    ¬ TPathLen G T m := Nat.findGreatest_is_greatest h hm

/-- The initial vertex of a longest path has no neighbour in `T` besides the second vertex. -/
lemma max_nbr_eq (hac : G.IsAcyclic) {a b v : V} {q : G.Walk b v} (h : G.Adj a b)
    (hp : IsTPath G T (Walk.cons h q)) (hlen : (Walk.cons h q).length = maxTPathLen G T)
    {y : V} (hyT : y ∈ T) (hadj : G.Adj a y) : y = b := by
  classical
  have hys : y ∈ (Walk.cons h q).support := by
    by_contra hys
    have hp2 : (Walk.cons hadj.symm (Walk.cons h q)).IsPath :=
      (Walk.cons_isPath_iff _ _).2 ⟨hp.1, hys⟩
    have hsup2 : ∀ x ∈ (Walk.cons hadj.symm (Walk.cons h q)).support, x ∈ T := by
      intro x hx
      rw [Walk.support_cons, List.mem_cons] at hx
      rcases hx with rfl | hx
      · exact hyT
      · exact hp.2 x hx
    have htp : IsTPath G T (Walk.cons hadj.symm (Walk.cons h q)) := ⟨hp2, hsup2⟩
    have hlt := htp.length_lt
    rw [Walk.length_cons] at hlt
    exact not_tPathLen_of_gt (m := (Walk.cons h q).length + 1) (by omega) (by omega)
      ⟨_, _, _, htp, by rw [Walk.length_cons]⟩
  rw [Walk.support_cons, List.mem_cons] at hys
  rcases hys with rfl | hyq
  · exact absurd rfl (G.ne_of_adj hadj)
  by_contra hyb
  obtain ⟨hq, haq⟩ := (Walk.cons_isPath_iff h q).1 hp.1
  have hr : (q.takeUntil y hyq).IsPath := hq.takeUntil hyq
  have hP : (Walk.cons h (q.takeUntil y hyq)).IsPath :=
    (Walk.cons_isPath_iff _ _).2 ⟨hr, fun hc => haq (Walk.support_takeUntil_subset q hyq hc)⟩
  have hedge : s(y, a) ∉ (Walk.cons h (q.takeUntil y hyq)).edges := by
    rw [Walk.edges_cons, List.mem_cons]
    rintro (heq | hmem)
    · rcases Sym2.eq_iff.1 heq with ⟨h1, h2⟩ | ⟨h1, h2⟩
      · exact h.ne h2
      · exact hyb h1
    · exact haq (Walk.support_takeUntil_subset q hyq
        (Walk.snd_mem_support_of_mem_edges _ hmem))
  exact hac _ (SimpleGraph.Path.cons_isCycle ⟨_, hP⟩ hadj.symm hedge)

/-- A neighbour of the second vertex of a longest path, other than the third vertex,
has no neighbour in `T` besides that second vertex. -/
lemma max_leaf (hac : G.IsAcyclic) {a b c v : V} {t : G.Walk c v} (h₁ : G.Adj a b) (h₂ : G.Adj b c)
    (hp : IsTPath G T (Walk.cons h₁ (Walk.cons h₂ t)))
    (hlen : (Walk.cons h₁ (Walk.cons h₂ t)).length = maxTPathLen G T)
    {x : V} (hxT : x ∈ T) (hbx : G.Adj b x) (hxc : x ≠ c) :
    ∀ y ∈ T, G.Adj x y → y = b := by
  classical
  obtain ⟨hr, har⟩ := (Walk.cons_isPath_iff h₁ (Walk.cons h₂ t)).1 hp.1
  obtain ⟨ht, hbt⟩ := (Walk.cons_isPath_iff h₂ t).1 hr
  -- `x` does not lie on the tail of the path
  have hxr : x ∉ (Walk.cons h₂ t).support := by
    rw [Walk.support_cons, List.mem_cons]
    rintro (rfl | hxt)
    · exact G.irrefl hbx
    · -- a shortcut from `b` to `x` would close a cycle
      have hP : (Walk.cons h₂ (t.takeUntil x hxt)).IsPath :=
        (Walk.cons_isPath_iff _ _).2
          ⟨ht.takeUntil hxt, fun hc => hbt (Walk.support_takeUntil_subset t hxt hc)⟩
      have hedge : s(x, b) ∉ (Walk.cons h₂ (t.takeUntil x hxt)).edges := by
        rw [Walk.edges_cons, List.mem_cons]
        rintro (heq | hmem)
        · rcases Sym2.eq_iff.1 heq with ⟨e1, e2⟩ | ⟨e1, e2⟩
          · exact G.irrefl (e1 ▸ hbx)
          · exact hxc e1
        · exact hbt (Walk.support_takeUntil_subset t hxt
            (Walk.snd_mem_support_of_mem_edges _ hmem))
      exact hac _ (SimpleGraph.Path.cons_isCycle ⟨_, hP⟩ hbx.symm hedge)
  -- replace the initial vertex by `x`
  have hP : IsTPath G T (Walk.cons hbx.symm (Walk.cons h₂ t)) := by
    refine ⟨(Walk.cons_isPath_iff _ _).2 ⟨hr, hxr⟩, ?_⟩
    intro z hz
    rw [Walk.support_cons, List.mem_cons] at hz
    rcases hz with rfl | hz
    · exact hxT
    · exact hp.2 z (by rw [Walk.support_cons, List.mem_cons]; exact Or.inr hz)
  have hlen' : (Walk.cons hbx.symm (Walk.cons h₂ t)).length = maxTPathLen G T := by
    rw [Walk.length_cons] at hlen ⊢
    exact hlen
  exact fun y hy hadj => max_nbr_eq hac hbx.symm hP hlen' hy hadj

/-- Every nonempty vertex set of a forest contains a vertex with at most one neighbour. -/
lemma exists_leaf (hac : G.IsAcyclic) (hT : T.Nonempty) :
    ∃ a ∈ T, ∃ b : V, ∀ y ∈ T, G.Adj a y → y = b := by
  obtain ⟨u, v, p, hp, hlen⟩ := tPathLen_max (G := G) hT
  cases p with
  | nil =>
    refine ⟨u, hp.2 u (by simp), u, fun y hy hadj => ?_⟩
    exfalso
    have hpath : IsTPath G T (Walk.cons hadj Walk.nil) := by
      refine ⟨(Walk.cons_isPath_iff _ _).2 ⟨Walk.IsPath.nil, by simp [hadj.ne]⟩, ?_⟩
      intro z hz
      rw [Walk.support_cons, Walk.support_nil, List.mem_cons, List.mem_singleton] at hz
      rcases hz with rfl | rfl
      · exact hp.2 z (by simp)
      · exact hy
    have hlt := hpath.length_lt
    rw [Walk.length_cons, Walk.length_nil] at hlt
    refine not_tPathLen_of_gt (m := 1) ?_ (by omega) ⟨_, _, _, hpath, by simp⟩
    rw [← hlen]
    simp
  | cons h q =>
    exact ⟨u, hp.2 u (by simp), _, fun y hy hadj => max_nbr_eq hac h hp hlen hy hadj⟩

/-- In a forest, every vertex set has independence number at least half its size. -/
lemma card_le_two_alpha (hac : G.IsAcyclic) :
    ∀ (n : ℕ) (T : Finset V), T.card ≤ n → T.card ≤ 2 * alpha G T := by
  intro n
  induction n with
  | zero => intro T hT; omega
  | succ n ih =>
    intro T hT
    rcases T.eq_empty_or_nonempty with rfl | hne
    · simp
    obtain ⟨a, haT, b, hb⟩ := exists_leaf hac hne
    have hsub : ((T.erase a).erase b) ⊆ T :=
      (Finset.erase_subset _ _).trans (Finset.erase_subset _ _)
    have hcard : T.card ≤ ((T.erase a).erase b).card + 2 := by
      have h1 := Finset.pred_card_le_card_erase (s := T) (a := a)
      have h2 := Finset.pred_card_le_card_erase (s := T.erase a) (a := b)
      have h3 : 1 ≤ T.card := Finset.card_pos.2 hne
      omega
    have hlt : ((T.erase a).erase b).card < T.card := by
      have h1 : ((T.erase a).erase b).card ≤ (T.erase a).card :=
        Finset.card_le_card (Finset.erase_subset _ _)
      have h2 : (T.erase a).card < T.card := Finset.card_erase_lt_of_mem haT
      omega
    have halpha : alpha G ((T.erase a).erase b) + 1 ≤ alpha G T := by
      obtain ⟨M, hM, hMc⟩ := exists_max_indFam G ((T.erase a).erase b)
      have hMsub : M ⊆ (T.erase a).erase b := (mem_indFam.1 hM).1
      have haM : a ∉ M := fun hh =>
        (Finset.mem_erase.1 (Finset.mem_of_mem_erase (hMsub hh))).1 rfl
      have hins : insert a M ∈ indFam G T := by
        refine mem_indFam.2 ⟨Finset.insert_subset haT (hMsub.trans hsub), ?_⟩
        intro y hy z hz hyz
        simp only [Finset.coe_insert, Set.mem_insert_iff, Finset.mem_coe] at hy hz
        rcases hy with rfl | hy <;> rcases hz with rfl | hz
        · exact absurd rfl hyz
        · exact fun hadj => (Finset.mem_erase.1 (hMsub hz)).1
            (hb z (hsub (hMsub hz)) hadj)
        · exact fun hadj => (Finset.mem_erase.1 (hMsub hy)).1
            (hb y (hsub (hMsub hy)) hadj.symm)
        · exact (mem_indFam.1 hM).2 (by exact_mod_cast hy) (by exact_mod_cast hz) hyz
      have hle := card_le_alpha hins
      rw [Finset.card_insert_of_notMem haM, hMc] at hle
      omega
    have := ih ((T.erase a).erase b) (by omega)
    omega

end FivefoldForest
