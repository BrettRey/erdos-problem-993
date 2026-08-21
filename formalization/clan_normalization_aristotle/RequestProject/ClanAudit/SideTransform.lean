import RequestProject.ClanAudit.Prefix

/-!
# The canonical local transformation on one side of the tree

The published spider transformation, transported from the abstract local model
`HubState` of `RequestProject/ClanAudit/LocalMap.lean` onto the clan maps of an actual
spider, *with the published deterministic choices*: the arm `A = i₀` is the one with the
shortest odd active prefix, ties broken by the arm order, and `B = i₁` is the next such
arm.  Both are read off from the minimum of the key `i + (pref s i) * k`, which orders
the arms lexicographically by (prefix length, index).

* `ActiveSide s` — the hub carries multiplicity one, the initial positive run of every
  arm consists of ones, and at least two arms have odd active prefix;
* `canonical_spec` — the deterministic choices are well defined;
* `transf_admissible` — an active side satisfies the hypotheses of the published local
  map, so `localMapP2_injective` and `localMapP2_preserves_total_order` apply;
* `transf_injective` — the transformation is injective on active sides;
* `sum_transf` — it preserves the total order of the clan map.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph

variable {k : ℕ} {len : Fin k → ℕ}

/-! ### The abstract local state of a side -/

/-- The abstract local state (`HubState`) of a clan map on a spider. -/
def toHubState (s : SpiderV k len → ℕ) : HubState where
  hub := s none
  nArms := k
  armLen := fun i => if h : i < k then len ⟨i, h⟩ else 0
  mult := fun i j => if h : i < k then armVal s ⟨i, h⟩ j else 0

theorem toHubState_mult (s : SpiderV k len → ℕ) (i : Fin k) (j : ℕ) :
    (toHubState s).mult i.val j = armVal s i j := by
  rw [toHubState]
  simp only [dif_pos i.isLt]

theorem toHubState_injective {s s' : SpiderV k len → ℕ} (h : toHubState s = toHubState s') :
    s = s' := by
  funext x
  match x with
  | none => exact congrArg HubState.hub h
  | some ⟨i, j⟩ =>
      have := congrFun (congrFun (congrArg HubState.mult h) i.val) j.val
      rw [toHubState_mult, toHubState_mult, armVal_of_lt, armVal_of_lt] at this
      exact this

/-! ### The deterministic choice of the two arms -/

/-- The key of arm `i`, ordering the arms lexicographically by (prefix length, index). -/
def armKey (s : SpiderV k len → ℕ) (i : Fin k) : ℕ := i.val + pref s i * k

theorem armKey_div (s : SpiderV k len → ℕ) (i : Fin k) : armKey s i / k = pref s i := by
  have hk : 0 < k := by have := i.isLt; omega
  rw [armKey, Nat.add_mul_div_right _ _ hk, Nat.div_eq_of_lt i.isLt, Nat.zero_add]

theorem armKey_mod (s : SpiderV k len → ℕ) (i : Fin k) : armKey s i % k = i.val := by
  rw [armKey, Nat.add_mul_mod_self_right, Nat.mod_eq_of_lt i.isLt]

theorem armKey_injective (s : SpiderV k len → ℕ) {i i' : Fin k} (h : armKey s i = armKey s i') :
    i = i' := by
  exact Fin.ext (by rw [← armKey_mod s i, ← armKey_mod s i', h])

/-- The set of keys of the arms with odd active prefix. -/
def keySet (s : SpiderV k len → ℕ) : Set ℕ :=
  {n | ∃ i : Fin k, pref s i % 2 = 1 ∧ armKey s i = n}

/-- The key of the arm `A`: the shortest odd active prefix, ties by arm order. -/
noncomputable def key0 (s : SpiderV k len → ℕ) : ℕ := sInf (keySet s)

/-- The key of the arm `B`: the next shortest odd active prefix. -/
noncomputable def key1 (s : SpiderV k len → ℕ) : ℕ := sInf (keySet s \ {key0 s})

/-- The index of the arm `A`. -/
noncomputable def idx0 (s : SpiderV k len → ℕ) : ℕ := key0 s % k

/-- The index of the arm `B`. -/
noncomputable def idx1 (s : SpiderV k len → ℕ) : ℕ := key1 s % k

/-- The length of the shortest odd active prefix. -/
noncomputable def plen (s : SpiderV k len → ℕ) : ℕ := key0 s / k

/-- Admissibility of the initial positive run of each arm: it consists of ones.  Any
clan map violating this has a `2` next to a positive vertex, hence weight zero. -/
def AdmRun (s : SpiderV k len → ℕ) : Prop :=
  ∀ (i : Fin k) (j : ℕ), (∀ j' ≤ j, 0 < armVal s i j') → armVal s i j = 1

theorem AdmRun.stop {s : SpiderV k len → ℕ} (h : AdmRun s) (i : Fin k)
    (hlt : pref s i < len i) : armVal s i (pref s i) = 0 := by
  by_contra hne
  refine pref_stop s i hlt (h i (pref s i) fun j' hj' => ?_)
  rcases lt_or_eq_of_le hj' with hlt' | rfl
  · rw [pref_ones s i hlt']; omega
  · omega

/-- **An active side.**  The hub is active, the initial positive run of every arm
consists of ones, and at least two arms have odd active prefix — exactly the situation
in which the published transformation applies. -/
structure ActiveSide (s : SpiderV k len → ℕ) : Prop where
  hub : s none = 1
  run : AdmRun s
  two_odd : 2 ≤ pNum s

theorem ActiveSide.pos_arms {s : SpiderV k len → ℕ} (h : ActiveSide s) : 0 < k := by
  by_contra hk
  have : k = 0 := by omega
  subst this
  have hz : pNum s = 0 := by
    rw [pNum]
    simp
  have := h.two_odd
  omega

/-- **The deterministic choice is well defined.** -/
theorem canonical_spec {s : SpiderV k len → ℕ} (h : ActiveSide s) :
    ∃ i₀ i₁ : Fin k, i₀.val = idx0 s ∧ i₁.val = idx1 s ∧ pref s i₀ = plen s ∧
      pref s i₀ % 2 = 1 ∧ pref s i₁ % 2 = 1 ∧ i₀ ≠ i₁ ∧ plen s ≤ pref s i₁ := by
  classical
  have h2 := h.two_odd
  rw [pNum] at h2
  have h2' : 1 < (Finset.univ.filter (fun i : Fin k => pref s i % 2 = 1)).card := by omega
  obtain ⟨a, ha, b, hb, hab⟩ := Finset.one_lt_card.mp h2'
  rw [Finset.mem_filter] at ha hb
  have hne : keySet s ≠ ∅ := by
    intro hc
    have : armKey s a ∈ keySet s := ⟨a, ha.2, rfl⟩
    rw [hc] at this
    exact this
  have hmem0 : key0 s ∈ keySet s := Nat.sInf_mem (Set.nonempty_iff_ne_empty.mpr hne)
  obtain ⟨i₀, hi₀odd, hi₀key⟩ := hmem0
  -- the second smallest key
  have hne2 : (keySet s \ {key0 s}).Nonempty := by
    by_cases hca : armKey s a = key0 s
    · refine ⟨armKey s b, ⟨b, hb.2, rfl⟩, ?_⟩
      simp only [Set.mem_singleton_iff]
      intro hcb
      exact hab (armKey_injective s (by rw [hca, hcb]))
    · exact ⟨armKey s a, ⟨a, ha.2, rfl⟩, hca⟩
  have hmem1 : key1 s ∈ keySet s \ {key0 s} := Nat.sInf_mem hne2
  obtain ⟨⟨i₁, hi₁odd, hi₁key⟩, hkne⟩ := hmem1
  have hle : key0 s ≤ key1 s := Nat.sInf_le ⟨i₁, hi₁odd, hi₁key⟩
  have hidx0 : i₀.val = idx0 s := by rw [idx0, ← hi₀key, armKey_mod]
  have hidx1 : i₁.val = idx1 s := by rw [idx1, ← hi₁key, armKey_mod]
  have hplen : pref s i₀ = plen s := by rw [plen, ← hi₀key, armKey_div]
  refine ⟨i₀, i₁, hidx0, hidx1, hplen, hi₀odd, hi₁odd, ?_, ?_⟩
  · intro hc
    exact hkne (by rw [Set.mem_singleton_iff, ← hi₁key, ← hc, hi₀key])
  · rw [← hplen]
    have e0 : armKey s i₀ = i₀.val + pref s i₀ * k := rfl
    have e1 : armKey s i₁ = i₁.val + pref s i₁ * k := rfl
    have h0 := i₀.isLt
    have h1 := i₁.isLt
    rw [← hi₀key, ← hi₁key, e0, e1] at hle
    by_contra hcon
    have : pref s i₁ + 1 ≤ pref s i₀ := by omega
    have : (pref s i₁ + 1) * k ≤ pref s i₀ * k := Nat.mul_le_mul_right k this
    have : pref s i₁ * k + k ≤ pref s i₀ * k := by rw [add_mul, one_mul] at this; omega
    omega

/-! ### The transformation -/

/-- The published local transformation, as a map on the clan maps of a spider: the hub
is switched off, the `L` vertices of the active prefix of the arm `A` are rewritten
`2,0,2,0,…,2`, and the first `L-1` vertices of the arm `B` are rewritten `2,0,…,2,0`. -/
def transfSide (s : SpiderV k len → ℕ) (i₀ i₁ L : ℕ) : SpiderV k len → ℕ := fun x =>
  match x with
  | none => 0
  | some p =>
      if p.1.val = i₀ ∧ p.2.val < L then altVal p.2.val
      else if p.1.val = i₁ ∧ p.2.val < L - 1 then altVal p.2.val
      else s (some p)

/-- The canonical transformation of an active side. -/
noncomputable def transf (s : SpiderV k len → ℕ) : SpiderV k len → ℕ :=
  transfSide s (idx0 s) (idx1 s) (plen s)

theorem transfSide_none (s : SpiderV k len → ℕ) (i₀ i₁ L : ℕ) :
    transfSide s i₀ i₁ L none = 0 := rfl

theorem armVal_transfSide (s : SpiderV k len → ℕ) (i₀ i₁ L : ℕ) (i : Fin k) (j : ℕ)
    (hj : j < len i) :
    armVal (transfSide s i₀ i₁ L) i j
      = if i.val = i₀ ∧ j < L then altVal j
        else if i.val = i₁ ∧ j < L - 1 then altVal j
        else armVal s i j := by
  rw [armVal, dif_pos hj, armVal, dif_pos hj]
  rfl

theorem armVal_transfSide_other (s : SpiderV k len → ℕ) {i₀ i₁ L : ℕ} {i : Fin k}
    (h0 : i.val ≠ i₀) (h1 : i.val ≠ i₁) (j : ℕ) :
    armVal (transfSide s i₀ i₁ L) i j = armVal s i j := by
  by_cases hj : j < len i
  · rw [armVal_transfSide s i₀ i₁ L i j hj, if_neg (by tauto), if_neg (by tauto)]
  · rw [armVal_of_ge _ _ (by omega), armVal_of_ge _ _ (by omega)]

/-- The transformation of a spider clan map is the published local map of its abstract
local state. -/
theorem toHubState_transfSide (s : SpiderV k len → ℕ) (i₀ i₁ : Fin k) (L : ℕ)
    (hL0 : L ≤ len i₀) (hL1 : L - 1 ≤ len i₁) :
    toHubState (transfSide s i₀.val i₁.val L) = localMapP2 (toHubState s) i₀.val i₁.val L := by
  refine HubState.ext rfl rfl rfl ?_
  funext i j
  rw [mult_localMapP2]
  by_cases hi : i < k
  · have hmi : (toHubState (transfSide s i₀ i₁ L)).mult i j
        = armVal (transfSide s i₀ i₁ L) ⟨i, hi⟩ j := toHubState_mult _ ⟨i, hi⟩ j
    have hms : (toHubState s).mult i j = armVal s ⟨i, hi⟩ j := toHubState_mult _ ⟨i, hi⟩ j
    rw [hmi, hms]
    by_cases hj : j < len ⟨i, hi⟩
    · rw [armVal_transfSide s i₀ i₁ L ⟨i, hi⟩ j hj]
    · have hz : armVal s ⟨i, hi⟩ j = 0 := armVal_of_ge _ _ (by omega)
      have hz' : armVal (transfSide s i₀ i₁ L) ⟨i, hi⟩ j = 0 := armVal_of_ge _ _ (by omega)
      rw [hz, hz']
      by_cases hc0 : i = i₀.val ∧ j < L
      · exfalso
        have hii : (⟨i, hi⟩ : Fin k) = i₀ := Fin.ext hc0.1
        rw [hii] at hj
        have := hc0.2
        omega
      · rw [if_neg hc0]
        by_cases hc1 : i = i₁.val ∧ j < L - 1
        · exfalso
          have hii : (⟨i, hi⟩ : Fin k) = i₁ := Fin.ext hc1.1
          rw [hii] at hj
          have := hc1.2
          omega
        · rw [if_neg hc1]
  · have hmi : (toHubState (transfSide s i₀ i₁ L)).mult i j = 0 := by
      rw [toHubState]; simp only [dif_neg hi]
    have hms : (toHubState s).mult i j = 0 := by
      rw [toHubState]; simp only [dif_neg hi]
    rw [hmi, hms, if_neg (fun hc => hi (by rw [hc.1]; exact i₀.isLt)),
      if_neg (fun hc => hi (by rw [hc.1]; exact i₁.isLt))]

/-- **An active side is admissible for the published local map.** -/
theorem transf_admissible {s : SpiderV k len → ℕ} (h : ActiveSide s) :
    Admissible (toHubState s) (idx0 s) (idx1 s) (plen s) := by
  obtain ⟨i₀, i₁, hidx0, hidx1, hplen, hodd0, hodd1, hne, hle⟩ := canonical_spec h
  have hL : plen s % 2 = 1 := by rw [← hplen]; exact hodd0
  have hhub : (toHubState s).hub = 1 := h.hub
  refine ⟨hhub, hL, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro i j hij
    by_cases hi : i < k
    · have harm : (toHubState s).armLen i = len ⟨i, hi⟩ := by
        rw [toHubState]; simp only [dif_pos hi]
      rw [harm] at hij
      exact (toHubState_mult s ⟨i, hi⟩ j).trans (armVal_of_ge _ _ hij)
    · rw [toHubState]; simp only [dif_neg hi]
  · intro i j hij
    by_cases hi : i < k
    · refine (toHubState_mult s ⟨i, hi⟩ j).trans (h.run ⟨i, hi⟩ j fun j' hj' => ?_)
      have := hij j' hj'
      rwa [toHubState_mult s ⟨i, hi⟩ j'] at this
    · exfalso
      have := hij 0 (by omega)
      rw [toHubState] at this
      simp only [dif_neg hi] at this
      omega
  · intro j hj
    rw [← hidx0, toHubState_mult]
    exact pref_ones s i₀ (by omega)
  · rw [← hidx0, toHubState_mult, ← hplen]
    by_cases hlt : pref s i₀ < len i₀
    · exact h.run.stop i₀ hlt
    · exact armVal_of_ge _ _ (by have := pref_le s i₀; omega)
  · intro j hj
    rw [← hidx1, toHubState_mult]
    exact pref_ones s i₁ (by omega)
  · rw [← hidx0, ← hidx1]
    intro hc
    exact hne (Fin.ext hc).symm
  · rw [← hidx0]; exact i₀.isLt
  · rw [← hidx1]; exact i₁.isLt

/-- **The transformation is injective on active sides.**  This is the global form of
`localMapP2_injective`. -/
theorem toHubState_transf {s : SpiderV k len → ℕ} (h : ActiveSide s) :
    toHubState (transf s) = localMapP2 (toHubState s) (idx0 s) (idx1 s) (plen s) := by
  obtain ⟨i₀, i₁, hidx0, hidx1, hplen, hodd0, hodd1, hne, hle⟩ := canonical_spec h
  rw [transf, ← hidx0, ← hidx1]
  exact toHubState_transfSide s i₀ i₁ (plen s) (by rw [← hplen]; exact pref_le s i₀)
    (by have := pref_le s i₁; omega)

theorem transf_injective {s s' : SpiderV k len → ℕ} (h : ActiveSide s) (h' : ActiveSide s')
    (heq : transf s = transf s') : s = s' := by
  have hmap : localMapP2 (toHubState s) (idx0 s) (idx1 s) (plen s)
      = localMapP2 (toHubState s') (idx0 s') (idx1 s') (plen s') := by
    rw [← toHubState_transf h, ← toHubState_transf h', heq]
  exact toHubState_injective
    (localMapP2_injective (transf_admissible h) (transf_admissible h') hmap)

/-! ### Preservation of the total order -/

theorem sum_side_eq_total (s : SpiderV k len → ℕ) :
    ∑ x : SpiderV k len, s x = (toHubState s).total := by
  classical
  rw [HubState.total, Fintype.sum_option, Fintype.sum_sigma]
  congr 1
  rw [toHubState]
  simp only
  rw [Finset.sum_range fun i => ∑ j ∈ Finset.range (if h : i < k then len ⟨i, h⟩ else 0),
    (if h : i < k then armVal s ⟨i, h⟩ j else 0)]
  refine Finset.sum_congr rfl fun i _ => ?_
  simp only [dif_pos i.isLt]
  have : (⟨i.val, i.isLt⟩ : Fin k) = i := Fin.ext rfl
  rw [this, Finset.sum_range fun j => armVal s i j]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [armVal_of_lt]

/-- **The transformation preserves the total order of the clan map.**  This is the
global form of `localMapP2_preserves_total_order`. -/
theorem sum_transf {s : SpiderV k len → ℕ} (h : ActiveSide s) :
    ∑ x : SpiderV k len, transf s x = ∑ x : SpiderV k len, s x := by
  rw [sum_side_eq_total, sum_side_eq_total, toHubState_transf h,
    localMapP2_preserves_total_order (transf_admissible h)]

end ClanAudit
