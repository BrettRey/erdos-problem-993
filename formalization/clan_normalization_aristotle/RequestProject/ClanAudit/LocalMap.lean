import Mathlib

/-!
# The published spider transformation at `p = 2`

This file formalizes the local pairing map of the published spider transformation and
proves the two purely combinatorial claims the request asks for at `p = 2`:

* `localMapP2_preserves_total_order` — the transformation preserves the total order of
  the clan graph (`∑ multiplicities`);
* `localMapP2_injective` — the transformation is injective on admissible local states.

## Model

A *local state at a hub* records the hub multiplicity, the number of pendant arms, the
length of each arm, and the multiplicity of the `j`-th vertex of the `i`-th arm.  This
is exactly the clan-map data the transformation reads and writes.

The transformation is parameterized by the two arm indices `i₀` (the arm `A`, whose odd
positive prefix has length `L`) and `i₁` (the arm `B`), together with `L`.  The
published recipe picks `i₀` as the first shortest odd prefix and `i₁` as the first other
odd prefix; here those choices are inputs constrained by `Admissible`, so that the
injectivity statement proved below is *stronger* than injectivity for one particular
tie-breaking rule: it shows that two admissible states with any admissible parameter
choices that produce the same output must be equal.

Note that the published case distinction `L = 1` versus `L ≥ 3` is unnecessary: with
`L = 1` the recipe for `B` rewrites the first `L - 1 = 0` vertices, i.e. does nothing,
and the recipe for `A` sets its unique prefix vertex to `2`, which is exactly the
published `L = 1` rule.
-/

namespace ClanAudit

open Finset

/-- The multiplicity written at position `j` of a rewritten run: `2, 0, 2, 0, …`. -/
def altVal (j : ℕ) : ℕ := if j % 2 = 0 then 2 else 0

theorem altVal_ne_one (j : ℕ) : altVal j ≠ 1 := by
  unfold altVal; split <;> omega

theorem altVal_zero : altVal 0 = 2 := rfl

theorem sum_altVal (n : ℕ) : ∑ j ∈ Finset.range n, altVal j = n + n % 2 := by
  induction n with
  | zero => simp
  | succ n ih =>
      rw [Finset.sum_range_succ, ih]
      unfold altVal
      rcases Nat.even_or_odd n with he | ho
      · have h1 : n % 2 = 0 := Nat.even_iff.mp he
        simp only [h1, if_pos]
        omega
      · have h1 : n % 2 = 1 := Nat.odd_iff.mp ho
        simp only [h1]
        rw [if_neg (by omega)]
        omega

/-- A local state at a hub: the hub multiplicity, the number of pendant arms, the length
of each arm, and the multiplicity of each arm vertex. -/
@[ext] structure HubState where
  hub : ℕ
  nArms : ℕ
  armLen : ℕ → ℕ
  mult : ℕ → ℕ → ℕ

/-- The total order of the clan graph of a local state. -/
def HubState.total (s : HubState) : ℕ :=
  s.hub + ∑ i ∈ Finset.range s.nArms, ∑ j ∈ Finset.range (s.armLen i), s.mult i j

/-- The published local transformation at `p = 2`: the hub multiplicity is set to zero,
the `L` vertices of the arm `A = i₀` are rewritten as `2,0,2,0,…,2`, and the first
`L-1` vertices of the arm `B = i₁` are rewritten as `2,0,2,0,…,2,0`. -/
def localMapP2 (s : HubState) (i₀ i₁ L : ℕ) : HubState where
  hub := 0
  nArms := s.nArms
  armLen := s.armLen
  mult i j :=
    if i = i₀ ∧ j < L then altVal j
    else if i = i₁ ∧ j < L - 1 then altVal j
    else s.mult i j

@[simp] theorem mult_localMapP2 (s : HubState) (i₀ i₁ L i j : ℕ) :
    (localMapP2 s i₀ i₁ L).mult i j
      = if i = i₀ ∧ j < L then altVal j
        else if i = i₁ ∧ j < L - 1 then altVal j
        else s.mult i j := rfl

/-- Admissibility of a local state together with a choice of the two arms and the length
of the shortest odd prefix.  These are exactly the hypotheses of the published
transformation: the hub is active with multiplicity one, every arm's positive prefix
consists of multiplicities one, the arm `A = i₀` has odd positive prefix exactly `L`, and
the arm `B = i₁` is a different arm whose positive prefix is at least as long. -/
structure Admissible (s : HubState) (i₀ i₁ L : ℕ) : Prop where
  hub_one : s.hub = 1
  L_odd : L % 2 = 1
  /-- Multiplicities vanish beyond the end of an arm. -/
  support : ∀ i j, s.armLen i ≤ j → s.mult i j = 0
  /-- Every arm's positive prefix consists of multiplicities one. -/
  ones : ∀ i j, (∀ j' ≤ j, 0 < s.mult i j') → s.mult i j = 1
  /-- The positive prefix of `A` has length exactly `L`. -/
  A_ones : ∀ j < L, s.mult i₀ j = 1
  A_stop : s.mult i₀ L = 0
  /-- The positive prefix of `B` has length at least `L`. -/
  B_ones : ∀ j < L, s.mult i₁ j = 1
  B_ne : i₁ ≠ i₀
  i₀_lt : i₀ < s.nArms
  i₁_lt : i₁ < s.nArms

theorem Admissible.L_pos {s : HubState} {i₀ i₁ L : ℕ} (h : Admissible s i₀ i₁ L) : 0 < L := by
  have := h.L_odd; omega

theorem Admissible.head_ne_two {s : HubState} {i₀ i₁ L : ℕ} (h : Admissible s i₀ i₁ L) (i : ℕ) :
    s.mult i 0 ≠ 2 := by
  intro hc
  have : s.mult i 0 = 1 := h.ones i 0 (fun j' hj' => by
    have : j' = 0 := by omega
    subst this; omega)
  omega

theorem Admissible.A_len {s : HubState} {i₀ i₁ L : ℕ} (h : Admissible s i₀ i₁ L) :
    L ≤ s.armLen i₀ := by
  by_contra hc
  have h1 : s.mult i₀ (L - 1) = 0 := h.support i₀ (L - 1) (by omega)
  have h2 : s.mult i₀ (L - 1) = 1 := h.A_ones (L - 1) (by have := h.L_pos; omega)
  omega

theorem Admissible.B_len {s : HubState} {i₀ i₁ L : ℕ} (h : Admissible s i₀ i₁ L) :
    L ≤ s.armLen i₁ := by
  by_contra hc
  have h1 : s.mult i₁ (L - 1) = 0 := h.support i₁ (L - 1) (by omega)
  have h2 : s.mult i₁ (L - 1) = 1 := h.B_ones (L - 1) (by have := h.L_pos; omega)
  omega

/-! ### The transformation preserves the total clan order -/

private theorem sum_arm_A {s : HubState} {i₀ i₁ L : ℕ} (h : Admissible s i₀ i₁ L) :
    ∑ j ∈ Finset.range (s.armLen i₀), (localMapP2 s i₀ i₁ L).mult i₀ j
      = (∑ j ∈ Finset.range (s.armLen i₀), s.mult i₀ j) + 1 := by
  have hL := h.A_len
  rw [← Finset.sum_range_add_sum_Ico _ hL, ← Finset.sum_range_add_sum_Ico _ hL]
  have e1 : ∑ j ∈ Finset.range L, (localMapP2 s i₀ i₁ L).mult i₀ j = L + 1 := by
    have hc : ∀ j ∈ Finset.range L, (localMapP2 s i₀ i₁ L).mult i₀ j = altVal j := by
      intro j hj
      simp only [Finset.mem_range] at hj
      simp only [mult_localMapP2]
      rw [if_pos (show _ ∧ j < _ from ⟨by trivial, hj⟩)]
    rw [Finset.sum_congr rfl hc, sum_altVal, h.L_odd]
  have e2 : ∑ j ∈ Finset.range L, s.mult i₀ j = L := by
    rw [Finset.sum_congr rfl (fun j hj => h.A_ones j (Finset.mem_range.mp hj))]
    simp
  have e3 : ∑ j ∈ Finset.Ico L (s.armLen i₀), (localMapP2 s i₀ i₁ L).mult i₀ j
      = ∑ j ∈ Finset.Ico L (s.armLen i₀), s.mult i₀ j := by
    refine Finset.sum_congr rfl fun j hj => ?_
    simp only [Finset.mem_Ico] at hj
    simp only [mult_localMapP2]
    rw [if_neg (by rintro ⟨-, hlt⟩; omega), if_neg (by rintro ⟨-, hlt⟩; omega)]
  rw [e1, e2, e3]
  ring

private theorem sum_arm_B {s : HubState} {i₀ i₁ L : ℕ} (h : Admissible s i₀ i₁ L) :
    ∑ j ∈ Finset.range (s.armLen i₁), (localMapP2 s i₀ i₁ L).mult i₁ j
      = ∑ j ∈ Finset.range (s.armLen i₁), s.mult i₁ j := by
  have hL : L - 1 ≤ s.armLen i₁ := by have := h.B_len; omega
  rw [← Finset.sum_range_add_sum_Ico _ hL, ← Finset.sum_range_add_sum_Ico _ hL]
  have e1 : ∑ j ∈ Finset.range (L - 1), (localMapP2 s i₀ i₁ L).mult i₁ j = L - 1 := by
    have hc : ∀ j ∈ Finset.range (L - 1), (localMapP2 s i₀ i₁ L).mult i₁ j = altVal j := by
      intro j hj
      simp only [Finset.mem_range] at hj
      simp only [mult_localMapP2]
      rw [if_neg (by rintro ⟨he, -⟩; exact h.B_ne he), if_pos (show _ ∧ j < _ from ⟨by trivial, hj⟩)]
    rw [Finset.sum_congr rfl hc, sum_altVal]
    have : (L - 1) % 2 = 0 := by have := h.L_odd; omega
    omega
  have e2 : ∑ j ∈ Finset.range (L - 1), s.mult i₁ j = L - 1 := by
    rw [Finset.sum_congr rfl (fun j hj => h.B_ones j (by
      have := Finset.mem_range.mp hj; omega))]
    simp
  have e3 : ∑ j ∈ Finset.Ico (L - 1) (s.armLen i₁), (localMapP2 s i₀ i₁ L).mult i₁ j
      = ∑ j ∈ Finset.Ico (L - 1) (s.armLen i₁), s.mult i₁ j := by
    refine Finset.sum_congr rfl fun j hj => ?_
    simp only [Finset.mem_Ico] at hj
    simp only [mult_localMapP2]
    rw [if_neg (by rintro ⟨he, -⟩; exact h.B_ne he), if_neg (by rintro ⟨-, hlt⟩; omega)]
  rw [e1, e2, e3]

/-- **The local transformation preserves the total clan-graph order.** -/
theorem localMapP2_preserves_total_order {s : HubState} {i₀ i₁ L : ℕ}
    (h : Admissible s i₀ i₁ L) : (localMapP2 s i₀ i₁ L).total = s.total := by
  classical
  have hmem₀ : i₀ ∈ Finset.range s.nArms := Finset.mem_range.mpr h.i₀_lt
  have hmem₁ : i₁ ∈ (Finset.range s.nArms).erase i₀ :=
    Finset.mem_erase.mpr ⟨h.B_ne, Finset.mem_range.mpr h.i₁_lt⟩
  have hrest : ∀ i ∈ ((Finset.range s.nArms).erase i₀).erase i₁,
      ∑ j ∈ Finset.range (s.armLen i), (localMapP2 s i₀ i₁ L).mult i j
        = ∑ j ∈ Finset.range (s.armLen i), s.mult i j := by
    intro i hi
    simp only [Finset.mem_erase] at hi
    refine Finset.sum_congr rfl fun j _ => ?_
    simp only [mult_localMapP2]
    rw [if_neg (by rintro ⟨he, -⟩; exact hi.2.1 he), if_neg (by rintro ⟨he, -⟩; exact hi.1 he)]
  have hn : (localMapP2 s i₀ i₁ L).nArms = s.nArms := rfl
  have hlen : (localMapP2 s i₀ i₁ L).armLen = s.armLen := rfl
  have hhub : (localMapP2 s i₀ i₁ L).hub = 0 := rfl
  rw [HubState.total, HubState.total, hn, hlen, hhub, h.hub_one,
    ← Finset.add_sum_erase _ _ hmem₀, ← Finset.add_sum_erase _ _ hmem₀,
    ← Finset.add_sum_erase _ _ hmem₁, ← Finset.add_sum_erase _ _ hmem₁,
    Finset.sum_congr rfl hrest, sum_arm_A h, sum_arm_B h]
  ring

/-! ### Injectivity -/

/-- On the arm `A`, the image is the alternating run below position `L`. -/
private theorem A_prefix {s : HubState} {i₀ i₁ L : ℕ} {j : ℕ} (hj : j < L) :
    (localMapP2 s i₀ i₁ L).mult i₀ j = altVal j := by
  simp only [mult_localMapP2]
  rw [if_pos (show _ ∧ j < _ from ⟨by trivial, hj⟩)]

/-- On the arm `B`, the image is the alternating run below position `L - 1`. -/
private theorem B_prefix {s : HubState} {i₀ i₁ L : ℕ} (h : Admissible s i₀ i₁ L) {j : ℕ}
    (hj : j < L - 1) : (localMapP2 s i₀ i₁ L).mult i₁ j = altVal j := by
  simp only [mult_localMapP2]
  rw [if_neg (by rintro ⟨he, -⟩; exact h.B_ne he), if_pos (show _ ∧ j < _ from ⟨by trivial, hj⟩)]

/-- The marker distinguishing the arm `B`: its image has value `1` at position `L - 1`. -/
private theorem B_marker {s : HubState} {i₀ i₁ L : ℕ} (h : Admissible s i₀ i₁ L) :
    (localMapP2 s i₀ i₁ L).mult i₁ (L - 1) = 1 := by
  simp only [mult_localMapP2]
  rw [if_neg (by rintro ⟨he, -⟩; exact h.B_ne he), if_neg (by rintro ⟨-, hlt⟩; omega)]
  exact h.B_ones (L - 1) (by have := h.L_pos; omega)

/-- The arms whose image starts with multiplicity `2` are exactly `A`, and `B` when
`L ≥ 3`. -/
private theorem head2_iff {s : HubState} {i₀ i₁ L : ℕ} (h : Admissible s i₀ i₁ L) (i : ℕ) :
    (localMapP2 s i₀ i₁ L).mult i 0 = 2 ↔ (i = i₀ ∨ (3 ≤ L ∧ i = i₁)) := by
  have hLpos := h.L_pos
  constructor
  · intro h2
    by_contra hc
    push_neg at hc
    obtain ⟨hc0, hc1⟩ := hc
    simp only [mult_localMapP2] at h2
    rw [if_neg (by rintro ⟨he, -⟩; exact hc0 he)] at h2
    rw [if_neg ?_] at h2
    · exact h.head_ne_two i h2
    · rintro ⟨he, hlt⟩
      have hL3 : 3 ≤ L := by have := h.L_odd; omega
      exact hc1 hL3 he
  · rintro (rfl | ⟨hL3, rfl⟩)
    · rw [A_prefix hLpos, altVal_zero]
    · rw [B_prefix h (by omega), altVal_zero]

/-- **The local pairing map is injective on admissible states**, including at `p = 2`.
Two admissible local states (with any admissible choices of the two arms and of the
prefix length) that have the same image under the transformation are equal. -/
theorem localMapP2_injective {s s' : HubState} {i₀ i₁ L i₀' i₁' L' : ℕ}
    (h : Admissible s i₀ i₁ L) (h' : Admissible s' i₀' i₁' L')
    (heq : localMapP2 s i₀ i₁ L = localMapP2 s' i₀' i₁' L') : s = s' := by
  have hmult : ∀ i j, (localMapP2 s i₀ i₁ L).mult i j = (localMapP2 s' i₀' i₁' L').mult i j :=
    fun i j => by rw [heq]
  have hnArms : s.nArms = s'.nArms := by
    have := congrArg HubState.nArms heq; simpa [localMapP2] using this
  have hlen : s.armLen = s'.armLen := by
    have := congrArg HubState.armLen heq; simpa [localMapP2] using this
  have hhead : ∀ i, (i = i₀ ∨ (3 ≤ L ∧ i = i₁)) ↔ (i = i₀' ∨ (3 ≤ L' ∧ i = i₁')) := by
    intro i
    rw [← head2_iff h i, ← head2_iff h' i, hmult i 0]
  have key : i₀ = i₀' ∧ L = L' ∧ (L = 1 ∨ i₁ = i₁') := by
    rcases Nat.lt_or_ge L 3 with hL | hL3
    · have hL1 : L = 1 := by have := h.L_odd; omega
      rcases Nat.lt_or_ge L' 3 with hL' | hL3'
      · have hL1' : L' = 1 := by have := h'.L_odd; omega
        refine ⟨?_, by omega, Or.inl hL1⟩
        rcases (hhead i₀).mp (Or.inl rfl) with h1 | ⟨h2, -⟩
        · exact h1
        · omega
      · exfalso
        have e1 : i₀' = i₀ := by
          rcases (hhead i₀').mpr (Or.inl rfl) with h1 | ⟨h2, -⟩
          · exact h1
          · omega
        have e2 : i₁' = i₀ := by
          rcases (hhead i₁').mpr (Or.inr ⟨hL3', rfl⟩) with h1 | ⟨h2, -⟩
          · exact h1
          · omega
        exact h'.B_ne (e2.trans e1.symm)
    · rcases Nat.lt_or_ge L' 3 with hL' | hL3'
      · exfalso
        have hL1' : L' = 1 := by have := h'.L_odd; omega
        have e1 : i₀ = i₀' := by
          rcases (hhead i₀).mp (Or.inl rfl) with h1 | ⟨h2, -⟩
          · exact h1
          · omega
        have e2 : i₁ = i₀' := by
          rcases (hhead i₁).mp (Or.inr ⟨hL3, rfl⟩) with h1 | ⟨h2, -⟩
          · exact h1
          · omega
        exact h.B_ne (e2.trans e1.symm)
      · -- both prefixes have length at least three
        have hi0 : i₀ = i₀' := by
          by_contra hne
          have e1 : i₀ = i₁' := by
            rcases (hhead i₀).mp (Or.inl rfl) with h1 | ⟨-, h2⟩
            · exact absurd h1 hne
            · exact h2
          have e2 : i₁ = i₀' := by
            rcases (hhead i₁).mp (Or.inr ⟨hL3, rfl⟩) with h1 | ⟨-, h2⟩
            · exact h1
            · exact absurd (h2.trans e1.symm) h.B_ne
          have m1 : (localMapP2 s i₀ i₁ L).mult i₀ (L' - 1) = 1 := by
            rw [hmult, e1]; exact B_marker h'
          have hge : L ≤ L' - 1 := by
            by_contra hc
            rw [A_prefix (show L' - 1 < L by omega)] at m1
            exact altVal_ne_one _ m1
          have m2 : (localMapP2 s' i₀' i₁' L').mult i₀' (L - 1) = 1 := by
            rw [← hmult i₀' (L - 1), ← e2]; exact B_marker h
          have hge' : L' ≤ L - 1 := by
            by_contra hc
            rw [A_prefix (show L - 1 < L' by omega)] at m2
            exact altVal_ne_one _ m2
          omega
        have hi1 : i₁ = i₁' := by
          rcases (hhead i₁).mp (Or.inr ⟨hL3, rfl⟩) with h1 | ⟨-, h2⟩
          · exact absurd (h1.trans hi0.symm) h.B_ne
          · exact h2
        refine ⟨hi0, ?_, Or.inr hi1⟩
        have m1 : (localMapP2 s i₀ i₁ L).mult i₁ (L - 1) = 1 := B_marker h
        have m2 : (localMapP2 s' i₀' i₁' L').mult i₁' (L' - 1) = 1 := B_marker h'
        have g1 : L' - 1 ≤ L - 1 := by
          by_contra hc
          rw [hmult, hi1, B_prefix h' (show L - 1 < L' - 1 by omega)] at m1
          exact altVal_ne_one _ m1
        have g2 : L - 1 ≤ L' - 1 := by
          by_contra hc
          rw [← hmult, ← hi1, B_prefix h (show L' - 1 < L - 1 by omega)] at m2
          exact altVal_ne_one _ m2
        omega
  obtain ⟨hi0, hLL, hi1⟩ := key
  subst hi0
  subst hLL
  have hmulteq : ∀ i j, s.mult i j = s'.mult i j := by
    intro i j
    have hij := hmult i j
    simp only [mult_localMapP2] at hij
    by_cases hA : i = i₀ ∧ j < L
    · obtain ⟨rfl, hj⟩ := hA
      rw [h.A_ones j hj, h'.A_ones j hj]
    · rw [if_neg hA, if_neg hA] at hij
      rcases hi1 with hL1 | rfl
      · rw [if_neg (by rintro ⟨-, hlt⟩; omega), if_neg (by rintro ⟨-, hlt⟩; omega)] at hij
        exact hij
      · by_cases hB : i = i₁ ∧ j < L - 1
        · obtain ⟨rfl, hj⟩ := hB
          rw [h.B_ones j (by omega), h'.B_ones j (by omega)]
        · rw [if_neg hB, if_neg hB] at hij
          exact hij
  refine HubState.ext ?_ hnArms hlen (funext fun i => funext fun j => hmulteq i j)
  rw [h.hub_one, h'.hub_one]

end ClanAudit
