import RequestProject.ClanAudit.BlockCU
import RequestProject.ClanAudit.LocalMap

/-!
# Active prefixes and the weight of one side of the tree

This file sets up the *global model* of a clan map on one side (one spider) of the
adjacent two-hub tree, and computes its normalized weight exactly.

* `runLen` — the length of the initial run of a boolean predicate; `pref s i` is the
  length of the initial run of multiplicity-one vertices ("the active prefix") of arm
  `i` in the clan map `s`;
* `prefW p` — the weight contributed by an active prefix of length `p` once the hub is
  switched off: `1` if `p = 0`, `z + z⁻¹` if `p` is odd, `2` if `p` is even and positive;
* `tailW s i` — the weight of the part of arm `i` strictly beyond its active prefix;
* `Wpoly_arm` — every arm factors as `prefW (pref s i) * tailW s i`;
* `Wpoly_side_inactive` — a side whose hub is *not* active is the product of its arms;
* `Wpoly_side_active` — a side whose hub is active is `(z^r + z^(-r)) * ∏ tails` with
  `r = |1 - p|`, `p` the number of arms with odd active prefix.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph

/-! ### Initial runs -/

/-- The length of the initial run of `f` inside `[0, n)`. -/
def runLen : (ℕ → Bool) → ℕ → ℕ
  | _, 0 => 0
  | f, (n + 1) => if f 0 then runLen (fun j => f (j + 1)) n + 1 else 0

theorem runLen_le : ∀ (f : ℕ → Bool) (n : ℕ), runLen f n ≤ n
  | _, 0 => le_rfl
  | f, (n + 1) => by
      rw [runLen]
      by_cases h : f 0
      · rw [if_pos h]
        have := runLen_le (fun j => f (j + 1)) n
        omega
      · rw [if_neg h]
        omega

theorem runLen_true : ∀ (f : ℕ → Bool) (n : ℕ) (j : ℕ), j < runLen f n → f j = true
  | _, 0, j, hj => by rw [runLen] at hj; omega
  | f, (n + 1), j, hj => by
      rw [runLen] at hj
      by_cases h : f 0
      · rw [if_pos h] at hj
        match j with
        | 0 => exact h
        | (j + 1) => exact runLen_true (fun j => f (j + 1)) n j (by omega)
      · rw [if_neg h] at hj; omega

theorem runLen_stop : ∀ (f : ℕ → Bool) (n : ℕ), runLen f n < n → f (runLen f n) = false
  | _, 0, h => by omega
  | f, (n + 1), h => by
      rw [runLen] at h ⊢
      by_cases hf : f 0
      · rw [if_pos hf] at h ⊢
        exact runLen_stop (fun j => f (j + 1)) n (by omega)
      · rw [if_neg hf]
        simpa using hf

/-- Characterisation of `runLen`. -/
theorem runLen_eq (f : ℕ → Bool) (n c : ℕ) (hc : c ≤ n) (h1 : ∀ j < c, f j = true)
    (h2 : c = n ∨ f c = false) : runLen f n = c := by
  rcases lt_trichotomy (runLen f n) c with h | h | h
  · have hlt : runLen f n < n := lt_of_lt_of_le h hc
    have := runLen_stop f n hlt
    have := h1 _ h
    simp_all
  · exact h
  · have hcn : c < n := lt_of_lt_of_le h (runLen_le f n)
    have hfc : f c = false := h2.resolve_left (by omega)
    have := runLen_true f n c h
    simp_all

/-! ### The clan map of one side -/

variable {k : ℕ} {len : Fin k → ℕ}

/-- The multiplicity of position `j` of arm `i`, extended by `0` past the end. -/
def armVal (s : SpiderV k len → ℕ) (i : Fin k) (j : ℕ) : ℕ :=
  if h : j < len i then s (some ⟨i, ⟨j, h⟩⟩) else 0

theorem armVal_of_lt (s : SpiderV k len → ℕ) (i : Fin k) (j : Fin (len i)) :
    armVal s i j.val = s (some ⟨i, j⟩) := by
  rw [armVal, dif_pos j.isLt]

theorem armVal_of_ge (s : SpiderV k len → ℕ) (i : Fin k) {j : ℕ} (hj : len i ≤ j) :
    armVal s i j = 0 := by
  rw [armVal, dif_neg (by omega)]

/-- The *active prefix* of arm `i`: the length of its initial run of multiplicity-one
vertices. -/
def pref (s : SpiderV k len → ℕ) (i : Fin k) : ℕ :=
  runLen (fun j => decide (armVal s i j = 1)) (len i)

theorem pref_le (s : SpiderV k len → ℕ) (i : Fin k) : pref s i ≤ len i :=
  runLen_le _ _

theorem pref_ones (s : SpiderV k len → ℕ) (i : Fin k) {j : ℕ} (hj : j < pref s i) :
    armVal s i j = 1 := by
  have := runLen_true _ _ _ hj
  simpa using this

theorem pref_stop (s : SpiderV k len → ℕ) (i : Fin k) (h : pref s i < len i) :
    armVal s i (pref s i) ≠ 1 := by
  have := runLen_stop (fun j => decide (armVal s i j = 1)) (len i) h
  simpa using this

theorem pref_eq (s : SpiderV k len → ℕ) (i : Fin k) (c : ℕ) (hc : c ≤ len i)
    (h1 : ∀ j < c, armVal s i j = 1) (h2 : c = len i ∨ armVal s i c ≠ 1) :
    pref s i = c := by
  refine runLen_eq _ _ _ hc (fun j hj => by simpa using h1 j hj) ?_
  rcases h2 with h2 | h2
  · exact Or.inl h2
  · exact Or.inr (by simpa using h2)

/-! ### Prefix and tail weights -/

/-- The weight of an active prefix of length `p` once the hub is switched off. -/
noncomputable def prefW (p : ℕ) : LaurentPolynomial ℚ := if p = 0 then 1 else Pw (p % 2)

theorem prefW_zero : prefW 0 = 1 := by rw [prefW, if_pos rfl]

theorem prefW_odd {p : ℕ} (h : p % 2 = 1) : prefW p = zz := by
  rw [prefW, if_neg (by omega), h, Pw_one]

theorem prefW_even {p : ℕ} (h : p % 2 = 0) (hp : p ≠ 0) : prefW p = (2 : ℚ) • (1 : LaurentPolynomial ℚ) := by
  rw [prefW, if_neg hp, h, Pw_zero]
  simp only [smul_eq_C_mul, mul_one]
  norm_num [C_two_eq]

theorem isGoodW_prefW (p : ℕ) : IsGoodW (prefW p) := by
  rcases Nat.eq_zero_or_pos p with rfl | hp
  · rw [prefW_zero]; exact IsGoodW.one
  rcases Nat.mod_two_eq_zero_or_one p with h | h
  · rw [prefW_even h (by omega)]
    exact Or.inr ⟨1, 0, by simp⟩
  · rw [prefW_odd h]
    exact Or.inr ⟨0, 1, by simp⟩

/-- The weight of the part of arm `i` strictly beyond its active prefix. -/
noncomputable def tailW (s : SpiderV k len → ℕ) (i : Fin k) : LaurentPolynomial ℚ :=
  Wpoly (pathGraph (len i - pref s i))
    (fun j : Fin (len i - pref s i) => armVal s i (pref s i + j.val))

theorem isGoodW_tailW (s : SpiderV k len → ℕ) (i : Fin k) : IsGoodW (tailW s i) :=
  Wpoly_pathGraph_isGood _ _

theorem Wpoly_prefix_ones (P : ℕ) (g : Fin P → ℕ) (hg : ∀ x, g x = 1) :
    Wpoly (pathGraph P) g = prefW P := by
  have hg' : g = fun _ => 1 := funext hg
  subst hg'
  rcases Nat.eq_zero_or_pos P with rfl | h0
  · rw [Wpoly_pathGraph_zero, prefW_zero]
  · rw [Wpoly_pathGraph_ones P h0, prefW, if_neg (by omega)]

/-! ### The weight of a single arm -/

/-- **Every arm factors as prefix times tail.** -/
theorem Wpoly_arm (s : SpiderV k len → ℕ) (i : Fin k)
    (hstop : pref s i < len i → armVal s i (pref s i) = 0) :
    Wpoly (pathGraph (len i)) (fun j : Fin (len i) => armVal s i j.val)
      = prefW (pref s i) * tailW s i := by
  set P := pref s i with hP
  set N := len i with hN
  have hPN : P ≤ N := pref_le s i
  have hsplit : N = P + (N - P) := by omega
  rw [Wpoly_pathGraph_split' P (N - P) hsplit _ ?_]
  · have hall : (fun x : Fin P => armVal s i (finCongr hsplit.symm (Fin.castAdd (N - P) x)).val)
        = fun x : Fin P => armVal s i x.val := rfl
    have hall2 : (fun y : Fin (N - P) => armVal s i (finCongr hsplit.symm (Fin.natAdd P y)).val)
        = fun y : Fin (N - P) => armVal s i (P + y.val) := rfl
    rw [hall, hall2, Wpoly_prefix_ones P _ (fun x => pref_ones s i x.isLt)]
    rfl
  · intro x y hxy
    right
    have hy : y.val = 0 := by have := x.isLt; omega
    have hidx : (finCongr hsplit.symm (Fin.natAdd P y)).val = P := by
      simp [hy]
    rw [hidx]
    exact hstop (by have := y.isLt; omega)

/-! ### The alternating pattern -/

/-- The rewritten runs `2,0,2,0,…` of the published transformation carry weight `1`:
each `2` is an isolated cloned `K₂`. -/
theorem Wpoly_pathGraph_altFrom : ∀ (N t : ℕ),
    Wpoly (pathGraph N) (fun j : Fin N => altVal (t + j.val)) = 1
  | 0, t => Wpoly_pathGraph_zero _
  | (N + 1), t => by
      have hsplit : N + 1 = 1 + N := by omega
      rw [Wpoly_pathGraph_split' 1 N hsplit _ ?_]
      · have h1 : (fun x : Fin 1 => altVal (t + (finCongr hsplit.symm (Fin.castAdd N x)).val))
            = fun _ : Fin 1 => altVal t := by
          funext x
          show altVal (t + x.val) = altVal t
          have hx : x.val = 0 := by omega
          rw [hx, Nat.add_zero]
        have h2 : (fun y : Fin N => altVal (t + (finCongr hsplit.symm (Fin.natAdd 1 y)).val))
            = fun y : Fin N => altVal ((t + 1) + y.val) := by
          funext y
          show altVal (t + (1 + y.val)) = altVal (t + 1 + y.val)
          congr 1
          omega
        rw [h1, h2, pathGraph_one_eq_bot, Wpoly_bot_fin, Fin.prod_univ_one,
          Wpoly_pathGraph_altFrom N (t + 1), mul_one]
        rcases Nat.mod_two_eq_zero_or_one t with h | h
        · rw [show altVal t = 2 from by rw [altVal, if_pos h]]
          rfl
        · rw [show altVal t = 0 from by rw [altVal, if_neg (by omega)]]
          rfl
      · intro x y hxy
        have hx : x.val = 0 := by omega
        have hy : y.val = 0 := by omega
        have hi1 : (finCongr hsplit.symm (Fin.castAdd N x)).val = 0 := hx
        have hi2 : (finCongr hsplit.symm (Fin.natAdd 1 y)).val = 1 := by
          show 1 + y.val = 1
          omega
        rw [hi1, hi2]
        unfold altVal
        rcases Nat.mod_two_eq_zero_or_one t with h | h
        · right
          rw [if_neg (by omega)]
        · left
          rw [if_neg (by omega)]

theorem Wpoly_pathGraph_alt (N : ℕ) :
    Wpoly (pathGraph N) (fun j : Fin N => altVal j.val) = 1 := by
  have := Wpoly_pathGraph_altFrom N 0
  simpa using this

/-! ### The weight of one side -/

/-- `Option X` splits off its base point as a single vertex. -/
def optionSumFin1 (X : Type*) : Option X ≃ X ⊕ Fin 1 where
  toFun o := o.elim (Sum.inr 0) Sum.inl
  invFun y := Sum.elim some (fun _ => none) y
  left_inv := by rintro (_ | x) <;> rfl
  right_inv := by
    rintro (x | j)
    · rfl
    · rw [Subsingleton.elim j 0]; rfl

/-- **A spider splits at its hub** whenever every hub edge has an inactive endpoint. -/
theorem Wpoly_spider_hub (s : SpiderV k len → ℕ)
    (hcross : ∀ (i : Fin k), 0 < len i → s none = 0 ∨ armVal s i 0 = 0) :
    Wpoly (spider k len) s
      = (∏ i : Fin k, Wpoly (pathGraph (len i)) (fun j : Fin (len i) => armVal s i j.val))
        * dotW (s none) := by
  have key := Wpoly_split_equiv (spider k len) s (ArmsG k len) (⊥ : SimpleGraph (Fin 1))
    (optionSumFin1 _) ?_ ?_ ?_
  · rw [key, Wpoly_ArmsG]
    congr 1
    · refine Finset.prod_congr rfl fun i _ => ?_
      refine congrArg _ ?_
      funext j
      rw [armVal_of_lt]
      rfl
    · have h1 : (fun y : Fin 1 => s ((optionSumFin1 _).symm (Sum.inr y)))
          = fun _ : Fin 1 => s none := by
        funext y
        rfl
      rw [h1, Wpoly_bot_one]
  · intro x y
    exact Iff.rfl
  · intro x y
    constructor
    · intro h; exact h.elim
    · intro h; exact h.elim
  · rintro ⟨i, j⟩ y hadj
    have hj : j.val = 0 := hadj
    have hlen : 0 < len i := by have := j.isLt; omega
    rcases hcross i hlen with h | h
    · right; exact h
    · left
      rw [armVal, dif_pos hlen] at h
      have : (⟨0, hlen⟩ : Fin (len i)) = j := Fin.ext (by simp [hj])
      rw [this] at h
      exact h

/-- **A side whose hub is inactive**: the weight is the product of the arm weights times
the single-vertex weight of the hub. -/
theorem Wpoly_side_inactive (s : SpiderV k len → ℕ)
    (hcross : ∀ (i : Fin k), 0 < len i → s none = 0 ∨ armVal s i 0 = 0)
    (hstop : ∀ i : Fin k, pref s i < len i → armVal s i (pref s i) = 0) :
    Wpoly (spider k len) s
      = (∏ i : Fin k, (prefW (pref s i) * tailW s i)) * dotW (s none) := by
  rw [Wpoly_spider_hub s hcross]
  congr 1
  exact Finset.prod_congr rfl fun i _ => Wpoly_arm s i (hstop i)

/-! ### The active hub component -/

/-- The number of arms with odd active prefix. -/
noncomputable def pNum (s : SpiderV k len → ℕ) : ℕ :=
  (Finset.univ.filter (fun i : Fin k => pref s i % 2 = 1)).card

/-- The number of arms with even *positive* active prefix; the derived scalar of the
published block is `c = 2 ^ eNum`. -/
noncomputable def eNum (s : SpiderV k len → ℕ) : ℕ :=
  (Finset.univ.filter (fun i : Fin k => pref s i % 2 = 0 ∧ pref s i ≠ 0)).card

/-- The product of the prefix weights over any set of arms is `2^e (z+z⁻¹)^p`, where `p`
counts the arms with odd active prefix and `e` those with even positive active prefix.
This is where the derived scalar `c = 2^e ≥ 1` of the published block comes from. -/
theorem prod_prefW (s : SpiderV k len → ℕ) (S : Finset (Fin k)) : ∃ e : ℕ,
    ∏ i ∈ S, prefW (pref s i)
      = ((2 : ℚ) ^ e) • zz ^ (S.filter (fun i => pref s i % 2 = 1)).card := by
  classical
  induction S using Finset.induction_on with
  | empty => exact ⟨0, by simp⟩
  | insert a S ha ih =>
      obtain ⟨e, he⟩ := ih
      rw [Finset.prod_insert ha, he, Finset.filter_insert]
      by_cases hodd : pref s a % 2 = 1
      · refine ⟨e, ?_⟩
        rw [if_pos hodd, Finset.card_insert_of_notMem (fun hc => ha (Finset.mem_filter.mp hc).1),
          prefW_odd hodd, pow_succ]
        rw [smul_eq_C_mul, smul_eq_C_mul]
        ring
      · rw [if_neg hodd]
        by_cases hz : pref s a = 0
        · exact ⟨e, by rw [hz, prefW_zero, one_mul]⟩
        · refine ⟨e + 1, ?_⟩
          rw [prefW_even (by omega) hz]
          rw [smul_eq_C_mul, smul_eq_C_mul, smul_eq_C_mul, mul_one, pow_succ]
          rw [map_mul]
          ring

/-- The all-ones spider on the active prefixes has weight `z^|1-p| + z^(-|1-p|)`. -/
theorem Wpoly_spider_ones (k : ℕ) (q : Fin k → ℕ) :
    Wpoly (spider k q) (fun _ => 1)
      = Pw ((1 - ((Finset.univ.filter (fun i : Fin k => q i % 2 = 1)).card : ℤ)).natAbs) := by
  classical
  set p := (Finset.univ.filter (fun i : Fin k => q i % 2 = 1)).card with hp
  rcases Nat.eq_zero_or_pos p with h0 | h0
  · have himb := spider_imbalance k q
    rw [← hp, h0] at himb
    have hgf : imbalanceGF (spider k q) = Pw 1 := by
      rw [imbalanceGF_connected (spider_connected k q) (spider_isStablePair k q), Pw]
      have e1 : ((spiderEven k q).card : ℤ) - (spiderOdd k q).card = 1 := by omega
      have e2 : ((spiderOdd k q).card : ℤ) - (spiderEven k q).card = -1 := by omega
      rw [e1, e2]
      norm_num
    rw [Wpoly, clanFactorial_one, imbalanceGF_iso (clanOneIso _), hgf, h0]
    simp
  · have := Wpoly_spider k q (p - 1) (by rw [← hp]; omega)
    rw [this]
    congr 1
    omega

/-- **A side whose hub is active**: the weight is `(z^r + z^(-r)) * ∏ tails` with
`r = |1 - p|`, `p` the number of arms with odd active prefix. -/
theorem Wpoly_side_active (s : SpiderV k len → ℕ) (hhub : s none = 1)
    (hstop : ∀ i : Fin k, pref s i < len i → armVal s i (pref s i) = 0) :
    Wpoly (spider k len) s
      = Pw ((1 - (pNum s : ℤ)).natAbs) * ∏ i : Fin k, tailW s i := by
  classical
  have hcut : ∀ (i : Fin k) (h : pref s i < len i), s (some ⟨i, ⟨pref s i, h⟩⟩) = 0 := by
    intro i h
    have := hstop i h
    rwa [armVal, dif_pos h] at this
  rw [Wpoly_spider_cut k len (pref s) (pref_le s) s hcut]
  congr 1
  · have hones : (fun x => s ((spiderCut k len (pref s) (pref_le s)).symm (Sum.inl x)))
        = fun _ : SpiderV k (pref s) => 1 := by
      funext x
      match x with
      | none => exact hhub
      | some ⟨i, j⟩ =>
          rw [spiderCut_symm_inl_some]
          have := pref_ones s i j.isLt
          rwa [armVal, dif_pos (show j.val < len i from by
            have := j.isLt; have := pref_le s i; omega)] at this
    rw [hones, Wpoly_spider_ones]
    rfl
  · rw [Wpoly_ArmsG]
    refine Finset.prod_congr rfl fun i _ => ?_
    rw [tailW]
    refine congrArg _ ?_
    funext j
    rw [spiderCut_symm_inr, armVal,
      dif_pos (show pref s i + j.val < len i from by have := j.isLt; omega)]

end ClanAudit
