import RequestProject.ClanAudit.Split

/-!
# Clan weights of paths

Away from the two hubs every component of a clan graph of the adjacent two-hub tree lives
inside a pendant path.  This file computes the normalized weight of the clan graph of a
path `pathGraph N` with an arbitrary multiplicity function:

* `Wpoly_bot_fin` — a graph with no edges: the weight is the product of the single-vertex
  weights `dotW` (`1` for multiplicity `0`, `z + z⁻¹` for `1`, `1` for `2` — the cloned
  `K₂` cancellation — and `0` from multiplicity `≥ 3`, which creates a triangle);
* `Wpoly_pathGraph_ones` — a path all of whose multiplicities are one contributes
  `z^(N mod 2) + z^(-(N mod 2))`;
* `Wpoly_pathGraph_split` — a path may be cut at a vertex of multiplicity zero;
* `Wpoly_pathGraph_isGood` — the weight of *any* path clan is either `0` or of the form
  `2^j (z+z⁻¹)^i`; in particular it is centrally unimodal.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph

/-! ### Good weights -/

/-- The shape of every clan weight occurring away from the hubs: either zero (a triangle
or an odd clique was present) or `2^j (z + z⁻¹)^i`. -/
def IsGoodW (f : LaurentPolynomial ℚ) : Prop :=
  f = 0 ∨ ∃ j i : ℕ, f = ((2 : ℚ) ^ j) • zz ^ i

theorem isCU_zero : IsCU (0 : LaurentPolynomial ℚ) where
  symm m := by simp [coeffL]
  nonneg m := by simp [coeffL]
  decr m _ := by simp [coeffL]

theorem IsGoodW.isCU {f : LaurentPolynomial ℚ} (h : IsGoodW f) : IsCU f := by
  rcases h with rfl | ⟨j, i, rfl⟩
  · exact isCU_zero
  · exact (isCU_zz_pow i).smul (by positivity)

theorem IsGoodW.mul {f g : LaurentPolynomial ℚ} (hf : IsGoodW f) (hg : IsGoodW g) :
    IsGoodW (f * g) := by
  rcases hf with rfl | ⟨j, i, rfl⟩
  · exact Or.inl (by ring)
  rcases hg with rfl | ⟨j', i', rfl⟩
  · exact Or.inl (by ring)
  refine Or.inr ⟨j + j', i + i', ?_⟩
  rw [smul_eq_C_mul, smul_eq_C_mul, smul_eq_C_mul]
  rw [map_pow, map_pow, map_pow, pow_add, pow_add]
  ring

theorem IsGoodW.one : IsGoodW (1 : LaurentPolynomial ℚ) :=
  Or.inr ⟨0, 0, by simp⟩

theorem IsGoodW.zero : IsGoodW (0 : LaurentPolynomial ℚ) := Or.inl rfl

theorem Pw_zero : Pw 0 = 2 := by
  rw [Pw]
  norm_num

theorem Pw_one : Pw 1 = zz := by rw [zz_eq, Pw]; norm_num

theorem IsGoodW.Pw_mod_two (t : ℕ) : IsGoodW (Pw (t % 2)) := by
  rcases Nat.mod_two_eq_zero_or_one t with h | h
  · rw [h, Pw_zero]
    exact Or.inr ⟨1, 0, by simp; rw [smul_eq_C_mul]; norm_num [C_two_eq]⟩
  · rw [h, Pw_one]
    exact Or.inr ⟨0, 1, by simp⟩

/-! ### Edgeless clan graphs -/

/-- The normalized weight of a single vertex of multiplicity `k`. -/
noncomputable def dotW : ℕ → LaurentPolynomial ℚ
  | 0 => 1
  | 1 => zz
  | 2 => 1
  | _ => 0

theorem isGoodW_dotW (k : ℕ) : IsGoodW (dotW k) := by
  match k with
  | 0 => exact IsGoodW.one
  | 1 => exact Or.inr ⟨0, 1, by simp [dotW]⟩
  | 2 => exact IsGoodW.one
  | (n + 3) => exact IsGoodW.zero

theorem Wpoly_bot_one (k : ℕ) : Wpoly (⊥ : SimpleGraph (Fin 1)) (fun _ => k) = dotW k := by
  classical
  match k with
  | 0 =>
      haveI : IsEmpty (Σ _ : Fin 1, Fin 0) := ⟨fun x => x.2.elim0⟩
      rw [Wpoly, imbalanceGF_isEmpty, clanFactorial]
      simp [dotW]
  | 1 =>
      have hconn : (⊥ : SimpleGraph (Fin 1)).Connected := by
        rw [SimpleGraph.connected_iff]
        exact ⟨fun x y => by rw [Subsingleton.elim x y], ⟨0⟩⟩
      have hst : IsStablePair (⊥ : SimpleGraph (Fin 1)) {0} ∅ := by
        refine ⟨fun x _ y _ h => h.elim, fun x hx => absurd hx (by simp), by simp, ?_⟩
        ext x
        simp [Subsingleton.elim x 0]
      rw [Wpoly, clanFactorial_one, imbalanceGF_iso (clanOneIso _),
        imbalanceGF_connected hconn hst]
      simp [dotW, zz_eq]
  | 2 => rw [Wpoly_isolated_two]; rfl
  | (n + 3) =>
      have h3 : imbalanceGF (clan (⊥ : SimpleGraph (Fin 1)) (fun _ => n + 3)) = 0 := by
        refine imbalanceGF_eq_zero_of_triangle
          (x := ⟨0, ⟨0, by omega⟩⟩) (y := ⟨0, ⟨1, by omega⟩⟩) (z := ⟨0, ⟨2, by omega⟩⟩)
          (Or.inl ⟨rfl, ?_⟩) (Or.inl ⟨rfl, ?_⟩) (Or.inl ⟨rfl, ?_⟩) <;>
        · simp only [ne_eq, Sigma.mk.injEq, heq_eq_eq, true_and]
          intro hcon
          exact absurd (congrArg Fin.val hcon) (by simp)
      rw [Wpoly, h3]
      simp [dotW]

/-- **A clan graph with no edges.**  Its normalized weight is the product of the
single-vertex weights. -/
theorem Wpoly_bot_fin : ∀ (N : ℕ) (alpha : Fin N → ℕ),
    Wpoly (⊥ : SimpleGraph (Fin N)) alpha = ∏ i, dotW (alpha i)
  | 0, alpha => by
      haveI : IsEmpty (Fin 0) := inferInstance
      rw [Wpoly_isEmpty]
      simp
  | (N + 1), alpha => by
      have hsplit := Wpoly_split_equiv (⊥ : SimpleGraph (Fin (N + 1))) alpha
        (⊥ : SimpleGraph (Fin N)) (⊥ : SimpleGraph (Fin 1))
        (finSumFinEquiv.symm) (by simp) (by simp) (by simp)
      have hl : ∀ x : Fin N, (finSumFinEquiv.symm.symm (Sum.inl x)) = Fin.castSucc x := by
        intro x; simp [Fin.castSucc, Fin.castAdd]
      have hr : (finSumFinEquiv.symm.symm (Sum.inr (0 : Fin 1))) = Fin.last N := by
        simp [Fin.last, Fin.natAdd]
      rw [hsplit]
      have h1 : Wpoly (⊥ : SimpleGraph (Fin N))
          (fun x => alpha (finSumFinEquiv.symm.symm (Sum.inl x)))
          = ∏ i : Fin N, dotW (alpha (Fin.castSucc i)) := by
        rw [Wpoly_bot_fin N _]
        exact Finset.prod_congr rfl fun i _ => by rw [hl i]
      have h2 : Wpoly (⊥ : SimpleGraph (Fin 1))
          (fun y => alpha (finSumFinEquiv.symm.symm (Sum.inr y)))
          = dotW (alpha (Fin.last N)) := by
        have : (fun y : Fin 1 => alpha (finSumFinEquiv.symm.symm (Sum.inr y)))
            = (fun _ : Fin 1 => alpha (Fin.last N)) := by
          funext y
          rw [Subsingleton.elim y 0, hr]
        rw [this, Wpoly_bot_one]
      rw [h1, h2, Fin.prod_univ_castSucc]

/-! ### Paths -/

theorem pathGraph_one_eq_bot : (pathGraph 1) = (⊥ : SimpleGraph (Fin 1)) := by
  ext x y
  simp only [pathGraph_adj, SimpleGraph.bot_adj]
  constructor
  · intro h
    have := x.isLt
    have := y.isLt
    omega
  · exact False.elim

theorem Wpoly_pathGraph_zero (alpha : Fin 0 → ℕ) : Wpoly (pathGraph 0) alpha = 1 :=
  Wpoly_isEmpty _ _

private theorem card_filter_par_range (m b : ℕ) (hb : b < 2) :
    ((Finset.range m).filter (fun j => j % 2 = b)).card = (m + 1 - b) / 2 := by
  induction m with
  | zero => simp; omega
  | succ m ih =>
      rw [Finset.range_add_one, Finset.filter_insert]
      by_cases h : m % 2 = b
      · rw [if_pos h, Finset.card_insert_of_notMem (by simp), ih]
        omega
      · rw [if_neg h, ih]
        omega

private theorem card_filter_par_fin (m b : ℕ) (hb : b < 2) :
    ((Finset.univ : Finset (Fin m)).filter (fun j => j.val % 2 = b)).card = (m + 1 - b) / 2 := by
  classical
  rw [← card_filter_par_range m b hb, Finset.card_filter, Finset.card_filter,
    Fin.sum_univ_eq_sum_range (fun j => if j % 2 = b then 1 else 0)]

/-- **A path all of whose multiplicities are one.**  Its clan graph is the path itself,
a connected bipartite graph whose colour classes differ by `N mod 2`. -/
theorem Wpoly_pathGraph_ones (N : ℕ) (hN : 0 < N) :
    Wpoly (pathGraph N) (fun _ => 1) = Pw (N % 2) := by
  classical
  obtain ⟨M, rfl⟩ : ∃ M, N = M + 1 := ⟨N - 1, by omega⟩
  set A : Finset (Fin (M + 1)) := Finset.univ.filter (fun i => i.val % 2 = 0) with hA
  set B : Finset (Fin (M + 1)) := Finset.univ.filter (fun i => i.val % 2 = 1) with hB
  have hst : IsStablePair (pathGraph (M + 1)) A B := by
    refine ⟨?_, ?_, ?_, ?_⟩
    · intro x hx y hy hadj
      rw [hA, Finset.mem_filter] at hx hy
      rw [pathGraph_adj] at hadj
      omega
    · intro x hx y hy hadj
      rw [hB, Finset.mem_filter] at hx hy
      rw [pathGraph_adj] at hadj
      omega
    · rw [Finset.disjoint_left]
      intro x hx hx'
      rw [hA, Finset.mem_filter] at hx
      rw [hB, Finset.mem_filter] at hx'
      omega
    · ext x
      simp only [hA, hB, Finset.mem_union, Finset.mem_filter, Finset.mem_univ, true_and,
        iff_true]
      omega
  have hcA : A.card = (M + 2) / 2 := by
    rw [hA, card_filter_par_fin (M + 1) 0 (by omega)]
    omega
  have hcB : B.card = (M + 1) / 2 := by
    rw [hB, card_filter_par_fin (M + 1) 1 (by omega)]
    omega
  have himb : imbalanceGF (pathGraph (M + 1)) = Pw ((M + 1) % 2) := by
    rw [imbalanceGF_connected (pathGraph_connected M) hst, hcA, hcB,
      show ((((M + 2) / 2 : ℕ)) : ℤ) - ((((M + 1) / 2 : ℕ)) : ℤ)
          = ((((M + 1) % 2 : ℕ)) : ℤ) from by omega,
      show ((((M + 1) / 2 : ℕ)) : ℤ) - ((((M + 2) / 2 : ℕ)) : ℤ)
          = -((((M + 1) % 2 : ℕ)) : ℤ) from by omega]
    rfl
  rw [Wpoly, clanFactorial_one, imbalanceGF_iso (clanOneIso _), himb]
  simp

/-- **Cutting a path.**  If the vertex at position `r` (or the one before it) has
multiplicity zero, the path clan splits. -/
theorem Wpoly_pathGraph_split (r s : ℕ) (alpha : Fin (r + s) → ℕ)
    (hcut : ∀ (x : Fin r) (y : Fin s), x.val + 1 = r + y.val →
      alpha (Fin.castAdd s x) = 0 ∨ alpha (Fin.natAdd r y) = 0) :
    Wpoly (pathGraph (r + s)) alpha
      = Wpoly (pathGraph r) (fun x => alpha (Fin.castAdd s x))
        * Wpoly (pathGraph s) (fun y => alpha (Fin.natAdd r y)) := by
  have key := Wpoly_split_equiv (pathGraph (r + s)) alpha (pathGraph r) (pathGraph s)
    finSumFinEquiv.symm ?_ ?_ ?_
  · simpa using key
  · intro x y
    simp only [Equiv.symm_symm, finSumFinEquiv_apply_left, pathGraph_adj, Fin.val_castAdd]
  · intro x y
    simp only [Equiv.symm_symm, finSumFinEquiv_apply_right, pathGraph_adj, Fin.val_natAdd]
    omega
  · intro x y hadj
    simp only [Equiv.symm_symm, finSumFinEquiv_apply_left, finSumFinEquiv_apply_right,
      pathGraph_adj, Fin.val_castAdd, Fin.val_natAdd] at hadj
    have hx := x.isLt
    have hy := y.isLt
    simp only [Equiv.symm_symm, finSumFinEquiv_apply_left, finSumFinEquiv_apply_right]
    exact hcut x y (by omega)

theorem Wpoly_pathGraph_congr {N N' : ℕ} (h : N = N') (alpha : Fin N → ℕ) :
    Wpoly (pathGraph N') (fun y => alpha (finCongr h.symm y)) = Wpoly (pathGraph N) alpha := by
  subst h
  simp

/-- The cutting lemma for a length that is only propositionally a sum. -/
theorem Wpoly_pathGraph_split' {N : ℕ} (r s : ℕ) (hN : N = r + s) (alpha : Fin N → ℕ)
    (hcut : ∀ (x : Fin r) (y : Fin s), x.val + 1 = r + y.val →
      alpha (finCongr hN.symm (Fin.castAdd s x)) = 0 ∨
      alpha (finCongr hN.symm (Fin.natAdd r y)) = 0) :
    Wpoly (pathGraph N) alpha
      = Wpoly (pathGraph r) (fun x => alpha (finCongr hN.symm (Fin.castAdd s x)))
        * Wpoly (pathGraph s) (fun y => alpha (finCongr hN.symm (Fin.natAdd r y))) := by
  rw [← Wpoly_pathGraph_congr hN alpha]
  exact Wpoly_pathGraph_split r s _ hcut

/-! ### An arbitrary path clan -/

theorem Wpoly_eq_zero_of_pathGraph_triangle {N : ℕ} (alpha : Fin N → ℕ) {u v : Fin N}
    (huv : (pathGraph N).Adj u v) (hu : 2 ≤ alpha u) (hv : 1 ≤ alpha v) :
    Wpoly (pathGraph N) alpha = 0 := by
  rw [Wpoly, clan_imbalanceGF_eq_zero_of_mult_two huv hu hv, smul_zero]

/-- **Every path clan weight is `0` or `2^j (z+z⁻¹)^i`.** -/
theorem Wpoly_pathGraph_isGood (N : ℕ) :
    ∀ alpha : Fin N → ℕ, IsGoodW (Wpoly (pathGraph N) alpha) := by
  induction N using Nat.strong_induction_on with
  | _ N ih =>
    intro alpha
    classical
    rcases Nat.eq_zero_or_pos N with rfl | hN
    · rw [Wpoly_pathGraph_zero]
      exact IsGoodW.one
    have hex : ∃ j : ℕ, N ≤ j ∨ ∃ h : j < N, alpha ⟨j, h⟩ ≠ 1 := ⟨N, Or.inl le_rfl⟩
    set rho := Nat.find hex with hrho
    have hspec := Nat.find_spec hex
    have hmin : ∀ j < rho, ¬ (N ≤ j ∨ ∃ h : j < N, alpha ⟨j, h⟩ ≠ 1) := fun j hj =>
      Nat.find_min hex hj
    have hones : ∀ (j : ℕ) (h : j < N), j < rho → alpha ⟨j, h⟩ = 1 := by
      intro j h hj
      by_contra hne
      exact hmin j hj (Or.inr ⟨h, hne⟩)
    have hrhoN : rho ≤ N := by
      by_contra hcon
      exact hmin N (by omega) (Or.inl le_rfl)
    by_cases hfull : rho = N
    · have hall : alpha = fun _ => 1 := by
        funext j
        exact hones j.val j.isLt (by omega)
      rw [hall, Wpoly_pathGraph_ones N hN]
      exact IsGoodW.Pw_mod_two N
    · have hlt : rho < N := by omega
      have hne1 : alpha ⟨rho, hlt⟩ ≠ 1 := by
        rcases hspec with h | ⟨h, h'⟩
        · omega
        · exact h'
      rcases Nat.eq_zero_or_pos rho with hz | hpos
      · -- the first vertex already has multiplicity `≠ 1`
        have hfin : (⟨rho, hlt⟩ : Fin N) = ⟨0, hN⟩ := Fin.ext (by omega)
        rw [hfin] at hne1
        by_cases hbig : 2 ≤ alpha ⟨0, hN⟩
        · by_cases hnext : ∃ h : 1 < N, 1 ≤ alpha ⟨1, h⟩
          · obtain ⟨h1, h1'⟩ := hnext
            refine Or.inl (Wpoly_eq_zero_of_pathGraph_triangle alpha
              (u := ⟨0, hN⟩) (v := ⟨1, h1⟩) ?_ hbig h1')
            rw [pathGraph_adj]
            exact Or.inl rfl
          · have hN1 : N = 1 + (N - 1) := by omega
            rw [Wpoly_pathGraph_split' 1 (N - 1) hN1 alpha ?_]
            · refine IsGoodW.mul ?_ (ih (N - 1) (by omega) _)
              rw [pathGraph_one_eq_bot, Wpoly_bot_fin, Fin.prod_univ_one]
              exact isGoodW_dotW _
            · intro x y hxy
              have hy0 : y.val = 0 := by have := x.isLt; omega
              have hidx : finCongr hN1.symm (Fin.natAdd 1 y) = ⟨1, by omega⟩ :=
                Fin.ext (by simp [hy0])
              right
              rw [hidx]
              by_contra hcon
              exact hnext ⟨by omega, by omega⟩
        · have hzero : alpha ⟨0, hN⟩ = 0 := by omega
          have hN1 : N = 1 + (N - 1) := by omega
          rw [Wpoly_pathGraph_split' 1 (N - 1) hN1 alpha ?_]
          · refine IsGoodW.mul ?_ (ih (N - 1) (by omega) _)
            rw [pathGraph_one_eq_bot, Wpoly_bot_fin, Fin.prod_univ_one]
            exact isGoodW_dotW _
          · intro x y _
            have hidx : finCongr hN1.symm (Fin.castAdd (N - 1) x) = ⟨0, hN⟩ :=
              Fin.ext (by simp)
            left
            rw [hidx]
            exact hzero
      · -- a genuine run of multiplicity-one vertices
        by_cases hbig : 2 ≤ alpha ⟨rho, hlt⟩
        · refine Or.inl (Wpoly_eq_zero_of_pathGraph_triangle alpha
            (u := ⟨rho, hlt⟩) (v := ⟨rho - 1, by omega⟩) ?_ hbig ?_)
          · rw [pathGraph_adj]
            exact Or.inr (by simp; omega)
          · rw [hones (rho - 1) (by omega) (by omega)]
        · have hzero : alpha ⟨rho, hlt⟩ = 0 := by omega
          have hN1 : N = rho + (N - rho) := by omega
          rw [Wpoly_pathGraph_split' rho (N - rho) hN1 alpha ?_]
          · refine IsGoodW.mul ?_ (ih (N - rho) (by omega) _)
            have hall : (fun x : Fin rho => alpha (finCongr hN1.symm (Fin.castAdd (N - rho) x)))
                = fun _ => 1 := by
              funext x
              have hidx : finCongr hN1.symm (Fin.castAdd (N - rho) x)
                  = ⟨x.val, by have := x.isLt; omega⟩ := Fin.ext (by simp)
              rw [hidx]
              exact hones x.val (by have := x.isLt; omega) x.isLt
            rw [hall, Wpoly_pathGraph_ones rho hpos]
            exact IsGoodW.Pw_mod_two rho
          · intro x y hxy
            have hy0 : y.val = 0 := by have := x.isLt; omega
            have hidx : finCongr hN1.symm (Fin.natAdd rho y) = ⟨rho, hlt⟩ :=
              Fin.ext (by simp [hy0])
            right
            rw [hidx]
            exact hzero

end ClanAudit
