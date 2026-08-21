import RequestProject.ClanAudit.P2Weight

/-!
# The normalized weight of a `p = 2` local block with arbitrary further even arms

This file closes the arbitrary-arm case of the `p = 2` weight theorem.  The hub carries

* two distinguished active arms whose positive prefixes have *odd* lengths `L ≤ M`
  (as in `RequestProject/ClanAudit/P2Weight.lean`), and
* an arbitrary finite indexed family of `e` further active arms, the `i`-th having a
  positive *even* prefix length `len i`.

Both clan states paired by `localMapP2` are modelled by their clan graphs:

* the **active** state: the hub component is the spider with `2 + e` arms of lengths
  `L, M, len 0, …, len (e-1)`, all multiplicities one.  Exactly two of these lengths are
  odd, so by `Wpoly_spider` its normalized weight is `z^r + z^(-r)` with
  `r = p - 1 = 1`;
* the **image** state: the hub is deactivated and the prefixes of the two odd arms are
  rewritten as `2,0,2,0,…`.  This produces `L` cloned `K₂` components (whose orientation
  counts cancel the new factors `2!`, `Wpoly_bot_two`), the untouched remainder of the
  arm `B` — a path with `M - L + 1` vertices, contributing `z + z⁻¹` — and, for each of
  the `e` even arms, a *separate* path component with `len i` vertices.  Each of those is
  balanced (it has an even number of vertices), so it contributes the factor `2`.

The disjoint union of the `e` even-arm components is built here as `armsGraph e len`
(vertex type `ArmsV e len`), and the reusable finite-family product theorem
`Wpoly_armsGraph` gives its normalized weight `2 ^ e`.  Consequently

```text
W(image state) = 2^e * (z + z⁻¹),
W(active state) + W(image state) = A(1, 2^e; z),
```

so the scalar of the block is *derived* to be `c = 2^e`; in particular `1 ≤ c`, which is
exactly the hypothesis the adjacent two-hub conclusion needs (`Fblock_decr`).
-/

namespace ClanAudit

open Finset LaurentPolynomial

/-! ### The empty clan graph -/

/-- A graph with no vertices has imbalance generating function `1`: its only stable
two-block partition is `(∅, ∅)`. -/
theorem imbalanceGF_isEmpty {W : Type*} [Fintype W] [DecidableEq W] [IsEmpty W]
    (H : SimpleGraph W) : imbalanceGF H = 1 := by
  classical
  have hpairs : stablePairs H = {(∅, ∅)} := by
    ext p
    simp only [Finset.mem_singleton, mem_stablePairs]
    constructor
    · intro _
      exact Prod.ext (Finset.eq_empty_of_isEmpty p.1) (Finset.eq_empty_of_isEmpty p.2)
    · rintro rfl
      exact ⟨by simp [IsIndep], by simp [IsIndep], by simp, by simp⟩
  rw [imbalanceGF, hpairs]
  simp

/-- The normalized weight of a clan graph over an empty vertex set is `1`. -/
theorem Wpoly_isEmpty {V : Type*} [Fintype V] [DecidableEq V] [IsEmpty V] (G : SimpleGraph V)
    (alpha : V → ℕ) : Wpoly G alpha = 1 := by
  have : IsEmpty (Σ v : V, Fin (alpha v)) := ⟨fun x => IsEmpty.elim ‹IsEmpty V› x.1⟩
  rw [Wpoly, imbalanceGF_isEmpty, clanFactorial]
  simp

/-! ### A finite family of even paths -/

/-- The vertex type of the disjoint union of `e` paths, the `i`-th having `len i`
vertices (modelled, as elsewhere in this development, as the spider with a single arm of
length `len i - 1`). -/
def ArmsV : (e : ℕ) → (Fin e → ℕ) → Type
  | 0, _ => Empty
  | (e + 1), len => SpiderV 1 (fun _ => len 0 - 1) ⊕ ArmsV e (fun i => len i.succ)

instance armsVFintype : ∀ (e : ℕ) (len : Fin e → ℕ), Fintype (ArmsV e len)
  | 0, _ => inferInstanceAs (Fintype Empty)
  | (e + 1), len =>
      letI := armsVFintype e (fun i => len i.succ)
      inferInstanceAs (Fintype (SpiderV 1 (fun _ => len 0 - 1) ⊕ ArmsV e (fun i => len i.succ)))

instance armsVDecEq : ∀ (e : ℕ) (len : Fin e → ℕ), DecidableEq (ArmsV e len)
  | 0, _ => inferInstanceAs (DecidableEq Empty)
  | (e + 1), len =>
      letI := armsVDecEq e (fun i => len i.succ)
      inferInstanceAs
        (DecidableEq (SpiderV 1 (fun _ => len 0 - 1) ⊕ ArmsV e (fun i => len i.succ)))

instance armsV0IsEmpty (len : Fin 0 → ℕ) : IsEmpty (ArmsV 0 len) :=
  inferInstanceAs (IsEmpty Empty)

/-- The disjoint union of `e` paths, the `i`-th having `len i` vertices. -/
def armsGraph : (e : ℕ) → (len : Fin e → ℕ) → SimpleGraph (ArmsV e len)
  | 0, _ => (⊥ : SimpleGraph Empty)
  | (e + 1), len => (spider 1 (fun _ => len 0 - 1)).sum (armsGraph e (fun i => len i.succ))

/-- A path with an even, positive number of vertices is a balanced bipartite component:
its normalized weight is `2`. -/
theorem Wpoly_path_even (k : ℕ) (hk : k % 2 = 0) (hpos : 0 < k) :
    Wpoly (spider 1 (fun _ => k - 1)) (fun _ => 1) = 2 := by
  classical
  have hcard : (Finset.univ.filter (fun i : Fin 1 => (k - 1) % 2 = 1)).card = 0 + 1 := by
    rw [Finset.filter_true_of_mem (fun i _ => by omega)]
    simp
  rw [Wpoly_spider 1 (fun _ => k - 1) 0 hcard, Pw]
  norm_num

/-- **Finite-family product theorem.**  The disjoint union of `e` paths with an even,
positive number of vertices each has normalized weight `2 ^ e`. -/
theorem Wpoly_armsGraph : ∀ (e : ℕ) (len : Fin e → ℕ), (∀ i, len i % 2 = 0) →
    (∀ i, 0 < len i) → Wpoly (armsGraph e len) (fun _ => 1) = 2 ^ e
  | 0, len, _, _ => by
      rw [pow_zero]
      exact Wpoly_isEmpty _ _
  | (e + 1), len, heven, hpos => by
      have hone : (fun _ : ArmsV (e + 1) len => (1 : ℕ))
          = Sum.elim (fun _ : SpiderV 1 (fun _ => len 0 - 1) => (1 : ℕ))
              (fun _ : ArmsV e (fun i => len i.succ) => (1 : ℕ)) := by
        funext x
        cases x <;> rfl
      have hrec := Wpoly_armsGraph e (fun i => len i.succ) (fun i => heven i.succ)
        (fun i => hpos i.succ)
      rw [hone]
      show Wpoly ((spider 1 (fun _ => len 0 - 1)).sum (armsGraph e (fun i => len i.succ)))
        (Sum.elim (fun _ => 1) (fun _ => 1)) = _
      rw [Wpoly_sum_elim, Wpoly_path_even (len 0) (heven 0) (hpos 0), hrec]
      ring

/-! ### The two normalized weights -/

/-- The hub of the active state has exactly two arms of odd prefix length. -/
theorem card_filter_odd_append {e : ℕ} {L M : ℕ} {len : Fin e → ℕ} (hL : L % 2 = 1)
    (hM : M % 2 = 1) (heven : ∀ i, len i % 2 = 0) :
    (Finset.univ.filter
        (fun i : Fin (2 + e) => Fin.append ![L, M] len i % 2 = 1)).card = 1 + 1 := by
  classical
  rw [Finset.card_filter, Fin.sum_univ_add]
  have h1 : ∑ i : Fin 2, (if Fin.append ![L, M] len (Fin.castAdd e i) % 2 = 1 then 1 else 0)
      = 1 + 1 := by
    rw [Fin.sum_univ_two]
    simp [Fin.append_left, hL, hM]
  have h2 : ∑ i : Fin e, (if Fin.append ![L, M] len (Fin.natAdd 2 i) % 2 = 1 then 1 else 0)
      = 0 := by
    refine Finset.sum_eq_zero fun i _ => ?_
    rw [Fin.append_right, if_neg (by have := heven i; omega)]
  rw [h1, h2, add_zero]

/-- **The active state.**  Its hub component is the spider with the two odd arms and the
`e` even arms; exactly two arms have odd prefix length, so its normalized weight is
`z + z⁻¹`, i.e. the term `z^r + z^(-r)` of `A(r,c;z)` with `r = p - 1 = 1`. -/
theorem Wpoly_active_even_arms {e L M : ℕ} {len : Fin e → ℕ} (hL : L % 2 = 1) (hM : M % 2 = 1)
    (heven : ∀ i, len i % 2 = 0) :
    Wpoly (spider (2 + e) (Fin.append ![L, M] len)) (fun _ => 1) = Pw 1 :=
  Wpoly_spider (2 + e) _ 1 (card_filter_odd_append hL hM heven)

/-- **The image state.**  The `L` cloned `K₂`s contribute `1`, the remainder of the arm
`B` contributes `z + z⁻¹`, and each of the `e` even arms contributes the balanced factor
`2`: the normalized weight is `2^e * (z + z⁻¹)`. -/
theorem Wpoly_image_even_arms {e L M : ℕ} {len : Fin e → ℕ} (hL : L % 2 = 1) (hM : M % 2 = 1)
    (hLM : L ≤ M) (heven : ∀ i, len i % 2 = 0) (hpos : ∀ i, 0 < len i) :
    Wpoly (((⊥ : SimpleGraph (Fin L)).sum (spider 1 (fun _ => M - L))).sum (armsGraph e len))
        (Sum.elim (Sum.elim (fun _ => 2) (fun _ => 1)) (fun _ => 1))
      = 2 ^ e * Pw 1 := by
  rw [Wpoly_sum_elim, Wpoly_image_two_arms hL hM hLM, Wpoly_armsGraph e len heven hpos,
    mul_comm]

/-! ### The block, with the derived scalar `c = 2^e` -/

theorem C_two_eq : (C (2 : ℚ) : LaurentPolynomial ℚ) = 2 :=
  LaurentPolynomial.ext_iff.mpr (congrFun rfl)

theorem two_pow_smul_zz (e : ℕ) : ((2 : ℚ) ^ e) • zz = (2 : LaurentPolynomial ℚ) ^ e * Pw 1 := by
  rw [smul_eq_C_mul, map_pow, C_two_eq, zz_eq, Pw]
  norm_num

/-- **The `p = 2` block with arbitrary further even arms has exactly the claimed
normalized weight**, with the scalar derived to be `c = 2^e`:

```text
W(active) + W(image) = A(1, 2^e; z).
```

The hub carries two active arms with odd positive prefixes of lengths `L ≤ M` together
with an arbitrary finite family of `e` further active arms with positive even prefix
lengths `len i`.  This is the general form of the request's
`localMapP2_has_claimed_normalized_two_row_weight` at `p = 2`; the scalar `c` is derived
from the normalized outside components, not postulated. -/
theorem localMapP2_normalized_weight_even_arms {e L M : ℕ} {len : Fin e → ℕ} (hL : L % 2 = 1)
    (hM : M % 2 = 1) (hLM : L ≤ M) (heven : ∀ i, len i % 2 = 0) (hpos : ∀ i, 0 < len i) :
    Wpoly (spider (2 + e) (Fin.append ![L, M] len)) (fun _ => 1)
      + Wpoly (((⊥ : SimpleGraph (Fin L)).sum (spider 1 (fun _ => M - L))).sum (armsGraph e len))
          (Sum.elim (Sum.elim (fun _ => 2) (fun _ => 1)) (fun _ => 1))
      = Ablock 1 ((2 : ℚ) ^ e) := by
  rw [Wpoly_active_even_arms hL hM heven, Wpoly_image_even_arms hL hM hLM heven hpos, Ablock,
    pow_one, two_pow_smul_zz]
  ring

/-- **The derived scalar is at least one.**  The scalar `c = 2^e` derived above satisfies
`1 ≤ c`, which is exactly the hypothesis under which the adjacent two-hub conclusion is
true (`Fblock_decr`; recall that for `c < 1` it is false, `Fblock_not_centrally_unimodal`). -/
theorem one_le_derived_scalar (e : ℕ) : (1 : ℚ) ≤ (2 : ℚ) ^ e :=
  one_le_pow₀ (by norm_num)

end ClanAudit
