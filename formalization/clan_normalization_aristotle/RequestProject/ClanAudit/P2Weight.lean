import RequestProject.ClanAudit.Spider

/-!
# The normalized weight of a `p = 2` local block (two active arms)

This file assembles the normalization kernel into a *derivation* of the proposed
one-hub block

```text
A(r,c;z) = c (z + z⁻¹)^r + z^r + z^(-r)
```

for the smallest published bad configuration outside the source's range, namely `p = 2`
(so `r = p - 1 = 1`), in the case where the hub carries exactly two active arms, whose
positive prefixes have odd lengths `L ≤ M`.

The two clan maps paired by the transformation `localMapP2` are modelled by their clan
graphs:

* the *active* state: the hub component is the spider `spider 2 ![L, M]` with all
  multiplicities one (`Wpoly = z + z⁻¹`, by `Wpoly_spider`, since `r = p - 1 = 1`);
* its *image*: the hub is deleted and the prefixes of the two arms are rewritten as
  `2,0,2,0,…`, which produces exactly `L` cloned `K₂` components (`(L+1)/2` from the arm
  `A` and `(L-1)/2` from the arm `B`) together with the untouched remainder of `B`, a
  path with `M - L + 1` vertices, modelled by `spider 1 (fun _ => M - L)` with all
  multiplicities one.

The conclusion `localMapP2_normalized_weight_two_arms` is that the *normalized* weight of
the two-state block is exactly `A(1, 1; z)`.  In particular the scalar `c` is here
*derived*, and equals `1`; it is not postulated.  Each cloned `K₂` contributes the
orientation count `2` and the factorial `2!`, which cancel (`Wpoly_bot_two`); this is the
cancellation the request asks to be proved rather than assumed.

Scope (partial, as the request requires such results to be labelled): only the two-arm
case is treated.  With `e` further arms of *even* prefix length the same computation
gives `c = 2^e`, but the `e`-fold disjoint union is not formalized here.
-/

namespace ClanAudit

open Finset LaurentPolynomial

/-! ### Invariance of `W` under isomorphisms of the base graph -/

variable {V₁ V₂ : Type*} [Fintype V₁] [DecidableEq V₁] [Fintype V₂] [DecidableEq V₂]

/-- An isomorphism of base graphs induces an isomorphism of clan graphs. -/
def clanIsoOfIso {G₁ : SimpleGraph V₁} {G₂ : SimpleGraph V₂} (e : G₁ ≃g G₂) (a : V₂ → ℕ) :
    clan G₁ (fun v => a (e v)) ≃g clan G₂ a where
  toEquiv := Equiv.sigmaCongrLeft (β := fun w => Fin (a w)) e.toEquiv
  map_rel_iff' := by
    rintro ⟨u, i⟩ ⟨v, j⟩
    show (clan G₂ a).Adj ⟨e u, i⟩ ⟨e v, j⟩ ↔ (clan G₁ (fun v => a (e v))).Adj ⟨u, i⟩ ⟨v, j⟩
    simp only [clan, ne_eq, Sigma.mk.injEq, EmbeddingLike.apply_eq_iff_eq, e.map_rel_iff]

theorem Wpoly_of_iso {G₁ : SimpleGraph V₁} {G₂ : SimpleGraph V₂} (e : G₁ ≃g G₂) (a : V₂ → ℕ) :
    Wpoly G₁ (fun v => a (e v)) = Wpoly G₂ a := by
  have hf : clanFactorial (fun v => a (e v)) = clanFactorial a :=
    Fintype.prod_equiv e.toEquiv _ _ (fun v => rfl)
  rw [Wpoly, Wpoly, hf, imbalanceGF_iso (clanIsoOfIso e a)]

/-- Any bijection of vertex sets is an isomorphism of edgeless graphs. -/
def botIso {A B : Type*} (e : A ≃ B) : (⊥ : SimpleGraph A) ≃g (⊥ : SimpleGraph B) where
  toEquiv := e
  map_rel_iff' := by simp

theorem sum_bot_bot {A B : Type*} :
    (⊥ : SimpleGraph A).sum (⊥ : SimpleGraph B) = (⊥ : SimpleGraph (A ⊕ B)) := by
  ext x y
  cases x <;> cases y <;> simp [SimpleGraph.sum]

/-! ### `L` cloned `K₂`s normalize to `1` -/

/-- **The cloned `K₂`s cancel their factorials.**  A clan graph consisting of `L + 1`
isolated vertices of multiplicity two — i.e. `L + 1` cloned `K₂` components — has
normalized imbalance polynomial `1`. -/
theorem Wpoly_bot_two (L : ℕ) : Wpoly (⊥ : SimpleGraph (Fin (L + 1))) (fun _ => 2) = 1 := by
  induction L with
  | zero => exact Wpoly_isolated_two
  | succ L ih =>
      have h1 : Wpoly (⊥ : SimpleGraph (Fin (L + 1 + 1))) (fun _ => 2)
          = Wpoly (⊥ : SimpleGraph (Fin (L + 1) ⊕ Fin 1)) (fun _ => 2) :=
        Wpoly_of_iso (botIso (finSumFinEquiv.symm :
          Fin (L + 1 + 1) ≃ Fin (L + 1) ⊕ Fin 1)) (fun _ => 2)
      have h2 : Wpoly (⊥ : SimpleGraph (Fin (L + 1) ⊕ Fin 1)) (fun _ => 2)
          = Wpoly ((⊥ : SimpleGraph (Fin (L + 1))).sum (⊥ : SimpleGraph (Fin 1)))
              (Sum.elim (fun _ => 2) (fun _ => 2)) := by
        rw [sum_bot_bot]
        congr 1
        funext x
        cases x <;> rfl
      rw [h1, h2, Wpoly_sum_elim, ih, Wpoly_isolated_two, mul_one]

/-! ### The two normalized weights -/

/-- A spider all of whose arms have even length has imbalance `1`: its hub component is
a path-like bipartite graph with one more vertex on the hub side. -/
theorem Wpoly_spider_zero (n : ℕ) (len : Fin n → ℕ)
    (hp : (Finset.univ.filter (fun i : Fin n => len i % 2 = 1)).card = 0) :
    Wpoly (spider n len) (fun _ => 1) = Pw 1 := by
  classical
  have himb := spider_imbalance n len
  rw [hp] at himb
  have h1 : ((spiderEven n len).card : ℤ) - (spiderOdd n len).card = 1 := by omega
  have h2 : ((spiderOdd n len).card : ℤ) - (spiderEven n len).card = -1 := by omega
  rw [Wpoly, clanFactorial_one, imbalanceGF_iso (clanOneIso (spider n len)),
    imbalanceGF_connected (spider_connected n len) (spider_isStablePair n len), h1, h2, Pw]
  simp

/-- The active state: hub plus two odd prefixes.  Its normalized weight is `z + z⁻¹`,
the term `z^r + z^(-r)` of `A(r,c;z)` for `r = p - 1 = 1`. -/
theorem Wpoly_active_two_arms {L M : ℕ} (hL : L % 2 = 1) (hM : M % 2 = 1) :
    Wpoly (spider 2 ![L, M]) (fun _ => 1) = Pw 1 := by
  classical
  refine Wpoly_spider 2 ![L, M] 1 ?_
  have hall : (Finset.univ.filter (fun i : Fin 2 => ![L, M] i % 2 = 1)) = Finset.univ := by
    refine Finset.filter_true_of_mem fun i _ => ?_
    fin_cases i
    · simpa using hL
    · simpa using hM
  rw [hall]
  simp

/-- The image state: `L` cloned `K₂`s and the remainder of the arm `B`.  Its normalized
weight is again `z + z⁻¹`; the cloned `K₂`s contribute the factor `1`. -/
theorem Wpoly_image_two_arms {L M : ℕ} (hL : L % 2 = 1) (hM : M % 2 = 1) (hLM : L ≤ M) :
    Wpoly ((⊥ : SimpleGraph (Fin L)).sum (spider 1 (fun _ => M - L)))
        (Sum.elim (fun _ => 2) (fun _ => 1)) = Pw 1 := by
  classical
  obtain ⟨K, rfl⟩ : ∃ K, L = K + 1 := ⟨L - 1, by omega⟩
  have hpath : Wpoly (spider 1 (fun _ => M - (K + 1))) (fun _ => 1) = Pw 1 := by
    refine Wpoly_spider_zero 1 _ ?_
    have hempty : (Finset.univ.filter
        (fun _ : Fin 1 => (M - (K + 1)) % 2 = 1)) = (∅ : Finset (Fin 1)) := by
      refine Finset.filter_false_of_mem fun i _ => ?_
      omega
    rw [hempty]
    simp
  rw [Wpoly_sum_elim, Wpoly_bot_two, hpath, one_mul]

/-- **The `p = 2` two-arm block has exactly the claimed normalized weight**, with the
scalar `c` derived to be `1`:

```text
W(active) + W(image) = A(1, 1; z) = (z + z⁻¹) + (z + z⁻¹).
```

This is the `p = 2` instance of the request's
`localMapP2_has_claimed_normalized_two_row_weight`, for a hub with exactly two active
arms (both prefixes odd, of lengths `L ≤ M`); see the file header for the scope. -/
theorem localMapP2_normalized_weight_two_arms {L M : ℕ} (hL : L % 2 = 1) (hM : M % 2 = 1)
    (hLM : L ≤ M) :
    Wpoly (spider 2 ![L, M]) (fun _ => 1)
      + Wpoly ((⊥ : SimpleGraph (Fin L)).sum (spider 1 (fun _ => M - L)))
          (Sum.elim (fun _ => 2) (fun _ => 1))
      = Ablock 1 1 := by
  rw [Wpoly_active_two_arms hL hM, Wpoly_image_two_arms hL hM hLM, Ablock, zz_eq, Pw]
  simp only [pow_one, one_smul, Nat.cast_one]

end ClanAudit
