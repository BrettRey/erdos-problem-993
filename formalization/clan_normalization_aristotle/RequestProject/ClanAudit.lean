import RequestProject.ClanAudit.Laurent
import RequestProject.ClanAudit.Blocks
import RequestProject.ClanAudit.Stable
import RequestProject.ClanAudit.Clan
import RequestProject.ClanAudit.Normalized
import RequestProject.ClanAudit.LocalMap
import RequestProject.ClanAudit.Weight
import RequestProject.ClanAudit.Spider
import RequestProject.ClanAudit.P2Weight

/-!
# `ClanAudit` — source-labelled summary of the audit

This file collects, in one place, the definitions and theorems of the audit, each
labelled with the item of the formalization request (and of Li–Li–Yang–Zhang,
arXiv:2501.04245) that it corresponds to.  Every statement here is a restatement of a
theorem proved in one of the imported files; no new mathematics is introduced.

Layout:

* `RequestProject/ClanAudit/Laurent.lean` — Laurent coefficient sequences and central
  unimodality;
* `RequestProject/ClanAudit/Blocks.lean` — the proposed blocks `A(r,c;z)` and `F(z)`;
* `RequestProject/ClanAudit/Stable.lean` — stable two-block partitions, the imbalance
  generating function, components and multiplicativity;
* `RequestProject/ClanAudit/Clan.lean` — clan graphs and the clan normalization;
* `RequestProject/ClanAudit/Normalized.lean` — the normalized polynomial `W` and the
  normalized two-row coefficient;
* `RequestProject/ClanAudit/LocalMap.lean` — the published spider transformation at
  `p = 2`;
* `RequestProject/ClanAudit/Weight.lean` — isomorphism invariance, multiplicativity of
  the normalized `W`, and the cloned `K₂` cancellation;
* `RequestProject/ClanAudit/Spider.lean` — the hub component of an active spider state
  and its imbalance `r = p - 1`;
* `RequestProject/ClanAudit/P2Weight.lean` — the normalized weight of a `p = 2` two-arm
  block.

See `README.md` for the traceability table and `RESULT.md` for the graded result note.
-/

namespace ClanAudit

open LaurentPolynomial

namespace Audit

variable {V : Type*} [Fintype V] [DecidableEq V]
variable {W₁ W₂ : Type*} [Fintype W₁] [DecidableEq W₁] [Fintype W₂] [DecidableEq W₂]

/-! ## Mandatory normalization kernel, item 1: clan normalization

Source: request, "Mandatory normalization kernel" 1; LLYZ Section 2 (clan graph and the
identity `Xchrom (clan G alpha) = (∏ (alpha v)!) * Xmulti G alpha`).

`clan G alpha` (`RequestProject/ClanAudit/Clan.lean`) is the clan graph, and
`IsMulticoloring` is a proper multicolouring of type `alpha`.  The two theorems below are
the two halves of the normalization identity at the level of finite colouring counts:
the fibres of the map "colouring of the clan graph ↦ its multicolouring" all have size
`∏ (alpha v)!`, and the monomial is preserved.  Summing over multicolourings therefore
gives `Xchrom (clan G alpha) = (∏ (alpha v)!) * Xmulti G alpha`. -/
omit [DecidableEq V] in
theorem clan_normalization_fibre_card {G : SimpleGraph V} {alpha : V → ℕ} {m : V → Finset ℕ}
    (hm : IsMulticoloring G alpha m) :
    Nat.card {f : (Σ v : V, Fin (alpha v)) → ℕ //
        IsProperColoring (clan G alpha) f ∧ toMulti alpha f = m}
      = ∏ v : V, Nat.factorial (alpha v) :=
  clan_fiber_card hm

omit [DecidableEq V] in
theorem clan_normalization_monomial {G : SimpleGraph V} {alpha : V → ℕ}
    {f : (Σ v : V, Fin (alpha v)) → ℕ} (hf : IsProperColoring (clan G alpha) f) (c : ℕ) :
    monomialOfColoring f c = monomialOfMulti (toMulti alpha f) c :=
  monomial_toMulti hf c

/-! ## Two-row coefficients

Source: request, "Source and notation"; LLYZ Proposition 2.4.  As the request permits,
`twoRowCoeff H k l = stableCount H k l - stableCount H (k+1) (l-1)` is taken as the
definition of `[s_(k,l)] Xchrom H`.  The theorem below expresses it through the
imbalance generating function `∑ z^(|A|-|B|)`. -/
theorem two_row_coeff_eq_imbalance_coeff {W : Type*} [Fintype W] [DecidableEq W]
    (H : SimpleGraph W) (k l : ℕ) (hl : 1 ≤ l) (h : k + l = Fintype.card W) :
    (twoRowCoeff H k l : ℚ)
      = coeffL (imbalanceGF H) ((k : ℤ) - l) - coeffL (imbalanceGF H) ((k : ℤ) - l + 2) :=
  twoRowCoeff_eq_coeff H k l hl h

/-! ## Mandatory normalization kernel, item 2: weighted bipartition formula

Source: request, "Mandatory normalization kernel" 2 and the candidate `W(G, alpha; z)`.

The candidate is *derived*, not assumed: `imbalanceGF_connected` computes the imbalance
polynomial of a connected bipartite graph, fixing all boundary conventions (a balanced
component contributes `2`, an isolated vertex contributes `z + z⁻¹`, a nonbipartite
graph contributes `0`), and `normalizedTwoRowCoeff_eq` is the interior formula
`normalizedTwoRowCoeff alpha k l = coeff W (k-l) - coeff W (k-l+2)`. -/
theorem component_orientation_weight {W : Type*} [Fintype W] [DecidableEq W]
    {H : SimpleGraph W} (hc : H.Connected) {A B : Finset W} (h : IsStablePair H A B) :
    imbalanceGF H = T ((A.card : ℤ) - B.card) + T ((B.card : ℤ) - A.card) :=
  imbalanceGF_connected hc h

theorem nonbipartite_weight_zero {W : Type*} [Fintype W] [DecidableEq W]
    {H : SimpleGraph W} (h : ¬ H.Colorable 2) : imbalanceGF H = 0 :=
  imbalanceGF_eq_zero_of_not_colorable h

theorem normalized_two_row_coeff_formula (G : SimpleGraph V) (alpha : V → ℕ) (k l : ℕ)
    (hl : 1 ≤ l) (h : k + l = ∑ v : V, alpha v) :
    normalizedTwoRowCoeff G alpha k l
      = coeffL (Wpoly G alpha) ((k : ℤ) - l) - coeffL (Wpoly G alpha) ((k : ℤ) - l + 2) :=
  normalizedTwoRowCoeff_eq G alpha k l hl h

/-! ## Mandatory normalization kernel, item 3: products over disjoint components

Source: request, "Mandatory normalization kernel" 3. -/
theorem imbalance_multiplicative (H₁ : SimpleGraph W₁) (H₂ : SimpleGraph W₂) :
    imbalanceGF (H₁.sum H₂) = imbalanceGF H₁ * imbalanceGF H₂ :=
  imbalanceGF_sum H₁ H₂

theorem imbalance_product_of_components {H₁ : SimpleGraph W₁} {H₂ : SimpleGraph W₂}
    (hc₁ : H₁.Connected) (hc₂ : H₂.Connected) {A₁ B₁ : Finset W₁} {A₂ B₂ : Finset W₂}
    (h₁ : IsStablePair H₁ A₁ B₁) (h₂ : IsStablePair H₂ A₂ B₂) :
    imbalanceGF (H₁.sum H₂)
      = (T ((A₁.card : ℤ) - B₁.card) + T ((B₁.card : ℤ) - A₁.card))
        * (T ((A₂.card : ℤ) - B₂.card) + T ((B₂.card : ℤ) - A₂.card)) :=
  imbalanceGF_sum_connected hc₁ hc₂ h₁ h₂

/-! ## Vanishing of the two-row contribution at a multiplicity-two neighbour

Source: request, "Published spider transformation", first paragraph ("If a neighbor of
the hub has multiplicity at least two, the clan graph contains a triangle and its
two-row contribution is zero"). -/
theorem two_row_zero_of_neighbour_multiplicity_two {G : SimpleGraph V} {alpha : V → ℕ}
    {u v : V} (huv : G.Adj u v) (hu : 2 ≤ alpha u) (hv : 1 ≤ alpha v) (k l : ℕ) (hl : 1 ≤ l)
    (h : k + l = ∑ v : V, alpha v) : normalizedTwoRowCoeff G alpha k l = 0 :=
  clan_normalizedTwoRowCoeff_eq_zero_of_mult_two huv hu hv k l hl h

/-! ## The published spider transformation at `p = 2`

Source: request, "Published spider transformation", `localMapP2_preserves_total_order`
and `localMapP2_injective`.  (`localMapP2_has_claimed_normalized_two_row_weight` is proved
below in the two-arm case only; see `RESULT.md`.) -/
theorem localMapP2_preserves_total_order' {s : HubState} {i₀ i₁ L : ℕ}
    (h : Admissible s i₀ i₁ L) : (localMapP2 s i₀ i₁ L).total = s.total :=
  localMapP2_preserves_total_order h

theorem localMapP2_injective' {s s' : HubState} {i₀ i₁ L i₀' i₁' L' : ℕ}
    (h : Admissible s i₀ i₁ L) (h' : Admissible s' i₀' i₁' L')
    (heq : localMapP2 s i₀ i₁ L = localMapP2 s' i₀' i₁' L') : s = s' :=
  localMapP2_injective h h' heq

/-! ## The normalized weight of the transformation at `p = 2`

Source: request, "Published spider transformation",
`localMapP2_has_claimed_normalized_two_row_weight`, and "Adjacent two-hub target" 4.

The cloned `K₂` cancellation ("the cloned `K2` orientation counts are supposed to cancel
the new factors of `2!`") is `cloned_K2_cancellation`; the term `z^r + z^(-r)` of the
proposed block `A(r,c;z)`, together with `r = p - 1`, is `active_hub_component_weight`;
and `p2_block_normalized_weight` computes the whole two-state block for a hub with
exactly two active arms, deriving `c = 1` there. -/
theorem cloned_K2_cancellation (L : ℕ) :
    Wpoly (⊥ : SimpleGraph (Fin (L + 1))) (fun _ => 2) = 1 :=
  Wpoly_bot_two L

theorem active_hub_component_weight (n : ℕ) (len : Fin n → ℕ) (r : ℕ)
    (hp : (Finset.univ.filter (fun i : Fin n => len i % 2 = 1)).card = r + 1) :
    Wpoly (spider n len) (fun _ => 1) = Pw r :=
  Wpoly_spider n len r hp

theorem p2_block_normalized_weight {L M : ℕ} (hL : L % 2 = 1) (hM : M % 2 = 1)
    (hLM : L ≤ M) :
    Wpoly (spider 2 ![L, M]) (fun _ => 1)
      + Wpoly ((⊥ : SimpleGraph (Fin L)).sum (spider 1 (fun _ => M - L)))
          (Sum.elim (fun _ => 2) (fun _ => 1))
      = Ablock 1 1 :=
  localMapP2_normalized_weight_two_arms hL hM hLM

/-! ## Adjacent two-hub target, item 5: the four-map identity

Source: request, "Adjacent two-hub target", identity for `F(z)`. -/
theorem four_map_identity {r s : ℕ} (h : s ≤ r) (c d : ℚ) :
    Fblock r s c d
      = (c * d) • zz ^ (r + s) + c • (zz ^ r * Pw s) + d • (zz ^ s * Pw r) + Pw (r - s) :=
  Fblock_expand h c d

/-! ## Adjacent two-hub target, item 6: central unimodality — REFUTED

Source: request, "Adjacent two-hub target" 6.  The claim fails for general `c, d > 0`;
the explicit witness is `r = 3, s = 1, c = 1/4, d = 1`, where `coeff F 0 = 3 < 4 =
coeff F 2`, i.e. the normalized two-row coefficient of the block is `-1`.  It already
fails for a single hub block `A(3, 1/4; z)`.  It does hold once `1 ≤ c` and `1 ≤ d`. -/
theorem central_unimodality_refuted :
    ∃ (r s : ℕ) (c d : ℚ) (m : ℤ), 1 ≤ r ∧ 1 ≤ s ∧ 0 < c ∧ 0 < d ∧ 0 ≤ m ∧
      coeffL (Fblock r s c d) m < coeffL (Fblock r s c d) (m + 2) :=
  Fblock_not_centrally_unimodal

theorem central_unimodality_witness :
    coeffL (Fblock 3 1 (1/4) 1) 0 = 3 ∧ coeffL (Fblock 3 1 (1/4) 1) 2 = 4 :=
  Fblock_coeffs_of_counterexample

theorem central_unimodality_of_one_le {r s : ℕ} (hr : 1 ≤ r) (hs : 1 ≤ s) {c d : ℚ}
    (hc : 1 ≤ c) (hd : 1 ≤ d) (m : ℤ) (hm : 0 ≤ m) :
    coeffL (Fblock r s c d) (m + 2) ≤ coeffL (Fblock r s c d) m :=
  Fblock_decr hr hs hc hd m hm

end Audit

end ClanAudit
