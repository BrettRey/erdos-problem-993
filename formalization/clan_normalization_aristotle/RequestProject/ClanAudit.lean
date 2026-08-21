import RequestProject.ClanAudit.Laurent
import RequestProject.ClanAudit.Blocks
import RequestProject.ClanAudit.Stable
import RequestProject.ClanAudit.Clan
import RequestProject.ClanAudit.Normalized
import RequestProject.ClanAudit.LocalMap
import RequestProject.ClanAudit.Weight
import RequestProject.ClanAudit.Spider
import RequestProject.ClanAudit.P2Weight
import RequestProject.ClanAudit.EvenArms
import RequestProject.ClanAudit.GlobalPartition
import RequestProject.ClanAudit.Collision
import RequestProject.ClanAudit.AdjacentLogConcave
import RequestProject.ClanAudit.TreeC2
import RequestProject.ClanAudit.ModelsAreTrees

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
  block;
* `RequestProject/ClanAudit/EvenArms.lean` — the same for a hub carrying, besides the two
  odd arms, an arbitrary finite family of further arms with even prefixes, deriving the
  scalar `c = 2^e`;
* `RequestProject/ClanAudit/PathClan.lean`, `Split.lean`, `Arms.lean` — the good-shape
  weight of an arbitrary path clan map and the vertex-splitting toolkit;
* `RequestProject/ClanAudit/DoubleSpider.lean` — the adjacent two-hub tree `dgraph` and
  its exact imbalance law `|p - q|`;
* `RequestProject/ClanAudit/Prefix.lean` — active prefixes, the exact weight of one side;
* `RequestProject/ClanAudit/SideTransform.lean` — the canonical transformation `transf` of
  a whole side, with the published shortest-odd-prefix choice; injectivity and order
  preservation;
* `RequestProject/ClanAudit/ImageWeight.lean` — the exact weight of a transformed side,
  with the scalar `c = 2^e ≥ 1` derived at arbitrary `p ≥ 2`;
* `RequestProject/ClanAudit/TwoHubWeight.lean` — the exact weight of a two-hub clan map;
* `RequestProject/ClanAudit/BlockCU.lean` — central unimodality of `A(r,c)` and
  `F(r,s;c,d)` on the derived range `c, d ≥ 1`;
* `RequestProject/ClanAudit/Vanishing.lean` — the triangle/nonbipartite vanishing and the
  good shape of every non-active side;
* `RequestProject/ClanAudit/BlockWeight.lean` — the exact normalized Laurent sum of each
  block of the partition;
* `RequestProject/ClanAudit/GlobalPartition.lean` — the explicit global partition and the
  final degreewise theorem;
* `RequestProject/ClanAudit/Collision.lean` — the explicit collision of the proof notes,
  resolved.

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
below, in the general form allowing arbitrarily many further even arms; see `RESULT.md`.) -/
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

/-! ## The normalized weight of the `p = 2` block with arbitrary further even arms

Source: request, "Published spider transformation",
`localMapP2_has_claimed_normalized_two_row_weight`, general (arbitrary-arm) case; and the
continuation task `FOLLOWUP_ARBITRARY_EVEN_ARMS_20260820.md`.

The hub carries the two arms with odd positive prefixes `L ≤ M` and an arbitrary finite
family of `e` further active arms with positive *even* prefixes `len i`.  The image clan
graph is decomposed explicitly: the `L` cloned `K₂`s (whose orientation counts cancel the
new factors `2!`), the untouched remainder of the arm `B`, and one balanced path component
per even arm.  `even_arm_family_weight` is the reusable finite-family product theorem for
those components, `image_state_normalized_weight` is the weight of the whole image state,
and `p2_block_normalized_weight_even_arms` computes the two-state block, deriving the
scalar `c = 2^e`; `p2_block_derived_scalar_one_le` records `1 ≤ c`, the hypothesis under
which the adjacent two-hub conclusion is true. -/
theorem even_arm_family_weight (e : ℕ) (len : Fin e → ℕ) (heven : ∀ i, len i % 2 = 0)
    (hpos : ∀ i, 0 < len i) : Wpoly (armsGraph e len) (fun _ => 1) = 2 ^ e :=
  Wpoly_armsGraph e len heven hpos

theorem image_state_normalized_weight {e L M : ℕ} {len : Fin e → ℕ} (hL : L % 2 = 1)
    (hM : M % 2 = 1) (hLM : L ≤ M) (heven : ∀ i, len i % 2 = 0) (hpos : ∀ i, 0 < len i) :
    Wpoly (((⊥ : SimpleGraph (Fin L)).sum (spider 1 (fun _ => M - L))).sum (armsGraph e len))
        (Sum.elim (Sum.elim (fun _ => 2) (fun _ => 1)) (fun _ => 1))
      = 2 ^ e * Pw 1 :=
  Wpoly_image_even_arms hL hM hLM heven hpos

theorem p2_block_normalized_weight_even_arms {e L M : ℕ} {len : Fin e → ℕ} (hL : L % 2 = 1)
    (hM : M % 2 = 1) (hLM : L ≤ M) (heven : ∀ i, len i % 2 = 0) (hpos : ∀ i, 0 < len i) :
    Wpoly (spider (2 + e) (Fin.append ![L, M] len)) (fun _ => 1)
      + Wpoly (((⊥ : SimpleGraph (Fin L)).sum (spider 1 (fun _ => M - L))).sum (armsGraph e len))
          (Sum.elim (Sum.elim (fun _ => 2) (fun _ => 1)) (fun _ => 1))
      = Ablock 1 ((2 : ℚ) ^ e) :=
  localMapP2_normalized_weight_even_arms hL hM hLM heven hpos

theorem p2_block_derived_scalar_one_le (e : ℕ) : (1 : ℚ) ≤ (2 : ℚ) ^ e :=
  one_le_derived_scalar e

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

/-! ## Adjacent two-hub target, checkpoint 1: the global model

Source: `FOLLOWUP_GLOBAL_PARTITION_20260820.md`, required checkpoint 1.

`dgraph m a n b` is the adjacent two-hub tree with arbitrary finite ordered families of
positive-length pendant paths at both hubs; `mapsOfOrder m a n b N` is the finite set of
all clan maps of total order `N`; `pref` is the active prefix of an arm, `ActiveSide` the
hub-active stratum, and `idx0`, `idx1`, `plen` the published deterministic choices
(shortest odd prefix, ties by arm order), proved well defined by `canonical_spec`. -/
theorem degree_map_space_membership {m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ} {N : ℕ}
    {alpha : DV m a n b → ℕ} : alpha ∈ mapsOfOrder m a n b N ↔ ∑ x, alpha x = N :=
  mem_mapsOfOrder

theorem canonical_choices_well_defined {k : ℕ} {len : Fin k → ℕ} {s : SpiderV k len → ℕ}
    (h : ActiveSide s) :
    ∃ i₀ i₁ : Fin k, i₀.val = idx0 s ∧ i₁.val = idx1 s ∧ pref s i₀ = plen s ∧
      pref s i₀ % 2 = 1 ∧ pref s i₁ % 2 = 1 ∧ i₀ ≠ i₁ ∧ plen s ≤ pref s i₁ :=
  canonical_spec h

/-! ## Adjacent two-hub target, checkpoint 2: classification

Source: `FOLLOWUP_GLOBAL_PARTITION_20260820.md`, required checkpoint 2.  A multiplicity
`≥ 2` next to a positive vertex makes the clan graph nonbipartite (it contains a
triangle), so the whole normalized weight — and hence every two-row coefficient —
vanishes.  Every side that is not active therefore has a weight of the good shape
`2^j (z + z⁻¹)^i`, possibly zero. -/
theorem two_row_zero_of_clan_triangle {V : Type*} [Fintype V] [DecidableEq V]
    {G : SimpleGraph V} {alpha : V → ℕ} {u v : V} (huv : G.Adj u v) (hu : 2 ≤ alpha u)
    (hv : 1 ≤ alpha v) : Wpoly G alpha = 0 :=
  Wpoly_eq_zero_of_mult_two huv hu hv

theorem inactive_side_good_shape {k : ℕ} {len : Fin k → ℕ} {s : SpiderV k len → ℕ}
    (h : ¬ ActiveSide s) : IsGoodW (Wpoly (spider k len) s) :=
  isGoodW_Wpoly_spider_of_not_active h

/-! ## Adjacent two-hub target, checkpoint 3: the collision audit

Source: `FOLLOWUP_GLOBAL_PARTITION_20260820.md`, required checkpoint 3, and
`two_hub_clan_cancellation_attack_2026-08-20.md`, "the global collision".

The source and image strata are disjoint because the transformation switches the hub off,
so a map whose side is an image is not the start of a second, competing block: it lies in
the *same* block as the untransformed map.  The explicit collision of the notes (three
unit arms at `u`, five at `v`) is formalized and shown to lie in a single four-element
block. -/
theorem source_and_image_strata_disjoint {k : ℕ} {len : Fin k → ℕ} {s : SpiderV k len → ℕ}
    (h : ImgSide s) : ¬ ActiveSide s :=
  not_activeSide_of_imgSide h

theorem single_and_double_hub_strata_do_not_compete {m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ}
    {u : SpiderV m a → ℕ} (h : ActiveSide u) (v : SpiderV n b → ℕ) :
    rep (Sum.elim (transf u) v) = rep (Sum.elim u v) :=
  collision_resolved_left h v

theorem notes_collision_in_one_block :
    collisionSource ∈ block 3 (unitArms 3) 5 (unitArms 5) 10 collisionSource
      ∧ collisionPartner ∈ block 3 (unitArms 3) 5 (unitArms 5) 10 collisionSource
      ∧ (block 3 (unitArms 3) 5 (unitArms 5) 10 collisionSource).card = 4 :=
  collision_block

/-! ## Adjacent two-hub target, checkpoints 2–4: the explicit global partition

Source: `FOLLOWUP_GLOBAL_PARTITION_20260820.md`, "Final target" 1–3 and required
checkpoint 4.  The blocks are the fibres of the explicit idempotent representative map
`rep`; each is exactly the set of maps obtained by keeping or transforming each side
independently, so it has `4`, `2` or `1` elements.  Disjointness is automatic for a fibre
decomposition, exhaustion is `blocks_cover`, and `sum_rep` is the preservation of the
total order. -/
theorem global_partition_blocks_explicit {m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ} {N : ℕ}
    {beta : DV m a n b → ℕ} (hb : rep beta = beta) (hmem : beta ∈ mapsOfOrder m a n b N) :
    block m a n b N beta
      = ((sideBlock (fun x => beta (Sum.inl x))) ×ˢ (sideBlock (fun y => beta (Sum.inr y)))).image
          (fun p => Sum.elim p.1 p.2) :=
  fiber_eq_image hb hmem

open scoped Classical in
theorem global_partition_block_sizes {m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ} {N : ℕ}
    {beta : DV m a n b → ℕ} (hb : rep beta = beta) (hmem : beta ∈ mapsOfOrder m a n b N) :
    (block m a n b N beta).card
      = if ActiveSide (fun x => beta (Sum.inl x)) then
          (if ActiveSide (fun y => beta (Sum.inr y)) then 4 else 2)
        else (if ActiveSide (fun y => beta (Sum.inr y)) then 2 else 1) :=
  card_block_cases hb hmem

theorem global_partition_disjoint {m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ} (N : ℕ)
    {beta gamma : DV m a n b → ℕ} (h : beta ≠ gamma) :
    Disjoint (block m a n b N beta) (block m a n b N gamma) :=
  blocks_pairwise_disjoint N h

theorem global_partition_exhaustive {m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ} {N : ℕ}
    {alpha : DV m a n b → ℕ} (h : alpha ∈ mapsOfOrder m a n b N) :
    alpha ∈ block m a n b N (rep alpha) ∧ rep alpha ∈ mapsOfOrder m a n b N :=
  blocks_cover h

theorem global_partition_preserves_total_order {m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ}
    (alpha : DV m a n b → ℕ) : ∑ x, rep alpha x = ∑ x, alpha x :=
  sum_rep alpha

/-! ## Adjacent two-hub target, checkpoint 5: the weight of each block

Source: `FOLLOWUP_GLOBAL_PARTITION_20260820.md`, required checkpoint 5.  The four-state
block sums to `F(r,s;c,d)` times the product of all outside components, with `r, s ≥ 1`
and the scalars `c, d = 2^e ≥ 1` *derived*; the two-state blocks sum to `A(r,c)` times
outside components; the singletons are products of good factors.  `Wpoly_side_image` is
the missing normalized-weight theorem at arbitrary `p ≥ 2` (in particular `p ≥ 3`). -/
theorem transformed_side_weight {k : ℕ} {len : Fin k → ℕ} {s : SpiderV k len → ℕ}
    (h : ActiveSide s) :
    ∃ c : ℚ, 1 ≤ c ∧
      Wpoly (spider k len) (transf s) = c • (zz ^ (pNum s - 1) * ∏ i : Fin k, tailW s i) :=
  Wpoly_side_image h

theorem four_state_block_weight {m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ}
    {u : SpiderV m a → ℕ} {v : SpiderV n b → ℕ} (hu : ActiveSide u) (hv : ActiveSide v) :
    IsCU (Wpoly (dgraph m a n b) (Sum.elim u v)
      + Wpoly (dgraph m a n b) (Sum.elim u (transf v))
      + (Wpoly (dgraph m a n b) (Sum.elim (transf u) v)
        + Wpoly (dgraph m a n b) (Sum.elim (transf u) (transf v)))) :=
  blockSum_both_active hu hv

theorem blockwise_central_unimodality {m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ} {N : ℕ}
    {beta : DV m a n b → ℕ} (hb : rep beta = beta) (hmem : beta ∈ mapsOfOrder m a n b N) :
    IsCU (∑ alpha ∈ block m a n b N beta, Wpoly (dgraph m a n b) alpha) :=
  isCU_sum_block hb hmem

/-! ## Adjacent two-hub target, checkpoint 6: the final degreewise theorem

Source: `FOLLOWUP_GLOBAL_PARTITION_20260820.md`, "the strongest desired final theorem"
and required checkpoint 6.  Summed over *all* clan maps of total order `N = k + l`, every
normalized two-row coefficient of the adjacent two-hub tree is nonnegative — for
arbitrary arm families and arbitrary `N`. -/
theorem total_degree_weight_centrally_unimodal (m : ℕ) (a : Fin m → ℕ) (n : ℕ)
    (b : Fin n → ℕ) (N : ℕ) : IsCU (totalW m a n b N) :=
  isCU_totalW m a n b N

theorem adjacent_two_hub_degreewise_nonneg (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ)
    (N k l : ℕ) (hl : 1 ≤ l) (hkl : l ≤ k) (hN : k + l = N) :
    0 ≤ ∑ alpha ∈ mapsOfOrder m a n b N, normalizedTwoRowCoeff (dgraph m a n b) alpha k l :=
  sum_normalizedTwoRowCoeff_nonneg m a n b N k l hl hkl hN

/-! ## Target A: the direct Li–Li–Yang–Zhang coefficient bridge

Source: the continuation request, "Target A: direct Li–Li–Yang–Zhang coefficient bridge".
The finite coefficient identity is *derived* from the clan/multicolouring definitions, and
turned into a reusable transfer theorem which is then instantiated for `dgraph`. -/

theorem coefficient_bridge {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    (k l N : ℕ) (hl : 1 ≤ l) (hN : k + l = N) :
    ∑ alpha ∈ multMaps V N, normalizedTwoRowCoeff G alpha k l
      = (indepCount G k : ℚ) * (indepCount G l : ℚ)
        - (indepCount G (k + 1) : ℚ) * (indepCount G (l - 1) : ℚ) :=
  sum_normalizedTwoRowCoeff_eq G k l N hl hN

theorem coefficient_bridge_diagonal {V : Type*} [Fintype V] [DecidableEq V]
    (G : SimpleGraph V) (k : ℕ) (hk : 1 ≤ k) :
    ∑ alpha ∈ multMaps V (2 * k), normalizedTwoRowCoeff G alpha k k
      = (indepCount G k : ℚ) ^ 2
        - (indepCount G (k - 1) : ℚ) * (indepCount G (k + 1) : ℚ) :=
  sum_normalizedTwoRowCoeff_diag G k hk

theorem coefficient_bridge_zero {V : Type*} [Fintype V] [DecidableEq V]
    (G : SimpleGraph V) :
    ∑ alpha ∈ multMaps V 0, normalizedTwoRowCoeff G alpha 0 0 = (indepCount G 0 : ℚ) ^ 2 :=
  sum_normalizedTwoRowCoeff_zero G

theorem degreewise_nonneg_implies_log_concavity {V : Type*} [Fintype V] [DecidableEq V]
    (G : SimpleGraph V)
    (h : ∀ N k l : ℕ, 1 ≤ l → l ≤ k → k + l = N →
      0 ≤ ∑ alpha ∈ multMaps V N, normalizedTwoRowCoeff G alpha k l) (j : ℕ) :
    indepCount G j * indepCount G (j + 2) ≤ indepCount G (j + 1) * indepCount G (j + 1) :=
  indepCount_logConcave_of_degreewise_nonneg G h j

theorem adjacent_two_hub_independence_polynomial_log_concave (m : ℕ) (a : Fin m → ℕ)
    (n : ℕ) (b : Fin n → ℕ) (j : ℕ) :
    (indepPoly (dgraph m a n b)).coeff j * (indepPoly (dgraph m a n b)).coeff (j + 2)
      ≤ (indepPoly (dgraph m a n b)).coeff (j + 1)
        * (indepPoly (dgraph m a n b)).coeff (j + 1) :=
  dgraph_indepPoly_logConcave m a n b j

/-! ### The `C₂` family: from the clan machinery to abstract trees -/

/-- Source: the connector two-hub tree, arbitrary connector length and arbitrary arms. -/
theorem connector_two_hub_independence_polynomial_log_concave (t m : ℕ) (a : Fin m → ℕ)
    (n : ℕ) (b : Fin n → ℕ) (j : ℕ) :
    (indepPoly (connGraph t m a n b)).coeff j
        * (indepPoly (connGraph t m a n b)).coeff (j + 2)
      ≤ (indepPoly (connGraph t m a n b)).coeff (j + 1)
        * (indepPoly (connGraph t m a n b)).coeff (j + 1) :=
  connGraph_indepPoly_logConcave t m a n b j

/-- Source: the one-hub case, i.e. spiders (and hence paths). -/
theorem spider_independence_polynomial_log_concave (k : ℕ) (len : Fin k → ℕ) (j : ℕ) :
    (indepPoly (spider k len)).coeff j * (indepPoly (spider k len)).coeff (j + 2)
      ≤ (indepPoly (spider k len)).coeff (j + 1) * (indepPoly (spider k len)).coeff (j + 1) :=
  spider_indepPoly_logConcave k len j

/-- Source: recognition of a tree with at most one branch vertex as a spider. -/
noncomputable def tree_one_branch_vertex_is_spider {V : Type*} [Fintype V] [DecidableEq V]
    {G : SimpleGraph V} [DecidableRel G.Adj] (hG : G.IsTree) (h : V)
    (hdeg : ∀ v : V, v ≠ h → G.degree v ≤ 2) :
    G ≃g spider (G.neighborFinset h).card (hubArmLen hG h) :=
  spiderIsoOfTree hG h hdeg

/-- Source: recognition of a tree with two branch vertices as a connector tree. -/
noncomputable def tree_two_branch_vertices_is_connector {V : Type*} [Fintype V] [DecidableEq V]
    {G : SimpleGraph V} [DecidableRel G.Adj] (hG : G.IsTree) {h1 h2 : V} (hne : h1 ≠ h2)
    (hdeg : ∀ w : V, w ≠ h1 → w ≠ h2 → G.degree w ≤ 2) :
    G ≃g connModel hG (h1 := h1) (h2 := h2) :=
  connIsoOfTree hG hne hdeg

/-- Source: the `C₂` conclusion for abstract finite trees. -/
theorem tree_with_two_branch_vertices_log_concave {V : Type*} [Fintype V] [DecidableEq V]
    {G : SimpleGraph V} [DecidableRel G.Adj] (hG : G.IsTree)
    (hb : (branchVerts G).card ≤ 2) (j : ℕ) :
    (indepPoly G).coeff j * (indepPoly G).coeff (j + 2)
      ≤ (indepPoly G).coeff (j + 1) * (indepPoly G).coeff (j + 1) :=
  indepPoly_logConcave_of_isTree_of_branchVerts_card_le_two hG hb j

/-- Source: the scope of the `C₂` family is exactly the finite trees with at most two
branch vertices. -/
theorem c2_family_is_exactly_two_branch_trees {V : Type*} [Fintype V] [DecidableEq V]
    {G : SimpleGraph V} [DecidableRel G.Adj] :
    IsC2Model G ↔ G.IsTree ∧ (branchVerts G).card ≤ 2 :=
  isC2Model_iff_isTree_branchVerts

end Audit

end ClanAudit
