# Result note

## Grade

**`C2_COMPLETE`** — the programme has been carried to its intended conclusion: every
finite tree with at most two vertices of degree at least three has a log-concave
independence polynomial (`ClanAudit.indepPoly_logConcave_of_isTree_of_branchVerts_card_le_two`).
This supersedes the earlier `ADJACENT_COMPLETE` grade, which recorded the explicit global
partition and the degreewise nonnegativity theorem for the adjacent two-hub tree (both
still proved; see "The explicit global partition and the final degreewise theorem"
below).

The earlier grades stand unchanged: the full **`NORMALIZATION_KERNEL`**, the `p = 2` local
package (**`P2_LOCAL`**) including the arbitrary-arm case of the weight theorem, and the
**`REFUTED`** verdict on item 6 of the original adjacent two-hub target — the proposed
central unimodality of `F` for arbitrary `c, d > 0` is false, and the global proof below
works precisely because the scalars it produces are *derived* to be powers of two.

* Items 1–3 of the mandatory normalization kernel are proved in full, from the
  definitions, with no axioms and no `sorry`.
* All three `p = 2` theorems asked for are proved in full: `localMapP2_injective`,
  `localMapP2_preserves_total_order`, and
  `localMapP2_has_claimed_normalized_two_row_weight` — the latter for a hub with the two
  arms of odd positive prefix (this is what `p = 2` means) together with an arbitrary
  finite family of `e` further active arms of even positive prefix, deriving the scalar
  `c = 2^e` (see "The normalized weight of the `p = 2` block" below).
* Item 5 of the adjacent two-hub target (the four-map identity) is proved.
* **Item 6 of the adjacent two-hub target is false as stated** and is refuted with an
  explicit finite witness.
* Items 2 and 3 of the adjacent two-hub target (the global exhaustive block partition)
  **are now formalized and proved**, together with the final degreewise theorem.

Nothing is hidden: the project contains no `sorry`, `admit`, `axiom`, or
`implemented_by`, and every theorem depends only on `propext`, `Classical.choice`,
`Quot.sound`.

## The refutation

The request proposes, for `A(r,c;z) = c (z + z⁻¹)^r + z^r + z^(-r)` with `c > 0` and

```text
F(z) = A(r,c;z) * A(s,d;z) - z^(r+s) - z^(-(r+s)),
```

that the Laurent coefficients of `F` are centrally unimodal, "so every normalized
two-row coefficient of the block is nonnegative".  This is false for general `c, d > 0`.

**Witness.**  `r = 3`, `s = 1`, `c = 1/4`, `d = 1`:

```text
coeff F 0 = 3,   coeff F 2 = 4,
```

so `coeff F 0 < coeff F 2`, and the corresponding normalized two-row coefficient
`coeff F 0 - coeff F 2 = -1` is *negative*.

* Lean: `ClanAudit.Fblock_coeffs_of_counterexample` (the two coefficients) and
  `ClanAudit.Fblock_not_centrally_unimodal` (the refutation), restated as
  `ClanAudit.Audit.central_unimodality_witness` and
  `ClanAudit.Audit.central_unimodality_refuted`.

The failure is not an artefact of taking a product of two blocks: it already occurs for a
single hub block, `A(3, 1/4; z)`, whose coefficients satisfy
`coeff 1 = 3/4 < 5/4 = coeff 3` (`ClanAudit.Ablock_not_centrally_unimodal`).

**Which proposed step fails.**  The four-map identity (item 5) is *correct*: for `s ≤ r`,

```text
F = c d (z+z⁻¹)^(r+s) + c (z+z⁻¹)^r (z^s+z^(-s)) + d (z+z⁻¹)^s (z^r+z^(-r))
      + z^(r-s) + z^(-(r-s))
```

(`ClanAudit.Fblock_expand`).  What fails is the step from this identity to central
unimodality.  Central unimodality is *not* preserved by the sum: the summand
`z^(r-s) + z^(-(r-s))` (and, for a single hub, `z^r + z^(-r)`) is a spike far from the
centre, and it is only dominated by the smooth summands `c (z+z⁻¹)^r …` when the scalar
`c` is large enough.  For small `c` the spike wins.  Hence the argument cannot be closed
without a lower bound on the derived scalars.

**Repaired statement (proved).**  As soon as `1 ≤ c` and `1 ≤ d`, the coefficients of `F`
*are* centrally unimodal, and hence every normalized two-row coefficient
`coeff F m - coeff F (m+2)`, `m ≥ 0`, of the block is nonnegative:
`ClanAudit.Fblock_decr` (= `ClanAudit.Audit.central_unimodality_of_one_le`).

So the request's insistence that "the scalar `c` must be derived from the normalized
outside components, it must not simply be postulated to be a positive integer" is exactly
the load-bearing point: the proposed conclusion needs `c, d ≥ 1`, which must come from a
derivation of `c`.

**On minimality of the witness.**  `c = 1/4, d = 1` was chosen for readability;
by `Fblock_decr` no counterexample exists with `c, d ≥ 1`, so the failure genuinely
requires a scalar below `1`.  For the exponents, a by-hand expansion of the remaining
pairs with `r + s ≤ 4` — `(1,1)`, `(2,1)`, `(2,2)` — shows no counterexample there, so
`(r,s) = (3,1)` is minimal in `r + s`; that side computation was done informally and is
**not** part of the Lean development, unlike the witness itself, which is fully verified.

## The normalized weight of the `p = 2` block

`localMapP2_has_claimed_normalized_two_row_weight` is proved in the following generality,
which is the full `p = 2` case: the hub has exactly two active arms with *odd* positive
prefixes, of lengths `L ≤ M` (so `p = 2` and `r = p - 1 = 1`), and an arbitrary finite
indexed family of `e ≥ 0` further active arms with *even* positive prefixes of lengths
`len i`.

```text
W(active state) = z + z⁻¹,
W(image state)  = 2^e * (z + z⁻¹),
W(active state) + W(image state) = A(1, 2^e; z).
```

* `ClanAudit.localMapP2_normalized_weight_even_arms`
  (= `ClanAudit.Audit.p2_block_normalized_weight_even_arms`) is the block identity, with
  the scalar **derived** to be `c = 2^e`; it is not postulated.
* `ClanAudit.one_le_derived_scalar` (= `ClanAudit.Audit.p2_block_derived_scalar_one_le`)
  records `1 ≤ c`, which by `Fblock_decr` is exactly what the adjacent two-hub conclusion
  needs — and, by the refutation above, what it cannot do without.
* `ClanAudit.localMapP2_normalized_weight_two_arms`
  (= `ClanAudit.Audit.p2_block_normalized_weight`) is the earlier special case `e = 0`,
  where `c = 1`.

The ingredients, all derived from the explicit graph models rather than assumed:

* `ClanAudit.Wpoly_spider`: the hub component of the active state is a connected
  bipartite spider — here with `2 + e` arms, of which exactly two have odd length — with
  colour-class imbalance exactly `p - 1`, so its normalized weight is `z^r + z^(-r)`.
  This derives both the `z^r + z^(-r)` term of `A(r,c;z)` and the identification
  `r = p - 1`.
* `ClanAudit.Wpoly_bot_two`: each cloned `K₂` created by the alternating run
  `2,0,2,0,…` contributes the orientation count `2`, which is exactly cancelled by the
  new factor `2!` in the denominator, so the `L` cloned `K₂`s contribute the factor `1`.
  This is the cancellation the request asks to be *proved*; it is proved, not assumed.
* `ClanAudit.Wpoly_armsGraph` (= `ClanAudit.Audit.even_arm_family_weight`): the reusable
  finite-family disjoint-union theorem.  After the hub is deactivated, each of the `e`
  even arms becomes a separate path component with an even number of vertices, hence
  balanced, hence of normalized weight `2`; the disjoint union of the family
  (`ClanAudit.armsGraph`) has normalized weight `2^e`.
* `ClanAudit.Wpoly_image_even_arms` (= `ClanAudit.Audit.image_state_normalized_weight`):
  the whole image state — `L` cloned `K₂`s, the untouched remainder of the arm `B` (a
  path with `M - L + 1` vertices, contributing `z + z⁻¹`), and the `e` balanced even-arm
  paths — has normalized weight `2^e (z + z⁻¹)`.

## The explicit global partition and the final degreewise theorem

This is the content of `FOLLOWUP_GLOBAL_PARTITION_20260820.md`, and it is now proved in
full, for the tree `dgraph m a n b` with **arbitrary** finite ordered families of
positive-length pendant paths at both hubs and **arbitrary** total order `N`.  Nothing is
bounded or enumerated.

### The final theorem

```text
ClanAudit.sum_normalizedTwoRowCoeff_nonneg
  (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) (N k l : ℕ)
  (hl : 1 ≤ l) (hkl : l ≤ k) (hN : k + l = N) :
  0 ≤ ∑ α ∈ mapsOfOrder m a n b N, normalizedTwoRowCoeff (dgraph m a n b) α k l
```

restated as `ClanAudit.Audit.adjacent_two_hub_degreewise_nonneg`.  It is deduced from

```text
ClanAudit.isCU_totalW : IsCU (∑ α ∈ mapsOfOrder m a n b N, Wpoly (dgraph m a n b) α)
```

through the interior formula `normalizedTwoRowCoeff_eq`.

### The partition

The blocks are the fibres of an explicit, idempotent, order-preserving representative map
`ClanAudit.rep`, computed side by side.  On one side:

* `ActiveSide s` — hub of multiplicity one, initial positive runs made of ones, at least
  two arms of odd active prefix;
* `ImgSide s` — `s = transf t` for an active `t`; `transf` is the published transformation
  with the published deterministic choices (shortest odd active prefix `idx0`, ties broken
  by arm order; second such arm `idx1`; prefix length `plen`), proved well defined by
  `canonical_spec`;
* `repSide s` is `s` if `s` is active, the (unique, by `transf_injective`) source `t` if
  `s` is an image, and `s` otherwise.

Then `rep α := Sum.elim (repSide (α ∘ inl)) (repSide (α ∘ inr))`, and
`block m a n b N β` is the fibre of `rep` over `β` inside `mapsOfOrder m a n b N`.

* **Explicitness** (`fiber_eq_image`): the block of `β` is exactly the image of
  `sideBlock (β ∘ inl) ×ˢ sideBlock (β ∘ inr)` under `Sum.elim`, where `sideBlock s` is
  `{s, transf s}` at an active side and `{s}` otherwise.  So the blocks really are the
  four-state, two-state and singleton blocks of the target.
* **Sizes** (`card_block_cases`): `4`, `2` or `1`.
* **Disjointness** (`blocks_pairwise_disjoint`) and **uniqueness of the block of a map**
  (`block_unique`) are automatic for a fibre decomposition.
* **Exhaustion** (`blocks_cover`): every map of order `N` is in the block of its
  representative, which is itself of order `N`.
* **Total order preserved** (`sum_rep`, from `sum_transf`, from
  `localMapP2_preserves_total_order`), so each block lies inside a single degree.

### The collision of `two_hub_clan_cancellation_attack_2026-08-20.md`

The collision is **separated, not fatal**.  The reason is structural: `transf` switches the
hub off, so an image can never be a source (`not_activeSide_of_imgSide`).  A map whose side
is an image therefore does not open a second, competing block; it lies in the *same* block
as its source (`collision_resolved_left`, `collision_resolved_right`).

The concrete collision of the notes is formalized in
`RequestProject/ClanAudit/Collision.lean`: three unit arms at `u`, five at `v`,

```text
collisionSource  = (hub 1; 1,1,1) ⊕ (hub 1; 1,1,1,1,1)
collisionPartner = (hub 0; 2,1,1) ⊕ (hub 1; 1,1,1,1,1)
```

both of total order `10`.  `transf_onesSide` proves that `collisionPartner` is exactly the
canonical image of the `u`-side of `collisionSource`, and `collision_block` proves that
both maps lie in one and the same **four-element** block.  What the notes called a
"single-hub source" is not a source at all.

### The weight of each block

All outside factors are carried explicitly; none is postulated.  Write
`Tu = ∏ tailW u`, `Tv = ∏ tailW v`, `r = p - 1`, `s = q - 1`, where `p`, `q` count the arms
of odd active prefix.

* Both hubs active (`blockSum_both_active`): the four weights are
  `Pw |r-s| · Tu·Tv`, `Pw r·Tu · d•(zz^s·Tv)`, `c•(zz^r·Tu) · Pw s·Tv`,
  `c•(zz^r·Tu) · d•(zz^s·Tv)`, and they sum to `Fblock r s c d · (Tu·Tv)` by the four-map
  identity `Fblock_expand`, with `r, s ≥ 1` and `c, d = 2^e ≥ 1` derived by
  `Wpoly_side_image`.  Hence `Fblock_isCU` applies — and, by the refutation above, it
  could not apply without `c, d ≥ 1`.
* Exactly one hub active (`blockSum_left_active`, `blockSum_right_active`): the block sums
  to `Ablock r c · (good factors)`, `Ablock (r+1) c · (…)` or `Ablock r (2c) · (…)`
  according to the state of the inactive side, or to `c•(zz^r · good)`, or to `0`.
  `Ablock_isCU` applies.
* Neither hub active (`blockSum_neither_active`): the singleton weight is a product of
  good factors `2^j (z+z⁻¹)^i`, or zero — in particular the two-hub component then has
  imbalance at most one.

### Classification and vanishing

`Wpoly_eq_zero_of_mult_two`: a multiplicity `≥ 2` next to a positive vertex creates a
triangle in the clan graph, so the whole normalized weight vanishes.  This handles every
non-admissible map with a positive hub, and every map with two positive adjacent hubs one
of which has multiplicity `≥ 2`.  `isGoodW_Wpoly_spider_of_not_active` then shows that
every side which is not active has a weight of the good shape.

### The `p ≥ 3` local transformation

`Wpoly_side_image` is proved for **every** active side, i.e. for every `p ≥ 2`, not only
`p = 2`; it produces the derived scalar `c = 2^e ≥ 1`, where `e` is the number of arms with
even positive active prefix.  So the missing normalized-weight theorem at `p ≥ 3` is
formalized rather than assumed.

## From degreewise nonnegativity to log-concavity

The gap recorded in earlier versions of this note — the passage from 2-Schur-positivity to
log-concavity of the independence polynomial — is now closed inside this development, and
closed by *derivation* rather than by citing an equivalence.

* `RequestProject/ClanAudit/Bridge.lean` defines independent sets, the independence counts
  `i(G, j)` and the independence polynomial `I(G; x)`, and proves the coefficient identity

  ```text
  ∑ α with ∑_v α v = k + l,  normalizedTwoRowCoeff G α k l  =  i_k i_l - i_(k+1) i_(l-1)
  ```

  directly from the clan / multicolouring definitions (`sum_normalizedTwoRowCoeff_eq`).
  Hence degreewise nonnegativity implies log-concavity of `i(G, ·)`, and of the
  coefficients of `I(G; x)` (`indepCount_logConcave_of_degreewise_nonneg`,
  `indepPoly_logConcave_of_degreewise_nonneg`).
* Applied to the adjacent two-hub tree this gives `dgraph_indepPoly_logConcave`; applied to
  the arbitrary-connector tree `connGraph t m a n b` — whose global partition, blockwise
  weights and degreewise theorem are proved in `ConnectorPartition.lean`,
  `ConnectorBlock.lean` and `ConnectorZero.lean` — it gives
  `connGraph_indepPoly_logConcave` and `connectorGraph_indepPoly_logConcave`; applied to a
  spider it gives `spider_indepPoly_logConcave`, and hence `pathGraph_indepPoly_logConcave`
  for every path.

## Trees with at most two branch vertices

Call a vertex of degree at least three a *branch vertex*.  `IsC2Model G` says that `G` is
isomorphic to a spider or to a connector tree, and every `C₂` model has log-concave
independence counts and independence polynomial
(`indepCount_logConcave_of_isC2Model`, `indepPoly_logConcave_of_isC2Model`, using the
isomorphism invariance of `IsoTransport.lean`).

The two *recognition* theorems make this intrinsic:

* `spiderIsoOfTree` — a finite tree in which every vertex other than a chosen `h` has
  degree at most two is isomorphic to `spider (G.neighborFinset h).card (hubArmLen hG h)`;
  the vertex map sends `v` to the pair (arm of `h` containing `v`, depth of `v`).
* `connIsoOfTree` — a finite tree with two distinct vertices `h1 ≠ h2` such that every
  other vertex has degree at most two is isomorphic to
  `connGraph (dist(h1,h2) - 1) …`, the connector model whose left/right arm data are the
  arms of `h1` and `h2` other than the ones facing the other hub.  The vertex map is the
  polar-coordinate map on each side, and injectivity, surjectivity and adjacency are all
  proved by the five-way classification `toConn_classify` of the vertices.

Combining them (`isC2Model_of_isTree_of_branchVerts_card_le_two`) with the `C₂`
log-concavity theorem yields, with no hypothesis other than being a finite tree with at
most two branch vertices,

```text
i(G, j) * i(G, j+2) ≤ i(G, j+1)^2      for all j
```

(`indepCount_logConcave_of_isTree_of_branchVerts_card_le_two`) and the same statement for
the coefficients of `I(G; x)`
(`indepPoly_logConcave_of_isTree_of_branchVerts_card_le_two`).

The converse is proved as well: the explicit models really are trees whose only possible
branch vertices are their hubs (`spider_isTree`, `connGraph_isTree`, via the rank/parent
criterion `isTree_of_rank_parent`, together with `spider_branchVerts_card_le_one` and
`connGraph_branchVerts_card_le_two`).  Hence, for a finite graph `G`,

```text
IsC2Model G  ↔  G.IsTree ∧ (branchVerts G).card ≤ 2
```

(`isC2Model_iff_isTree_branchVerts`): the class handled by this development is exactly the
class of finite trees with at most two branch vertices.

## Remaining scope

The refutation above is unaffected: central unimodality of `F` for arbitrary `c, d > 0` is
false, and the argument goes through only because the scalars it produces are derived to
be powers of two.  Trees with three or more branch vertices, and graphs that are not
trees, are outside the scope of this development.
