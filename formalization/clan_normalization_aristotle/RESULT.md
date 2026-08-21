# Result note

## Grade

**`REFUTED`** — with the full **`NORMALIZATION_KERNEL`** also completed, and the
`p = 2` local package (**`P2_LOCAL`**) now completed, including the arbitrary-arm case of
the weight theorem.

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
  are not formalized.

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

## First exact unfilled mathematical obligation

The first obligation that is stated in the request and is *not* discharged here is items
2 and 3 of the adjacent two-hub target: the **exhaustive global clan-map block
partition** — that the global clan maps partition disjointly and exhaustively into
four-state blocks, one-active two-state blocks, and individually nonnegative or
two-row-zero maps, and that no partner map is used twice.  (Injectivity of the local
pairing maps, which is the "no partner map is reused at a single hub" half of item 3, *is*
proved: `localMapP2_injective`.)

Note that this remaining obligation concerns only the global bookkeeping: the local
`p = 2` package it would consume — the weight of each block, with its derived scalar
`c = 2^e ≥ 1` — is complete.
