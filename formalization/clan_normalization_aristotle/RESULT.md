# Result note

## Grade

**`REFUTED`** — with the full **`NORMALIZATION_KERNEL`** also completed, and the
`p = 2` local package (**`P2_LOCAL`**) completed except for the arbitrary-arm case of the
weight theorem.

* Items 1–3 of the mandatory normalization kernel are proved in full, from the
  definitions, with no axioms and no `sorry`.
* Of the three `p = 2` theorems asked for, `localMapP2_injective` and
  `localMapP2_preserves_total_order` are proved in full;
  `localMapP2_has_claimed_normalized_two_row_weight` is proved for a hub with exactly two
  active arms (see "Partial item" below).
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

## Partial item: the normalized weight of the `p = 2` block

`localMapP2_has_claimed_normalized_two_row_weight` is proved in the following case, and
is labelled partial as the request requires.

For a hub whose active arms are exactly two, with odd positive prefixes of lengths
`L ≤ M` (so `p = 2` and `r = p - 1 = 1`):

```text
W(active state) + W(image state) = A(1, 1; z),
```

`ClanAudit.localMapP2_normalized_weight_two_arms`
(= `ClanAudit.Audit.p2_block_normalized_weight`).  In particular the scalar `c` is here
**derived** and equals `1` — it is not postulated.  The two ingredients are:

* `ClanAudit.Wpoly_spider`: the hub component of the active state is a connected
  bipartite spider with colour-class imbalance exactly `p - 1`, so its normalized weight
  is `z^r + z^(-r)`.  This derives both the `z^r + z^(-r)` term of `A(r,c;z)` and the
  identification `r = p - 1`.
* `ClanAudit.Wpoly_bot_two`: each cloned `K₂` created by the alternating run
  `2,0,2,0,…` contributes the orientation count `2`, which is exactly cancelled by the
  new factor `2!` in the denominator, so `L` cloned `K₂`s contribute the factor `1`.
  This is the cancellation the request asks to be *proved*; it is proved, not assumed.

## First exact unfilled mathematical obligation

The first obligation that is stated in the request and is *not* discharged here is:

> the normalized weight of the `p = 2` local block at a hub carrying, besides the two
> arms with odd positive prefixes, an arbitrary finite number `e ≥ 1` of further active
> arms with **even** positive prefixes.

Concretely, what is missing is the decomposition of the image clan graph into its
connected components in that generality: after the hub is deactivated, each even arm
becomes a separate path component of even order, contributing the balanced factor `2`,
so the same computation gives `c = 2^e` (hence `c ≥ 1`, which by `Fblock_decr` is exactly
what the adjacent two-hub conclusion needs).  Formalizing it requires an `e`-fold
disjoint-union statement (`imbalanceGF` of a finite family of components), whereas the
project currently provides the binary version `imbalanceGF_sum` / `Wpoly_sum_elim` and
uses it for the two-arm case.

Beyond that, items 2 and 3 of the adjacent two-hub target — that the global clan maps
partition disjointly and exhaustively into four-state blocks, one-active two-state
blocks, and individually nonnegative or two-row-zero maps, and that no partner map is
used twice — are not formalized at all.  (Injectivity of the local pairing maps, which
is the "no partner map is reused at a single hub" half of item 3, *is* proved:
`localMapP2_injective`.)
