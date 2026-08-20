This project was edited by [Aristotle](https://aristotle.harmonic.fun).

To cite Aristotle:
- Tag @Aristotle-Harmonic on GitHub PRs/issues
- Add as co-author to commits:
```
Co-authored-by: Aristotle (Harmonic) <aristotle-harmonic@harmonic.fun>
```

# Normalized clan cancellation for adjacent two-hub trees — Lean 4 audit

An independent Lean 4 / Mathlib formalization and audit of the normalized clan-graph
bookkeeping described in `clan_normalization_aristotle_input_20260820.md`
(Li–Li–Yang–Zhang, *Log-concavity of independence polynomials of some families of
graphs*, arXiv:2501.04245, Section 2 and Theorem 3.1).

* Aristotle project: `76e84c3b-9ce0-4fc6-9074-45ea38dbb576`.
* Independent local replay: `lake build` succeeded on 2026-08-20 (8,037 jobs).
* Local `#print axioms` checks on the twelve principal declarations found only
  `propext`, `Classical.choice`, and `Quot.sound`.

* Build: `lake build` (Lean 4.28.0, Mathlib v4.28.0).
* The project contains **no** `sorry`, `admit`, `axiom`, or `implemented_by`.
  All theorems depend only on `propext`, `Classical.choice`, `Quot.sound`.
* Graded result note: [`RESULT.md`](RESULT.md).

## Files

| File | Contents |
| --- | --- |
| `RequestProject/ClanAudit/Laurent.lean` | Laurent coefficient sequences over `ℚ`, central unimodality (`IsCU`) and its closure properties |
| `RequestProject/ClanAudit/Blocks.lean` | the proposed blocks `A(r,c;z)`, `F(z)`, the four-map identity, the refutation, and the repaired statement |
| `RequestProject/ClanAudit/Stable.lean` | stable two-block partitions, `stableCount`, the imbalance generating function, connected components, multiplicativity |
| `RequestProject/ClanAudit/Clan.lean` | clan graphs and the clan normalization identity |
| `RequestProject/ClanAudit/Normalized.lean` | `clanFactorial`, `Wpoly`, `normalizedTwoRowCoeff` and the interior formula |
| `RequestProject/ClanAudit/LocalMap.lean` | the published spider transformation at `p = 2` |
| `RequestProject/ClanAudit/Weight.lean` | isomorphism invariance, multiplicativity of `W`, the cloned `K₂` cancellation |
| `RequestProject/ClanAudit/Spider.lean` | the hub component of an active spider state, its bipartition and imbalance `r = p - 1` |
| `RequestProject/ClanAudit/P2Weight.lean` | the normalized weight of a `p = 2` block with two active arms |
| `RequestProject/ClanAudit.lean` | source-labelled summary of every audited item |

## Traceability table

Every entry of the request is listed, with the Lean declaration that formalizes it.  The
`ClanAudit.Audit.*` names are the source-labelled restatements in `RequestProject/ClanAudit.lean`;
the third column gives the file where the mathematics is proved.

### Definitions

| Mathematics (request) | Lean declaration | File |
| --- | --- | --- |
| clan graph `clan G alpha` | `ClanAudit.clan` | `Clan.lean` |
| proper multicoloring of type `alpha` | `ClanAudit.IsMulticoloring` | `Clan.lean` |
| monomial of a coloring / multicoloring | `ClanAudit.monomialOfColoring`, `ClanAudit.monomialOfMulti` | `Clan.lean` |
| stable (independent) set, stable two-block partition | `ClanAudit.IsIndep`, `ClanAudit.IsStablePair` | `Stable.lean` |
| `stableCount H k l` (semi-ordered two-block count; empty blocks allowed, ordered pairs) | `ClanAudit.stableCount` | `Stable.lean` |
| two-row coefficient `[s_(k,l)] Xchrom H` (taken as a definition, as the request permits) | `ClanAudit.twoRowCoeff` | `Stable.lean` |
| imbalance generating function `∑ z^(|A|-|B|)` | `ClanAudit.imbalanceGF` | `Stable.lean` |
| normalization denominator `∏_v (alpha v)!` | `ClanAudit.clanFactorial` | `Normalized.lean` |
| candidate normalized imbalance polynomial `W(G, alpha; z)` | `ClanAudit.Wpoly` | `Normalized.lean` |
| normalized two-row coefficient | `ClanAudit.normalizedTwoRowCoeff` | `Normalized.lean` |
| central unimodality of a Laurent coefficient sequence | `ClanAudit.IsCU` | `Laurent.lean` |
| `z^t + z^(-t)`, `(z + z⁻¹)^r` | `ClanAudit.Pw`, `ClanAudit.zz` | `Blocks.lean` |
| the blocks `A(r,c;z)` and `F(z)` | `ClanAudit.Ablock`, `ClanAudit.Fblock` | `Blocks.lean` |
| local state at a hub (clan-map data read by the transformation) | `ClanAudit.HubState`, `ClanAudit.Admissible` | `LocalMap.lean` |
| the `p = 2` local transformation | `ClanAudit.localMapP2` | `LocalMap.lean` |
| hub component of an active spider state | `ClanAudit.spider` | `Spider.lean` |

### Theorems

| Mathematics (request) | Lean declaration | Status |
| --- | --- | --- |
| Kernel 1: clan normalization — fibres of size `∏ (alpha v)!` | `ClanAudit.Audit.clan_normalization_fibre_card` (`clan_fiber_card`) | proved |
| Kernel 1: clan normalization — equal monomials | `ClanAudit.Audit.clan_normalization_monomial` (`monomial_toMulti`) | proved |
| LLYZ Prop. 2.4 through the imbalance polynomial | `ClanAudit.Audit.two_row_coeff_eq_imbalance_coeff` (`twoRowCoeff_eq_coeff`) | proved |
| Kernel 2: weighted bipartition formula for a component (`z^delta + z^(-delta)`; balanced ⇒ `2`; isolated vertex ⇒ `z + z⁻¹`) | `ClanAudit.Audit.component_orientation_weight` (`imbalanceGF_connected`) | proved |
| Kernel 2: nonbipartite ⇒ weight `0` | `ClanAudit.Audit.nonbipartite_weight_zero` (`imbalanceGF_eq_zero_of_not_colorable`) | proved |
| Kernel 2: `normalizedTwoRowCoeff alpha k l = coeff W (k-l) - coeff W (k-l+2)` | `ClanAudit.Audit.normalized_two_row_coeff_formula` (`normalizedTwoRowCoeff_eq`) | proved |
| Kernel 3: multiplicativity over disjoint unions | `ClanAudit.Audit.imbalance_multiplicative` (`imbalanceGF_sum`) | proved |
| Kernel 3: product form for two components | `ClanAudit.Audit.imbalance_product_of_components` (`imbalanceGF_sum_connected`) | proved |
| Kernel 3: multiplicativity of the *normalized* `W` | `ClanAudit.Wpoly_sum_elim` | proved |
| triangle vanishing (neighbour of multiplicity `≥ 2`) | `ClanAudit.Audit.two_row_zero_of_neighbour_multiplicity_two` | proved |
| `localMapP2_preserves_total_order` | `ClanAudit.Audit.localMapP2_preserves_total_order'` (`localMapP2_preserves_total_order`) | proved |
| `localMapP2_injective` | `ClanAudit.Audit.localMapP2_injective'` (`localMapP2_injective`) | proved |
| cloned `K₂` orientation counts cancel the new factors of `2!` | `ClanAudit.Audit.cloned_K2_cancellation` (`Wpoly_bot_two`, `Wpoly_add_cloned_K2`) | proved |
| the `z^r + z^(-r)` term of `A(r,c;z)`, with `r = p - 1` derived | `ClanAudit.Audit.active_hub_component_weight` (`Wpoly_spider`) | proved |
| `localMapP2_has_claimed_normalized_two_row_weight` | `ClanAudit.Audit.p2_block_normalized_weight` (`localMapP2_normalized_weight_two_arms`) | proved **for a hub with exactly two active arms** (`c = 1` derived); general arms open |
| Adjacent target 1: local pairing maps injective, including `p = 2` | `ClanAudit.Audit.localMapP2_injective'` | proved |
| Adjacent target 2, 3: global exhaustive disjoint block partition | — | not formalized |
| Adjacent target 4: each block has the claimed normalized contribution | `ClanAudit.Audit.p2_block_normalized_weight` | partial (two-arm `p = 2` block) |
| Adjacent target 5: the four-map identity | `ClanAudit.Audit.four_map_identity` (`Fblock_expand`) | proved |
| Adjacent target 6: coefficients of `F` centrally unimodal for `c, d > 0` | `ClanAudit.Audit.central_unimodality_refuted` (`Fblock_not_centrally_unimodal`), `ClanAudit.Audit.central_unimodality_witness` | **refuted** |
| repaired form of target 6 (`1 ≤ c`, `1 ≤ d`) | `ClanAudit.Audit.central_unimodality_of_one_le` (`Fblock_decr`) | proved |
| the same failure already for one hub, `A(3, 1/4; z)` | `ClanAudit.Ablock_not_centrally_unimodal` | proved |
