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

* Build: `lake build` (Lean 4.28.0, Mathlib v4.28.0).
* The project contains **no** `sorry`, `admit`, `axiom`, or `implemented_by`.
  All theorems depend only on `propext`, `Classical.choice`, `Quot.sound`.
* Graded result note: [`RESULT.md`](RESULT.md).
* Follow-up target `FOLLOWUP_GLOBAL_PARTITION_20260820.md` (the explicit global partition
  and the final degreewise theorem) is **proved**; see the last traceability table below.
* Independent local replay on 2026-08-21: `lake build` completed 8,051 jobs. Direct
  `#print axioms` checks on the final degreewise theorem, partition declarations, and
  collision theorem reported only `propext`, `Classical.choice`, and `Quot.sound`.
* Final continuation: `FOLLOWUP_FULL_C2_ASSEMBLY_20260821.md`, Aristotle task
  `147581ec-6d3a-423f-b359-85cababab583`.

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
| `RequestProject/ClanAudit/EvenArms.lean` | the normalized weight of a `p = 2` block with two odd arms and an arbitrary finite family of further even arms, deriving the scalar `c = 2^e` |
| `RequestProject/ClanAudit/PathClan.lean` | the clan weight of an arbitrary path map: always of the good shape `2^j (z+z⁻¹)^i` or zero |
| `RequestProject/ClanAudit/Split.lean`, `Arms.lean` | the vertex-splitting toolkit and the arm/spider cut equivalences |
| `RequestProject/ClanAudit/DoubleSpider.lean` | the adjacent two-hub tree `dgraph`, its connectivity, bipartition, and the exact imbalance law `\|p - q\|` |
| `RequestProject/ClanAudit/BlockCU.lean` | central unimodality of `A(r,c)` and `F(r,s;c,d)` on the derived range `c, d ≥ 1` |
| `RequestProject/ClanAudit/Prefix.lean` | active prefixes `pref`, prefix/tail weights, the exact weight of one side of the tree |
| `RequestProject/ClanAudit/SideTransform.lean` | the canonical transformation `transf` of a whole side with the published shortest-odd-prefix choice; injectivity and total-order preservation |
| `RequestProject/ClanAudit/ImageWeight.lean` | the exact weight of a transformed side at arbitrary `p ≥ 2`, with the scalar `c = 2^e ≥ 1` derived |
| `RequestProject/ClanAudit/TwoHubWeight.lean` | the exact normalized weight of a clan map of the two-hub tree |
| `RequestProject/ClanAudit/Vanishing.lean` | triangle/nonbipartite vanishing; every non-active side has a good-shape weight |
| `RequestProject/ClanAudit/BlockWeight.lean` | the exact normalized Laurent sum of each block of the partition |
| `RequestProject/ClanAudit/GlobalPartition.lean` | the explicit global partition of the degree-`N` map space and the final degreewise theorem |
| `RequestProject/ClanAudit/Collision.lean` | the explicit collision of the proof notes, formalized and resolved |
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
| disjoint union of the `e` even-arm components of the image state | `ClanAudit.ArmsV`, `ClanAudit.armsGraph` | `EvenArms.lean` |
| the adjacent two-hub tree `D(a,b)` | `ClanAudit.DV`, `ClanAudit.dgraph` | `DoubleSpider.lean` |
| active prefix of an arm; number of odd / even-positive prefixes | `ClanAudit.pref`, `ClanAudit.pNum`, `ClanAudit.eNum` | `Prefix.lean` |
| prefix and tail weights of an arm | `ClanAudit.prefW`, `ClanAudit.tailW` | `Prefix.lean` |
| the hub-active stratum of one side | `ClanAudit.AdmRun`, `ClanAudit.ActiveSide` | `SideTransform.lean` |
| the published deterministic choices (shortest odd prefix, ties by arm order) | `ClanAudit.armKey`, `ClanAudit.idx0`, `ClanAudit.idx1`, `ClanAudit.plen` | `SideTransform.lean` |
| the canonical local transformation of a whole side | `ClanAudit.transfSide`, `ClanAudit.transf` | `SideTransform.lean` |
| the degree-`N` clan-map space | `ClanAudit.mapsOfOrder` | `GlobalPartition.lean` |
| block representative and blocks of the global partition | `ClanAudit.ImgSide`, `ClanAudit.repSide`, `ClanAudit.sideBlock`, `ClanAudit.rep`, `ClanAudit.block` | `GlobalPartition.lean` |
| total normalized weight of a degree | `ClanAudit.totalW` | `GlobalPartition.lean` |

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
| `localMapP2_has_claimed_normalized_two_row_weight`, two active arms | `ClanAudit.Audit.p2_block_normalized_weight` (`localMapP2_normalized_weight_two_arms`) | proved (`c = 1` derived) |
| `localMapP2_has_claimed_normalized_two_row_weight`, two odd arms and `e` further even arms | `ClanAudit.Audit.p2_block_normalized_weight_even_arms` (`localMapP2_normalized_weight_even_arms`) | proved (`c = 2^e` derived) |
| finite family of even-arm components: `W = 2^e` | `ClanAudit.Audit.even_arm_family_weight` (`Wpoly_armsGraph`) | proved |
| normalized weight of the whole image state, `2^e (z + z⁻¹)` | `ClanAudit.Audit.image_state_normalized_weight` (`Wpoly_image_even_arms`) | proved |
| the derived scalar satisfies `1 ≤ c` | `ClanAudit.Audit.p2_block_derived_scalar_one_le` (`one_le_derived_scalar`) | proved |
| Adjacent target 1: local pairing maps injective, including `p = 2` | `ClanAudit.Audit.localMapP2_injective'` | proved |
| Adjacent target 2, 3: global exhaustive disjoint block partition | `ClanAudit.Audit.global_partition_blocks_explicit`, `global_partition_block_sizes`, `global_partition_disjoint`, `global_partition_exhaustive`, `global_partition_preserves_total_order` | proved |
| Adjacent target 4: each block has the claimed normalized contribution | `ClanAudit.Audit.blockwise_central_unimodality`, `four_state_block_weight`, `p2_block_normalized_weight_even_arms` | proved (all blocks of the partition, arbitrary arms) |
| Adjacent target 5: the four-map identity | `ClanAudit.Audit.four_map_identity` (`Fblock_expand`) | proved |
| Adjacent target 6: coefficients of `F` centrally unimodal for `c, d > 0` | `ClanAudit.Audit.central_unimodality_refuted` (`Fblock_not_centrally_unimodal`), `ClanAudit.Audit.central_unimodality_witness` | **refuted** |
| repaired form of target 6 (`1 ≤ c`, `1 ≤ d`) | `ClanAudit.Audit.central_unimodality_of_one_le` (`Fblock_decr`) | proved |
| the same failure already for one hub, `A(3, 1/4; z)` | `ClanAudit.Ablock_not_centrally_unimodal` | proved |

### Follow-up target: the explicit global partition (`FOLLOWUP_GLOBAL_PARTITION_20260820.md`)

| Required checkpoint | Lean declaration | Status |
| --- | --- | --- |
| 1. Global model: tree, arm lengths, clan maps, total order, active prefixes, deterministic local choices, degree-`N` map space | `ClanAudit.dgraph`, `ClanAudit.pref`, `ClanAudit.ActiveSide`, `ClanAudit.idx0`/`idx1`/`plen`, `ClanAudit.Audit.canonical_choices_well_defined`, `ClanAudit.Audit.degree_map_space_membership` | proved |
| 2. Classification: triangle / nonbipartite cases are two-row-zero; every non-active side has a good-shape weight | `ClanAudit.Audit.two_row_zero_of_clan_triangle`, `ClanAudit.Audit.inactive_side_good_shape` | proved |
| 3. Collision audit: single-hub and double-hub strata cannot reuse a partner; the recorded collision is separated | `ClanAudit.Audit.source_and_image_strata_disjoint`, `single_and_double_hub_strata_do_not_compete`, `notes_collision_in_one_block` | proved (separated, not refuted) |
| 4. Partition: disjointness, exhaustion, preservation of total order | `ClanAudit.Audit.global_partition_blocks_explicit`, `global_partition_block_sizes`, `global_partition_disjoint`, `global_partition_exhaustive`, `global_partition_preserves_total_order` | proved |
| 5. Weight: exact normalized Laurent sum for each block, with all outside factors derived | `ClanAudit.Audit.transformed_side_weight`, `four_state_block_weight`, `blockwise_central_unimodality` | proved |
| the published `p ≥ 3` local transformation and its missing normalized-weight theorem | `ClanAudit.Audit.transformed_side_weight` (`Wpoly_side_image`, any `p ≥ 2`) | proved |
| local pairing map injective / order preserving on a whole side | `ClanAudit.transf_injective`, `ClanAudit.sum_transf` | proved |
| 6. Final sum: degreewise normalized two-row nonnegativity | `ClanAudit.Audit.adjacent_two_hub_degreewise_nonneg`, `ClanAudit.Audit.total_degree_weight_centrally_unimodal` | proved |
