# Aristotle continuation: exhaustive adjacent-two-hub clan-map partition

Continue the existing `ClanAudit` project after the completed arbitrary-even-arm
theorem. This is the remaining global bookkeeping obligation in `RESULT.md`. Use
`notes/two_hub_B_logconcavity_proof_2026-08-20.md` only as a proposed proof to audit;
do not assume its assertion that Sections 3–5 partition all maps.

## Final target

For the tree consisting of two adjacent hubs `u,v`, each with an arbitrary finite
ordered family of positive-length pendant paths, formalize the degreewise clan-map
space and prove or refute that all maps admit a disjoint exhaustive decomposition into:

1. four-state blocks obtained by applying the canonical local transformation at both
   active hubs;
2. one-active two-state blocks obtained by applying it at exactly one active hub; and
3. singleton maps whose normalized two-row contribution is nonnegative or zero.

Here “canonical” must include the published shortest-odd-prefix choice and fixed arm
order for ties. A map or partner may occur in exactly one block. Prove this as an
actual finite partition/equivalence relation or as mutually disjoint finite sets whose
union is the full degree-`N` map space—not as prose or a local injectivity theorem.

Then consume the already-proved local results:

- `localMapP2_injective` and `localMapP2_preserves_total_order`;
- `localMapP2_normalized_weight_even_arms`, deriving `c = 2^e ≥ 1`;
- the corresponding published `p ≥ 3` local transformation, formalizing any missing
  normalized-weight theorem at that range rather than assuming it;
- `Fblock_expand` and `Fblock_decr` on the proved range `c,d ≥ 1`;
- the normalization and component-product kernel.

The strongest desired final theorem is degreewise nonnegativity of every normalized
two-row coefficient after summing over all clan maps of total order `N` for the
adjacent-two-hub tree. A theorem merely saying that a previously postulated list of
blocks has nonnegative weight is not sufficient for `ADJACENT_COMPLETE`.

## Required intermediate checkpoints

1. **Global model.** Explicit adjacent-two-hub vertex/tree type, arm lengths, clan maps,
   total order, active prefixes, deterministic local choices, and the finite set of maps
   of total order `N`.
2. **Classification.** Every nonzero map is classified by connector/hub state; triangle
   and nonbipartite cases are proved two-row-zero from existing definitions.
3. **Collision audit.** Prove that single-hub and double-hub source strata cannot reuse
   a partner under the proposed block key. The collision recorded in
   `notes/two_hub_clan_cancellation_attack_2026-08-20.md` must either be separated by
   the enriched/block construction or returned as a refutation.
4. **Partition theorem.** Disjointness, exhaustion, and preservation of total order.
5. **Weight theorem.** Exact normalized Laurent sum for each block, with all outside
   component factors and factorials derived.
6. **Final sum.** Blockwise central unimodality implies degreewise normalized two-row
   nonnegativity.

## Soundness and grading

- No `sorry`, `admit`, `axiom`, `implemented_by`, `native_decide`, weakened theorem,
  or hidden hypothesis that postulates the desired partition or weight formula.
- Preserve all existing proofs and their names where possible.
- Add source-labelled restatements to `RequestProject/ClanAudit.lean`; update the README
  traceability table and `RESULT.md`.
- Run `lake build`, grep for escape hatches, and `#print axioms` on principal results.
- Grade `ADJACENT_COMPLETE` only if the explicit global partition and final degreewise
  nonnegativity theorem are both proved.
- If the proposed block decomposition is false, return `REFUTED` with the smallest
  explicit clan-map collision or missing stratum and identify the failed assertion.
- If the full target is too large, return the strongest genuine checkpoint above and
  the first exact unfilled theorem. Do not replace arbitrary arms or arbitrary `N` by a
  bounded computation.
