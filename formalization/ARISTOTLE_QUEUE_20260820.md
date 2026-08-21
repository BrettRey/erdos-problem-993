# Aristotle formalization queue

Date: 2026-08-20.

This queue sends only frozen, definition-complete results. Every returned project must
be replayed locally with `lake build`, searched for proof escape hatches, checked with
`#print axioms`, and compared declaration-by-declaration with its mathematical packet.

## Active

1. **Exhaustive adjacent-two-hub global partition.** Define the full degreewise clan-map
   space and prove disjointness, exhaustion, collision-freedom, exact block weights,
   and final normalized two-row nonnegativity. Packet:
   `formalization/clan_normalization_aristotle/FOLLOWUP_GLOBAL_PARTITION_20260820.md`.
   Project `76e84c3b-9ce0-4fc6-9074-45ea38dbb576`.
2. **Poset/Pascal bridge.** Transport the universal smoothing theorem to the arbitrary
   finite ground type of the poset project and prove the extendable-profile
   inequalities (7)/(7a) end to end. Packet:
   `formalization/matching_bag_poset_aristotle/FOLLOWUP_PASCAL_BRIDGE_20260820.md`.
   Project `2a2d691b-91bb-4204-bb00-ad47a51861ba`.
3. **Even-connector Laurent block.** Prove central unimodality of the uniform
   `c,d≥1` even-connector block without bounded computation. Packet:
   `formalization/c2_even_connector_aristotle_input_20260820.md`. Project
   `f9cd8c48-bd2e-4aa8-bc67-3a9316657b61`.

## Completed and locally replayed

1. **Clan normalization continuation: arbitrary even arms.** The explicit finite-family
   component model proves `W(image)=2^e(z+z⁻¹)`, the full `p=2` block identity
   `A(1,2^e;z)`, and `2^e≥1`. Project
   `76e84c3b-9ce0-4fc6-9074-45ea38dbb576`, task
   `f2b3e260-0547-4694-b2a5-b103ab16d201`; local build: 8,038 jobs.
2. **Universal Pascal smoothing.** Formalized
   `notes/pascal_smoothing_shadow_lemma_2026-08-20.md`, including the downward-closed
   incidence inequality, the `8/9` Pascal-transform bound, and its three stated
   consequences. Project `6519b82f-e80b-4c41-b633-990b72c7e3b9`, task
   `e6301287-ecb0-4d75-8b65-b6711945038f`.
3. **Order-ideal projection and graph-transform identities.** Formalize the frozen
   poset layer in `notes/matching_bag_poset_reduction_2026-08-20.md`, especially
   equations (3), (5), (8), (9), and (10), independently of the unresolved blocked
   correction. Project `2a2d691b-91bb-4204-bb00-ad47a51861ba`.

## Ready after the clan partition

1. The even-connector Laurent/binomial theorem in
   `notes/c2_connector_clan_reduction_2026-08-20.md`.
2. The complete adjacent and arbitrary-connector `C₂` theorem, assembled only after
   the global partition and connector algebra have both replayed.
3. The tree-to-matching-bag contraction and maximum-cover-to-forest-poset bridge in
   `notes/matching_bag_csp_attack_2026-08-20.md`; the abstract order-ideal layer is
   already active above.

## Not ready

- the unresolved blocked-correction terms `b₂,b₃,b₄`;
- broad LAW-W or head-strip conjectures;
- Conjecture A or the full Erdős 993 conjecture;
- any quarantined broad STP2 closure statement.
