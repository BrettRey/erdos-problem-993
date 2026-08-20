# Aristotle formalization queue

Date: 2026-08-20.

This queue sends only frozen, definition-complete results. Every returned project must
be replayed locally with `lake build`, searched for proof escape hatches, checked with
`#print axioms`, and compared declaration-by-declaration with its mathematical packet.

## Active

1. **Clan normalization continuation: arbitrary even arms.** Extend the completed
   `p=2` two-arm theorem to `e` additional even-prefix arms and derive `c=2^e`.
   Packet:
   `formalization/clan_normalization_aristotle/FOLLOWUP_ARBITRARY_EVEN_ARMS_20260820.md`.
   Project `76e84c3b-9ce0-4fc6-9074-45ea38dbb576`, task
   `f2b3e260-0547-4694-b2a5-b103ab16d201`.
2. **Universal Pascal smoothing.** Formalize
   `notes/pascal_smoothing_shadow_lemma_2026-08-20.md`, including the downward-closed
   incidence inequality, the `8/9` Pascal-transform bound, and its three stated
   consequences. Project `6519b82f-e80b-4c41-b633-990b72c7e3b9`, task
   `e6301287-ecb0-4d75-8b65-b6711945038f`.
3. **Order-ideal projection and graph-transform identities.** Formalize the frozen
   poset layer in `notes/matching_bag_poset_reduction_2026-08-20.md`, especially
   equations (3), (5), (8), (9), and (10), independently of the unresolved blocked
   correction. Project `2a2d691b-91bb-4204-bb00-ad47a51861ba`.

## Next continuation

After item 1 is replayed, continue the same clan project on the exhaustive global
block partition. The target is a disjoint and exhaustive classification of clan maps
into four-state blocks, one-active two-state blocks, and individually nonnegative or
two-row-zero states, with no partner reused. This packet should build directly on the
derived `c=2^e` theorem rather than re-encode normalized weights.

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
