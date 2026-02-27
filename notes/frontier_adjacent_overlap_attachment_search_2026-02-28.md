# Frontier Adjacent Overlap Attachment Search (2026-02-28)

## Context
Goal was to lift known adjacent same-`(m,lambda)` pair to same-`(m,lambda,rho)` by bounded gadget attachments, while preserving odd/even adjacency.

Pair used:
- A: `O??????_A?C?E?@_WG@j?` (`N=13`)
- B: `P????A?OD?E?E?B??o?E?OO?` (`N=14`)

## In-repo script
- `scripts/frontier_adjacent_overlap_attachment_search.py`

Search setup in this local reproduction:
- rooted gadgets: all unique `(size,F,G)` rooted d_leaf<=1 types with `2 <= |V| <= 8`
- `gadgets_count = 57`
- attachment patterns: multiset attachments with
  - `max_gadgets = 4`
  - `max_added_size = 24`
  - even added size only
- mode gate: `m >= 4`

Command:
- `nice -n 18 python3 scripts/frontier_adjacent_overlap_attachment_search.py --gadget-max-n 8 --max-gadgets 4 --max-added-size 24 --m-min 4 --out results/frontier_adjacent_overlap_attachment_search_k8_g4_s24_local.json`

## Artifact
- `results/frontier_adjacent_overlap_attachment_search_k8_g4_s24_local.json`

## Result
- `split_found = false`
- `collision_count = 0`
- `common_size_m_lambda_keys = [{deltaN:0, m:5, lambda:[7,9]}]` only

Closest overlap metric:
- key `(deltaN,m,lambda) = (0,5,7/9)`
- `|rhoA - rhoB| = 10519070037889952 / 59484897518860809`

This matches C's reported strongest-gap fraction exactly.

## Note on pattern count
- Local multiset-pattern count is `44444`.
- C reported `68196` patterns under their generation method.
- Despite pattern-generation differences, both analyses agree on the decisive output:
  no nontrivial common `(deltaN,m,lambda)` outcome beyond the base point, and no lifted same-`(m,lambda,rho)` witness in tested bounds.
