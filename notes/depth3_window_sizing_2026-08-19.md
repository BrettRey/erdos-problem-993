# Depth-3 window sizing: the dist-3 frontier is 531B trees, and the ladder ends there
<!-- SUMMARY: Exact Otter-DP counts of trees by matching number (alpha = n - mu, Konig) size the dist-3 frontier: rung 3 (n in [33,38], alpha in {17,18,19}) is exactly 531,192,598,181 trees inside 74.2T generated; rung 4 is 129.7T inside 34.1 quadrillion, so enumeration ends at rung 3. Filter engine built and validated; measured 2.0-4.2M trees/s/pipeline prices rung 3 at ~5,300 pipeline-hours (Modal-shaped; OUTPROC fusion cuts ~8x) · status: sized, engine validated, run not launched · updated: 2026-08-19 -->

Lane decision input for consequence 4 of
`notes/break_depth_lane_2026-08-14.md`: "sizing that is the next lane
decision." Sized exactly, without enumerating a tree.

## Method: count, don't enumerate

Trees are bipartite, so by König **alpha = n − mu** (mu = matching
number). Counting unlabeled trees by mu therefore counts them by alpha.
An Otter-dissymmetry DP does this exactly in milliseconds:

- Rooted trees carry the classical two-value state s = 1 iff some maximum
  matching leaves the root exposed. Recurrence: if every child has
  s = 0 the root stays exposable (mu = sum, s = 1); if some child has
  s = 1 the root is matched in every maximum matching (mu = sum + 1,
  s = 0).
- Unrooted counts via Otter t = r − (r² − r(x²))/2 with each pairing term
  carrying the edge-join weight mu = mu_A + mu_B + [s_A = s_B = 1].

Script: `scripts/count_trees_by_matching_20260819.py`. Results:
`results/depth3_window_sizing_20260819.json`.

**Validation (all passing):**
1. Totals sum_mu t[n][mu] match all 20 grounded orders: the A000055
   table for n ≤ 26 and the six 2026-08-14 census totals (n = 27..32),
   each of which had matched an independent `gentreeg -u` count.
2. Exhaustive brute force n = 4..14: full (n, mu) histograms from
   `gentreeg -p`, mu by exact greedy leaf-matching, alpha independently
   as deg I(T) via the tree DP, with alpha = n − mu asserted per tree.
3. The C filter engine (below) reproduces the DP's per-alpha counts
   exactly on full enumerations at n = 10, 20, 26.
4. **OEIS, discovered after the fact:** the triangle is A339829
   ("number of unlabeled trees on n vertices with independence number
   k", b-file through n = 50), and the perfect-matching diagonal is
   A005751 ("matched trees with 2n nodes"). All 990 b-file terms with
   n ≤ 44 match this DP exactly — a full independent external check by
   a different method. Candidly: the window counts could have been read
   off A339829 directly (the inward-search lesson, OEIS edition); the
   DP still earns its keep as independent verification, as the source
   of the mu-state machinery the filter engine reuses, and as the rung
   analysis. Conversely, nothing here is OEIS-new below n = 51.

## Rung 3: exactly 531,192,598,181 trees

Window = dist-3 ladder rung d = 3: alpha in {17,18,19}, n in [33,38]
(n ≤ 32 census-exhausted; alpha ≥ ceil(n/2) automatic).

| n  | window trees        | of total              | fraction | by alpha (17 / 18 / 19) |
|----|---------------------|-----------------------|----------|--------------------------|
| 33 | 90,059,916,288      | 300,628,862,480       | 3.0e-1   | 2.11B / 21.9B / 66.0B    |
| 34 | 125,405,065,758     | 823,779,631,721       | 1.5e-1   | 0.58B / 19.5B / 105.3B   |
| 35 | 136,719,894,645     | 2,262,366,343,746     | 6.0e-2   | — / 10.9B / 125.9B       |
| 36 | 108,850,200,415     | 6,226,306,037,178     | 1.7e-2   | — / 2.9B / 106.0B        |
| 37 | 56,113,092,079      | 17,169,677,490,714    | 3.3e-3   | — / — / 56.1B            |
| 38 | 14,044,428,996      | 47,436,313,524,262    | 3.0e-4   | — / — / 14.0B            |
| Σ  | **531,192,598,181** | 74,219,071,890,101    | 7.2e-3   |                          |

The alpha filter is worthless at n = 33 (30% pass) and decisive at
n = 38 (3×10⁻⁴). The window is ~3.1× the entire n ≤ 32 census (173.4B).
In matching-deficiency terms the window is delta = 2·alpha − n ≤ 38 − n:
n = 38 is exactly the perfect-matching class, n = 37 the near-perfect
class — the classes that never fail LC anywhere in the n ≤ 32 census.

## Rung 4 kills the enumeration ladder

Rung d = 4 (alpha in {20,21,22}, n in [33,44]): window
**129,653,876,445,655** trees (244× rung 3) inside
**34,055,248,568,141,887** generated (34.1 quadrillion, 459× rung 3's
74.2T). Even a perfect generator emitting only window trees would face
130T full DPs. So:

- **Rung 3 is the last enumeration-reachable point of the dist-3
  question.** Beyond it, only proof.
- **If rung 3 comes back empty:** dist ≥ 4 becomes unconditional for
  *every tree with alpha ≤ 19* (alpha ≤ 16 forces n ≤ 32 = census;
  17–19 at n ≥ 33 = this window; the run checks all breaks in every
  window tree, not just depth-3 ones). Any dist ≤ 3 break anywhere would
  then need depth ≥ 4 and alpha ≥ 20. Clean quotable boundary,
  erdosproblems-comment-shaped.
- **If it hits:** first break of depth ≥ 2 below alpha = 19, first
  dist-3 break ever, one step from the Levit-Kadrawi prefix. Either
  outcome pays.

## Engine: built and validated, not launched

`scripts/lc_census_alpha.c` extends the census engine to n ≤ 38:

- O(n) exact greedy-matching prefilter (children-before-parents order;
  alpha = n − mu), full DP + LC/unimodality checks only on passers;
- every DP'd tree cross-checks greedy alpha against deg I(T) and aborts
  on mismatch (exit 4), so a production run validates its own filter on
  every passed tree;
- PASS_ALPHA histogram lines let every shard's pass counts be summed and
  reconciled against the exact DP counts above — a completeness check
  the unfiltered census couldn't have;
- Turán products go through __uint128_t (coefficients reach
  C(38,19) ≈ 2³⁴ at n = 38, so uint64 products would overflow; the DP
  itself cannot overflow: every partial product counts a subset of
  independent sets and is bounded by a final coefficient).

Validated, all passing:

- full n = 10 and n = 20 enumerations: per-alpha pass counts equal the
  Otter DP values exactly;
- full n = 26 run at alpha_max = 19 (279,793,450 trees, 3m32s on one
  pipeline): passed = 278,830,502 — the DP's prediction exactly — with
  all seven per-alpha counts matching (alpha 13..19: 1,133,602 /
  23,091,910 / 75,054,840 / 93,224,539 / 59,027,099 / 22,023,518 /
  5,274,994); found exactly 2 LC failures, both alpha = 14 with the
  lone break at k = 13 (dist 4), and both polynomials byte-match the
  stored ground truth in `results/analysis_n26.json`;
- zero greedy-vs-DP alpha mismatches across all 278.8M DP'd trees.

## Measured throughput and the run decision

Measured on this laptop (M-series, single pipeline = gentreeg -p
producer + filter consumer, ~1.05 cores):

- gentreeg counting mode (`-u`): ~53M trees/s/core — the generation
  ceiling.
- gentreeg `-p` through a pipe: ~9.7M trees/s — formatting + pipe costs
  5× the generation itself.
- Filtered pipeline, n = 33 (39% of shard passed to DP): **1.96M
  trees/s**.
- Filtered pipeline, n = 38 (~0 pass): **4.2M trees/s** — producer-bound.

Whole-rung price at measured rates: ≈ 19M pipeline-seconds ≈ **5,300
pipeline-hours** (n = 38 alone is ~3,200: 47.4T trees at 4.2M/s). That
is ~28 days on this laptop at 8 pipelines — not a laptop job as-is.

Options, in recommended order:

1. **Modal, as-is** (precedent: the 2026-03 n = 29 sweep, 1024 workers).
   At 1024 workers and ~half-laptop per-core speed: roughly half a day
   wall. Order-of-magnitude ~10k core-hours of credits — check the
   balance before launching.
2. **Fuse the filter into gentreeg first (recommended before any
   launch).** The heavy orders are producer-bound: the tree exists in
   gentreeg's memory and we pay 5× to print and re-parse it. nauty's
   geng-family supports compile-time output/prune plugins; if
   `gentreeg.c` (source not on this machine — download nauty, check)
   exposes its OUTPROC-style hook, the greedy filter runs in-process at
   near counting-mode rate and only passers get piped out (or DP'd
   in-process). Expected ~8× on the filter-bound orders → ~700–800
   core-hours total: 3–4 days on the laptop, or ~1 hour on Modal at
   1024 workers. A day of engineering that turns a $currency-sized
   Modal job into a small one, or the laptop into a viable fallback.
3. **Deficiency-bounded structural generation** (generate only
   delta ≤ 38 − n trees, e.g. expanding the contracted matching
   skeleton): up to 3,400× fewer trees generated at n = 38, but needs
   new isomorphism-free machinery (canonical augmentation). Not worth
   it for rung 3 given option 2; it is the only route anyone could ever
   take toward rung-4-scale questions, and rung 4 is out of reach even
   so.

## Provenance

Session 2026-08-19 (Fable, token-burn directive, autonomous lane choice:
the named "sizing" next step from consequence 4 of the 2026-08-14 lane).
Exact counts from `scripts/count_trees_by_matching_20260819.py`; engine
`scripts/lc_census_alpha.c`; timings above measured this session on this
machine. No manuscript or public surface touched. Run not launched:
Modal spend and the fuse-first choice are Brett's call.
