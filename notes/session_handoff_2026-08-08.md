# Erdős #993 research handoff — 2026-08-08

## Bottom line

Erdős Problem #993 remains unresolved. This session produced neither a
counterexample to #993 nor a proof of #993. It did, however, prove
log-concavity for every double star whose hub connector is subdivided an
arbitrary number of times, refute one plausible auxiliary induction statement
exactly, correct a misleading search score, check the saved ratio-order paper
and its corrigendum, and close several bounded search lanes with clean negative
evidence. The connector theorem is a strict subfamily result; no novelty claim
is made without a fuller prior-art search.

## Exact result: the EADS auxiliary statement is false

EADS was the proposed statement that every tree `T` has a vertex `v` such
that the modal intervals of `I(T-v)` and `x I(T-N[v])` are at distance at
most one. The tree in `results/eads_counterexample_n60_20260808.json` has
order 60 and distance exactly two at every vertex (60/60 vertices). This is
an exact counterexample to EADS, **not** to Erdős #993: its full independence
polynomial is unimodal and log-concave.

Independent verification compares the directed-message tree DP with direct
subset summation over the smaller side of each induced bipartite forest:

```text
python3 verify_eads_counterexample_20260808.py
order = 60
minimum split distance = 2
distance counts = {2: 60}
orientation counts = {left: 11, right: 49}
full polynomial unimodal = true
full polynomial log-concave = true
near-miss ratio = 0.890136102957537 at index 25
```

Certificate SHA-256:

```text
2693c5f0ab9062fd28e60f940079abb3753b8bb7e5fcec7e60263b35ec2126bd
```

The project test suite also passed: 55 tests.

## Search-score correction

`results/rooted_product_valley_20260808.json` records a ratio
`0.9999999769395256` at `k=84`. Exact replay with the corrected `near_miss`
definition returns no rebound: the ratios decrease monotonically after the
first descent, and `k=84` is the first descent itself. The stored number is a
historical scorer/plateau artifact, not a genuine post-descent recovery.
The same exact replay removes the apparent signals in
`rooted_product_valley_smoke` and `rooted_product_later_slope_smoke`.

Operational rule: rank only a rise occurring *after* a genuine earlier
descent. Do not use proximity of the first descent ratio to 1 as evidence of
a valley.

Attaching a leaf-rooted `m`-star to every vertex of the strongest order-28
log-concavity-failure tree was also screened numerically for `1 <= m <= 500`.
No candidate appeared. At `m=500`, the first descent ratio was
`0.999992110934132`, followed by `0.9997066959325532`: flattening, not a
rebound. This is numerical negative evidence only.

## Completed and interrupted negative searches

- Exact rooted-attachment formula search: no hit among all distinct rooted
  tree states through rooted order 10 with `t <= 100` (1,422,000 cases), nor
  orders 11–12 with `t <= 50` (3,832,200 cases).
- Random spanning-tree descendants of the Bhattacharyya–Kahn graph: 10,000
  sampled spanning-tree cores, no counterexample. The exploratory scorer's
  best value was `0.9676983209794777`; it is not a certificate.
- Full pattern-family screen: interrupted at 1,530,000 parameter tuples,
  with no candidate.
- Seeded EADS leaf/local-constraint optimizers were intentionally stopped.
  The leaf run reached generation 110 / 63,350 evaluations; best diagnostics
  remained `balanced=[1,25,1,24]`, `best_leaf=[2,0,18]`,
  `best_edge=[1,24,115]`. The local-constraint run reached generation 340 /
  197,300 evaluations; best objective remained 1 at order 162. These partial
  runs wrote no final JSON reports.

These are bounded computational exclusions, not theorems about unsearched
orders or families.

## Saved-paper audit

The supplied PDF
`/Users/brettreynolds/pdf-inbox/1-s2.0-S0012365X18303042-mainext.pdf` is a
15-page bundle containing the 2019 article and its 2023 corrigendum. All
article pages and corrigendum pages were rendered and inspected. The
corrigendum adds the support condition `delta(p) <= partial(q) + 1` to
Definition 2.3; the paper states that the original results remain valid with
the corrections. The repository's strict-core experiment had already used
the corrected partial ratio order. Its simple synchronization invariants
fail, so this paper sharpens the algebraic language but does not close #993.

Source hashes:

```text
PDF a106b4c612d4d59bf3bef93412c9ec300c4fe93239b55599639129d796677952
MD  d342d602ddde449b891c62651dec76375901c0c672a9405c4d272e746d9c3723
```

## Machine and process state

The Mac had not restarted: system uptime exceeded 27 days when checked. The
three exploratory searches and all `caffeinate` processes were stopped at
the user's request, and the process table was rechecked afterward. Research
browser tabs were closed. No source PDF was changed.

## Sensible restart point

Do not revive EADS as stated, and do not optimize the rooted-product
first-descent plateau artifact. A bounded continuation after the initial
wrap proved an explicit effective leaf threshold for every fixed multi-arm
vector, and four tested vectors have exact all-count certificates. This is
only a quantitative refinement: every such tree is a spider, and the known
Li--Li--Yang--Zhang theorem already gives strong log-concavity for all
spiders. See
`notes/fixed_arm_effective_unimodality_2026-08-08.md` and
`results/fixed_arm_unimodality_threshold_20260808.json`.

If the underlying problem is reopened, do not attack “no cross-gap rebound”
for arbitrary tree-DP states as if it were a smaller lemma; that is #993 in
new language. The next genuine frontier is multiple interacting hubs. The
most direct corrected-ratio-order proof and the synchronization proof via
subdivision already fail on exact double-star witnesses; see
`notes/two_hub_ratio_order_killtest_2026-08-08.md`. A bounded aggregate
factorization gives a short reproof of the unsubdivided double-star base case,
which was already covered by Galvin--Hilyard. Product-binomial-basis
inequalities now prove the required partial-synchronization base relations for
all leaf counts. Consequently, every double star with an arbitrarily
subdivided hub connector is log-concave; see
`notes/connector_partial_sync_route_2026-08-08.md`. The next bounded two-hub
target is pendant-arm subdivision. Require exact witness replay before
allocating a long run.

One bounded pendant-arm stage is now complete. Exact enumeration of all
163,523 unordered pendant-arm multiset pairs with total pendant weight at
most 24 verifies the recurrence-compatible base invariant. The resulting
computer-assisted theorem covers every connector length and therefore
unbounded total order for these cores. It does not cover pendant weight above
24 or prove full `C_2`; see
`notes/c2_bounded_pendant_core_2026-08-08.md` and
`results/c2_bounded_pendant_core_20260808.json`. The next proof target is the
same base invariant without the pendant-weight bound; ordinary
synchronization is already exactly false for arm vectors `()` and `(1,3)`.

The sharpest observed boundary is the family with pendant-arm vectors
`(1,1)` and `(m,m)`. Exact computation verifies both full base
partial-synchronization relations for `1 <= m <= 500`. For every odd
`m=2h+1`, the top-degree diagonal has the exact all-`h` margin `h+2>0`;
relative slack nevertheless tends to zero, so this family is the natural
place to test any proposed unbounded proof. The all-`h` statement concerns
only that top diagonal, while the full relations through 500 remain finite
computational evidence. Replay with
`verify_c2_bottleneck_family_20260808.py`; its certificate is
`results/c2_bottleneck_family_20260808.json`.
