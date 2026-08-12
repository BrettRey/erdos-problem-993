# Adversarial kill-test of the universal C2 base facts beyond pendant weight 24

Replay status: the harness was written by a Sonnet subagent (validated
against the c2 note's known (1,1) x (9,9) margin of 6 before running); the
full sweep then ran twice, once by the subagent (16,531 pairs) and once as
an independent replay in the main session (16,525 pairs; the counts differ
only through wall-clock cutoffs in the time-budgeted phases). Both runs
found zero violations and the same global minima to displayed precision.
The canonical JSON is the main-session replay.

## Question

`notes/c2_bounded_pendant_core_2026-08-08.md` proves that for a fixed
pendant-arm pair `(a, b)` on two hubs, four base facts imply log-concavity of
the two-hub tree for every connector length: `B` log-concave, `C`
log-concave, `C ~_p B`, `xC ~_p B`. Every pair with total pendant weight at
most 24 was checked exhaustively. This note asks whether any pair with
total pendant weight strictly above 24 breaks one of the three genuinely
open facts (`B` log-concave, `C ~_p B`, `xC ~_p B`). `C` log-concavity is a
theorem, not a target: `C` is the independence polynomial of the spider
formed by merging both arm multisets at one hub, and Li, Li, Yang, and Zhang
prove every spider strongly log-concave (arXiv:2501.04245).

## Method

Script: `scripts/killtest_c2_universal_20260811.py`. All arithmetic exact
(Python integers and `fractions.Fraction`); floats appear only for display.
For every checked pair `(a, b)` the script builds `B` and `C` from the exact
recurrence in the bounded-core note (`P_{-2} = 0`, `P_{-1} = P_0 = 1`,
`P_n = P_{n-1} + xP_{n-2}`), then computes, for each of `B` log-concavity,
`C ~_p B`, and `xC ~_p B`, the exact minimum normalized margin
`(required - bound) / max(1, |bound|)` over either the full index grid or a
banded subset for very large degrees (see coverage below). A negative
margin is a violation; every candidate violation is independently confirmed
by rebuilding the actual tree (the `t = 0` two-hub tree for `B`, the merged
spider for `C`) and rerunning a from-scratch rooted-tree independence-set
DP, matching the coefficient lists exactly.

Search plan, per the dispatch brief:
1. Structured grid: `(1^r) x (l^s)`, `(k,k) x (m,m)`, `(1,1) x (m,m)` out to
   large `m`, `(1,1) x (2h,2h)` and `(1,1) x (2h+1,2h+3)` neighbor families,
   single long arm `(1^r) x (M)`, and a handful of three-arm shapes against
   a single long arm.
2. Random search: arm-multiset pairs at total weight 25-150, with skewed
   length distributions (many 1s plus one giant arm; a balanced partition;
   a geometric-length partition).
3. Hill-climbing from the tightest pairs found in 1-2 (weight-capped so
   each step stays cheap), mutating by adding/removing/tweaking an arm or
   moving one across hubs, accepting a move only if it lowers the minimum
   normalized margin.
4. Independent DP confirmation of any candidate violation before it is
   reported.

### Coverage tiers and why

Computing `B` and `C` exactly for the doubled-arm family `(1,1) x (m,m)`
costs roughly `O(m^{3.5-3.9})` empirically (build plus a full two-index
partial-synchronization scan): 0.29s at `m = 500`, 3.3s at `m = 1000`, 46s
at `m = 2000`, 221s at `m = 3000` (measured directly, not extrapolated). A
dense per-integer full scan out to `m = 3000` would cost on the order of a
day of wall time. The search therefore uses three tiers per large family:
a dense full scan (every integer) up to a moderate bound where the cost is
small; a sparse full scan at specific larger parameter values, largest
first so a time-budget cutoff still reaches the frontier; and a cheap sweep
across the full requested range that computes the log-concavity margin
exactly (cheap, linear in degree) plus a banded partial-synchronization
scan restricted to the top and bottom `~80` coefficients and a `~15`-wide
near-diagonal stripe, rather than the full `O(d^2)` grid. The banded check
is exact wherever it looks; it is not a claim of full-grid coverage at
every parameter value in that tier.

## Coverage achieved

16,531 pairs evaluated in total, 2,951 seconds of wall time. Zero violations.
(The table below reflects the subagent's run. The canonical JSON is the
main-session replay, which the time budgets cut at 16,525 pairs, dense tier
to m = 635; its global minima agree with the numbers below to displayed
precision and its verdict is identical.)
Per phase (full = exact full `O(d^2)` grid on both sync facts; banded = exact
log-concavity plus the top/bottom-80-coefficient and 15-wide near-diagonal
band described above):

| phase | pairs | coverage | stopped on budget |
|---|---|---|---|
| `(1^r) x (l^s)`, r,s<=12, l<=40 | 5,760 | full | no |
| `(k,k) x (m,m)`, k<=5, m<=200 | 1,000 | full | no |
| `(1,1) x (m,m)`, m=1..641 | 641 | full | yes (budget) |
| `(1,1) x (m,m)`, m in {3000,2500,2000} | 3 | full | yes (budget; 6 of 9 requested points not reached) |
| `(1,1) x (m,m)`, m=1,11,...,2921 (step 10) | 293 | banded | yes (7 of 300 points short of m=3000) |
| `(1,1) x (2h,2h)`, h=1..300 | 300 | full | no |
| `(1,1) x (2h,2h)`, h=1,6,...,796 (step 5) | 160 | banded | no (full requested range reached) |
| `(1,1) x (2h+1,2h+3)`, h=0..300 | 301 | full | no |
| `(1,1) x (2h+1,2h+3)`, h=0,5,...,800 (step 5) | 161 | banded | no (full requested range reached) |
| `(1^r) x (M)`, r=1..12, M=1..300 | 3,600 | full | no |
| `(1^r) x (M)`, r in {1,6,12}, M=1,11,...,991 (step 10) | 300 | banded | no |
| three-arm shapes x single/scaled arm | 1,812 | full | no |
| random pairs, weight 25-150, 3 length-distribution styles | 2,000 | full | no |
| hill-climb, 20 seeds x 10 restarts x up to 200 steps, weight capped at 220 | 200 trajectories | full (each step) | no |

The `(1,1) x (m,m)` family is the only one where a full `O(d^2)` scan at
every integer parameter was infeasible: a single full check costs 0.3s at
m=500, 3.3s at m=1000, 46s at m=2000, 221s at m=3000 (measured directly).
The banded tier still gives an exact answer at every point it touches; it
just does not certify the interior of the grid away from the top/bottom
corners and the near-diagonal at points between the sampled ones.

## Result

No violation of `B` log-concavity, `C ~_p B`, or `xC ~_p B` at any total
pendant weight above 24, anywhere in this search. The global minimum
normalized margin for each fact, exact:

- **B log-concave:** `0.0017314`, at `(1,1)` vs `(3000,3000)` (weight 6002,
  degree 3003), interior index `k = 1269` (about 42% of the way through
  the degree, not at either end). Found by the full exhaustive scan at
  `m = 3000`.
- **C ~_p B:** `0.0017314`, same pair, same index (`lower = upper = 1269`).
  Also full exhaustive. The near-equality of these two minima at the same
  `(a,b,k)` looks structural rather than coincidental, but that was not
  investigated further here.
- **xC ~_p B:** `0.00034235`, at `(1,1)` vs `(2921,2921)` (weight 5844,
  degree 2924), top corner `lower = upper = 2924`, margin exactly 1462.
  Found by the banded sweep, `m = 2921` is odd (`h = 1460`), and 1462 is
  exactly `h + 2`, matching the closed-form formula in
  `c2_bounded_pendant_core_2026-08-08.md` on the nose. `m = 3000` is even
  and falls outside that closed form; the three full-exhaustive sparse
  points (2000, 2500, 3000) happen to all be even, so the sparse tier never
  directly re-hit the odd-`m` bottleneck family, but the banded sweep (which
  walks every 10th integer starting from the odd number 1, so it lands on
  odd `m` throughout) does, independently reconfirming the theorem's closed
  form at many points up to `h = 1460`.

## Frontier structure

Every one of the tightest 20 pairs found, across the whole search, is an
instance of `(1,1) x (m,m)` at the largest sampled `m`, all tight in
`xC ~_p B`, all at the top-corner index. The margin shrinks monotonically
as `m` grows (consistent with the proved `h + 2` numerator against a
denominator growing like `h^2`), so nothing here suggests the family stops
being the bottleneck at higher weight; if anything the search confirms the
asymptotic approach to zero continues smoothly through `h = 1460` with no
irregularity.

Hill-climbing from the 20 tightest pairs (weight capped at 220 so each
step stays cheap) never found a configuration competitive with the large-m
family: its best endpoint, a near-caterpillar shape with arms
`(1,1,1,1,1,2,3,3)` on one hub (weight 13) against 95 arms on the other
hub (74 ones, 6 twos, 13 threes, one 4, and one 79; weight 208; total
weight 221), reached normalized margin `0.0283` in the `xC ~_p B` fact,
still tight in the same fact as the global bottleneck but two orders of
magnitude looser than its normalized margin. 68 of 200 proposed steps were
accepted in that run; other restarts plateaued similarly. So within the
explored weight-capped neighborhood of structurally varied shapes, nothing
rivals the symmetric-double-arm family for tightness, and the interesting
frontier stays exactly where the closed-form analysis already said it was.

## Relation to `scan_two_branch_lc.py`

`scan_two_branch_lc.py` exhaustively enumerates C2 trees (at most two branch
vertices) by total vertex order `n` (default up to 24), which bounds
pendant weight plus hubs plus connector vertices together. It therefore
cannot reach pendant weight above 24 combined with an unbounded connector,
and does not cover the search space explored here. It was not rerun.

## Relation to concurrent work

`notes/c2_single_index_laws_2026-08-11.md` (same day) reformulates the two
open synchronization facts as single-index laws `U`, `V`, `W` plus a proved
assembly lemma, and reports the same qualitative structure found here: the
laws hold by aggregate cancellation across bilinear generator-pair terms,
not termwise, and `U` has isolated exceptions confined to a shallow
reflected-depth head strip. That note's search is independent of this one
(different code, different decomposition). Where both notes report a
number for the same object, they should be read as two independently
computed data points, not as one confirming the other's arithmetic.

## Verdict

No violation. 16,531 pairs at total pendant weight above 24 (up to weight
6,002), spanning structured grids, random search, and hill-climbing, all
satisfy `B` log-concave, `C ~_p B`, and `xC ~_p B` (`C` log-concavity also
checked throughout as a harness self-check and never failed, consistent
with its being a cited theorem). The `(1,1) x (2h+1,2h+1)` family stays
the global bottleneck at every scale reached, with margins tracking the
proved `h + 2` formula exactly at points checked as far as `h = 1460`. This
is finite-corpus evidence for the universal C2 target, not a proof: the
banded and sparse tiers do not certify every parameter value in
`[25, 3000]`, and nothing here rules out a violation at weights or shapes
outside this search.

## Limitations

The banded and sparse tiers are not full-grid proofs at every parameter
value they touch; they are exact at the indices actually checked. No tier
here is a proof for unchecked parameter values, and nothing here proves or
disproves the universal target, the bounded-pendant-core theorem's scope,
or Erdos Problem 993.

Replay:

```text
python3 scripts/killtest_c2_universal_20260811.py --out results/c2_universal_killtest_20260811.json
```
