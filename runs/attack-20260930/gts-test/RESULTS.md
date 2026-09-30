# Tree-shift test: does Csikvári's GTS push the central margin toward the star?
<!-- SUMMARY: Exact test of a transformation route to central log-concavity: the one-step and two-step forms are false; a bounded-step form survives small trees but needs up to 5 shifts in hub constructions (n=225) with no uniform bound in sight; route graded weak, not pursued · status: done · updated: 2026-09-30 -->

Run on 30 September 2026 at Brett's direction ("proceed"). This was the
extremal-graph-theory idea from the discussion of idea-generation techniques.
Everything is exact, with `Fraction` margins. Scripts are in this directory.

## Idea

The census (`../census/`) found the star to be the unique tree with the
smallest central margin
`m(T) = min_{k in [ceil(n/4), q]} (1 - i_{k-1} i_{k+1}/i_k^2)` at every
`n <= 27`. Csikvári's generalized tree shift (GTS) orders the trees on `n`
vertices with the star on top.

- **Source:** P. Csikvári, "On a poset of trees", Def. 2.1 (p. 2) and
  Thm. 2.4 (p. 3). The author's copy is
  <https://csikvarip.web.elte.hu/csikvarin2002uj.pdf>.
- **The shift:** take `x`, `y` such that every interior vertex of the `x-y`
  path has degree 2, and let `z` be `y`'s path neighbour. Move every edge
  `y-w` with `w != z` to `x-w`. The shift is *proper* when neither `x` nor
  `y` is a leaf, and it then adds exactly one leaf.

**Why it would matter.** Suppose every non-star tree had a proper shift that
does not raise `m`. Chaining such shifts up to the star would give
`m(T) >= m(star) = 4/(n+2) > 0` for every tree. That would prove
log-concavity on the window, and so #993.

## Results

| Form tested | Verdict | Evidence |
|---|---|---|
| **Strong:** every proper shift has `m(T') <= m(T)` | false | about 14% of proper shifts raise `m`, by up to 16%; all trees `n <= 14` (`gts_margin_test.py`, `gts_margin_n4_14.json`) |
| **Weak, 1 step:** some proper shift has `m(T') <= m(T)` | false | stuck trees at every `n` from 7 to 16 (1, 1, 3, 6, 5, 5, 11, 14, 20, 38). At n = 9, 10, 11 and 13 the non-star tree with the smallest margin (within 2–3% of the star's; two hubs joined through a middle vertex) is itself stuck (`wfail_profile.py`) |
| **Weak, 2 steps** | false | all trees `n <= 16` repaired except the star of stars SoS(3,3), n = 13. SoS(m,3) needs 3 steps at m = 3 to 6, n up to 25 (`two_step.py`, `star_of_stars_steps.py`; SoS(4,3) rechecked without deduplication in `verify_sos43.py`) |
| **Weak, 3 steps** | survives small trees, fails in hub constructions | holds for all trees `n <= 16`, stars of stars up to n = 31, and doubly nested stars up to n = 31 (`nested_steps.py`). The R1 counterexample (n = 237) needs 4 steps (`r1_killer_steps.py`). The family MSH(c; 9, 2) needs 1, 1, 2, **5**, 4 steps for c = 2, 4, 6, 8, 10 (n up to 281) (`msh_scan.py`) |

One repair is circular: SoS(3,3) at n = 13 is exactly 3 shifts from the
star, and reaching the star is its only repair. The others are not. The star
is 84 shifts away for the n = 237 tree, and SoS(m,3) for m >= 4 repairs in 3
shifts while the star is m shifts away. All the stuck hub trees have margins
1.6 to 1.9 times the star's.

## Verdict: weak, not pursued

The central margin is not monotone along the GTS order. Its valleys sit next
to the star (two-hub trees, stars of stars), and in hub constructions they
are deep enough to need at least 5 shifts. A lemma stated for a bounded
number of shifts would need `k >= 5` already. Nothing here suggests a
uniform bound, and an inequality for a composition of five or more shifts is
not a workable proof device. Letting `k` grow without bound makes the
statement circular: the star lies above every tree, so "some tree above has
a smaller margin" just restates `m(T) >= m(star)`.

A route by transformations would need either a different potential that is
monotone along the shifts and implies the window margin, or a different
transformation. Neither was tried.
