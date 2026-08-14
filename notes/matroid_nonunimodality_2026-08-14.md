# Rota's unimodality conjecture is false: what arXiv:2608.07342 does and does not mean for #993
<!-- SUMMARY: Read of Divoux-Larson-Lowen-Wang (Rota flat-count disproof) + exact transfer probe; trees already have steps 1-2 of their recipe (severe multi-break LC failure, ratio 2.4e8) but confined to the LK tail; missing step is the tilt · status: assessed, no decision taken · updated: 2026-08-14 -->

**Paper:** Divoux, Larson, Lowen, Wang, *Matroid flat counts can have many
peaks*, arXiv:2608.07342v1 (math.CO, 7 Aug 2026). On disk:
`notes/literature/arxiv_2608_07342.pdf` / `.txt`; abstract note in the
research vault (`260807342-matroid-flat-counts-can-have-many-peaks`).

**Probe artifacts (this note's numeric claims, exact arithmetic):**
`scripts/probe_matroid_recipe_transfer_20260814.py`,
`results/matroid_recipe_transfer_20260814.json`. Family builders reused from
`scripts/probe_absorption_margin_20260811.py` (validated 2026-08-11 against
Bautista-Ramos's exact break sets for TG_{4,6} and TG_{5,6}).

## What the paper proves

Rota's 1970/71 ICM conjecture: for any matroid of rank r, the Whitney numbers
of the second kind (flat counts by rank) (W_0, ..., W_r) are unimodal
(Conjecture 1.1, p. 1). Theorem 1.2 (p. 2): "For each positive integer m,
there is a matroid whose sequence of Whitney numbers of the second kind has
exactly m peaks." This also refutes Mason's log-concavity conjecture for
flats. Scale (p. 2): smallest non-unimodal example their method gives has
ground set 4957; smallest non-log-concave, ground set 78.

Boundary they state themselves (p. 2): "Our techniques seem incapable of
producing a counterexample to the points-lines-planes conjecture." So the
bottom of the rank sequence (W_2^2 >= W_1 W_3) survives their method; every
peak they build lives near the top (coranks m+2s+1, proof of Thm 1.2, p. 7).

Provenance note (p. 3): "In producing this work, ChatGPT played an important
role in finding some initial examples which the authors were able to
generalize." The seed-by-machine-search, amplify-by-hand pattern is the same
one this campaign runs.

## Anatomy of the construction (three separable steps)

1. **Severe LC failure in a parametric family.** G_t = graphic matroid of the
   generalized theta graph: two endpoints, four internally disjoint paths,
   three of length t and one a single edge. Lemma 2.1 (by corank eta):
   W^eta = Theta(t^{eta+2}) for eta in {1,2} but Theta(t^{eta+3}) for
   eta >= 3. The exponent increment jumps from 2 to 3, so the violation
   ratio W^1 W^3 / (W^2)^2 "grows linearly in t" (p. 2). Severity is
   unbounded, and it sits at corank 1-3, i.e. near the top, on the
   decreasing end of the sequence.

2. **Replication by direct sums.** Flat counts of a direct sum convolve.
   Lemma 2.4: the exponent function of G_t^{oplus m} is piecewise linear in
   three regimes; in the middle band consecutive ratios alternate Theta(t)
   and Theta(t^2), giving ~m spaced LC failures. Note G_t^{oplus m} is still
   unimodal; direct sums alone never produce a valley.

3. **Geometric tilt, realized in-category by the Whittle q-lift.** The
   elementary observation (p. 2): a nonnegative sequence is "log-concave with
   no internal zeroes if and only if, for every c > 0, the sequence
   (a_0, c a_1, c^2 a_2, ..., c^m a_m) is unimodal." LC is tilt-invariant;
   unimodality is not; LC is exactly unimodality under every tilt. The q-lift
   realizes the tilt inside the matroid category: W_i(M~) = q^i W_i(M) +
   W_{i-1}(M) (Bonin-Qin, Lemma 3.1, p. 7). Choosing q ~ t^{3/2}, between the
   two alternating ratio scales, flips every other middle-band step downhill:
   m+1 peaks. Truncation then prescribes the exact peak count (p. 8).

The tilt is the step that converts. Steps 1-2 manufacture tilt-exposable LC
damage; step 3 slides the mode into the damage zone with zero smoothing.

## Transfer probe: where trees stand on each step (measured today)

**Step 1 transfers, and more strongly than we had on record.** The known tree
LC failures at n = 26 are marginal (ratios 1.1453 and 1.0295,
`results/analysis_n26.json`). But along Bautista-Ramos's alternating-break
family TG_{m,t} (arXiv:2511.00334) the violation ratio is unbounded, growing
exponentially in t: at m = 4 it runs 6.5 (t=8), 3.6e3 (t=17), 8.0e7 (t=32);
grid max 2.39e8 within m <= 12, t <= 32 (n = 786 tree). Severe LC failure is
not a matroid luxury; trees have it.

**Step 2 transfers, inside trees, and is already in this family.** TG_{m,t}
grafts m copies of T_{3,t} at a hub; for t large it has exactly m breaks at
the odd reflected offsets 1, 3, ..., 2m-1. That is the direct-sum replication
step realized in-category (with the pendant correction x*I(S_{2,t})^{3m} as
the price of staying a tree). Powers of the marginal n = 26 failure die
instead (I(T)^p fully LC for every p in 2..12): replication preserves only
severe damage, which is exactly why severity mattered in step 1.

**But all of it sits in the immunized zone.** Every break in the entire grid
lies above the Levit-Kadrawi threshold ceil((2*alpha - 1)/3) (their
Lemma 2.19 / Cor 2.20 reduction: prefix LC below that index implies
unimodality for trees). Structurally, within this family the break band can
never reach the prefix: break depth is <= 2m - 1 while the threshold sits at
reflected depth ~ alpha/3 = (t+1)m + O(1), and 2m - 1 < (t+1)m for every
t >= 1. Same geometry as the matroid side: their peaks live at corank O(m)
near the top; their method cannot touch points-lines-planes at the bottom.
Pathology is a head phenomenon in both categories. The difference: trees have
a theorem at the head (the decreasing tail), matroids only had top-heaviness,
which peaks respect.

**Step 3 is what trees lack, and the lack is quantified.** Per break k the
tilt window c in (i_k/i_{k+1}, i_{k-1}/i_k) is nonempty iff LC fails there,
so a tilted tree sequence goes non-unimodal trivially; for TG_{8,8} a single
c = 4764 exposes 5 peaks (windows overlap consecutively; their global
intersection is empty, so one tilt cannot expose all m at once without the
three-regime exponent engineering the theta construction supplies). No tree
operation realizes the tilt: every mode-raising move trees admit is
binomial-type (pendant attachment convolves with (1+x)-factors, moving the
mode by ~s/2 while moving alpha by s, so the head recedes faster than the
mode advances), and binomial convolution smooths precisely the marginal dips
it approaches. That is the content of the D17/D20/D23 obstructions and the
manuscript's mean bounds, now nameable in one sentence: **trees pin the
sequence to the c = 1 slice of the tilt orbit, and everything this campaign
has measured says the c = 1 slice cannot reach the damage zone.**

## Answers to the two questions posed

**"What is the corresponding object, and does the exhaustive search cover the
analogous shapes?"** The corresponding object is not missing: it is the
pattern-graph / TG family already on disk, which is a theta analogue in the
only sense trees permit (m parallel channels from a hub instead of between
two hubs; the two-hub versions are the connectors and double brooms the
campaign has hammered, and the one-hub theta limit, spiders, is fully
log-concave by Li-Li-Yang-Zhang). The exhaustive record (n <= 32, external)
does not cover these shapes; they start near n = 100. The targeted lanes do
cover them, and today's probe shows they carry the theta phenomenon
quantitatively. What trees cannot express is not the shape but the operation:
there is no q-lift for independence polynomials.

**"Should the prior on #993 move?"** Yes, but the useful movement is a
redirection, not a straight discount:

- The generic Bayesian point is real. A 55-year unimodality conjecture with
  the same evidence profile ("natural sequence, small cases, Hodge-adjacent
  confidence") fell to an explicit construction, found partly by machine
  search. Exhaustive small-case verification is weak evidence; their
  counterexample needs ~5000 elements, far beyond exhaustive matroid
  enumeration. Our n <= 32 wall should be weighted accordingly (the
  structured-family searches at n up to ~11k are the right hedge and already
  exist).
- The specific mechanism, though, is one whose tree-side absence this
  campaign has spent months measuring, mostly without naming it. After this
  paper the ledger reads: severity available (unbounded, today's probe),
  replication available (TG), tilt unavailable, damage confined to the
  theorem-protected tail. #993's entire content is now two statements: the
  prefix-LC conjecture (no LC break below the LK threshold) and mode
  pinning. Rota's fall is a warning against trusting the first on small-case
  evidence alone.
- It also cuts one way that is mildly reassuring: the failure their method
  cannot produce (points-lines-planes, the bottom of the sequence) is the
  exact analogue of the zone #993 needs protected (the prefix). The
  cross-category pattern is bottom rigid, top loose.

## Actionable follow-ups (not started; Brett's call)

1. **New search objective: break depth, not severity or rebound.** No lane
   has run "maximize reflected depth of the deepest LC break relative to
   alpha" as an explicit objective (lanes to date optimized rebound deficit,
   root angle, near-miss ratio, absorption margin). All known LC breaks sit
   at reflected depth O(m) = o(alpha); the sync-violation band, a weaker
   invariant, already penetrates to 0.54*alpha and escapes the decreasing
   zone at m >= 21 (F4/F5, `notes/partial_sync_obstruction_2026-08-11.md`).
   The gap between those two facts is where a counterexample would have to
   live. An evolutionary lane with depth as fitness is cheap and is the
   honest kill-test the prefix-LC conjecture now owes us.
   **[RUN same day at Brett's direction: see
   `notes/break_depth_lane_2026-08-14.md`. Outcome: dist = k - thr >= 4
   everywhere (58.5M exact evaluations, complete small mutation balls);
   dist = 3 reduced to a finite window (depth-2 break, alpha in {15,16},
   n in [29,32]).]**
2. **Manuscript hook, for a revision round only.** If E-JC comes back, the
   "why trees resist" discussion can now cite 2608.07342 as the
   cross-category contrast: the recipe that killed Rota stalls for trees at
   the tilt step, and the mean-bound results are the measurement of that
   stall. Do not touch `paper/main_v2.tex` now (submitted 2026-04-25, OJS
   15526).
3. **erdosproblems.com/993 comment.** The 2026-08-12 comments predate this
   paper. A short pointer ("Rota unimodality fell 2026-08-07; the
   construction's LC-severity + replication steps have tree analogues in the
   Bautista-Ramos family, confined to the Levit-Kadrawi tail; the missing
   ingredient is a tilt-realizing operation") would be on-topic. Needs
   Brett's sign-off; nothing posted.

## Provenance

Session 2026-08-14, prompted by Brett supplying the abstract. Paper read in
full (9 pp). All numeric claims above replay from
`scripts/probe_matroid_recipe_transfer_20260814.py` (exact integer/rational
arithmetic; runtime ~1 min). No manuscript, STATUS, or public surface
touched. No decision recorded; this is an assessment note.
