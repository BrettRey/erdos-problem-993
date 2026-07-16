# Why trees resist: two empirical barrier laws from the July 2026 disproof campaign
<!-- SUMMARY: Quantitative failure laws from the 2026-07-15/16 valley-first disproof campaign: rebound deficit >= ~11/n in coefficient space, dominance-angle tradeoff in root space · status: private research note, empirical · updated: 2026-07-16 -->

Status: private research note. Everything here is exact computation plus
measured scaling, not theorem. No manuscript change; the E-JC submission
is untouched. Scripts and result files are listed at the end; every
number is replayable.

## 1. Setup

For a tree `T` with independent-set counts `i_0, ..., i_alpha`, define
the valley margin

    V(T) = max_b min( max_{a<b} i_a , max_{c>b} i_c ) / i_b .

`V > 1` (exact integer comparison) is precisely a unimodality
counterexample. Two background facts frame everything:

- The sequence is weakly decreasing from `ceil((2 alpha - 1)/3)` on
  (Levit--Mandrescu), so the rebound index `c` of any witness must sit
  before that threshold.
- All trees on `n <= 29` are unimodal (exhaustive, 8,691,747,673 trees).

The campaign (2026-07-15/16) scored trees directly by `V` and, after a
correction described in section 3, by the descent-thresholded rebound

    R_gap(theta) = max over b with prefixmax(b) >= (1+theta) i_b
                   of suffixmax(b) / i_b ,

reported with the rise distance `c - b` of the achieving pair.

## 2. Barrier law 1: rebound deficit is Theta(1/n) with a growing constant

Roughly 160,000 exact evaluations across every composition mechanism we
could construct or mutate into:

| mechanism | scale reached | best genuine rebound |
|---|---|---|
| two-phase spider bouquets (incl. perturbed `T_{3,M,N}`) | n = 7,299 | `n(1-R) ~ 13.3` at n ~ 320 |
| two-type hybrids, optimized mix, phase gap t = 6..90 | n = 5,073 | constant flat in t; grows with n (20.5 at n = 2,537) |
| dumbbells (product realizations) | n = 6,563 | 15.8--16.8, growing |
| depth-3 phase stacking (4,298 configs) | n = 2,600 | 10.9 at n = 1,378 |
| depth-4 / connector architectures (independent Codex search) | n = 3,000 | ~4--5 on plateau-scale descents only |
| subtractive edge-join carving (1,406 joins, all vertex roles) | n ~ 600 | no carving; ~0.2% sensitivity to join vertex |
| exact-pathology interference (21 LC-failure trees, 540,492 joins) | n ~ 100--200 | 5.3 at n ~ 104; dilutes under scaling |
| free mutation (263k + 612k evals) | n <= 600 | never escapes the envelope |

Three uniform findings:

1. Under a genuine prior descent (`theta = 0.001`), the rebound deficit
   never beat `~11/n`, and the constant grows with n and with theta.
2. Every achieving rebound in the entire campaign has rise distance
   `c - b = 1`: locally slowed decay, never a cross-gap recovery. A
   counterexample IS a cross-gap recovery.
3. The `T_{3,M,N}`-type late shoulder (the mechanism behind every known
   log-concavity failure, including all 19 exact n = 28 failures) slides
   out of the legal window as it scales: its rise index converges to the
   decreasing-tail threshold from below and crosses it. That mechanism
   cannot produce a witness at any size.

## 3. A cautionary artifact: mode flatness and the lattice offset

The raw margin `V` is contaminated near the mode: for any smooth hump,
`1 - V ~ (local ratio slope, ~C/n) x (fractional offset of the r = 1
crossing from the integer lattice)`. The offset is quasi-random in
(0,1) along a family ray, which produces spurious "deficit
acceleration" and non-monotone constants. One depth-3 ray looked like a
new `1/sqrt(n)` decay class until component decomposition showed its
scoring index was the mode itself, with the root-included phase nowhere
dominant. Plateau-approach can push `V` arbitrarily close to 1 on
perfectly unimodal families; only descent-thresholded rebounds carry
disproof information.

## 4. Barrier law 2: the dominance-angle tradeoff in root space

A rebound at distance `c - b >> 1` corresponds spectrally to a
conjugate root pair at substantial angle from the negative real axis
with near-dominant modulus (subdominant pairs' oscillations are damped
like `(z_min/|z|)^k`). Measured with certified Arb root isolation
(python-flint; a claw-free control exposed float64 root-finding as
unusable ~-- it reported angle 0.685 for a path, certified value 0):

- Certified frontier of known objects: smooth families 0; the n = 28
  LC-failure trees 0.062--0.105; `T_{6,6,1}` 0.256 (angle tracks known
  convexity pathology exactly).
- Certified evolution (611,897 evals): the angle at modulus ratio <= 2
  rose to 0.444 in seconds (the old ceiling was family bias), then
  plateaued for 49 minutes with zero improvements, no witness, and
  `R <= 0.90`, `c - b = 1` on every champion.
- The champion's own spectrum shows the invariant: angle falls
  monotonically as modulus dominance rises (ratio 2.0 -> 0.44 rad;
  2.5 -> 0.26; 2.7 -> 0.13; 2.76 -> 0.03), heading to 0 at amplitude
  parity.

A counterexample needs the empty corner of that frontier: near-dominant
modulus AND substantial angle. This is the spectral form of law 1.

## 5. What would falsify these barriers

- A tree family whose genuine-descent rebound deficit decays faster
  than `c/n` with a shrinking constant, or any single tree with
  `c - b > 1` on a thresholded descent.
- A certified conjugate pair at modulus ratio < 1.5 with angle > 0.3.
- Either would justify reopening the disproof program at scale; nothing
  in ~1.4 million exact evaluations came near.

## 6. Caveats

- Empirical, one session's compute, one mutation kernel; the evolution
  plateaus could be kernel artifacts (though 610k evals without
  improvement is severe).
- Root-space search bounded at n <= 150 and dominance factor 2.0.
- The conjectural positive reading (some universal `V <= 1 - c/n`, which
  would prove Erdos #993) is consistent with all data but is a separate
  program; see the 2026-07-11 posture in STATUS.md.

## 7. Replay

Tools (repo `scripts/`): `valley_search.py` (grammar sweeps + self-test
against the generic DP), `product_valley_search.py` (forest/product lane
with dumbbell realization), `valley_scaling_probe.py`
(Kronecker-packed exact convolution for n in the thousands).

Dated probes (repo root): `scratch_valley_freeform_20260715.py`,
`scratch_deficit_constant_map_20260715.py`,
`scratch_depth3_valley_20260715.py`,
`scratch_depth3_scaling_20260715.py`,
`scratch_depth3_mix_reopt_20260715.py`,
`scratch_rgap_rescore_20260715.py`,
`scratch_join_carving_20260716.py`,
`scratch_bump_interference_20260716.py`,
`scratch_root_herding_20260716.py`,
`scratch_certified_root_evolution_20260716.py`.

Results: `results/root_herding_best_20260716.json`,
`results/certified_root_evolution_20260716.json`. The n = 26/28
failure-tree blocks come from `results/analysis_n26.json` and
`results/analysis_n28_modal_lc_nm.json`. Certified root isolation needs
`python-flint` (0.9.0); the session venv lived in the scratchpad, so
install fresh: `python3 -m venv venv && venv/bin/pip install
python-flint networkx`. Full decision trail: DECISIONS.md, 2026-07-15
and 2026-07-16.
