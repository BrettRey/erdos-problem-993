# Effective central log-concavity attack, 30 September 2026
<!-- SUMMARY: Multi-agent attack after Fang et al.: N0 audited (10^(1.68e7) as formalised, ~10^2070 certified-sharpened, floor >= 10^66-10^100); binomial-smoothing hybrid verified sound but only ~10^495; exact central-margin census n<=27 (no window LC failure, star extremal); deterministic lemma R1 exhaustive to n=24 but refuted at n=237 · status: complete · updated: 2026-09-30 -->

Brett authorized the attack ("attack whatever you think is worth attacking")
with expiring Claude credits. Research only: no manuscript touched. Each lane's
full structured report is `RETURN.json` in its directory (the harness blocked
subagent `.md` reports, so the returns were saved verbatim by the parent);
scripts, data and `manifest.yaml` sit beside it. Wave-1 workflow run
`wf_f7bb4547-c54`; project HEAD at launch `e94d47f`.

## The target

Fang et al. (arXiv:2609.20961) make both ends deterministic: `i_k` is
nondecreasing through `ceil(n/4)` (Prop. 8.2) and nonincreasing from
`q = ceil((2 alpha - 1)/3)` (Levit–Mandrescu). So unimodality of every tree
follows from log-concavity on the window `[ceil(n/4), q]`. Their proof of
that step needs `n >= N0`, with N0 never computed.

## Wave 1 findings

**N0 audit (two independent agents, `audit-a`, `audit-b`).** They agree on
everything that matters:

- N0 is computable. No step is existential in the compactness sense.
- Instantiating the Lean proof's own witnesses (commit `b2a1d3e`) gives
  `log10 N0 = 1.68 x 10^7`. Both audits got this figure independently.
- The paper route with its own best parameters gives `log10 N0` of
  2.4 x 10^5 (b) or 3.6 x 10^5 (a); the difference comes from parameter
  choices.
- **Bottleneck:** the centroid-decomposition mean-shift term `b^(a-1)` in
  (5.6). Lemma 4.2 pins the root-moment exponent within 0.0016 of 1, because
  `lambda_max = 12` sits just below `e^(5/2)`.
- **Cheapest fix, found independently by both:** replace Lemma 4.2's
  interpolation with a certified one-dimensional supremum. Audit a:
  `F_12(1.8969) <= 0.99518`. Audit b: `p = 19/10`, `rho <= 199/200`. Either
  gives `1 - a` of about 0.053 and `log10 N0` of about 2 x 10^3 to 8 x 10^3.
- **Architecture floor:** even with idealized constants, `N0 >= 10^80` to
  `10^540`. The local-limit step needs accuracy of order `R^-3`.
- Heuristic float probes put the root-moment growth exponent at about 0.8 at
  high fugacity, so `a` near 1 is not only loose constants.
- Audit a also finds `lambda_max = 3` refuted for forests (a 17-vertex tree
  has `E_3 X/alpha = 0.6542`).

**Central-margin census (`census`).** Exact and exhaustive over
1,198,737,961 trees with `10 <= n <= 27`, with counts matching A000055 and a
Python cross-check for `n <= 14`:

- There is no log-concavity failure in the window `[ceil(n/5), q]`.
- The star is the unique minimizer at every n, with minimum normalized margin
  `4/(n+2)` (n even) or `4n/(n+1)^2` (n odd), so `n * min delta -> 4`.
- No member of 473 family trees up to `n` of about 300 goes below the star.
- The scale-free margin `V_k * delta_k` (tilted variance times margin) is at
  least 0.51 everywhere seen and tends to 1. This is a float diagnostic.
- Newton's inequality fails on the window (minimum ratio about 0.92), so a
  real-rootedness route cannot work as stated.
- Observed on everything seen, not proved:
  `delta_k >= 0.9 (alpha+1)/((k+1)(alpha-k+1))` on the window.

**Routes.**

- `route-lclt` (graded *plausible*): an exact bipartite binomial-smoothing
  identity. It is the conditioning inside Fang's Lemma 6.1. A certificate
  then reduces the local limit to a **fixed-accuracy** Kolmogorov bound for a
  conditional count `M'`. What remains unproved is a uniform, effective,
  fixed-accuracy CLT for `M'`, whose fugacities reach about 2.5 to 4.3, above
  bounded-degree uniqueness thresholds.
- `route-freecount` (graded *weak*, but with the only deterministic target):
  - LC at `k` holds iff `sum_v D_v(k) <= 0`, where
    `D_v(k) = k i_{k-1}(T) i_k(T-N[v]) - (k+1) i_k(T) i_{k-1}(T-N[v])`.
  - The pointwise lemma `D_v <= 0` is **false**: the hub-star H(9,2), n = 28,
    violates it by a constant factor, up to about 1.15.
  - The neighbourhood average **R1**, `sum_{v in N[u]} D_v/(deg v + 1) <= 0`
    on the window for every `u`, survived every wave-1 test. It would prove
    #993 with no size threshold.
- `route-pbmix` (graded *weak*): the Poisson-binomial theorem's constant 1/4
  is too small for densely branching trees. Its bipartition-class variant
  turns into the lclt route.

## Wave 2: PAUSED at 14:57 (Brett closing the laptop)

Four agents were launched at 14:38 (run `wf_0dbf2e94-a82`) and stopped
after about 19 minutes, before any returned. Their partial scripts and data
are in `r1-exhaustive/`, `r1-adversarial/`, `r1-analytic/` and
`lclt-check/`. There are no `RETURN.json` files for wave 2. Three orphaned
`sa_r1.py` searches were killed.

**R1 is very likely refuted.** Two agents, working separately, found
violations in the same family: a centre joined to 8 or 9 hub-stars H(m,2)
with m = 9 or 10.

- Smallest recorded case: `MSH(4x9, 4x10)`, n = 237, alpha = 160, window
  [60, 107].
- At k = 107 = q, `L_centre` is positive in exact arithmetic. LC still holds
  there, with relative margin 0.029, so R1 fails while its conclusion is
  true.
- More cases: `S(8,10,2)`, n = 249, k = 111 and 112; `S(9,9,2)`, n = 253.

Records: `r1-adversarial/data/counterexamples_summary.json` and
`r1-analytic/data/recheck_MSH_4x9_4x10_generic.json`. The second holds the
explicit parent array, rechecked with the wave-1 generic DP (`pv_lib`).
**Before recording R1 as dead, confirm with a third, independent
implementation.**

`r1-exhaustive` finished `19 <= n <= 24` (log `run_19_26.log`, raw counts
in `raw/`); its summary was not written. `lclt-check` has partial exact tests
of the smoothing identity and kernel constants.

### Resumed 15:49

- **R1 refuted, confirmed three ways.** I wrote
  `r1_third_check_20260930.py` fresh; it shares no code with `pv_lib`,
  `r1lib` or the closed forms, and reads only the edge list.
  - `networkx.is_tree` holds, n = 237, alpha = 160, window [60, 107].
  - It reproduces the recorded polynomial and the two double-count
    identities at every k.
  - It finds exactly one violation: `u` = the degree-8 centre, `k = 107 = q`,
    `L_u > 0`, normalized load `3.3226e-5` (matching the record), with LC
    holding at that k.
  - Construction: a degree-8 centre joined to 8 hubs; each hub has 9 or 10
    further neighbours, each carrying 2 leaves.
- **R1 exhaustive, 19 <= n <= 24** (`r1-exhaustive/results.json`, merged
  and summarized after the pause):
  - 63,037,252 trees, every count equal to A000055.
  - 6,830,511,974 window `(u, k)` checks: 0 R1 violations, 0 pointwise (PV)
    violations.
  - The largest PV ratio rises from 0.911 at n = 19 to 0.971 at n = 24,
    always at a hub-star. So the smallest PV counterexample has
    `25 <= n <= 28`; H(9,2) at n = 28 is known.
  - The C code agreed exactly with `pv_lib` on all trees with n <= 13
    (`crosscheck_n2_13.json`).
- So R1 holds through n = 24 and fails in a constructed tree at n = 237.
  The deterministic free-count route has no live target. A radius-2 average
  would likely be defeated by the same construction nested one level deeper;
  that has not been tested.
- The smoothing-certificate check (`lclt-check`) was relaunched as a single
  agent. It finished at about 16:20; its report is `lclt-check/REPORT.md`.

### lclt-check result: sound, marginal gain

- **Claims A–C are correct.** The binomial-smoothing identity, the Abel
  bound and the certificate were re-derived and checked exactly on all 5,445
  trees with n <= 14 plus 7 larger trees. The certificate's tail term can be
  sharpened by at least a factor of 2, which removes every small-n
  certificate failure route-lclt reported.
- **`c_* = 0.26420` holds only as a limit.** As a uniform bound,
  `TV(h_B) v^(3/2) <= int|phi'''|` is exactly refuted (B = 22 at lambda = 12
  gives about 1.58).
- **The hybrid removes Fang's `R^-3` floor, but not the root-moment
  exponent.**
  - Like for like, the floor falls from 10^240 to 10^66.
  - At the certified constants, the explicit threshold falls from 10^2070
    to about 10^495.
  - The centroid mean shift is still paid at exponent `2/(1-a)`. With `a`
    of about 0.8 to 0.88 (the audits' growth probes), N1 stays above 10^100
    even with unit constants.
  - Using only constants the audits can prove, N1 is about 10^1218.
- **Nothing reaches enumeration range.** That would need an unproved
  Kolmogorov rate `K n^(-1/2)` for `M'`, which is exactly the missing
  content of Lemma L (conditions V, K, T).
- **Aristotle candidate: low value, not sent.** The exact finite
  certificate, with sub-targets the bipartite identity
  `I_T(x) = sum_J x^|J| (1+x)^(B_J)` and the Abel inequality.

## Where the afternoon leaves #993

Unimodality of every tree now reduces to log-concavity on
`[ceil(n/4), q]`. That holds for every tree with n <= 27, exhaustively, with
the star as the worst case at margin about 4/n. No known method proves it
below a threshold of about 10^495:

- Fang's architecture as sharpened: about 10^2070.
- The binomial-smoothing hybrid: about 10^495.
- Idealized constants still leave at least 10^66, and the true root-moment
  exponent probably forces at least 10^100.

The one deterministic candidate, R1, is refuted at n = 237. A proof that
reaches the enumeration record needs a genuinely different idea: an
effective CLT whose error does not pass through the root-moment exponent,
or a deterministic argument that survives hub-star constructions.

**Original resume plan (done except lclt-check):** confirm the MSH counterexample independently; run
`r1-exhaustive/summarize.py` or read `raw/` for the n <= 24 counts; rerun the
certificate check (`lclt-check`) if the fixed-accuracy route is still
wanted. Restarting the wave-2 workflow (`resumeFromRunId`) would replay all
four agents from scratch, because none had completed. If R1 is dead, the
deterministic route has no live target, and the remaining live item is the
fixed-accuracy smoothing route in `route-lclt`.

## Follow-up: tree-shift test (same day)

After the attack, Brett asked whether mathematicians have established
techniques for generating ideas. One extremal-graph-theory idea came out of
that discussion: push each tree toward the star with Csikvári's generalized
tree shift, and show the central margin never rises on the way. It was tested
exactly and is **weak, not pursued**:

- the one-step and two-step forms are false;
- some trees are stuck even though they come within 2% of the star's margin;
- hub constructions need up to 5 shifts, with no uniform bound in sight.

Details: `gts-test/RESULTS.md`.

## Follow-up: LP search for discharging certificates (same day)

This came from the second idea-generation technique, certificate search.
Details: `lp-discharge/RESULTS.md`.

- **Constant rules fit everything tested, at first.** An LP over 1.8M
  constraints found simple one-hop rules certifying window LC on all trees
  `n <= 16` and on every hub family, including the n = 237 R1 killer. With
  the constant rule "keep 2/3, give 1/3 equally to neighbours", every tree
  with `n <= 27` passes: about 1.2e9 trees exhaustively, 0 violations. The
  pointwise lemma also holds through n = 27, so its smallest counterexample
  is exactly n = 28.
- **Then the adversary refuted every constant rule.**
  - The 2/3 rule fails at the hub-star H(75,5), n = 451. I confirmed this
    independently.
  - A single tree with n = 1293 excludes every constant.
- **Degree-based rules survive the hub families only with zero slack.** The
  float LP optimum fails an exact check at the `1e-20` scale.
- **Verdict.** One-hop discharging is exhausted as a route. Small trees
  predict nothing here: every candidate rule holds through `n <= 27` and
  fails only in constructions with hundreds to thousands of vertices.

## Follow-up: radius check (same day)

Before any larger multi-hop search, a cheap exact test asked whether deeper
nesting pushes positive defect farther from negative defect. It doesn't:

- in 13,444 layered trees of depth 2–6, positive defect sits only at
  isolated depths next to negative ones;
- every positive vertex, including those in all the killers, can be paid by
  its own children and grandchildren with at least 5x surplus.

So one-hop rules fail through uniform shares, not distance. A future search
should use short-range transfers with structure-aware shares, not more hops.
Details: `radius-check/RESULTS.md`, which also records a float-overflow trap
in networkx max-flow.
