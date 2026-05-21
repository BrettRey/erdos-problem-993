# SG3 Route-2 Packet: Mode-Tie Bridge for `d_leaf <= 1` Trees

Use this packet for the next serious model run. It replaces the broad
"prove mode <= ceil(mu)" request with a narrower proof target that matches the
strongest existing diagnostics.

## Objective

Prove the focused tie-fugacity inequality for every tree `T` with
`d_leaf(v) <= 1`:

```text
mu_T(lambda_m) >= m - 1,
```

where:

- `m` is the leftmost mode of `I_T(x)` at fugacity `1`,
- `lambda_m = i_{m-1}(T) / i_m(T)`,
- `mu_T(lambda)` is the expected independent-set size under weights
  proportional to `i_k(T) lambda^k`.

This proves `mode(I(T)) <= ceil(mu_T(1))`. Combined with the proved mean bound
`mu_T(1) < n/3` for `d_leaf <= 1` trees, it proves Conjecture A in the reduced
regime.

## Known Bad Targets

Do not try to prove either of these:

- `mu(T) < n/3` for all trees.
- `mode(I(T)) <= floor(n/3)+1` for all trees.

Stars `K_{1,d}` disprove both for large `d`; run
`python3 gpt_attack/sanity_check.py` for exact examples.

## Decomposition

Choose a leaf `l` whose support vertex `s` has degree `2`, and let `u` be the
other neighbor of `s`.

Set:

```text
A = T - l
B = T - {l, s}
lambda = lambda_m(T)
Phi_q(F; lambda) = lambda I'_F(lambda) - (q-1) I_F(lambda)
```

Then:

```text
Phi_m(T; lambda) = Phi_m(A; lambda) + lambda Phi_{m-1}(B; lambda).
```

Therefore the focused tie inequality is proved if both bridge terms are
nonnegative:

```text
Phi_m(A; lambda) >= 0
Phi_{m-1}(B; lambda) >= 0.
```

## Existence Lemma

Every `d_leaf <= 1` tree with at least one leaf-support pair has a support
vertex of degree `2`.

Counting sketch:

- Let `L` be leaves and `S` be their supports.
- In a `d_leaf <= 1` tree, `|L| = |S|`.
- Let `C = V \ (L union S)`.
- If all supports had degree at least `3`, then
  `2(n-1) = sum_v deg(v) >= |L| + 3|S| + 2|C| = 2n`, contradiction.

This gives a canonical induction leaf.

## Route-2 Target

For the degree-2 support leaf above, write:

```text
I(B) = sum_k b_k x^k
tau = b_{m-2} / b_{m-1}
```

The empirical Route-2 target is:

```text
mu_B(lambda) >= m - 3/2.
```

A stronger target also held in scans:

```text
mu_B(lambda) >= m - 1 - lambda/(1+lambda).
```

## Exact Compensation Form

Define:

```text
deficit_tau = (m - 3/2) - mu_B(tau)
gain = mu_B(lambda) - mu_B(tau)
```

Route-2 is exactly:

```text
gain >= deficit_tau.
```

Since:

```text
d mu_B(t) / dt = Var_B(X_t) / t,
```

we have:

```text
gain = integral_tau^lambda Var_B(X_t)/t dt.
```

If `mu_B` is concave on `[tau, lambda]`, then:

```text
gain >= (lambda - tau) Var_B(X_lambda) / lambda.
```

Thus a sufficient endpoint inequality is:

```text
((lambda - tau) / lambda) * Var_B(X_lambda) >= deficit_tau.    (C1)
```

## Tie-Gap Formula

Let `P = dp_B[u][0] = sum_k p_k x^k`. The bridge coefficients give:

```text
i_{m-1}(T) = b_{m-1} + b_{m-2} + p_{m-2}
i_m(T)     = b_m     + b_{m-1} + p_{m-1}
```

So:

```text
lambda - tau =
[
  (b_{m-1}^2 - b_m b_{m-2})
  + (p_{m-2} b_{m-1} - p_{m-1} b_{m-2})
]
/
[
  b_{m-1}(b_m + b_{m-1} + p_{m-1})
].
```

The numerator is a log-concavity surplus plus a mismatch determinant.

## Verified Evidence

From `notes/mode_tie_leaf_bridge_2026-02-18.md` and
`notes/route2_compensation_algebraic_2026-02-19.md`, with an additional
all-degree-2-support endpoint scan on 2026-05-20:

- Degree-2 support bridge check through the full `d_leaf <= 1`, `n <= 23`
  frontier: `4,543,370` degree-2-support leaves, `0` failures.
- Focused mode-tie inequality through `n <= 23`: `931,596` trees, `0`
  failures, minimum margin `0.36225903819956606`, attained by `S(2^10)`.
- Route-2 canonical scan through `n <= 23`: `931,595` checked leaves, `0`
  failures.
- Endpoint slope check through `n <= 23`: minimum `M_lam =
  0.17929045776096725`.
- Endpoint slope check for every degree-2-support leaf through `n <= 23`:
  `4,543,368` checked leaves, `0` route-2 failures, `0` stronger-threshold
  failures, `0` nonpositive tie gaps in deficit cases. The minimum endpoint
  margin is `M_lam = 0.1757370039751993`, attained by the spider
  `S(1,2^7,4)`.
- Structured spider-lane scan of `S(1^j,2^a,r)` for `j <= 1`, `a <= 500`,
  `r <= 80`: minimum endpoint margin `0.164846737171075` at `S(2^12,8)`,
  removing a length-2 arm. Large-`a` route-2 slacks in this lane appear to
  decrease toward `1/6` from above.
- Fixed-`r` probe of `S(2^a,r)` through `a <= 1000` for
  `r in {2,3,4,5,6,7,8,9,10,12,16,20,40,80}`: the endpoint margin stays
  positive, and the route-2 slack in every tested lane approaches the same
  `1/6` binomial-spine asymptote from above. The finite endpoint dip remains
  `S(2^12,8)`.
- Pure-spider closed-form analysis for `S(2^a)`: exact verification for
  `2 <= a <= 250` gives no mode-formula failures and minimum Route-2 slack
  `0.169322709163`; stable float samples through `a=5000` match the
  binomial-spine prediction
  `Delta_a^0 = 3/2 - 2m_a/(a+1) -> 1/6` from above.
- The pure-spider mode formula is now proved in
  `notes/route2_pure_spider_asymptotic_2026-05-20.md` by the sign of
  `c_k(a)-c_{k-1}(a)`.
- Route-2 itself is proved for the pure spider lane `S(2^a)`, removing a
  length-2 arm: exact finite check for `2 <= a <= 29`, and for `a >= 30` the
  elementary envelope
  `|Delta_a-Delta_a^0| <= 4a^2/2^((2a-5)/3) + 4a(3/4)^(a-1) < 1/6`.
- Fixed-`r` mode probe for `S(2^a,r)`, `r <= 80`, `200 <= a <= 260`: every
  tested `r` has a stabilized residue-class shift
  `3m_{a,r}-2a = D_{r,a mod 3}`. A wider `50 <= a <= 260` probe found only
  three delayed stabilizations (`r=26,55,61`), all settled by `a >= 200`.
- First fixed-`r` symbolic target `S(2^a,8)`: exact formulas now exist for
  the stabilized boundary differences and for the hub-off reserve
  `Delta^0_{a,8}-1/6`. In all three residue classes the reserve is a rational
  function in `t` with positive numerator/denominator coefficients. Exact
  full Route-2 check through `a <= 80`: `0` failures, minimum slack
  `0.1710112823235669` at `a=78`.
- `S(2^a,8)` domination probe through `a <= 250`: proposed stabilized mode
  formula fails only in tiny cases (`a=6,9`), and all domination diagnostics
  pass from `a=30` onward. At `a=30`, the hub-on/hub-off tail ratio after the
  `G_a` peak is `6.355e-06` versus hub-off increasing margin `2.781e-03`.
- Conservative `S(2^a,8)` proof split: exact Route-2 finite check through
  `a <= 200` has `0` failures, min slack `0.168459434752` at `a=198`; for
  `a >= 200`, symbolic residue-class certificates prove hub-off margins
  `>=1/(10a)`, hub-off reserve `>=1/(4a)`, and `3/4<=lambda_0<=2`.
  The hub-on peak satisfies `G_k/F_m <= 918(2/3)^(a/3)`, and the two-sided
  coefficient comparison uses `|G_m-G_k| <= 2G_max`, so
  `2*918(2/3)^(a/3)<1/(10a)` gives mode localization. The remaining
  perturbation terms are bounded by
  `10000a^8/2^(a/2) + 10000a^3(3/4)^(a-1) < 1/(4a)`.
- Parameterized hub-off fixed-`r` certificate now works for selected lanes
  `r in {4,8,12,16,20,24,32,40,60,80}` at threshold `a >= 200`. The script now
  uses normalized coefficient ratios and shifted-first reserve polynomials;
  the old `r=20` symbolic-expansion bottleneck is gone. The certificate has
  also been strengthened to prove `3/4<=lambda_0<=2`.
- Parameterized fixed-`r` hub-on certificates now work for the same selected
  lanes through `r=80`: `fixed_r_hubon_mode_certificate.py` proves full-mode
  localization from `2 max G/F_m < 1/(Da)`, and
  `fixed_r_hubon_route2_perturbation.py` proves the crude full Route-2
  perturbation is below the hub-off reserve for `a >= 200`.
- Exact finite Route-2 checker `fixed_r_finite_route2_check.py` verifies the
  same selected lanes for `a <= 199` by cleared-denominator integer signs:
  `1990` records, `0` Route-2 failures, `0` stronger-threshold failures.
  The global finite minimum is again in the `r=8` lane at `a=198`, with
  Route-2 slack about `0.1684594347526992`.
- Consolidated runner `fixed_r_sampled_certificate_suite.py` reruns finite,
  hub-off, hub-on mode, and hub-on perturbation certificates for the sampled
  lanes and writes `results/fixed_r_sampled_certificate_suite.json`. Latest
  run: `all_ok=True`.
- Generalization probe `fixed_r_hubon_threshold_probe.py` shows the hub-on
  threshold must grow with `r`: in a coarse grid, first passing thresholds were
  `r=80 -> 200`, `r=100 -> 200`, `r=120 -> 250`, `r=160 -> 300`,
  `r=200 -> 350`, and then `r=240 -> 400`, `r=280 -> 450`,
  `r=320 -> 500`, `r=400 -> 650`. See
  `notes/fixed_r_generalization_probe_status.md`.
- Hub-off reserve construction was optimized using the path recurrence for
  `P_r(lambda)` and `lambda P'_r(lambda)`. With this optimization, the hub-off
  symbolic certificate now passes at `r=100, threshold=200` and
  `r=120, threshold=250`. A GMP-backed reserve helper
  `fixed_r_huboff_cpp_certificate.py` then certifies
  `r=160, threshold=300` and `r=200, threshold=350` in all three residue
  classes. After adding direct `mpz_addmul` accumulation and optional threaded
  multiplication, the helper also certifies `r=240, threshold=400` and
  `r=280, threshold=450` in all three residue classes. The next unrun
  stress point is `r=320, threshold=500`.
- The finite exact checker now evaluates
  `B(lambda)` and `lambda B'(lambda)` directly from
  `B=P_r(1+2x)^(a-1)+xP_{r-1}(1+x)^(a-1)` after clearing one common
  denominator. The larger finite batch
  `r in {100,120,160,200}`, `a <= 349` now completes: `1396` records,
  `0` Route-2 failures, `0` stronger-threshold failures, global minimum at
  `r=100, a=349` with Route-2 slack about `0.17438116131632425`.
  The next lane finite checks also pass: `r=240`, `a <= 399`, `399`
  records, global minimum slack about `0.18036271095442716`; `r=280`,
  `a <= 449`, `449` records, global minimum slack about `0.1808975747567671`;
  and `r=320`, `a <= 499`, `499` records, global minimum slack about
  `0.18125313899706172`. All have `0` Route-2 failures and `0`
  stronger-threshold failures.
- Sampled concavity of `mu_B` on the relevant intervals: no failures in tested
  deficit cases.
- The fixed-`r` work has been repackaged as a theorem-plan/risk-register note
  after reading the published unit-distance CoT:
  `notes/fixed_r_theorem_plan_from_unit_distance_cot_2026-05-21.md`.  Use it
  for the next proof-writing pass; it separates complete certified lanes from
  the still-missing arbitrary fixed-`r` theorem.
- The formal certificate implication has been drafted in
  `notes/fixed_r_certificate_lemma_2026-05-21.md`.  It packages the finite,
  mode, hub-off, and hub-on certificates into one Route-2 conclusion and marks
  the fugacity-shift sublemma as the next proof-writing target.
- The fugacity-shift sublemma has now been drafted in
  `notes/fixed_r_fugacity_shift_sublemma_2026-05-21.md`; it justifies the
  `2(a+r)^2T` and `(a+r)C_r(3/4)^(a-1)` perturbation terms used by
  `fixed_r_hubon_route2_perturbation.py`.  The next proof-writing gap is the
  global `F`-mode margin from checked adjacent margins.
- The global `F`-mode margin sublemma has now also been drafted in
  `notes/fixed_r_global_f_mode_margin_sublemma_2026-05-21.md`.  It uses
  real-rootedness/unimodality of `P_r(x)(1+2x)^a` to show that the adjacent
  margin checks imply the global margin needed for mode domination.
- The shifted-coefficient positivity principle and hub-off reserve recurrence
  certificate have now been drafted in
  `notes/fixed_r_shifted_coefficient_positivity_2026-05-21.md` and
  `notes/fixed_r_huboff_reserve_recurrence_certificate_2026-05-21.md`.  These
  close the remaining local proof-writing justifications for the current
  certificate-composition scaffold.
- Hub-off reserve asymptotics are now extracted by
  `gpt_attack/fixed_r_huboff_asymptotic.py` and recorded in
  `notes/fixed_r_huboff_asymptotic_reserve_2026-05-21.md`.  The expansion is
  `reserve = C_{r,q}/a + O_r(a^-2)`, and all tested constants through `r=400`
  are positive.  The coefficient is now reduced to a closed moment formula in
  the first three moments of the fixed path polynomial and the stabilized
  shift `D_{r,q}`.  The next analytic target is proving this expression is
  positive for the stabilized shifts and making the remainder effective.
- The stabilized-shift rule has been separated in
  `notes/fixed_r_shift_rule_from_path_moments_2026-05-21.md`.  The first-order
  F-mode inequalities give the candidate interval
  `3M_1-1 <= D <= 3M_1+2`, with `D == q mod 3`; the current script matches
  those candidates to the extracted shifts through `r=400`, and a boundary
  probe through `r=2000` found only `r=2,3`.  The note now proves this
  boundary classification from the path moment formula
  `P'_r(1)=(rL_{r+1}+2F_r)/5`.
- Positivity of the hub-off reserve constant has been reduced in
  `notes/fixed_r_c_positivity_reduction_2026-05-21.md`.  With
  `alpha=M_1-D_{r,q}/3`, the constant is
  `(9(V+3K_3)+24alpha+16)/12`; the shift interval makes the `alpha` term
  positive for `r>=4`, so the clean remaining target is the path inequality
  `V+3K_3>=0` or its equivalent root/trigonometric sum.  That note now proves
  the inequality by exact Fibonacci/Lucas formulas for the first three raw
  path moments.  Leading-constant positivity is no longer the bottleneck; the
  next bottleneck is an effective `O_r(a^-2)` remainder.
- The upgraded arbitrary fixed-`r` theorem schema is now in
  `notes/fixed_r_eventual_route2_theorem_schema_2026-05-21.md`.  It explains
  how the first-order shift rule, positive reserve constant, rational-function
  root isolation, and exponential hub-on domination combine to give a
  computable eventual threshold `A(r)` for each fixed `r`, leaving a finite
  exact check below `A(r)`.
- Effective threshold probing is implemented in
  `gpt_attack/fixed_r_effective_threshold_probe.py`.  Exact root isolation now
  gives compact eventual thresholds for the hub-off/mode/lambda pieces:
  `r=4 -> A=15`, `r=8 -> A=27`, `r=20 -> A=97`.  The script uses the cleared
  path recurrence for the reserve numerator.  Scaling note: exact isolation
  was slow by `r=24`; Cauchy mode completed `r=40` but produced a huge
  232-digit threshold, so it is a computability fallback rather than a useful
  certificate.
- The effective threshold probe now also has a shifted-coefficient mode, which
  is the practical route.  It gives hub-off/mode/lambda thresholds
  `r=4 -> 12`, `r=8 -> 24`, `r=20 -> 94`, `r=24 -> 22`, `r=40 -> 26`,
  `r=80 -> 44`, and `r=120 -> 62`.  The `r=160` shifted run was stopped before
  first output.
- The old hub-on perturbation monotonicity check has been corrected: because
  the reserve target is proportional to `1/a`, the scripts now check decrease
  of `a * perturbation`, not just decrease of the perturbation.  The old
  coarse thresholds still pass under the corrected check.
- Effective hub-on probing is implemented in
  `gpt_attack/fixed_r_hubon_effective_threshold_probe.py`.  Using half of the
  actual mode-margin constants and half of `C_{r,q}`, the first passing
  thresholds in a step-5 grid are `r=80 -> 145`, `r=120 -> 200`,
  `r=160 -> 255`, `r=200 -> 305`, `r=240 -> 360`, `r=280 -> 415`,
  `r=320 -> 470`, and `r=400 -> 575`.  Hub-on perturbation is now the practical
  threshold driver.

## What To Prove

Best next theorem:

```text
For every d_leaf <= 1 tree T, every degree-2-support leaf l, and the notation
above, endpoint condition (C1) holds.
```

Fallback theorem:

```text
For every such bridge, mu_B is concave on [tau, lambda], and lambda >= tau.
```

Fallback output if proof fails:

- Produce an explicit symbolic obstruction showing why `(C1)` cannot follow
  from only log-concavity of `B` and the determinant split.
- Or identify the exact additional tree-DP invariant needed to make `(C1)`
  true.

## Files To Read First

- `notes/mode_tie_leaf_bridge_2026-02-18.md`
- `notes/route2_compensation_algebraic_2026-02-19.md`
- `conjecture_a_mode_tie_deg2_support_scan.py`
- `verify_route2_slope_compensation.py`
- `gpt_attack/route2_exact_witness.py`
- `gpt_attack/route2_spider_lane_scan.py`
- `gpt_attack/pure_spider_route2_asymptotic.py`
- `gpt_attack/pure_spider_route2_bound.py`
- `gpt_attack/fixed_r_spider_mode_probe.py`
- `gpt_attack/fixed_r8_symbolic_certificate.py`
- `gpt_attack/fixed_r8_symbolic_inequality_certificate.py`
- `gpt_attack/fixed_r8_domination_probe.py`
- `gpt_attack/fixed_r8_asymptotic_bounds.py`
- `gpt_attack/fixed_r_huboff_certificate.py`
- `results/route2_slope_compensation_n23_merged.json`
- `results/route2_slope_compensation_all_deg2_n23.json`
- `results/route2_spider_lane_j0_1_a200_r80_float.json`
- `results/pure_spider_route2_asymptotic.json`
- `results/fixed_r_spider_mode_shifts_r80_a200_260.json`
- `notes/route2_spider_lane_2026-05-20.md`
- `notes/route2_pure_spider_asymptotic_2026-05-20.md`
- `notes/route2_fixed_r_spider_modes_2026-05-20.md`
