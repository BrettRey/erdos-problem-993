# Fixed-`r` Generalization Probe Status

## Purpose

After certifying the sampled lanes

```text
r in {4,8,12,16,20,24,32,40,60,80},
```

the next question is whether the same proof template plausibly extends to
arbitrary fixed `r`.

This note records the first scaling probe.  It is not an arbitrary-`r` proof.

## Hub-On Threshold Probe

Script:

```bash
python3 gpt_attack/fixed_r_hubon_threshold_probe.py \
  --r-values 80,100,120,160,200 \
  --start 200 \
  --stop 700 \
  --step 50 \
  --out results/fixed_r_hubon_threshold_probe_r80_200.json
```

Output:

```text
r=80: first_ok threshold=200 max_perturb=8.368e-10 min_reserve=4.950e-06
r=100: first_ok threshold=200 max_perturb=1.837e-06 min_reserve=4.950e-06
r=120: first_ok threshold=250 max_perturb=2.628e-09 min_reserve=3.968e-06
r=160: first_ok threshold=300 max_perturb=7.763e-09 min_reserve=3.311e-06
r=200: first_ok threshold=350 max_perturb=2.206e-08 min_reserve=2.841e-06
```

The next grid probe was then run:

```bash
python3 gpt_attack/fixed_r_hubon_threshold_probe.py \
  --r-values 240,280,320,400 \
  --start 350 \
  --stop 700 \
  --step 50 \
  --out results/fixed_r_hubon_threshold_probe_r240_400.json
```

Output:

```text
r=240: first_ok threshold=400 max_perturb=6.100e-08 min_reserve=2.488e-06
r=280: first_ok threshold=450 max_perturb=1.653e-07 min_reserve=2.212e-06
r=320: first_ok threshold=500 max_perturb=4.413e-07 min_reserve=1.992e-06
r=400: first_ok threshold=650 max_perturb=1.807e-12 min_reserve=1.534e-06
```

Interpretation:

```text
There is no universal threshold a=200 for all fixed r.
The hub-on side behaves as expected: the needed threshold grows with r.
The growth is mild in this probe, roughly linear at these values.
```

The failed check that confirms this:

```bash
python3 gpt_attack/fixed_r_hubon_route2_perturbation.py \
  --r 120 --threshold 200 --reserve-denom 1000
```

failed only because the hub-on mixture bound exceeded the weak reserve:

```text
max mixture_bound ~= 4.013e-03
reserve ~= 5.000e-06
```

The same lane passed at threshold `a=300`, and the threshold-grid probe found
the coarser first passing grid point `a=250`.

Update, 2026-05-21: the monotonicity condition in the hub-on perturbation
certificate has been corrected.  Since the comparison target is proportional
to `1/a`, the script now verifies decrease of `a` times each perturbation
term.  The old coarse thresholds above still pass under the corrected check:

```text
r=80 -> 200, r=120 -> 250, r=160 -> 300, r=200 -> 350,
r=240 -> 400, r=280 -> 450, r=320 -> 500, r=400 -> 650.
```

There is now also an effective-target probe using half of the actual
first-order mode-margin constant and half of the proved reserve constant
`C_{r,q}`:

```bash
python3 gpt_attack/fixed_r_hubon_effective_threshold_probe.py \
  --r-values 80,120,160,200,240,280,320,400 \
  --start 100 --stop 620 --step 5 \
  --out results/fixed_r_hubon_effective_threshold_probe_r80_400_step5.json
```

Output:

```text
r=80:  first_ok threshold=145
r=120: first_ok threshold=200
r=160: first_ok threshold=255
r=200: first_ok threshold=305
r=240: first_ok threshold=360
r=280: first_ok threshold=415
r=320: first_ok threshold=470
r=400: first_ok threshold=575
```

This confirms that the large old thresholds were partly caused by the weak
reserve target `1/(1000a)`, but hub-on perturbation remains the practical
threshold driver.

## Hub-Off Reserve Scaling

The hub-on certificates still pass for `r=100` at threshold `a=200`:

```text
mode domination max 2G/F <= 9.493e-34
Route-2 perturbation max ~= 1.837e-06
reserve ~= 5e-06
```

The original hub-off symbolic certificate stalled at `r=100`, because it
expanded `P_r(lambda_0)` and `lambda_0 P'_r(lambda_0)` by explicit powers.
This has been replaced by the path recurrences for:

```text
E_n = Z^deg(P_n) P_n(L/Z),
N_n = Z^deg(P_n) (L/Z) P'_n(L/Z).
```

After that optimization, the hub-off certificate passes:

```bash
python3 gpt_attack/fixed_r_huboff_certificate.py \
  --r 100 --threshold 200 --margin-denom 1000 --reserve-denom 1000
```

and also:

```bash
python3 gpt_attack/fixed_r_huboff_certificate.py \
  --r 120 --threshold 250 --margin-denom 1000 --reserve-denom 1000
```

Both certify the hub-off margins, `3/4 <= lambda_0 <= 2`, and the weak
hub-off reserve `>=1/(1000a)` in every residue class.

The original Python reserve recurrence still stalled at the next stress point,

```bash
python3 gpt_attack/fixed_r_huboff_certificate.py \
  --r 160 --threshold 300 --margin-denom 1000 --reserve-denom 1000
```

after passing the first residue's cheap margin/fugacity checks.  The reserve
coefficient arithmetic was then moved to the GMP helper:

```bash
python3 gpt_attack/fixed_r_huboff_cpp_certificate.py \
  --r 160 --threshold 300 --margin-denom 1000 --reserve-denom 1000
```

Output summary:

```text
stabilized shifts D_q=3m-2a: {0: 135, 1: 133, 2: 134}
all three residue classes:
  margins true
  3/4 <= lambda_0 <= 2 true
  reserve >= 1/(1000a) true
reserve polynomial degrees: 4779, 4779, 4698
```

The next stress point,

```bash
python3 gpt_attack/fixed_r_huboff_cpp_certificate.py \
  --r 200 --threshold 350 --margin-denom 1000 --reserve-denom 1000
```

also passes.

Output summary:

```text
stabilized shifts D_q=3m-2a: {0: 168, 1: 166, 2: 167}
all three residue classes:
  margins true
  3/4 <= lambda_0 <= 2 true
  reserve >= 1/(1000a) true
reserve polynomial degrees: 7373, 7373, 7373
```

The next attempted stress point,

```bash
python3 gpt_attack/fixed_r_huboff_cpp_certificate.py \
  --r 240 --threshold 400 --margin-denom 1000 --reserve-denom 1000
```

originally passed the first residue's cheap margin/fugacity checks but did not
finish the first reserve positivity check after about 17 minutes.  The GMP
helper was then changed to use direct `mpz_addmul` accumulation and optional
threaded multiplication:

```bash
python3 gpt_attack/fixed_r_huboff_cpp_certificate.py \
  --r 240 --threshold 400 --margin-denom 1000 --reserve-denom 1000 \
  --reserve-threads 8
```

Output summary:

```text
stabilized shifts D_q=3m-2a: {0: 201, 1: 199, 2: 200}
all three residue classes:
  margins true
  3/4 <= lambda_0 <= 2 true
  reserve >= 1/(1000a) true
reserve polynomial degrees: 10648, 10648, 10527
```

The next stress point,

```bash
python3 gpt_attack/fixed_r_huboff_cpp_certificate.py \
  --r 280 --threshold 450 --margin-denom 1000 --reserve-denom 1000 \
  --reserve-threads 8
```

also passes.

Output summary:

```text
stabilized shifts D_q=3m-2a: {0: 234, 1: 232, 2: 233}
all three residue classes:
  margins true
  3/4 <= lambda_0 <= 2 true
  reserve >= 1/(1000a) true
reserve polynomial degrees: 14382, 14382, 14382
```

Thus the hub-off certificate now reaches `r=280` at the hub-on threshold grid
points.  The next unrun stress point is `r=320, threshold=500`.

## Finite Exact Scaling

The finite checker was first optimized to use the direct coefficient formula

```text
I(S(2^a,r)) = P_r(x)(1+2x)^a + xP_{r-1}(x)(1+x)^a
```

instead of generic spider polynomial exponentiation.

It was then optimized again to evaluate the post-removal polynomial

```text
B(x) = P_r(x)(1+2x)^(a-1) + xP_{r-1}(x)(1+x)^(a-1)
```

and `xB'(x)` directly at `lambda=L/M`, clearing one common denominator.  This
replaces the old degree-by-degree cleared-denominator summation.

The selected-lane check was rerun successfully after this change:

```bash
python3 gpt_attack/fixed_r_finite_route2_check.py \
  --r-values selected \
  --a-max 199 \
  --out results/fixed_r_finite_route2_selected_a199.json
```

and reproduced:

```text
checked records: 1990
Route-2 failures: 0
stronger-threshold failures: 0
global min: r=8, a=198, Route-2 slack ~= 0.1684594347526992
```

The larger exact finite batch that previously stalled now completes:

```bash
python3 gpt_attack/fixed_r_finite_route2_check.py \
  --r-values 100,120,160,200 \
  --a-max 349 \
  --out results/fixed_r_finite_route2_r100_200_a349.json
```

Output:

```text
checked records: 1396
Route-2 failures: 0
stronger-threshold failures: 0
global min: r=100, a=349, Route-2 slack ~= 0.17438116131632425
```

This covers the finite side up to the current hub-on threshold grid point for
`r=200` (`a=350`).  The finite checker is no longer the immediate bottleneck
for the `r <= 200` probe.

For the next lanes, the finite side also passes:

```bash
python3 gpt_attack/fixed_r_finite_route2_check.py \
  --r-values 240 \
  --a-max 399 \
  --out results/fixed_r_finite_route2_r240_a399.json
```

Output:

```text
checked records: 399
Route-2 failures: 0
stronger-threshold failures: 0
global min: r=240, a=399, Route-2 slack ~= 0.18036271095442716
```

and:

```bash
python3 gpt_attack/fixed_r_finite_route2_check.py \
  --r-values 280 \
  --a-max 449 \
  --out results/fixed_r_finite_route2_r280_a449.json
```

Output:

```text
checked records: 449
Route-2 failures: 0
stronger-threshold failures: 0
global min: r=280, a=447, Route-2 slack ~= 0.1808975747567671
```

and:

```bash
python3 gpt_attack/fixed_r_finite_route2_check.py \
  --r-values 320 \
  --a-max 499 \
  --out results/fixed_r_finite_route2_r320_a499.json
```

Output:

```text
checked records: 499
Route-2 failures: 0
stronger-threshold failures: 0
global min: r=320, a=498, Route-2 slack ~= 0.18125313899706172
```

## Current Status

What now looks solid:

```text
For each fixed r, the hub-on part should be controllable once a is large
enough, with a threshold growing with r.
```

The hub-off reserve also continues beyond the original sampled packet: it is
now certified through `r=280` at the hub-on threshold grid points.  The
`r=200`, `r=240`, and `r=280` lanes are now complete under this certificate
split.

The next lane, `r=320`, has finite exact coverage through `a <= 499` and
hub-on certificates for `a >= 500`; the missing piece is the hub-off reserve at
`a >= 500`.

What remains to make this a genuine arbitrary fixed-`r` theorem:

1. A scalable hub-off reserve certificate beyond the current `r=280` frontier.
2. A theorem or certified algorithm that produces stabilized residue shifts
   `D_{r,q}` and a threshold `A(r)`.
3. Larger finite exact coverage once `A(r)` is known for larger `r`; the
   direct-evaluation checker is adequate for the current `r <= 200` probe.
