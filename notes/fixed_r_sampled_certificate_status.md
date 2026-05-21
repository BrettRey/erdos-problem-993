# Sampled Fixed-`r` Route-2 Certificate Status

## Scope

This note records the current certificate status for the sampled lanes

```text
S(2^a,r),    r in {4,8,12,16,20,24,32,40,60,80},
```

with a length-2 arm removed in the Route-2 bridge.

It does **not** prove the arbitrary fixed-`r` theorem.  It gives a fully
rerunnable certificate suite for the sampled lanes.

## Certificate Runner

Main driver:

```bash
python3 gpt_attack/fixed_r_sampled_certificate_suite.py \
  --out results/fixed_r_sampled_certificate_suite.json
```

Run result:

```text
all_ok=True
```

The output summary is saved in:

```text
results/fixed_r_sampled_certificate_suite.json
```

After correcting the hub-on perturbation monotonicity check to verify decrease
of `a * perturbation`, the suite was rerun:

```bash
PYTHONPATH=gpt_attack python3 -u gpt_attack/fixed_r_sampled_certificate_suite.py \
  --out results/fixed_r_sampled_certificate_suite_corrected_monotonicity.json
```

Result:

```text
all_ok=True
```

## What Is Certified

For each sampled lane, the proof split is:

```text
finite exact check: a <= 199
asymptotic certificate: a >= 200
```

The finite side uses exact cleared-denominator integer signs from:

```text
gpt_attack/fixed_r_finite_route2_check.py
```

The asymptotic side uses three certificates:

```text
gpt_attack/fixed_r_huboff_certificate.py
gpt_attack/fixed_r_hubon_mode_certificate.py
gpt_attack/fixed_r_hubon_route2_perturbation.py
```

For `r=8`, the strong denominators are used:

```text
hub-off margin >= 1/(10a)
hub-off reserve >= 1/(4a)
```

For the other sampled lanes, weak but sufficient denominators are used:

```text
hub-off margin >= 1/(1000a)
hub-off reserve >= 1/(1000a)
```

All sampled lanes also certify:

```text
3/4 <= lambda_0 <= 2
```

for the hub-off fugacity in the asymptotic range.

## Finite Results

Exact finite check:

```bash
python3 gpt_attack/fixed_r_finite_route2_check.py \
  --r-values selected \
  --a-max 199 \
  --out results/fixed_r_finite_route2_selected_a199.json
```

Result:

```text
checked records: 1990
skipped records: 0
Route-2 failures: 0
stronger-threshold failures: 0
```

Global finite minimum:

```text
r=8, a=198, n=405, m=135
Route-2 slack ~= 0.1684594347526992
stronger-threshold slack ~= 0.1682996084886373
lambda ~= 0.9993608992377847
```

## Asymptotic Results

For every sampled `r`, the suite reports:

```text
huboff_ok=True
hubon_mode_ok=True
hubon_route2_perturbation_ok=True
asymptotic_ok=True
```

The outer sampled lane `r=80` is the useful stress case.  Its asymptotic
certificates include:

```text
mode domination:
  max reported 2G/F <= 8.580e-34
  target ~= 5.000e-06

Route-2 perturbation:
  max reported mixture_bound ~= 8.368e-10
  reserve ~= 5.000e-06
```

So the sampled-lane asymptotic estimates have large numerical slack even at
the weakest reserve level.

## Theorem-Ready Statement

The current theorem-level claim supported by the scripts is:

```text
For r in {4,8,12,16,20,24,32,40,60,80} and every a >= 1,
the Route-2 inequality holds for S(2^a,r) after removing a length-2 arm:

    mu_B(lambda) >= m - 3/2,

where m is the leftmost mode of S(2^a,r) and
lambda = i_{m-1}(S(2^a,r)) / i_m(S(2^a,r)).
```

The finite checker also verifies the stronger empirical threshold

```text
mu_B(lambda) >= m - 1 - lambda/(1+lambda)
```

for `a <= 199` in the same sampled lanes.

## Remaining Gaps

1. This is a sampled fixed-`r` theorem, not a proof for all fixed `r`.
2. The asymptotic certificates use lane-specific stabilized shifts computed
   at threshold `a=200`; the next general theorem needs a symbolic or
   algorithmic statement that produces these shifts for arbitrary fixed `r`.
3. The proof should be written as a small family of lemmas:
   hub-off margin/reserve, hub-on mode domination, hub-on Route-2 perturbation,
   and finite verification.
4. The current manuscript has not been edited; this is still a certificate
   packet, not paper text.
