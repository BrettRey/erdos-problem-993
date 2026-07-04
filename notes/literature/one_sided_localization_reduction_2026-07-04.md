# One-Sided Localization Reduction
Date: 2026-07-04

## Purpose

This note narrows the one-sided effective-drop route. It does not prove the signed reserve lemma. It reduces the remaining one-sided `1/4` target to a small elementary symmetric-polynomial inequality; that local-mode mean inequality is proved in the later proof note cited below.

## Starting Point

For a one-sided low-probability Poisson-binomial law, the previous note gives:

```text
Delta_eff >= 1/(D+1).
```

Thus `Delta_eff >= 1/(4V)` follows from

```text
D + 1 <= 4V.
```

The crude Darroch/Newton chain already gives

```text
D+1 <= 2V+3.
```

Therefore, if `V >= 3/2`, then

```text
D+1 <= 2V+3 <= 4V.
```

So the only possible gap is

```text
1 <= V < 3/2.
```

In that band, the crude chain gives

```text
D+1 < 6,
```

hence `D <= 4`. If `D <= 3`, then `D+1 <= 4 <= 4V`. Therefore the only remaining one-sided case is:

```text
D = 4.
```

## Elementary D=4 Target

Write odds

```text
w_i = p_i/(1-p_i),     0 <= w_i <= 1.
```

The coefficient sequence is proportional to the elementary symmetric sequence

```text
e_k(w_1,...,w_n).
```

If the first strict descent is at `D=4`, then

```text
e_1 >= e_0 = 1,
e_2 >= e_1,
e_3 >= e_2,
e_4 <  e_3.
```

Thus the one-sided localization would follow from the following elementary inequality:

> If `0 <= w_i <= 1` and
>
> ```text
> e_1 >= e_0 = 1,
> e_2 >= e_1,
> e_3 >= e_2,
> ```
>
> then
>
> ```text
> sum_i w_i/(1+w_i)^2 >= 5/4.
> ```

The right-hand side is exactly `V >= (D+1)/4` for `D=4`.

This is slightly stronger than the exact need, since it does not use the condition `e_4 < e_3`. The narrower D=4 statement is:

> If `0 <= w_i <= 1`, the first strict descent of the elementary symmetric sequence occurs at `D=4`, and the support continues far enough for the post-descent ratio to be relevant, then `V >= 5/4`.

The tempting shortcut

```text
e_3 >= e_2  =>  V >= 5/4
```

is too weak as stated: it admits irrelevant low-support equality cases such as one nonzero weight, where `e_2=e_3=0`. The actual reduction needs the pre-D4 guard, not merely the necessary condition `e_3 >= e_2`.

However, after excluding the degenerate case `e_2=0`, the pre-D4 side conditions are automatic. Let `m` be the number of nonzero weights and `S=e_1`. If `e_3 >= e_2 > 0`, then

```text
e_3 = (1/3) sum_{|A|=2} w_A sum_{j notin A} w_j <= (S/3)e_2,
```

so `S >= 3`, hence `e_1 >= 1`. Also `e_3/e_2 <= (m-2)/3`, so `m >= 5`. Newton's inequalities for the nonzero weights give

```text
e_3/e_2 <= [2(m-2)/(3(m-1))] e_2/e_1.
```

Since `e_3/e_2 >= 1`, this implies

```text
e_2/e_1 >= 3(m-1)/(2(m-2)) > 1.
```

Thus the remaining mean lemma can be stated more cleanly as:

> If `0 <= w_i <= 1`, `e_2 > 0`, and `e_3 >= e_2`, then `mu=sum_i w_i/(1+w_i) >= 5/2`.

Equivalently, for a Poisson-binomial law `S=sum_i Bernoulli(p_i)` with all `p_i <= 1/2`, the target is:

> If `P(S=3) >= P(S=2) > 0`, then `E S >= 5/2`.

## Probe

I added:

```bash
python3 scripts/probe_d4_variance_bound.py \
  --out results/d4_variance_bound_probe_2026-07-04.json
```

The probe samples random odds vectors and the structured family consisting of `m` ones plus one remainder. It checks the D=4 implication above and separately records failures of the over-broad `e_3 >= e_2` shortcut. This is only a falsification probe.

Default run:

```text
processed = 159080
weak_guard_e3_ge_e2_rows = 138120
weak_guard_failures = 1001
weak_positive_guard_cases = 137119
weak_positive_guard_failures = 0
pre_d4_guard_cases = 137119
pre_d4_guard_failures = 0
d4_cases = 7569
d4_failures = 0
```

The `1001` weak-guard failures are the expected irrelevant low-support rows with `e_2=e_3=0`. They show why the raw shortcut should not be used as the proof target. The nondegenerate weak-positive guard and the pre-D4 guard had identical counts in this corpus, no sampled failures, and the same five-fair-weight boundary below.

The same boundary appears if one minimizes the mean

```text
mu = sum_i w_i/(1+w_i)
```

under the nondegenerate local-mode guard:

```text
best sampled mean = 5/2.
```

This gives a plausible proof subroute. Since `w_i <= 1` is equivalent to `p_i <= 1/2`,

```text
V = sum_i p_i(1-p_i) >= (1/2) sum_i p_i = mu/2.
```

Therefore the elementary variance target would follow from the mean target:

> If `0 <= w_i <= 1` and `e_1 >= 1`, `e_2 >= e_1`, `e_3 >= e_2`, then `mu >= 5/2`.

Equivalently, using the automatic-side-condition reduction above:

> If `0 <= w_i <= 1`, `e_2 > 0`, and `e_3 >= e_2`, then `mu >= 5/2`.

This mean inequality is proved in:

```text
notes/literature/local_mode_mean_bound_proof_2026-07-04.md
```

The proof uses a deletion lemma plus a size-biased random-pair argument. It closes this elementary target, not the signed reserve target.

I also added an adversarial optimizer:

```bash
python3 scripts/optimize_pre_d4_mean_bound.py \
  --min-n 3 --max-n 12 --restarts 4 --maxiter 600 --popsize 12 \
  --out results/pre_d4_mean_optimizer_2026-07-04.json
```

Default run:

```text
processed = 48
feasible = 36
failures = 0
```

The optimizer uses a penalty formulation, so it is only a diagnostic. Exact structured boundary rows were included separately. The best rows are the boundary `w_1=...=w_5=1` plus zeros, with `mu=5/2` and `V=5/4`; optimizer-only rows converge to the same boundary from nearby feasible points.

The expected boundary example is:

```text
w_1 = ... = w_5 = 1,
```

corresponding to `Bin(5,1/2)`. Here

```text
e_3 = e_2 = 10,
mu = 5/2,
V = 5/4.
```

## Current Status

The one-sided proof route is now:

1. The elementary mean lemma is proved in `local_mode_mean_bound_proof_2026-07-04.md`.
2. Therefore the remaining `D=4` case has `V >= 5/4`, hence `D+1=5 <= 4V`.
3. The crude Darroch/Newton chain already handled `V >= 3/2` and the `D <= 3` cases.
4. Thus the sufficient one-sided localization `D+1 <= 4V` is proved for the route considered here.
5. Combining with Newton gives the one-sided effective-drop bound `Delta_eff >= 1/(4V)`.
6. The signed case remains open: conditional averaging and boundary terms still have to be handled.

This closes the one-sided localization gap for this proof route. It does not prove the full signed reserve needed for issue #5.
