# Signed Ratio-Drop Breaker Pass
Date: 2026-07-04

## Purpose

This is a breaker-track note for issue #5. It records an adversarial search
against the signed Poisson-binomial sufficient diagnostic

```text
Delta_eff = 1 - (c_{D+1}/c_D)/(c_D/c_{D-1})
```

and against the raw reserve

```text
1 - c_{D+1}/c_D
```

for signed laws `Z = X - Y`, where `X` and `Y` are independent
low-probability Poisson-binomial sums.

This is computational evidence only. It is not a proof of the signed reserve,
the effective-drop diagnostic, or the hub-bouquet reserve.

## Command

```bash
python3 scripts/signed_ratio_drop_breaker.py \
  --random-samples 2000 \
  --max-side-n 1600 \
  --de-maxiter 16 \
  --de-popsize 5 \
  --de-restarts 1 \
  --out results/signed_ratio_drop_breaker_extended_2026-07-04.json
```

The harness combines:

- finite Skellam/Poisson approximants near the sparse one-sided boundary;
- asymmetric near-one-sided perturbations;
- boundary-heavy components with `p=1/2` blocks;
- larger random grouped supports;
- differential-evolution plus local polishing on grouped patterns, optimizing
  both `V * Delta_eff` and `V * reserve`.

## Result

The extended run processed:

```text
attempted = 39747
analyzed = 39233
feasible with V >= 1 = 38417
effective-drop failures below 1/4 = 0
raw-reserve failures below 1/4 = 0
```

The best effective-drop row was:

```text
V = 1.0000000005267604
D = 2
c_D/c_{D-1} = 0.5002912510712524
c_{D+1}/c_D = 0.33325534793272427
Delta_eff = 0.3338773220216449
V * Delta_eff = 0.33387732219751826
V * reserve = 0.6667446524184905
side variance fraction = 0.00021623535870586368
X = Bin(1200, 0.0008338484408622789)
Y = Bin(20, 0.000010811884837842142)
```

This improves the previous optimizer minimum `0.336426` and pushes the
effective-drop boundary closer to the one-sided Poisson value `1/3`.

The best raw-reserve row was:

```text
V = 1.0000001546920794
D = 1
c_D/c_{D-1} = 0.9999995363624967
c_{D+1}/c_D = 0.49944385537668246
reserve = 0.5005561446233175
V * reserve = 0.5005562220553884
V * Delta_eff = 0.5005559904943431
side variance fraction = 0.002663617496161419
X = Bin(500, 0.0019986677463279664)
Y = Bin(500, 0.0000053272641961477115)
```

Balanced-side pressure remained larger. The best observed `V * Delta_eff`
under side-variance fraction floors was:

| floor | best `V * Delta_eff` |
|---:|---:|
| `0` | `0.3338773222` |
| `0.02` | `0.3424558781` |
| `0.05` | `0.3571608070` |
| `0.10` | `0.3849490902` |
| `0.25` | `0.4550086785` |

## Interpretation

No counterexample to the working `1/4` effective-drop diagnostic or raw reserve
was found. The best rows are still near-one-sided sparse Poisson approximants,
not genuinely two-sided examples.

The breaker pass does change calibration:

1. The effective-drop route should not be stated with a constant near `1/3`.
   The observed best row is already at `0.333877`, and the limiting one-sided
   Poisson calculation gives `1/3`.
2. The raw-reserve route retains more slack; its observed boundary remains near
   `1/2`.
3. The two-regime proof split is still the right target: handle the nearly
   one-sided sparse boundary by perturbation, and use smoothing/local-ratio
   arguments only for the side-balanced regime.

Issue #5 remains open.
