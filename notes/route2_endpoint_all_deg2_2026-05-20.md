# Route-2 Endpoint Scan for All Degree-2-Support Leaves (2026-05-20)

## Purpose

The previous Route-2 slope artifact checked the endpoint sufficient condition
for one canonical degree-2-support leaf per `d_leaf <= 1` tree through `n <= 23`.
The SG3 proof packet proposes a stronger theorem for every degree-2-support
leaf. This note records the matching all-leaf scan.

## Command

```bash
python3 verify_route2_slope_compensation.py \
  --max-n 23 \
  --all-deg2-leaves \
  --out /tmp/route2_all_deg2_n23.json \
  --no-resume

cp /tmp/route2_all_deg2_n23.json \
  results/route2_slope_compensation_all_deg2_n23.json
```

## Result

Scope: all `d_leaf <= 1` trees through `n <= 23`, all leaves whose support has
degree `2`.

- Trees seen: `23,942,356`
- `d_leaf <= 1` trees considered: `931,596`
- Checked trees: `931,595`
- Checked degree-2-support leaves: `4,543,368`
- Route-2 failures: `0`
- Stronger threshold failures: `0`
- Deficit-at-`tau` cases: `2,987,519`
- Nonpositive tie gaps in deficit cases: `0`
- Minimum route-2 slack: `0.20033413665157518`
- Minimum stronger-threshold slack: `0.1913484628930444`
- Minimum endpoint margin at `lambda`: `0.1757370039751993`
- Maximum deficit at `tau`: `0.10009595509631719`

Artifact:

- `results/route2_slope_compensation_all_deg2_n23.json`

## Extremal Witnesses

Minimum all-leaf endpoint margin:

- `n = 20`
- graph6: `S???????C?G?G?C?@??G??_?@?C@?B~_?`
- spider signature: `S(1,2,2,2,2,2,2,2,4) = S(1,2^7,4)`
- `m = 7`
- `lambda_m = 8729/8857`
- `tau = 2782/3255`
- `M_lam = 0.1757370039751993`

Maximum deficit at `tau`:

- `n = 20`
- graph6: `S???????C?G?G?C?@??G??_?@??@?F~_?`
- spider signature: `S(1,2^9)`
- `m = 7`
- `lambda_m = 1589/1678`
- `tau = 213/260`
- `deficit_tau = 0.10009595509631719`
- `M_lam = 0.17929045776096725`

Minimum stronger-threshold slack:

- `n = 21`
- graph6: `T???????C?G?G?C?@??G??_?@??@???_B~o?`
- spider signature: `S(2^10)`
- `m = 7`
- `lambda_m = 2282/2595`
- `tau = 99/131`
- exact slack: `0.1913484628930444`

Exact rational values for these witnesses are printed by:

```bash
python3 gpt_attack/route2_exact_witness.py
```

## Interpretation

The endpoint inequality is now supported on the same full `d_leaf <= 1`,
`n <= 23` frontier and the same all-degree-2-support-leaf scope as the bridge
term checks. The tight cases are spiders or near-balanced spiders, so a proof
should first settle the spider endpoint inequality explicitly, then explain why
non-spider branching has more slack.

The remaining proof obligations are:

1. Prove concavity of `mu_B(lambda)` on the relevant interval `[tau, lambda_m(T)]`
   or replace it with a weaker monotonicity of `Var_B(lambda)/lambda`.
2. Prove the endpoint inequality
   `((lambda - tau)/lambda) Var_B(X_lambda) >= (m - 3/2) - mu_B(tau)`.
3. Connect the Route-2 endpoint inequality to the two nonnegative bridge terms
   in the focused tie-fugacity proof.
