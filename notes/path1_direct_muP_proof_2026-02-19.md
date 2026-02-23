# Path 1 Direct `\mu_P` Proof via Sum-of-Means (2026-02-19)

## Goal

For the canonical degree-2 bridge decomposition

- pick leaf `l` by minimum support degree (tie: largest leaf id),
- support `s` with `deg(s)=2`,
- `u` the other neighbor of `s`,
- `B = T - {l,s}`,
- `P = dp_B[u][0]`,
- `m = mode(I(T))`,
- `lambda = lambda_m(T) = i_{m-1}(T)/i_m(T)`,

prove:

`mu_P(lambda) >= m - 2`.

Equivalent sum-of-means form (root `B` at `u`, children `c`):

`P = prod_c I(T_c)`, so

`mu_P(lambda) = sum_c mu_{T_c}(lambda)`.

Hence target is:

`sum_c mu_{T_c}(lambda_m(T)) >= m - 2`.

---

## 1) Exact algebraic bridge for Approach 3

Define:

- `a := lambda/(1+lambda)`,
- `D := mu_B(lambda) - mu_P(lambda)`,
- `exact_slack_B := mu_B(lambda) - (m - 1 - a)`,
- `exact_excess_D := D - (1 - a)`.

Then exactly:

`mu_P - (m-2) = exact_slack_B - exact_excess_D`.

Proof:

`mu_P - (m-2) = mu_B - D - (m-2)`

`= [mu_B - (m-1-a)] - [D - (1-a)]`

`= exact_slack_B - exact_excess_D`.

So a sufficient template is:

If `exact_slack_B >= sigma` and `exact_excess_D <= kappa`, then

`mu_P - (m-2) >= sigma - kappa`.

---

## 2) Verdict on the proposed chain

Question: can

- Route-2 exact: `mu_B >= m-1-a`,
- Transfer bound: `D <= 1-a+0.006`,

imply `mu_P >= m-2`?

Answer:

- **No**, not by itself. It only gives `mu_P >= m-2-0.006`.
- **Yes**, if Route-2 is strengthened by additive slack `sigma >= 0.006`.

From full `n<=23` data, we have much more:

- `min exact_slack_B = 0.1913484628930444`,
- `max exact_excess_D = 0.005107340588772491`,

thus

`mu_P - (m-2) >= 0.1913484628930444 - 0.005107340588772491`

`= 0.18624112230427192 > 0`.

So Approach 3 closes decisively on this frontier once the observed Route-2 slack is used.

---

## 3) Direct sum-of-means verification (new unified scan)

New script:

- `verify_muP_sum_of_means_2026_02_19.py`

Run:

```bash
python3 verify_muP_sum_of_means_2026_02_19.py \
  --min-n 4 --max-n 23 \
  --out results/verify_muP_sum_of_means_n23_2026_02_19.json \
  --no-resume
```

Result (d_leaf<=1, canonical leaf, full `n<=23`):

- seen: `23,942,356`
- considered/checked: `931,596`
- `mu_P >= m-2` failures: `0`
- sum-of-means identity failures (`mu_P` vs `sum_c mu_{T_c}`): `0`
- chain-identity failures: `0`
- partition-size failures (`sum_c |T_c| = n-3`): `0`
- transfer-cap violations (`exact_excess_D > 0.006`): `0`

Extrema:

- `min(mu_P-(m-2)) = 0.38345129375473075`
- `min exact_slack_B = 0.1913484628930444`
- `max exact_excess_D = 0.005107340588772491`
- `min chain gap = 0.38345129375473075`

Identity accuracy:

- `max |mu_P - sum_c mu_{T_c}| = 6.217248937900877e-15`
- `max |(mu_P-(m-2)) - (exact_slack_B-exact_excess_D)| = 4.440892098500626e-16`

Structural split:

- one-child-at-`u`: `590,779` trees, min gap `0.38345129375473075`
- multi-child-at-`u`: `340,817` trees, min gap `0.4882403086634657`

So the direct target is strongest in one-child cases, but still has large positive margin.

---

## 4) Artifact-only chain calibration script (new)

New script:

- `verify_path1_chain_from_existing_results_2026_02_19.py`

Run:

```bash
python3 verify_path1_chain_from_existing_results_2026_02_19.py \
  --out results/path1_direct_chain_from_existing_results_2026_02_19.json
```

This combines existing artifacts:

- transfer: `results/whnc_route1_transfer_scan_n23_merged.json`,
- pendant canonical: `results/whnc_pendant_bonus_scan_n20_canonical.json` + `results/whnc_pendant_bonus_scan_n23_canonical_tail.json`,
- direct `mu_P` cross-check: `results/whnc_phiP_scan_n23_merged.json`.

Derived bound from those artifacts:

- `lower_bound_from_chain = 0.18624112230427192`.

---

## 5) Status

- The requested inequality `mu_P(lambda_m(T)) >= m-2` is now verified in two independent ways on full `n<=23` frontier:
  1. direct `mu_P`/sum-of-means scan;
  2. transfer-chain calibration with positive net margin.
- Algebraically, Approach 3 is valid in the form:
  `mu_P-(m-2) = exact_slack_B - exact_excess_D`.
- The only missing piece for a fully symbolic universal proof is a non-computational bound ensuring `exact_slack_B >= exact_excess_D` (or directly `sum_c mu_{T_c}(lambda_m(T)) >= m-2`) for all sizes.
