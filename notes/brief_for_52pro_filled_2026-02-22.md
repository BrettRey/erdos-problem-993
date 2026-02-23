# Filled Proof Packet for 5.2 Pro (2026-02-22)

Use this as the complete replacement for the five placeholder blocks in the prior prompt.

## [PASTE: Setup + notation]
- Tree polynomial: `I(T;x) = sum_k i_k(T) x^k`.
- Mode index: `m = mode(I(T))` (leftmost mode convention in current scripts).
- Canonical bridge setup:
  - choose a leaf `l` whose support `s` has minimum degree;
  - require `deg(s)=2`; let `u` be the other neighbor of `s`;
  - define `B = T - {l,s}`;
  - rooted-at-`u` polynomials in `B`: `P = dp_B[u][0]`, `Q = dp_B[u][1]`, `I(B)=P+Q`.
- Coefficient notation around the target index:
  - `p0 = p_{m-2}`, `p1 = p_{m-1}`, `pm = p_m`;
  - `q0 = q_{m-2}`, `q1 = q_{m-1}`, `qm = q_m`;
  - `b0 = b_{m-2}`, `b1 = b_{m-1}`, `b2 = b_m`.
- Mode fugacity:
  - `lambda = lambda_m(T) = i_{m-1}(T)/i_m(T)`,
  - `a = lambda/(1+lambda)`.
- Means:
  - `mu_P(lambda)`, `mu_B(lambda)` are coefficient-distribution means at fugacity `lambda`.
- Route-1 transfer quantities:
  - `D := mu_B(lambda) - mu_P(lambda)`,
  - `exact_slack_B := mu_B(lambda) - (m - 1 - a)`,
  - `exact_excess_D := D - (1 - a)`.
- STRONG C2 split quantities:
  - `mismatch := p0*b1 - p1*b0`,
  - `lc_surplus := b1^2 - b2*b0`,
  - `combined := lc_surplus + mismatch`.
- Hard-regime symbols:
  - `db := b1-b0`,
  - `dq := q1-q0`,
  - `neg := -mismatch`,
  - `rise := b1*(b1-b0)`.

## [PASTE: Open targets]
1. Route 1 target:
   - `mu_P(lambda_m(T)) >= m-2`
   - equivalent calibrated form:
     - `exact_slack_B >= exact_excess_D`.
2. STRONG C2 target:
   - `a_{m-1} b_{m-1} - a_m b_{m-2} >= 0`
   - equivalent split form:
     - `lc_surplus + mismatch >= 0`.
3. Hard-regime reduction target (when `dq<0`):
   - `(-dq/db) <= p1/b1`.

## [PASTE: Proved identities/reductions]
Treat these as trusted givens:

1. Route-1 exact chain identity:
   - `mu_P - (m-2) = exact_slack_B - exact_excess_D`.
2. STRONG C2 determinant split:
   - `a_{m-1} b_{m-1} - a_m b_{m-2} = lc_surplus + mismatch`.
3. Rise-compensation identity:
   - `rise - neg = p1*db + b1*dq`.
4. Combined decomposition:
   - `lc_surplus + mismatch = (rise-neg) + b0*(b1-b2)`.
5. Hard-regime ratio equivalence:
   - if `dq<0` and `db>0`, then `rise-neg >= 0` iff `(-dq/db) <= p1/b1`.
6. Alternative LC decomposition:
   - `combined = lc_P + lc_Q + R`,
   - `R = cross + mismatch`,
   - `cross = 2*p1*q1 - pm*q0 - p0*qm`,
   - `mismatch = p0*q1 - p1*q0`.

## [PASTE: Computational certificates + min margins]
All on canonical `d_leaf<=1` frontier through `n<=23` unless stated.

1. Route 1 (`results/verify_muP_sum_of_means_n23_2026_02_19.json`):
   - checked trees: `931,596`
   - failures:
     - `muP_fail = 0`
     - `sum_identity_fail = 0`
     - `chain_identity_fail = 0`
     - `transfer_cap_violation = 0`
   - tight margins:
     - `min(mu_P-(m-2)) = 0.38345129375473075`
     - `min exact_slack_B = 0.1913484628930444`
     - `max exact_excess_D = 0.005107340588772491`
     - `min chain gap = 0.38345129375473075`
2. STRONG C2 split/hard regime (`results/verify_strong_c2_rise_identity_2026_02_19.json`):
   - checked trees: `931,596`
   - failures:
     - `combined_neg = 0`
     - `rise_fail = 0`
     - `hard_ratio_fail = 0`
   - profile:
     - `mismatch_neg = 129`
     - `hard_regime = 2`
   - tight margins:
     - `min combined = 9`
     - `min hard ratio gap = 0.9016058072627923`
     - `max hard transfer ratio = 0.026770775237032907`
     - `min hard need ratio = 0.9283765824998251`
3. Route-B decomposition (`results/verify_strong_c2_route_b_pair_bounds_2026_02_19_n23.json`):
   - checked trees: `931,595`
   - key facts:
     - `R_neg = 0`, `R_min = 4`
     - `cross_neg = 0`, `cross_min = 3`
     - `cross_minus_abs_mismatch_neg = 0`, `cross_minus_abs_mismatch_min = 2`
   - known non-sign-definite subterms:
     - `D1_neg = 1` with `D1_min = -72576`
     - `D2_neg = 14283` with `D2_min = -835354`
     - `D1_plus_D3_neg = 7` with `D1_plus_D3_min = -628992`
4. Unified STRONG C2 certificate (`results/prove_strong_c2_unified_n23.json`):
   - `combined_neg = 0`
   - `R_neg = 0`
   - `rise_neg_fail = 0`
   - `hard_ratio_fail = 0`
   - `min_p_minus_q = 1`
5. Route-1 exact-transfer failure pack (`results/route1_transfer_failures_all4_2026_02_19.json`):
   - exact-transfer counterexamples count: `4`
   - worst excess: `max_exact_excess = 0.005107340588772491`
6. Route-C support artifact (`results/prove_strong_c2_p_dominance_n20.json`, `n<=20`):
   - `p_ge_q = 77141/77141`
   - `min_gap = 1`

## [PASTE: Dead ends / invalid routes]
Do not propose any route that depends on these:

1. `mode(P') >= m-1` as a universal claim (false; explicit counterexample in Route-C note).
2. TP2 multiplicative closure for `>=3` factors (fails; do not use as proof engine).
3. Exact transfer inequality `D <= 1/(1+lambda)` as universally true (false; 4 witness failures with small positive excess).
4. Claim `D1 = p1*q1 - pm*q0` is always nonnegative (false; one explicit witness).
5. Route suggestions that are only brute-force reruns of already completed full-frontier scans.

