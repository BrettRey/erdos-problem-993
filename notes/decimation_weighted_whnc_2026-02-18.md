# Decimation + Weighted WHNC (2026-02-18)

This note records the new “stat-phys decimation” analogy and what it buys.

## 1) Exact decimation reduction (proved computationally, exact algebra)

For `d_leaf<=1` tree `T`:

- `L` = leaves,
- `A` = support vertices adjacent to leaves,
- `C = V(T)\L` (leaf-stripped core).

Define inhomogeneous hard-core model on `C`:

- `lambda_v = 1/2` for `v in A`,
- `lambda_v = 1` for `v in C\A`.

At fugacity 1:

- `Z_T = 2^{|L|} Z_C^(lambda)`,
- `P_T(v)=P_C(v)` for `v in C`,
- if leaf `l` is attached to support `s`: `P_T(l) = (1-P_C(s))/2`.

Mean rewrite:

`mu(T) = |A|/2 + sum_{v in C\A} P_C(v) + (1/2) sum_{v in A} P_C(v).`

Equivalent gap form:

`n/3 - mu(T)`
`= (|C|/3 - |A|/6) - [sum_{v in C\A} P_C(v) + (1/2) sum_{v in A} P_C(v)].`

Verification artifact:

- `results/whnc_decimation_core_n21.json`
- Through `n<=21`: `z_fail=0`, `core_prob_fail=0`, `leaf_prob_fail=0`,
  `gap_identity_fail=0`.

## 2) Decimated weighted-WHNC (new candidate inequality)

On core `C`, define heavy set only on non-support vertices:

`H = {h in C\A : P_C(h) > 1/3}.`

Define weighted supply:

- `s(u)=1/3 - P_C(u)` for `u in C\A`,
- `s(u)= (1/2)(2/3 - P_C(u))` for `u in A`.

Candidate:

`sum_{h in H}(P_C(h)-1/3) <= sum_{u in N(H)} s(u).`

This is the decimated analogue of WHNC, with support vertices carrying half-weight.

Verification artifact:

- `results/whnc_decimated_whnc_n23.json`
- Through full `n<=23` (`931,596` d_leaf<=1 trees):
  - `global_fail=0`.
  - minimum global margin: `0.27615847319174214`.

So this global weighted inequality is numerically robust at full current frontier.

## 3) Local weighted overlap is close but false

Stronger local condition:

`sum_{h~u, h in H}(P_C(h)-1/3) <= s(u)` for each `u in N(H)`.

Through `n<=23`:

- `local_fail=54`.
- worst local margin: `-0.10158506746109458`.

Hence the right target is global Hall-type compensation, not per-vertex domination.

## 4) Decimated marginal obstruction disappears on current frontier

Define decimated removal marginal for non-empty `S ⊆ H`:

`M(h,S) = s(N_priv(h,S)) - (P_C(h)-1/3),`

`N_priv(h,S)=N(h)\N(S\{h}).`

A peeling proof would follow if every non-empty `S` had some `h` with
`M(h,S) >= 0`.
Equivalently, there should be no all-negative subsets (`M(h,S)<0` for all `h∈S`).

Verification artifact:

- `results/whnc_decimated_marginal_scan_n23.json`
- Through full `n<=23` (`931,596` d_leaf<=1 trees):
  - `trees_with_h = 674,393`,
  - `trees_with_allneg = 0`,
  - `allneg_subsets = 0`.

So in the decimated weighted model, the exact obstruction that appears in the
original WHNC scan (all-negative subsets) is absent on the full frontier.

An analytic proof attempt is now written in:

- `notes/decimated_peeling_proof_attempt_2026-02-18.md`

It gives a clean peeling theorem (private-neighbor existence + strict edge
surplus) proving decimated singleton-argmin/Hall directly.

## 5) Why this route is different and useful

- It integrates out leaves exactly (renormalization step), instead of handling
  leaf/support overlap directly in original graph.
- It moves the problem to a single inhomogeneous hard-core model on the core.
- It keeps a large positive global margin computationally, while isolating only
  a small nonlocal overlap obstruction.

This now looks like the strongest proof lane:

1. formalize decimation identities;
2. prove no-all-negative-marginal subsets in the decimated core model (peeling);
3. pull back to `mu(T)<n/3`.

## 6) Steiner peeling with gap-formula weights (updated exact scan)

Gap-formula weights:

- `s(u) = (1/2)(1/3-P(u))` for `u ∈ A`,
- `s(u) = 1/3-P(u)` for `u ∈ C\A`.

Exact full scan script:

- `conjecture_a_steiner_gap_scan.py`

Run:

```bash
python3 conjecture_a_steiner_gap_scan.py \
  --min-n 3 --max-n 23 \
  --out results/whnc_steiner_gap_scan_n23.json
```

Result through full `n<=23` (`931,596` d_leaf<=1 trees):

- `with_H=674,393`,
- non-empty subsets checked: `2,881,985`,
- non-singleton subsets checked for Steiner peeling: `1,580,936`,
- `fgap_fail=0` (`F_gap(S) >= 0` for every non-empty `S ⊆ H_core`),
- `steiner_fail=0` (every non-singleton `S` has a Steiner leaf with `M_gap>=0`),
- `fgap_zero=0`, `steiner_zero=0` (strict positivity on this frontier),
- minimum `F_gap(S)`: `0.1380801129345332`,
- minimum best-Steiner marginal: `0.0031308091184815007`.

This substantially upgrades the earlier partial BP check (`n<=18`, `|H_core|<=8`)
to an exact full-frontier verification.
