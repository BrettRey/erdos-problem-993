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

## 4) Why this route is different and useful

- It integrates out leaves exactly (renormalization step), instead of handling
  leaf/support overlap directly in original graph.
- It moves the problem to a single inhomogeneous hard-core model on the core.
- It keeps a large positive global margin computationally, while isolating only
  a small nonlocal overlap obstruction.

This looks like a viable proof lane:

1. formalize decimation identities;
2. prove weighted Hall inequality on the core for `H ⊆ C\A`;
3. pull back to `mu(T)<n/3`.
