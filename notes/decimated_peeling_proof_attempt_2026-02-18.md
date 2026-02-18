# Decimated Peeling Proof Attempt (2026-02-18)

This note records a concrete analytic attempt in the decimated core model.

## Setup

For a `d_leaf<=1` tree `T`, decimate leaves:

- `L` = leaves,
- `A` = supports (each support has exactly one leaf),
- `C = V(T)\L` (core tree),
- activities on `C`: `lambda_v=1` on `C\A`, `lambda_v=1/2` on `A`.

Let `P(v)` be core marginals under this weighted model.
Define:

- `H = {h in C\A : P(h) > 1/3}`,
- demand `d(h) = P(h) - 1/3` for `h in H`,
- supply
  - `s(u) = 1/3 - P(u)` for `u in C\A`,
  - `s(u) = 1/3 - P(u)/2` for `u in A`.

For non-empty `S ⊆ H`:

- `F(S) = s(N(S)) - d(S)`,
- `M(h,S) = F(S)-F(S\{h}) = s(N_priv(h,S)) - d(h)`.

## Lemma 1 (Heavy-core structure)

1. `H` is independent in `C`.
2. Every `h in H` has `deg_C(h) >= 2`.

Reason:

- If `xy` is any edge in `C`, it is also an edge in `T`, and Theorem 6 gives
  `P_T(x)+P_T(y)<2/3`. Since `P_T=P_C` on core vertices, two endpoints cannot
  both exceed `1/3`.
- `h in C` means `h` is not an original leaf, so `deg_T(h)>=2`. Also `h notin A`
  means `h` has no deleted leaf neighbor, hence `deg_C(h)=deg_T(h)>=2`.

## Lemma 2 (Private-neighbor existence on every non-empty subset)

For every non-empty `S ⊆ H`, some `h in S` has a private neighbor:
`N_priv(h,S) != empty`.

Proof sketch:

- Build bipartite graph `B_S` with vertex set `S ∪ N(S)` and only the
  `S`-to-`N(S)` edges from `C`.
- `B_S` is a forest (subgraph of a tree).
- By Lemma 1 and independence, each `h in S` has all its neighbors in `N(S)`,
  and `deg_{B_S}(h)=deg_C(h)>=2`.
- If no `h` had a private neighbor, then every `u in N(S)` would satisfy
  `deg_{B_S}(u)>=2` (a degree-1 `u` is private to its unique `h`).
- So every vertex of the finite forest `B_S` would have degree at least 2,
  impossible.

Hence a private-neighbor vertex always exists.

## Lemma 3 (Private-neighbor marginal is strictly positive)

For any `h in H` and any neighbor `u~h` in `C`:

`s(u) - d(h) > 0`.

Indeed:

- If `u in C\A`, then
  `s(u)-d(h) = 2/3 - (P(u)+P(h)) > 0` (Theorem 6 on edge `uh`).
- If `u in A`, then
  `s(u)-d(h) = 2/3 - P(h) - P(u)/2 > 2/3 - (P(h)+P(u)) > 0`.

## Theorem (Decimated singleton/Hall by peeling)

For every non-empty `S ⊆ H`, there is `h in S` with `M(h,S) > 0`.
Consequently:

1. `F(S) > min_{x in S} F({x})`,
2. `min_{S!=empty} F(S) = min_{h in H} F({h})`,
3. `F(S) > 0` for every non-empty `S`,
4. in particular `F(H) > 0` (strict decimated weighted-WHNC).

Proof:

- Lemma 2 gives `h` with a private neighbor `u`.
- Then `M(h,S) >= s(u)-d(h) > 0` by Lemma 3.
- Iterative peeling strictly decreases set size and strictly decreases `F`,
  ending at a singleton.
- Singleton positivity: for any `h in H`,
  `F({h}) = s(N(h)) - d(h) >= s(u)-d(h) > 0` for any `u~h`.

## Bridge to `n/3 - mu(T)` (partial)

Exact identity (verified numerically):

`n/3 - mu(T)`
`= sum_{v in (C\A)\(H U N(H))} (1/3 - P(v))`
`+ sum_{u in A\N(H)} (1/6 - P(u)/2)`
`+ [F(H) - |A cap N(H)|/6].`

So the remaining analytic gap is not decimated Hall itself, but the
`-|A cap N(H)|/6` support penalty term.

## New computational checks from this pass

All checks below were run in-session through full `n<=23` (`931,596` `d_leaf<=1` trees):

1. Private-neighbor property on all non-empty `S ⊆ H`:
   - `trees_with_H = 674,393`
   - `subsets_checked = 2,881,985`
   - failures: `0`.
2. Strict edge marginal for decimated supplies:
   - evaluated heavy-core edge incidences: `2,603,599`
   - minimum `s(u)-d(h) = 0.021903464487004864` (strictly positive).
3. Supports outside `N(H)`:
   - checked support vertices not adjacent to heavy set: `5,631,569`
   - maximum `P(u) = 0.3307522458485829 < 1/3`.
4. Bridge lower-bound diagnostics:
   - `g := n/3-mu(T) >= R := F(H)-|A cap N(H)|/6` had `0` violations.
   - for `H != empty`, minimum observed `R`:
     `0.1380801129345332` (strictly positive).

Open point: convert these bridge diagnostics into a full analytic inequality.

## Reproducible packaged scan

Script:

- `conjecture_a_decimated_bridge_diagnostics.py`

Run:

```bash
python3 conjecture_a_decimated_bridge_diagnostics.py \
  --min-n 3 --max-n 23 --verify-original-mu \
  --out results/whnc_decimated_bridge_diag_n23.json
```

Summary from artifact `results/whnc_decimated_bridge_diag_n23.json`:

- `seen=23,942,357`, `considered=931,596`, `with_H=674,393`,
- `subsets_checked=2,881,985`, `private_fail=0`,
- `edge_checks=2,603,599`, `edge_margin_fail=0`,
- `support_outside_checks=5,631,569`, `support_outside_fail=0`,
- `bridge_viol=0`,
- `minimum_edge_margin=0.021903464487004864`,
- `maximum_support_outside_p=0.3307522458485829`,
- `minimum_r_nonempty_h=0.1380801129345332`,
- `minimum_non_support_neighbors_at_heavy=0` (witness confirms the stronger
  “heavy has non-support private neighbor” shortcut is false),
- `decomposition_max_abs_err=1.6653345369377348e-15`,
- `mu_core_max_abs_err=1.7763568394002505e-15`.
