# Next Prompt Set (A/B/C) — 2026-02-27

## Prompt A (Proof Lane: separation theorem)
You are working in the min-`u` canonical gated class (`d_leaf<=1`).

Known:
- `(m,lambda)->N` is proved for `m<=3`.
- For `m>=4`, we have rigorous lower bounds and a finite upper bound on `N` from `(m,lambda,rho)`:
  - `N >= ceil(m - 4 + m/lambda)`
  - `N >= 2*r_min - 1`, where `r_min` is least integer with `(1+lambda)^r >= lambda/rho`
  - `(1-lambda)*((1+2lambda)+(1+lambda)rho)*(1+lambda)^{ceil(N/2)} <= binom(N+3,m)`
- Exhaustive scans: no split for `(m,lambda,rho)` through `n<=25` with `m>=4` (min-`u`).

Task:
Prove or refute the following finite-window separation claim:
For fixed `(m,lambda)` with `m>=4`, there exists an explicit `K=K(m,lambda)` such that
`R_{m,lambda}(N) ∩ R_{m,lambda}(N') = empty` whenever `0 < |N-N'| <= K`.

Requirements:
1. Use only theorem-grade implications (no heuristics).
2. If proving fails, return the strongest exact obstruction in lemma form.
3. If refuting, give an explicit canonical-tree witness pair with same `(m,lambda,rho)` and different `N`.
4. Keep every hidden assumption explicit (well-definedness at `lambda`, canonical-triplet dependence, etc.).

Deliverable:
A clean theorem/lemma chain with full hypotheses and either:
- a proof of the finite-window claim with explicit `K`, or
- an explicit counterexample.

## Prompt B (Audit Lane: close logic gaps in current proof stack)
Audit the current `(m,lambda,rho)->N` route as a hostile referee.

Scope:
- Message recursion (`t_v = dp1[v]/dp0[v`) and root equation for `rho`.
- The `N >= 2r-1` step from `d_leaf<=1` (leaf-child transfer from `T` to `B`).
- Coefficient-ratio inequality `i_{k-1}/i_k >= k/(n-k+1)` and its use at `k=m`.
- Upper-bound proposition using `(1-lambda)I(lambda) <= i_m <= binom(N+3,m)`.

Task:
1. Identify every hidden assumption and every implication that needs a condition.
2. Rewrite each result with minimal corrected hypotheses.
3. State whether any corrected lemma becomes false; if yes, give counterexample.
4. Produce a final “safe” theorem list (only statements that survive audit).

Deliverable:
A patched theorem stack suitable for direct insertion into notes/paper, with exact assumptions and no silent leaps.

## Prompt C (Construction/Kernel Lane: heterogeneous-root collisions)
Work in the two-type and three-type canonical-root template:
`P = prod_i F_i^{a_i}`, `Q = x*prod_i G_i^{a_i}` at fixed `lambda`.

Known exact criterion:
Equal `(lambda,rho)` between exponent vectors is equivalent to lattice-kernel condition
`prod_i (G_i(lambda)/F_i(lambda))^{Delta a_i} = 1`.

Task:
1. For rational `lambda in (0,1)`, derive necessary and sufficient arithmetic conditions for nontrivial kernel directions.
2. Characterize when kernel directions necessarily preserve `N` (`Delta N=0`) vs allow `Delta N!=0`.
3. Test the smallest realizable canonical families beyond the proven two-type no-go:
   - include at least one 3-type rooted family realizable in gated min-`u` trees,
   - either produce same `(lambda,rho)` with different `N`, or prove no-go for that family.
4. If you claim realizability, give explicit labelled construction and prove min-`u` canonical triplet selection.

Deliverable:
Either a concrete lifted split witness in a realizable finite-type family, or a theorem-grade no-go criterion that strictly generalizes the current two-type no-go.
