# Finite-Window Blocker + 3-Type Family Check (2026-02-27)

## 1) Finite-window separation reduction (fixed `(m,lambda)`)

Let
- `R_{m,lambda}(N)` be feasible `rho=Q(lambda)/P(lambda)` values at size `N=[x]P`.

A nontrivial window claim
- `exists K>=1` such that `R_{m,lambda}(N) ∩ R_{m,lambda}(N') = empty` for `0<|N-N'|<=K`

is equivalent to adjacent disjointness as the minimal requirement:
- necessarily `R_{m,lambda}(N) ∩ R_{m,lambda}(N+1) = empty` for all `N`.
- conversely, adjacent disjointness implies the claim with `K=1`.

So the exact blocker is an adjacent-separation theorem (or explicit adjacent counterexample).

## 2) Why current envelope bounds do not settle adjacency

Available universal envelopes at fixed `(m,lambda,rho)` provide:
- lower bounds on `N` from incidence and root-message arity,
- finite upper bound on `N` from exponential-vs-binomial inequality.

But for adjacent pair `(N,N+1)=(2k-1,2k)`, both parity exponents used by those bounds agree:
- `floor((N+1)/2)=k`, `ceil(N/2)=k`, and same for `N+1=2k`.

Therefore those envelopes are structurally blind to the critical adjacent pair and cannot by themselves prove any nontrivial `K>=1` separation.

## 3) Realizable 3-type family no-go (exact bounded stress check)

Family (canonical-root child types):
- `E`: rooted edge endpoint,
- `P`: rooted path3 endpoint,
- `F`: rooted depth-2 fork root.

At fixed `lambda`:
- `rho = lambda * R_E^e * R_P^p * R_F^f`,
- `N = 2e + 3p + 5f`.

Added script:
- `scripts/three_type_family_no_go_check.py`

Runs:

```bash
python3 scripts/three_type_family_no_go_check.py \
  --max-e 10 --max-p 10 --max-f 10 --max-den 20 \
  --out results/three_type_family_no_go_check_den20_e10p10f10.json
```
- lambdas: 127
- triples/lambda: 1331
- split_collisions: 0

```bash
python3 scripts/three_type_family_no_go_check.py \
  --max-e 20 --max-p 20 --max-f 20 --max-den 30 \
  --out results/three_type_family_no_go_check_den30_e20p20f20.json
```
- lambdas: 277
- triples/lambda: 9261
- split_collisions: 0

Interpretation:
- No same-`(lambda,rho)` / different-`N` split found in this realizable 3-type family on a large exact grid.
- This is consistent with (but does not prove) the claimed 3-type no-go theorem.

## 4) Practical next target

Highest-impact theorem target remains unchanged:
- prove or refute adjacent separation
  `R_{m,lambda}(N) ∩ R_{m,lambda}(N+1) = empty`.

Without this, no nontrivial finite-window separation result can be concluded.
