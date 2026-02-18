# Sub-Claim B (Compensation Form) Progress (2026-02-18)

Target Sub-claim B in tie-fugacity leaf decomposition:

`p * slack_c2 >= (1-p) * deficit`

when condition (1) fails (`mu(T-l,lambda_m^T) < m-1`).

Here:
- `p = P(l,lambda_m^T)`,
- `slack_c2 = mu(T-{l,s},lambda_m^T) - (m-2)`,
- `deficit = (m-1) - mu(T-l,lambda_m^T)`.

Equivalent form:

`beta*c2 >= alpha*deficit`,

with `alpha = I(T-l,lambda)/I(T,lambda)`, `beta=lambda I(T-{l,s},lambda)/I(T,lambda)`,
`alpha+beta=1`, `c1 = mu(T-l,lambda)-(m-1)`, `c2 = mu(T-{l,s},lambda)-(m-2)`.

If `c1<0`, then `deficit=-c1` and

`beta*c2 - alpha*deficit = beta*c2 + alpha*c1 = margin(T)`.

So in the failing-`c1` regime, Sub-claim B is exactly `margin(T) >= 0`.

---

## Fully closed lane: mixed-spider unit-leaf branch with `j=1`

Take `T = S(2^k,1)`, and choose the unique unit leaf `l` (adjacent to hub):
- `A = T-l = S(2^k)`,
- `B = T-{l,s}` has polynomial `(1+2x)^k`.

At `lambda=lambda_m(T)`, Sub-claim B compensation gap is

`gap_comp := beta*c2 - alpha*deficit`.

When `c1<0` this equals `margin(T)`.
But branch-2 theorem already proved algebraically:

`margin(S(2^k,1),lambda_m) > 0` for all `k>=1`

(see `prove_mixed_spider_j0_j1_branches.py`, branch `j=1` proof template with exact finite bases).

Therefore:

**Sub-claim B compensation holds for all `S(2^k,1)` (all `k>=1`) algebraically.**

---

## What remains open

For general mixed spiders `S(2^k,1^j)` with `j>=2`, the same identity shows Sub-claim B
is equivalent to margin positivity in the `c1<0` regime.
So a full algebraic proof of Sub-claim B on that family reduces to proving
`margin(k,j) >= 0` for all `(k,j)` in that regime (or proving the minimizer-reduction
claims that reduce worst cases to `j in {0,1}`).

Current status:
- `j=0` and `j=1` branch positivity: algebraically proved.
- residue comparison lanes:
  - Sub-claim B (mod-1): proved in `notes/subclaim_B_mod1_algebra_2026-02-18.md`,
  - Sub-claim C (mod-0/2): algebraic script lane in `prove_subclaim_c_algebra.py`.
- parity-tail Sub-claim A (`margin(k,j+2) > margin(k,j)` for large `k`) remains the
  main missing reduction step.
