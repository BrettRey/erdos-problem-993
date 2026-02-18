# Tie-Fugacity Condition: Leaf Decomposition Analysis (2026-02-18)

Goal: prove `mu(lambda_m) >= m-1` for all trees, where `m = mode(I_T)` and
`lambda_m = i_{m-1}/i_m`.

## Key identity: mean decomposition at any fugacity lambda

For leaf `l` with support `s` in tree `T`, and any `lambda > 0`:

```
mu(T, lambda) = P(l, lambda) * (1 + mu(T-{l,s}, lambda))
              + (1 - P(l, lambda)) * mu(T-l, lambda)
```

This is the law of total expectation: the IS distribution splits by whether `l`
is included or excluded.

Rearranging (with `p = P(l, lambda) > 0`):

```
mu(T, lambda) - (m-1)
  = p * (mu(T-{l,s}, lambda) - (m-2))
  + (1-p) * (mu(T-l, lambda) - (m-1))
```

This is positive whenever BOTH of the following hold at `lambda = lambda_m^T`:

1.  `mu(T-l,  lambda_m^T) >= m-1`   (same threshold as T)
2.  `mu(T-{l,s}, lambda_m^T) >= m-2` (threshold lowered by 1)

## Why the induction closes when delta(T, l) = 0 or +1

If `mode(T-l, 1) = m` (mode unchanged by leaf removal, `delta = 0`), then:
- Tie-fugacity condition for T-l at its own tie point `lambda_m^{T-l}` gives
  `mu(T-l, lambda_m^{T-l}) >= m-1`.
- If additionally `lambda_m^T <= lambda_m^{T-l}`, then by monotonicity of `mu`:
  `mu(T-l, lambda_m^T) >= mu(T-l, lambda_m^{T-l}) >= m-1`. Condition (1) holds.

If `mode(T-l, 1) = m+1` (`delta = +1`, rare: 0.09% of cases), condition (1) is
even easier to satisfy since T-l has a higher mode at the same fugacity.

## The obstruction: delta = -1 case (mode drops when removing leaf l)

If `mode(T-l, 1) = m-1` (`delta = -1`, occurs in 0.53% of cases), then:
- Tie-fugacity for T-l gives `mu(T-l, lambda_{m-1}^{T-l}) >= m-2`.
- At `lambda_m^T` (different fugacity): we only know `mu(T-l, lambda_m^T) >= m-2`
  if `lambda_m^T >= lambda_{m-1}^{T-l}` (monotonicity).
- But we need `mu(T-l, lambda_m^T) >= m-1`, one unit more.
- The deficit `-(1-p) * 1` from condition (1) must be covered by the first term.

For condition (2) when `mode(T-{l,s}, 1) >= m-1`:
- Tie-fugacity gives `mu(T-{l,s}, lambda_{m-1}^{T-{l,s}}) >= m-2`.
- Again by monotonicity (if `lambda_m^T >= lambda_{m-1}^{T-{l,s}}`):
  `mu(T-{l,s}, lambda_m^T) >= m-2`. Condition (2) holds.

In the delta=-1 case, the decomposition becomes:
```
mu(T, lambda_m^T) - (m-1)
  = p * (mu(T-{l,s}, lambda_m^T) - (m-2))
  - (1-p) * 1    (worst case: mu(T-l, lambda_m^T) = m-2, one below threshold)
  + (1-p) * (mu(T-l, lambda_m^T) - (m-2))
```

For this to be non-negative, we need:
```
p * (mu(T-{l,s}, lambda_m^T) - (m-2)) >= (1-p) * (m-1 - mu(T-l, lambda_m^T))
```

I.e., the contribution of the leaf (its inclusion probability times the surplus
of T-{l,s} above m-2) must cover the deficit from T-l being below m-1.

## What additional structure could close the gap

Option A: Show that when delta(T, l) = -1, we have a lower bound
`P(l, lambda_m^T) >= c` for some `c > 1/2`, which would force `mu(T-{l,s})` to
contribute enough.

Option B: Show that `lambda_m^T >= lambda_{m-1}^{T-l}` always holds, giving
`mu(T-l, lambda_m^T) >= mu(T-l, lambda_{m-1}^{T-l}) >= m-2` -- but this is
already built into the estimate above and only gives m-2, not m-1.

Option C: Show that when delta = -1, the "cause" of the mode drop is a specific
structural feature of l (e.g., l is the unique leaf contributing to the top IS
coefficient), and exploit this structure to bound `mu(T-{l,s})` from below by m-1
(one unit above the threshold). Then condition (1) can be violated but condition
(2) provides extra cover.

Option D: Find a different inductive structure (e.g., induct on max IS size, or
on a structural parameter like number of leaves, instead of vertex count).

## Relation to Phi_m recursion

Note: `phi_m(T, lambda) = I(T,lambda) * (mu(T,lambda) - (m-1))` is the same as
the Phi_m functional studied earlier. The mean decomposition identity above is
equivalent to the Phi_m recursion:

```
Phi_m(T, lambda) = Phi_{m-1}(T-l, lambda) + lambda * Phi_{m-2}(T-{l,s}, lambda)
                 - I(T-l, lambda)
```

The -I(T-l, lambda) obstacle in the Phi recursion corresponds exactly to the
"deficit from condition (1) not holding" in the mean decomposition.

## Computational check: is lambda_m^T vs lambda_{m-1}^{T-l} monotone?

A useful experiment: for each n<=18 tree T and each leaf l:
- Compute lambda_m^T (tie fugacity for T)
- Compute lambda_{m-1}^{T-l,*}: the fugacity at which mu(T-l, lambda) = m-1
  (or the tie-fugacity lambda_{mode(T-l)-1}^{T-l} if that's the mode for T-l)
- Check whether lambda_m^T >= lambda_{m-1}^{T-l,*}

If this holds universally, it would establish condition (1) for delta=0 case
(and give m-2 for delta=-1 case, leaving a one-unit gap).

If it fails, the proof must go through condition (2) covering the deficit.
