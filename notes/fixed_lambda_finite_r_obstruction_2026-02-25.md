# Fixed-lambda Finite-r Obstruction (2026-02-25)

## Scope

This note proves the following statement in the **unconstrained fixed-lambda product model**:

- For any fixed finite jet order `r >= 1`, the tuple `(d, lambda, mu_1, ..., mu_r)` does not determine `N = [x]P`.

This is **not** yet a canonical-tree theorem. The canonical lift requires additional bridges listed at the end.

## Model

Fix `r >= 1` and rational `lambda > 0`. Let `a := 1 + lambda`.

Allowed child factors are star polynomials:

```text
P_t(x) := I(K_{1,t}; x) = (1 + x)^t + x,   t >= 1.
```

Consider products

```text
P(x) = prod_j P_{t_j}(x),
```

with integer `t_j >= 1`.

Define invariants at `lambda`:

```text
d := deg(P),
mu_k := lambda^k P^{(k)}(lambda)/P(lambda),   1 <= k <= r,
N := [x]P.
```

## Theorem

For every fixed `r`, there exist two such products `P^+` and `P^-` with identical

```text
(d, lambda, mu_1, ..., mu_r)
```

but different `N`.

Equivalently: no fixed finite `r` forces `N` in this model.

## Proof

### 1) Additive coordinates

Define factorial cumulants:

```text
kappa_k(P) := lambda^k (log P)^{(k)}(lambda),   k >= 1.
```

Under products, cumulants add:

```text
kappa_k(prod_j P_j) = sum_j kappa_k(P_j).
```

`(mu_1, ..., mu_r)` and `(kappa_1, ..., kappa_r)` are equivalent via universal triangular transforms (Bell-polynomial relations). So matching cumulants up to `r` is equivalent to matching moments up to `r`.

### 2) Expansion for star factors

For `P_t(x) = (1+x)^t + x`, write

```text
P_t(x) = (1+x)^t (1 + g_t(x)),   g_t(x) = x(1+x)^(-t).
```

Then

```text
log P_t(x) = t log(1+x) + log(1+g_t(x)).
```

For each fixed `k >= 1`, at `x = lambda`:

```text
kappa_k(P_t) = alpha_k t + a^(-t) R_k(t) + O(a^(-2t) t^(2k)),
alpha_k = (-1)^(k-1) (k-1)! (lambda/a)^k,
```

where `R_k(t)` is a polynomial of degree `k` with nonzero leading coefficient.

A concrete closed form for the first correction term is:

```text
lambda^k g_t^{(k)}(lambda) = a^(-t) R_k(t),
R_k(t) = (-1)^(k-1) (lambda^k/a^k) t^{overline{k-1}} (k + lambda - lambda t).
```

So the degree-`k` leading term is explicit and nonzero.

### 3) Linear independence of profile functions

Define functions on integers `t >= 1`:

```text
f_0(t)=1,
f_1(t)=t,
f_{k+1}(t)=kappa_k(P_t),   1 <= k <= r.
```

Claim: `{f_0, ..., f_{r+1}}` are linearly independent.

Reason: if

```text
u_0 + u_1 t + sum_{k=1}^r v_k kappa_k(P_t) = 0
```

for all `t`, use the expansion above and let `t -> infinity`.

- Constant and linear terms force `u_0=0` and `u_1 + sum v_k alpha_k = 0`.
- Multiplying residual by `a^t` yields a polynomial identity

```text
sum v_k R_k(t) = 0
```

for all large `t`, hence identically.
- Since `deg R_k = k` with distinct degrees and nonzero top coefficients, all `v_k=0`, then `u_1=0`.

So only the trivial combination vanishes.

### 4) Choose `r+2` sizes with full augmented rank

Pick `t_1 < ... < t_{r+2}` such that the evaluation matrix

```text
M = [ f_i(t_j) ]_{0<=i<=r+1, 1<=j<=r+2}
```

has nonzero determinant.

Because `lambda` is rational, entries are rational; `det(M)` is nonzero rational.

Define augmented vectors

```text
wtilde_{t_j} := (1, t_j, kappa_1(P_{t_j}), ..., kappa_r(P_{t_j})) in Q^(r+2).
```

These are linearly independent.

Now drop the first coordinate:

```text
w_{t_j} := (t_j, kappa_1(P_{t_j}), ..., kappa_r(P_{t_j})) in Q^(r+1).
```

Since there are `r+2` vectors in dimension `r+1`, there is a nontrivial rational relation

```text
sum_j c_j w_{t_j} = 0.
```

Scale to integer `c_j`.

Key point: `sum_j c_j != 0`.

If `sum_j c_j = 0`, then

```text
sum_j c_j wtilde_{t_j} = 0,
```

contradicting independence of the augmented vectors.

So this integer relation simultaneously satisfies:

```text
sum_j c_j t_j = 0,
sum_j c_j kappa_k(P_{t_j}) = 0   (1 <= k <= r),
sum_j c_j != 0.
```

### 5) Build the two products

Split `c = c^+ - c^-` with nonnegative integer parts. Define:

```text
P^+(x) = prod_j P_{t_j}(x)^(c_j^+),
P^-(x) = prod_j P_{t_j}(x)^(c_j^-).
```

Then:

- Degrees match: `deg(P^+) - deg(P^-) = sum c_j t_j = 0`.
- Cumulants match up to `r`: `kappa_k(P^+) - kappa_k(P^-) = sum c_j kappa_k(P_{t_j}) = 0`.
- Therefore `mu_1, ..., mu_r` also match.

Now compare `N = [x]P`.

For stars:

```text
[x]P_t = t + 1.
```

For products with constant term `1`:

```text
[x] prod_j F_j = sum_j [x]F_j.
```

Hence

```text
N(P^+) - N(P^-)
= sum_j (c_j^+ - c_j^-) (t_j + 1)
= sum_j c_j t_j + sum_j c_j
= 0 + sum_j c_j
!= 0.
```

So `N` differs while `(d, lambda, mu_1, ..., mu_r)` agrees.

QED.

## What this does not yet prove

This theorem is fixed-lambda and unconstrained. It does not, by itself, prove the same no-go on the canonical class.

To lift, two additional bridges are needed:

1. Canonical admissibility bridge:
   replace the analytic star family with an infinite family of factors that is actually realized by canonical child components under gate checks.

2. Lambda-coupling bridge:
   preserve canonical derived `(m, lambda)` from

```text
I(T;x) = (1+2x)P(x) + (1+x)Q(x)
```

while performing the `P`-side kernel swap, with canonical gates still satisfied.

Without these bridges, the result remains a model-level obstruction, not a canonical theorem.
