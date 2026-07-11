# D21 matching-block drift: exact reduction and local obstruction

## Result

Fix a maximum-matching partition into `alpha` singleton/double blocks and,
for an independent `r`-set `A`, write

```text
e(A)=|V(T-N[A])|,
b(A)=2(alpha-r)-e(A),
h(A)=|E(T-N[A])|.
```

The proposed averaged drift

```text
E_(r+1)b >= E_r b - 3
```

is **exactly equivalent to prefix GSB**, rather than a strictly simpler
matching-block invariant.  Moreover, the smallest prefix star `K_(1,4)`
rules out a proof by a pointwise containment-edge map or involution with loss
at most three.  Thus this route must use genuinely aggregate compensation;
the deterministic one-edge drift cannot by itself absorb the size bias.

## Exact algebra

Put `mu=E_r e` and `eta=E_r h`.  Sampling a uniform rank-`r+1` set by first
size-biasing a rank-`r` set by `e(A)` and then choosing a uniform addable
vertex gives

```text
b(A+v)-b(A)=deg_(T-N[A])(v)-1,

E_(r+1)b
  = E_r[e b]/mu - 1 + 2 eta/mu.
```

At fixed rank, `b=2(alpha-r)-e`, hence

```text
Cov_r(b,e)=-Var_r(e),
E_r[e b]/mu=E_r b-Var_r(e)/mu,

E_(r+1)b-E_r b
  = -1 + (2 eta-Var_r(e))/mu.                 (1)
```

Therefore

```text
E_(r+1)b >= E_r b-3
iff Var_r(e) <= 2 mu+2 eta
iff mu_(r+1) <= mu_r+1
iff (r+2)i_r i_(r+2)
       <= (r+1)i_(r+1)^2+i_r i_(r+1).
```

The matching partition disappears completely from (1), because `b` is an
affine transform of `e`.  Any covariance lemma strong enough to prove the
requested drift is consequently the original GSB inequality in equivalent
form unless it adds a new decomposition of `Var(e)-2 eta`.

## Smallest pointwise obstruction

Take `T=K_(1,4)`.  Then `alpha=4`, and

```text
ceil((2 alpha-1)/3)-2=1,
```

so rank `r=1` is in the required prefix.  Use the maximum matching consisting
of the center and one leaf; the other three leaves are singleton blocks.

At rank one:

```text
A={center}:  e(A)=0, b(A)=6,
A={leaf}:    e(A)=3, b(A)=3  (four choices).
```

Every rank-two independent set is a pair of leaves and has

```text
e(B)=2, b(B)=2.
```

Thus no map from all uniform rank-one sets to rank-two sets can satisfy
`b(B)>=b(A)-3`: the center set would require `b(B)>=3`, while every possible
target has `b(B)=2`.  In particular, the center has no outgoing containment
edge at all.  This is the smallest star where both rank one lies in the
prefix and the pointwise loss exceeds three (`K_(1,m)` has loss `m`).

The aggregate inequality nevertheless has ample slack:

```text
mu=12/5,       Var(e)=36/25,       eta=0,
E_1 b=18/5,    E_e-size-biased b=3,
E_2 b=2,       E_2 b-E_1 b=-8/5 >= -3.
```

So the exact compensation is fractional/aggregate: the exceptional center
set has zero size-biased mass, while its uniform mass is only `1/5`.

## Reproduction

```bash
python3 scratch_d21_block_drift_obstruction_certificate_20260711.py
```

The certificate enumerates both independent-set levels, verifies the
displayed maximum-matching blocks, checks the local residual-degree identity
on all 12 containment edges, and checks every fraction above exactly.

## Route decision

Freeze the bare `b`-drift/involution route.  A viable reopening condition is
an additional tree-structural decomposition of

```text
Var_r(e)-2 E_r h
```

whose positive pieces admit aggregate cancellation in the prefix.  Merely
renaming the size-bias term as `Cov(b,e)` does not create such a reduction.
