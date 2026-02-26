# Two-Type Root Mixture Algebra at Fixed `(lambda,rho)`

## Template
At canonical root `u`, suppose a two-type factorization:

- `P(x) = F1(x)^a * F2(x)^b`
- `Q(x) = x * G1(x)^a * G2(x)^b`

with `a,b in Z_{\ge 0}`.

At a fixed `lambda > 0`, define

- `r1 := G1(lambda)/F1(lambda) > 0`
- `r2 := G2(lambda)/F2(lambda) > 0`

Then

`rho = Q(lambda)/P(lambda) = lambda * r1^a * r2^b`.  (E1)

For two exponent pairs `(a,b)` and `(a',b')` with same `(lambda,rho)`, we get:

`r1^(a-a') * r2^(b-b') = 1`.  (E2)

Equivalently, with `Delta a := a-a'`, `Delta b := b-b'`,

`(Delta a, Delta b)` lies in

`H(r1,r2) := {(x,y) in Z^2 : r1^x r2^y = 1}`.  (E3)

## Exact uniqueness criterion
Equal `(lambda,rho)` forces `(a,b)=(a',b')` iff

`H(r1,r2) = {(0,0)}`.

That is, iff no nonzero integer pair `(x,y)` satisfies `r1^x r2^y = 1`.

## Implication for `N`
If type sizes are `n1, n2`, then

`N = a*n1 + b*n2`.

So along a kernel move `(Delta a, Delta b) in H(r1,r2)`:

`Delta N = n1*Delta a + n2*Delta b`.  (Ndiff)

Therefore even when `(a,b)` is not unique, `(lambda,rho)` still determines `N`
iff `Ndiff=0` for every kernel direction.

## Rank-1 kernel form
In the dependent case, `H(r1,r2)=Z*(p,q)` with primitive `(p,q)` satisfying

`r1^p r2^q = 1`.

Then all collisions satisfy

`(a',b') = (a,b) + k*(p,q)` (subject to nonnegativity),
and
`Delta N = k*(n1*p + n2*q)`.

Hence:

- if `n1*p + n2*q != 0`, equal `(lambda,rho)` can move `N`;
- if `n1*p + n2*q = 0`, equal `(lambda,rho)` keeps `N` fixed along the kernel.

## Use in current proof lane
This gives an exact obstruction map for heterogeneous-root attempts:
to build a split at fixed `(lambda,rho)`, one needs a nontrivial kernel direction
with nonzero `Delta N`; to prove injectivity in a restricted family, prove either
trivial kernel or kernel orthogonality to `(n1,n2)`.
