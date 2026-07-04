# Local-Mode Mean Bound
Date: 2026-07-04

## Purpose

This note attacks the issue #6 proof target:

> If `0 <= w_i <= 1`, `e_2(w) > 0`, and `e_3(w) >= e_2(w)`, then
>
> ```text
> mu = sum_i w_i/(1+w_i) >= 5/2.
> ```

Equivalently, for `S=sum_i Bernoulli(p_i)` with all `p_i <= 1/2`, if

```text
P(S=3) >= P(S=2) > 0,
```

then

```text
E S >= 5/2.
```

The equivalence uses the finite odds transform `w_i=p_i/(1-p_i)`;
the assumption `p_i <= 1/2` ensures `p_i<1` and hence

```text
P(S=k) = prod_i(1-p_i) e_k(w).
```

The result below proves this target. It does not prove the signed hub-bouquet reserve and does not close the signed conditional/boundary-term part of issue #5.

## Notation

Delete all zero weights; they do not affect `e_2`, `e_3`, or `mu`. Let `m` be the number of remaining positive weights. Write

```text
W = e_1 = sum_i w_i,
mu = sum_i p_i,          p_i = w_i/(1+w_i),
q_i = 1/(1+w_i).
```

For a set of weights with elementary symmetric functions `e_k`, Newton's inequalities imply the normalized sequence

```text
e_k / binom(r,k)
```

is log-concave, where `r` is the number of positive weights in that set.

The only use of Newton here is the following small consequence.

**Lemma 1.** For positive weights, if `0 < e_2 < e_1`, then `e_3 < e_2`.

**Proof.** If there are fewer than three positive weights, then `e_3=0`, and the claim is immediate from `e_2>0`. Otherwise, with `r >= 3` positive weights, Newton's inequalities give

```text
e_3/e_2 <= [2(r-2)/(3(r-1))] e_2/e_1.
```

The prefactor is less than `1`, and `e_2/e_1 < 1`, so `e_3/e_2 < 1`. Thus `e_3 < e_2`. QED.

## Deletion Lemma

Assume now that the full weight vector satisfies

```text
e_3 >= e_2 > 0.
```

For a fixed index `i`, let `a=w_i`, and let

```text
A = e_1^{(-i)},     B = e_2^{(-i)},     C = e_3^{(-i)}
```

be the elementary symmetric functions after deleting `i`. Then

```text
e_2 = B + aA,
e_3 = C + aB.
```

I claim that

```text
B >= A
```

for every positive coordinate `i`.

Indeed, suppose instead that `B < A`. If `B=0`, then `C=0`; since the full `e_2=B+aA` is positive and `a>0`, also `A>0`. Thus

```text
e_3 = 0 < aA = e_2,
```

contradicting `e_3 >= e_2`. If `B>0`, Lemma 1 applied to the deleted weight vector gives `C < B`. Since `a>0` and `B<A`,

```text
e_3 = C + aB < B + aB <= B + aA = e_2,
```

contradicting `e_3 >= e_2`. Therefore `B >= A`.

Consequently,

```text
e_2 = B + w_i A >= (1+w_i)A.
```

## Size-Biased Pair Argument

Choose a random unordered pair `Q={i,j}` from the positive support with probability

```text
Pr(Q={i,j}) = w_i w_j / e_2.
```

This is well-defined because `e_2>0`.

For a fixed `i`,

```text
Pr(i in Q) = w_i e_1^{(-i)} / e_2.
```

By the deletion lemma,

```text
e_2 >= (1+w_i)e_1^{(-i)}.
```

Hence

```text
Pr(i in Q) <= w_i/(1+w_i) = p_i,
```

and therefore

```text
Pr(i notin Q) >= q_i = 1/(1+w_i).
```

The condition `e_3 >= e_2` has a direct interpretation under this pair measure. Since each triple is counted once for each of its three pairs,

```text
E[ sum_{i notin Q} w_i ] = 3e_3/e_2 >= 3.
```

Now write `u_i=1-w_i`. For every pair `Q`,

```text
sum_{i notin Q} u_i = (m-2) - sum_{i notin Q} w_i.
```

Taking expectations and using the previous display,

```text
E[ sum_{i notin Q} u_i ] <= m-5.
```

In particular, the assumptions force `m >= 5`; otherwise the nonnegative
left-hand side could not be at most `m-5`.

On the other hand, the one-coordinate lower bound gives

```text
E[ sum_{i notin Q} u_i ]
  = sum_i u_i Pr(i notin Q)
  >= sum_i u_i q_i.
```

But

```text
u_i q_i = (1-w_i)/(1+w_i)
        = 1 - 2 w_i/(1+w_i)
        = 1 - 2p_i.
```

Therefore

```text
sum_i u_i q_i = m - 2mu.
```

Combining the upper and lower bounds,

```text
m - 2mu <= m - 5,
```

so

```text
mu >= 5/2.
```

This proves the local-mode mean bound.

The restriction `w_i <= 1` is used in the last comparison through
`u_i=1-w_i >= 0`; Lemma 1 and the deletion lemma themselves only need
positive weights.

## Boundary Case

The example

```text
w_1 = ... = w_5 = 1
```

with all other weights zero has

```text
e_2 = e_3 = 10,
mu = 5/2.
```

Thus the constant `5/2` is sharp. In probability language this is `Bin(5,1/2)`, where `P(S=2)=P(S=3)`.

The proof above only needs sharpness, not a uniqueness classification for equality cases.

## Consequence for the One-Sided Route

The one-sided localization reduction had isolated the following live target:

```text
0 <= w_i <= 1, e_2>0, e_3>=e_2
    => sum_i w_i/(1+w_i) >= 5/2.
```

That target is now proved here. Since `w_i<=1` implies

```text
w_i/(1+w_i)^2 >= (1/2) w_i/(1+w_i),
```

the variance bound

```text
V = sum_i w_i/(1+w_i)^2 >= 5/4
```

also follows.

This closes the remaining `D=4` gap in the one-sided localization route. It does not by itself prove the signed reserve needed for the hub-bouquet argument.
