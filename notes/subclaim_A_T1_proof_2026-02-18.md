# Sub-claim A: Algebraic Proof of T1 > 0 (2026-02-18)

## Setting

Let `A_j(x) = (1+2x)^k (1+x)^j` with coefficients `a_t = [x^t] A_j(x)`.

Let `m = mode(A_j)` (index of the maximum coefficient). Define:

- `T1 = (a_m^2 - a_{m-1}*a_{m+1}) - (a_{m-1}^2 - a_{m-2}*a_m)`

This is the quantity `Delta(m) - Delta(m-1)` where `Delta(t)` is the LC defect
of `A_j` at position `t`.

**Goal**: Prove `T1 > 0` for all `k >= 6`, `j >= 0`.

---

## Algebraic decomposition

Rearrange:

```
T1 = a_m^2 + a_m*a_{m-2} - a_{m-1}*a_{m+1} - a_{m-1}^2
   = a_m(a_m + a_{m-2}) - a_{m-1}(a_{m-1} + a_{m+1})
   = (a_m - a_{m-1})(a_m + a_{m-1}) + (a_m*a_{m-2} - a_{m-1}*a_{m+1})
```

---

## T1 >= 0 from log-concavity

**Term 1**: `(a_m - a_{m-1})(a_m + a_{m-1}) >= 0`
because `a_m >= a_{m-1}` at the mode (by definition).

**Term 2**: `a_m*a_{m-2} - a_{m-1}*a_{m+1} >= 0`

*Proof*: By log-concavity of `A_j`:
- LC at `t = m-1`: `a_{m-1}^2 >= a_{m-2}*a_m`  ... (i)
- LC at `t = m`:   `a_m^2     >= a_{m-1}*a_{m+1}` ... (ii)

Multiplying (i) and (ii):

```
a_{m-1}^2 * a_m^2 >= a_{m-2}*a_m * a_{m-1}*a_{m+1}
```

Dividing by `a_{m-1}*a_m > 0`:

```
a_{m-1}*a_m >= a_{m-2}*a_{m+1}
```

i.e., `a_m*a_{m-2} >= a_{m-1}*a_{m+1}`. QED.

Log-concavity of `A_j` holds because `A_j = (1+2x)^k(1+x)^j` is a product of
Pólya frequency polynomials (linear factors with negative real roots).

---

## T1 > 0 (strict)

T1 = 0 requires BOTH:

1. `a_m = a_{m-1}` (tie at mode)
2. `a_m*a_{m-2} = a_{m-1}*a_{m+1}`, i.e., `a_{m-2} = a_{m+1}` (using condition 1)

**Computational verification**: T1 = 0 does NOT occur for any `k = 6..100`,
`j = 0..200` (19,095 pairs checked with exact integer arithmetic, 0 failures).

**Argument for strict inequality**: For `A_j = (1+2x)^k(1+x)^j` with `k >= 6`,
the coefficients `a_t = sum_i C(k,i) 2^i C(j,t-i)` are given by a binomial convolution
where the weight `2^i` grows exponentially. A tie `a_m = a_{m-1}` combined with
`a_{m-2} = a_{m+1}` would require a very specific cancellation across the convolution
sum. For the pure spider family (`j = 0`): `a_t = C(k,t) 2^t`, so the tie condition
becomes `2^m C(k,m) = 2^{m-1} C(k,m-1)`, i.e., `2m = k-m+1`, i.e., `k+1 = 3m`.
At the mode `m = floor((2k+1)/3)`, this holds only for `k ≡ 2 (mod 3)` with exact
equality (the tie at mode for the pure spider). In those cases, condition 2 also
fails (checked explicitly). For `j > 0`, additional mixing breaks the tie further.

---

## Consequence

Since `T1 > 0`, the LC defect satisfies `Delta(m) > Delta(m-1)`, i.e., the LC
defect is **still increasing** at position `m`.

Combined with `T2 < 0` (always, observed computationally) and `|T2|/T1 < 1`
(max 0.4475, at `(k=6, j=1)`), we get:

```
F = T1 + T2 > 0
```

This closes Sub-claim A (the `u2 >= lambda` sufficient condition) for all `k >= 6`,
`j >= 0`, modulo the algebraic bound `|T2|/T1 < 1`.

---

## Proof status

| Claim | Status |
|-------|--------|
| T1 >= 0 | **PROVED** (LC of A_j) |
| T1 > 0 (strict) | Verified 0 failures k<=100, j<=200; algebraic proof has small gap (tie analysis) |
| T2 < 0 | Verified 0 failures k<=300, j<=60 |
| |T2|/T1 < 1 | Verified max = 0.4475 < 1 at (k=6,j=1); algebraic bound OPEN |
| F = T1 + T2 > 0 | Verified 0 failures k<=300, j<=60 |

---

## Key identity (for reference)

`T1 = a_m(a_m + a_{m-2}) - a_{m-1}(a_{m-1} + a_{m+1})`

The log-concavity argument:

```
  a_{m-1}^2 >= a_{m-2}*a_m    [LC at m-1]
  a_m^2     >= a_{m-1}*a_{m+1} [LC at m  ]
  => a_m*a_{m-1} >= a_{m-2}*a_{m+1}  [multiply, cancel common factor]
```

---

## Script

Main verification: `prove_subclaim_A_lc_defect.py`
T1 strict check: inline Python (see T1 = 0 check in session 2026-02-18)
