This project was edited by [Aristotle](https://aristotle.harmonic.fun).

To cite Aristotle:
- Tag @Aristotle-Harmonic on GitHub PRs/issues
- Add as co-author to commits:
```
Co-authored-by: Aristotle (Harmonic) <aristotle-harmonic@harmonic.fun>
```

# The even-connector Laurent block

Formalization of the request in
[`c2_even_connector_aristotle_input_20260820.md`](c2_even_connector_aristotle_input_20260820.md).

- Aristotle project: `f9cd8c48-bd2e-4aa8-bc67-3a9316657b61`.
- Independent local replay succeeded on 2026-08-20 (8,029 jobs).
- Local `#print axioms` checks on the decomposition, main theorem, and coefficient
  corollary found only `propext`, `Classical.choice`, and `Quot.sound`.

**Result: the requested theorem is true and is proved in full generality**
(all natural numbers `r, s ≥ 1` and all rationals `c, d ≥ 1`), with no `sorry`,
`admit`, `axiom`, `implemented_by`, `native_decide`, bounded replacement or
weakening of the statement.

## Building

```
lake build
```

Axiom check for the main theorem (`#print axioms EvenConnector.Gpoly_centrallyUnimodal`):

```
'EvenConnector.Gpoly_centrallyUnimodal' depends on axioms:
  [propext, Classical.choice, Quot.sound]
```

i.e. only the three standard Lean axioms.

## Files

| file | contents |
| --- | --- |
| `RequestProject/Binomial.lean` | the combinatorial core: binomial-row monotonicity, the Catalan bound on binomial descents, and the two master inequalities |
| `RequestProject/Coeff.lean` | coefficient calculus for `ℚ[T;T⁻¹]`, binomial coefficients with integer lower index, `X`, `H` and the coefficients of `X ^ n` |
| `RequestProject/Main.lean` | the block `G(r,s,c,d)`, its coefficients, central unimodality (main theorem) |

## Traceability table

| item of the request | formal counterpart | file |
| --- | --- | --- |
| Laurent polynomials with rational coefficients | `LaurentPolynomial ℚ` (`ℚ[T;T⁻¹]`), coefficient map `lcoeff` | `Coeff.lean` |
| `X = z + z⁻¹` | `EvenConnector.XL` | `Coeff.lean` |
| `H n = z^n + z^(-n)` | `EvenConnector.HL` | `Coeff.lean` |
| `G(r,s,c,d) = c*d*X^(r+s+1) + d*X^s*H(r+1) + c*X^r*H(s+1) + H(r+s+1)` | `EvenConnector.Gpoly` | `Main.lean` |
| centrally unimodal coefficients | `EvenConnector.CentrallyUnimodal` (symmetry under `m ↦ -m`, nonnegativity, and `coeff f (m+2) ≤ coeff f m` for `m ≥ 0`) | `Main.lean` |
| **Main theorem**: `G(r,s,c,d)` is centrally unimodal for `r,s ≥ 1`, `c,d ≥ 1` | `EvenConnector.Gpoly_centrallyUnimodal` | `Main.lean` |
| Consequence: `coeff (G r s c d) m - coeff (G r s c d) (m+2) ≥ 0` for `m ≥ 0` | `EvenConnector.Gpoly_coeff_sub_nonneg` | `Main.lean` |
| "First verify the exact identity `c*d*A + d*B + c*C + D = (c-1)*(d-1)*A + (c-1)*(A+C) + (d-1)*(A+B) + (A+B+C+D)`" | `EvenConnector.Gpoly_decomposition` | `Main.lean` |
| the ordinary-polynomial coefficient translation | `EvenConnector.coefG` together with `EvenConnector.lcoeff_Gpoly` and `EvenConnector.lcoeff_Gpoly_odd` | `Main.lean` |
| binomial-difference bounds used in the audited proof | `EvenConnector.key_two`, `EvenConnector.key_three` | `Binomial.lean` |

## The proof

Write `N = r + s + 1`. Since `X = z⁻¹(1+z²)` and `H k = z^(-k)(1+z^(2k))`, the block
is `z^(-N) P(z²)` with

```
P(y) = c*d*(1+y)^N + d*(1+y)^s*(1+y^(r+1)) + c*(1+y)^r*(1+y^(s+1)) + 1 + y^N .
```

Formally this is the content of `lcoeff_Gpoly`: the coefficient of `G` at the
exponent `2i - N` is `coefG r s c d i`, and all coefficients at exponents of the
other parity vanish (`lcoeff_Gpoly_odd`). Central unimodality of `G` is therefore
equivalent to: `P` is symmetric of degree `N` with nonnegative coefficients, and
its coefficient sequence is nonincreasing from the centre outwards.

Substituting `i = N - 1 - k` and using the symmetry `C(n, n-j) = C(n, j)`, the
required coefficient inequality collapses (see `coefG_step`) to the following
three statements about the natural numbers `k` with `2k + 2 ≤ N`, plus two
monotonicity facts on the decreasing halves of the rows `r` and `s`:

* `C(N,k) ≤ C(N,k+1)`;
* `C(N,k) + C(r,k) ≤ C(N,k+1) + C(r,k+1)`   (`key_two`, and its mirror in `s`);
* `C(N,k) + C(r,k) + C(s,k) + [k = 0] ≤ C(N,k+1) + C(r,k+1) + C(s,k+1)`
  (`key_three`).

These are exactly the four summands of the suggested decomposition
`(c-1)(d-1)A + (c-1)(A+C) + (d-1)(A+B) + (A+B+C+D)`; the `[k = 0]` term is the
contribution of `D = H(N)`, and this boundary term is where the naive version of
the split fails. The audit of the source argument is closed as follows.

The `j = s` / `s < j < N/2` split of the source proof is replaced by a uniform
Catalan estimate, which also repairs the boundary and parity cases:

1. `choose_sub_succ_le_catalan`: for **all** `n` and `k`,
   `C(n,k) - C(n,k+1) ≤ catalan k`. Indeed `n ↦ C(n,k) - C(n,k+1)` increases up
   to `n = 2k` and decreases afterwards, and its value at `n = 2k` is
   `catalan k` (`central_diff_eq_catalan`).
2. `catalan_succ_le_choose_diff`: for `2k + 2 ≤ n`,
   `catalan (k+1) ≤ C(n,k+1) - C(n,k)`. Indeed the ascent `C(n,k+1) - C(n,k)` is
   nondecreasing in `n` once `2k + 1 ≤ n` (`ascent_mono`), and at `n = 2k+2` it
   equals `catalan (k+1)`.
3. `two_catalan_le_catalan_succ`: `2 * catalan k ≤ catalan (k+1)` for `k ≥ 1`,
   from `(k+2) * catalan (k+1) = 2 * (2k+1) * catalan k`.

For `k ≥ 1` these give
`C(N,k+1) - C(N,k) ≥ catalan (k+1) ≥ 2 * catalan k ≥ (C(r,k) - C(r,k+1)) + (C(s,k) - C(s,k+1))`,
which is `key_three`; the case `k = 0` (the boundary case carrying the extra `+1`
coming from `H(N)`) is checked directly and uses `r, s ≥ 1`.
