This project was edited by [Aristotle](https://aristotle.harmonic.fun).

To cite Aristotle:
- Tag @Aristotle-Harmonic on GitHub PRs/issues
- Add as co-author to commits:
```
Co-authored-by: Aristotle (Harmonic) <aristotle-harmonic@harmonic.fun>
```

# Log-concavity under convolution, and the π/3 root-sector criterion

Formalization of the request in
`lc_convolution_sector_aristotle_input_20260826.md`.

**Outcome: COMPLETE — G0, G1 and G2 are all proved**, in the stated generality,
with no `sorry`, `admit`, `axiom`, `@[implemented_by]` or `native_decide`.
Every named result depends only on `propext`, `Classical.choice`, `Quot.sound`.

## Representation choice

Everything is stated in the **sequence form** `f : ℕ → ℝ`, exactly as in the
packet.  Polynomial corollaries are derived from the sequence statements through
`Polynomial.coeff`.

The `G1` statements are written with `Fin n`-indexed vectors `![w, 1]` and
`![ρ², 2ρ cos φ, 1]`; these are read as sequences `ℕ → ℝ` through the padding
map `LC.ofVec` (the only device added purely to make those statements
typecheck).

## Files

| File | Contents |
| --- | --- |
| `RequestProject/Defs.lean` | the four requested definitions, plus `ofVec` |
| `RequestProject/ZBasic.lean` | zero-extension to `ℤ`, and the key "inner product dominates outer product" lemma |
| `RequestProject/CauchyBinet.lean` | the `2 × 2` Cauchy–Binet identity and its positivity form |
| `RequestProject/Conv.lean` | **G0** (sequence and polynomial forms) |
| `RequestProject/Factors.lean` | **G1** (both factor shapes, with the sharp iff) |
| `RequestProject/Sector.lean` | **G2** (the π/3 root-sector criterion) |
| `RequestProject/Sharpness.lean` | explicit witness that `NoInternalZeros` cannot be dropped from G0 |

## Traceability

`source = request` means the statement was given verbatim in the packet;
`source = added` means it was introduced here.

| Goal | Lean name | File | Source | Status |
| --- | --- | --- | --- | --- |
| `Nonneg` | `LC.Nonneg` | `Defs.lean` | request | definition |
| `LogConcave` | `LC.LogConcave` | `Defs.lean` | request | definition |
| `NoInternalZeros` | `LC.NoInternalZeros` | `Defs.lean` | request | definition |
| `Conv` | `LC.Conv` | `Defs.lean` | request | definition |
| padding of `![…]` to `ℕ → ℝ` | `LC.ofVec` | `Defs.lean` | added | definition |
| **G0** | `LC.logConcave_conv` | `Conv.lean` | request | proved |
| G0, polynomial corollary | `LC.logConcave_coeff_mul` | `Conv.lean` | request | proved |
| G0, log-concavity part | `LC.logConcave_of_conv` | `Conv.lean` | added | proved |
| G0, support part | `LC.noInternalZeros_of_conv` | `Conv.lean` | added | proved |
| `2 × 2` Cauchy–Binet | `LC.cauchyBinet_two` | `CauchyBinet.lean` | added | proved |
| TP2 closure under products | `LC.cauchyBinet_two_nonneg` | `CauchyBinet.lean` | added | proved |
| inner ≥ outer for log-concave sequences | `LC.inner_mul_ge` | `ZBasic.lean` | added | proved |
| **G1**, linear factor | `LC.logConcave_linear` | `Factors.lean` | request | proved |
| **G1**, quadratic factor (iff) | `LC.logConcave_quadratic` | `Factors.lean` | request | proved |
| G1, `\|φ\| ≤ π/3` consequence | `LC.logConcave_quadratic_of_abs_le_pi_div_three` | `Factors.lean` | request | proved |
| **G2** | `LC.logConcave_of_roots_in_sector` | `Sector.lean` | request | proved |
| `NoInternalZeros` is necessary | `LC.conv_not_logConcave_without_noInternalZeros` | `Sharpness.lean` | added | proved |

## Notes on the statements

* `LC.logConcave_conv` is stated exactly as requested, including the finiteness
  hypotheses `hfin`, `hgin`.  They are kept because they were asked for, but the
  proof does not use them: each value of `Conv f g` is already a finite sum, and
  the argument is purely local.
* `LC.logConcave_quadratic` is proved as an **iff**, so the sharpness of the
  `π/3` threshold (`cos²φ ≥ 1/4 ⟺ |cos φ| ≥ 1/2`, and `cos(π/3) = 1/2`) is on
  the record, not only the direction the application needs.
* `NoInternalZeros` is not weakened anywhere.  That it cannot be dropped is
  witnessed formally by `LC.conv_not_logConcave_without_noInternalZeros`:
  with `f = (1,0,0,1)` (coefficients of `1 + x³`, nonnegative and log-concave,
  but with an internal zero) and `g = (1,1)` (nonnegative, log-concave,
  internal-zero-free), the convolution `(1,1,0,1,1)` fails log-concavity at
  `k = 1`, where `c₂² = 0 < 1 = c₁c₃`.

## Proof routes

* **G0** follows route 1 of the packet: a nonnegative sequence with no internal
  zeros and log-concave is TP₂, in the concrete form
  `LC.inner_mul_ge` (`a ≤ b ≤ c ≤ d`, `a + d = b + c` ⟹ `f b · f c ≥ f a · f d`),
  proved by integer induction from the one-step shift inequality.  The
  convolution is the relevant entry of the product of the two Toeplitz matrices,
  and the `2 × 2` minor of the product is expanded by `LC.cauchyBinet_two` into
  a sum of products of `2 × 2` minors of the factors, each nonnegative.  All the
  index bookkeeping is done on `ℤ`-indexed zero-extensions, so no truncated
  natural subtraction appears.
* **G2** is by strong induction on the degree: extract a complex root by the
  fundamental theorem of algebra, split into a real root (giving a factor
  `X + w`, `w > 0`) and a non-real root (giving a real quadratic factor
  `X² + 2ρ cos φ · X + ρ²`), apply G1 to the factor and the induction hypothesis
  to the cofactor, and combine with G0.  Along the way one also gets that all
  coefficients up to the degree stay strictly positive
  (`LC.PosCoeffs`), which is what keeps the induction going.

## Building

```
lake build
```
