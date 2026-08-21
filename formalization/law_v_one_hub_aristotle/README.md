This project was edited by [Aristotle](https://aristotle.harmonic.fun).

To cite Aristotle:
- Tag @Aristotle-Harmonic on GitHub PRs/issues
- Add as co-author to commits:
```
Co-authored-by: Aristotle (Harmonic) <aristotle-harmonic@harmonic.fun>
```

# LAW-V when one hub is trivial — Lean 4 formalization

Formalization of the request in `law_v_one_hub_aristotle_input_20260820.md`.

Aristotle project `059ba8f5-39e6-43f0-9944-8d1268d90740`, completed task
`9a6c6951-03eb-497b-bdb9-ed9f9013a1b7`. Independent local replay on 2026-08-21
completed all 8,032 jobs; the final axiom audit reported only `propext`,
`Classical.choice`, and `Quot.sound`.

Everything is developed from scratch on top of Mathlib (`Polynomial ℤ`), in the namespace
`LawV`.  There is **no** `sorry`, `admit`, `axiom`, `@[implemented_by]`, `native_decide`,
and no bounded/finite substitute for any general statement.

```
lake build                      -- builds every module (default target `RequestProject`)
lake build RequestProject.Main  -- final assembly + axiom audit (prints `#print axioms`)
```

## Status summary

| Layer | Content | Status |
|---|---|---|
| 1 | Path coefficient formula and adjacent likelihood-ratio inequality | **proved** |
| 2 | Common-convolution preservation of the likelihood-ratio order (Cauchy–Binet / TP2) | **proved** |
| 3 | Products of path polynomials are log-concave with no internal zeros; `R ≤lr Q` | **proved** |
| 4 | Log-concavity of the spider polynomial `C = Q + x·R` for arbitrary arms | **open** (see below) |
| 5 | Shifted cross term `c k * q (k-1) − c (k+1) * q (k-2) ≥ 0`, and `V = Turán(C) + cross` | **proved** |

Consequently:

* the LAW-V inequality `0 ≤ V k` is proved **unconditionally** for every family of at most
  two arms — in particular for the empty family and for the one-arm family
  (`V_nonneg_of_length_le_two`, `V_nonneg_nil`, `V_nonneg_single`); a spider with at most
  two arms *is* a path, and its polynomial is a path polynomial;
* it is also proved **unconditionally** for every star `K_{1,m}`, i.e. for arbitrarily many
  arms of length `1` (`Cof_logConcave_star`, `V_nonneg_star`);
* for an arbitrary family of arms the LAW-V inequality is proved **conditionally**, with
  layer 4 entering as an explicit hypothesis and never as an axiom
  (`V_nonneg_conditional`, and the sharper `V_nonneg_of_K`).

Grade: **not COMPLETE** — layer 4 is the only missing item.  This is the case the request
allows to be returned as "a clearly labelled conditional core and the exact theorem
signature still required".

### The exact theorem signature still required

```lean
theorem Cof_logConcave (arms : List ℕ) : LogConcave (Cof arms)
```

Adding this one theorem (in `RequestProject/Main.lean`, where `Cof` and `LogConcave` are in
scope) turns `V_nonneg_conditional` into the unconditional final theorem, with no other
change to the development.

The gap has been sharpened: it suffices to prove the single likelihood-ratio inequality

```lean
theorem Qof_LR_X_Cof (arms : List ℕ) : LR (Qof arms) (X * Cof arms)
```

(`Cof_logConcave_of_K` derives layer 4 from it, `V_nonneg_of_K` derives LAW-V from it).
Unfolding `C = Q + x·R`, this hypothesis is exactly the coefficient inequality

```
q (k+1) * r (k-2)  ≤  q k * r (k-1) + (q k ^ 2 − q (k+1) * q (k-1))     for all k,
```

a statement about the two products of path polynomials `Q = ∏ P(aᵢ)` and `R = ∏ P(aᵢ-1)`
alone: the Turán gap of `Q` has to absorb the failure of the "synchronisation" inequality
`q (k+1) * r (k-2) ≤ q k * r (k-1)`.

## Traceability table

Every declaration below is user-facing.  "Source" records where the statement comes from:
`request` = explicitly demanded by the request document, `agent` = auxiliary API introduced
here.

| Declaration | File | Layer | Source | Meaning |
|---|---|---|---|---|
| `cf` | `Basic.lean` | — | agent | coefficient at an **integer** index (`0` outside the natural support) |
| `LR` | `Basic.lean` | — | agent | likelihood-ratio order `p ≤lr q`: `cf p (k+1) * cf q k ≤ cf p k * cf q (k+1)` |
| `LogConcave` | `Basic.lean` | — | agent | `cf p (k+1) * cf p (k-1) ≤ cf p k * cf p k` |
| `NiceD`, `Nice` | `Basic.lean` | — | agent | positive on `0,…,d` and zero afterwards ("log-concave with no internal zeros" support data) |
| `NiceD.mul`, `NiceD.add_X_mul` | `Basic.lean` | — | agent | closure of the support shape under `*` and `p + x*q` |
| `W` | `PathPoly.lean` | 1 | request | path independence polynomials, `W n = P (n-1)`; an arm `a` contributes `W (a+1)` to `Q` and `W a` to `R` |
| `W_coeff` | `PathPoly.lean` | 1 | request | `[x^j] W n = choose (n-j) j`, i.e. `[x^j] P n = choose (n+1-j) j` |
| `path_choose_LR` | `PathPoly.lean` | 1 | request | `p(n,j) * p(n-1,j-1) ≥ p(n,j-1) * p(n-1,j)` in binomial form, all boundary cases included |
| `W_LR_succ`, `W_LR_X`, `W_LR_X_X` | `PathPoly.lean` | 1 | agent | `W n ≤lr W (n+1)`, `W (n+1) ≤lr x*W n`, `W (n+1) ≤lr x²*W n` |
| `W_niceD`, `W_logConcave` | `PathPoly.lean` | 1 | agent | `W n` is positive exactly on `0,…,⌊n/2⌋`, and log-concave |
| `W_add` | `PathPoly.lean` | 1 | agent | splitting identity `W (a+b+2) = W (a+1) W (b+1) + x W a W b` |
| `cf_mul_eq_sum` | `LR.lean` | 2 | agent | convolution formula with a fixed finite summation range |
| `LR.gen`, `LogConcave.gen`, `LRgen` | `LR.lean` | 2 | agent | from adjacent minors to *all* `2×2` minors |
| `LR.mul_left` | `LR.lean` | 2 | request | **common-convolution preservation of `≤lr`**: `h` log-concave with no internal zeros and `p ≤lr q` give `h*p ≤lr h*q` (Cauchy–Binet) |
| `LogConcave.mul` | `LR.lean` | 2 | agent | products of log-concave polynomials without internal zeros are log-concave |
| `LR.trans'`, `LR.trans_X` | `LR.lean` | 2 | agent | transitivity, and the shifted form `p ≤lr x*u`, `u ≤lr v` ⟹ `p ≤lr x*v` |
| `Qof`, `Rof` | `Products.lean` | 3 | request | `Q = ∏ P(aᵢ) = ∏ W (aᵢ+1)`, `R = ∏ P(aᵢ-1) = ∏ W aᵢ` |
| `Qof_niceD`, `Rof_niceD` | `Products.lean` | 3 | request | log-concave **with no internal zeros**: positive on `0,…,deg` |
| `Qof_logConcave`, `Rof_logConcave` | `Products.lean` | 3 | request | both products are log-concave |
| `Rof_LR_Qof` | `Products.lean` | 3 | request | **`R ≤lr Q`**, obtained by replacing the factors one at a time |
| `Cof`, `Bof`, `Vco` | `Spider.lean` | — | request | `C = Q + x*R`, `B = C + x*Q`, `V k = c k * b k − c (k+1) * b (k-1)` |
| `Rof_LR_Cof`, `Rof_LR_X_Qof` | `Spider.lean` | 5 | agent | `R ≤lr C` and `R ≤lr x*Q` |
| `Cof_LR_X_X_Qof` | `Spider.lean` | 5 | agent | `C ≤lr x²*Q`, the polynomial form of the shifted cross term |
| `cross_term_nonneg` | `Spider.lean` | 5 | request | **`c k * q (k-1) − c (k+1) * q (k-2) ≥ 0`** |
| `Vco_nonneg_iff` | `Spider.lean` | 5 | agent | `(∀ k, 0 ≤ V k) ↔ C ≤lr x*B` |
| `V_nonneg_of_logConcave` | `Spider.lean` | 5 | request | **`V k ≥ 0` from the Turán gap of `C` plus the cross term** |
| `Cof_nil_eq`, `Cof_single_eq`, `Cof_pair_eq` | `Main.lean` | 4 | agent | a spider with `0`, `1`, `2` arms is a path: `C = W 2`, `W (a+2)`, `W (a+b+2)` |
| `Cof_logConcave_of_length_le_two` | `Main.lean` | 4 | agent | **layer 4, unconditional for at most two arms** |
| `Qof_star`, `Rof_star`, `Cof_star` | `Main.lean` | 4 | agent | the star `K_{1,m}`: `Q = (1+x)^m`, `R = 1`, `C = (1+x)^m + x` |
| `star_choose_key` | `Main.lean` | 4 | agent | `choose m 3 * (m+1) ≤ (choose m 2)^2`, the one binomial inequality behind the star case |
| `Cof_logConcave_star` | `Main.lean` | 4 | agent | **layer 4, unconditional for every star `K_{1,m}`** |
| `V_nonneg_star` | `Main.lean` | final | agent | **LAW-V, unconditional, for every star `K_{1,m}`** |
| `V_nonneg_of_length_le_two`, `V_nonneg_nil`, `V_nonneg_single` | `Main.lean` | final | request | **LAW-V, unconditional, including the empty and one-arm boundary cases** |
| `Cof_logConcave_of_K`, `V_nonneg_of_K` | `Main.lean` | 4 | agent | reduction of the remaining gap to `Q ≤lr x*C` |
| `V_nonneg_conditional` | `Main.lean` | final | request | **LAW-V for arbitrary arms, with layer 4 as an explicit hypothesis** |
| `two_hub_cross_term_counterexample` | `Main.lean` | — | request | the optional adversarial check |

## Adversarial check (optional item of the request)

`two_hub_cross_term_counterexample` proves, by exact coefficient computation, that for the
two spiders with arms `(1,1,5)` and `(1,1)`,

```
c⁽¹¹⁵⁾ 3 * c⁽¹¹⁾ 2 − c⁽¹¹⁵⁾ 4 * c⁽¹¹⁾ 1 = −3 < 0,
```

so the naive general two-hub cross term fails.  (`C (1,1) = 1 + 3x + x²` and
`C (1,1,5) = 1 + 8x + 21x² + 21x³ + 8x⁴ + x⁵`.)  This does not affect the one-trivial-hub
theorem, where the second factor is the arm product `Q`; that case is exactly layer 5,
`cross_term_nonneg`, which is proved.

## Axiom audit

`RequestProject/Main.lean` ends with a `#print axioms` audit.  `lake build
RequestProject.Main` reports, for each of

`path_choose_LR`, `LR.mul_left`, `Rof_LR_Qof`, `cross_term_nonneg`,
`V_nonneg_of_logConcave`, `V_nonneg_of_length_le_two`, `V_nonneg_star`,
`V_nonneg_conditional`, `V_nonneg_of_K`, `two_hub_cross_term_counterexample`:

```
depends on axioms: [propext, Classical.choice, Quot.sound]
```

i.e. only the three standard Lean axioms.

## Notes on layer 4

Layer 4 is the Li–Li–Yang–Zhang spider theorem.  The pinned Mathlib contains no
independence-polynomial API (a text search of the pinned Mathlib sources finds no
occurrence of "independence polynomial"), and layer 4 is not derivable from the layer 1–3
machinery by the usual likelihood-ratio manipulations.  For the record, the following
reductions are proved in the development, and the listed obstructions are the reason the
induction does not close:

* `LogConcave C` follows from `K : Q ≤lr x*C` together with the proved `R ≤lr C`
  (`Cof_logConcave_of_K`); and `0 ≤ V k` follows from `K` alone, since the other three
  cross conditions (`Q ≤lr x²*Q`, `R ≤lr C`, `R ≤lr x*Q`) are proved.
* `K` unfolds to `q(k+1) r(k-2) ≤ q k · r(k-1) + Turán_Q(k)`.  The "synchronisation"
  inequality `q(k+1) r(k-2) ≤ q k · r(k-1)` (equivalently `Q ≤lr x²·R`) is **false** in
  general — it already fails for arms `(1,1,1)` — so the Turán gap of `Q` must be used, and
  every splitting of the sum `C = Q + x·R` into pieces comparable by `LR.mul_left`
  reproduces one of the false statements `Q ≤lr C`, `Q ≤lr x·R`, `Q ≤lr x²·R`.
* The star case `K_{1,m}` shows what has to replace those splittings: there `R = 1`, the
  splitting statement `Q ≤lr x²·R` fails from `m = 3` on, and log-concavity of `C` comes
  from the Turán gap of `Q = (1+x)^m` at the single index `k = 2`, i.e. from the binomial
  inequality `choose m 3 * (m+1) ≤ (choose m 2)^2`.

An unverified numerical exploration (not part of the Lean development) over all arm
families with at most 4 arms of length at most 7 supports the truth of layer 4 and of `K`,
and also of the auxiliary statements `C ≤lr x·Q`, `Q ≤lr B`, `B ≤lr x²·C`; it likewise
confirms the falsity of `Q ≤lr C`, `Q ≤lr x·R`, `Q ≤lr x²·R` and of the two-hub cross term
(the last of which *is* verified in Lean).

## File layout

```
RequestProject/Basic.lean      -- integer-indexed coefficients, ≤lr, log-concavity, support shape
RequestProject/LR.lean         -- layer 2: minors, transitivity, Cauchy–Binet convolution
RequestProject/PathPoly.lean   -- layer 1: path polynomials W n, coefficients, ≤lr facts
RequestProject/Products.lean   -- layer 3: Q, R, their support/log-concavity, R ≤lr Q
RequestProject/Spider.lean     -- C, B, V; layer 5 cross term; V ≥ 0 from log-concavity of C
RequestProject/Main.lean       -- final assembly, boundary cases, conditional core, axiom audit
```
