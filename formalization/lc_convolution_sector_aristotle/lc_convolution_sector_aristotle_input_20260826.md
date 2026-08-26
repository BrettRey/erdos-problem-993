# Formalization request: log-concavity is closed under convolution, and the π/3 root-sector criterion

## Objective

Two results, one reusable and one that this project used today and should not
have to keep on trust.

The **core** is that log-concavity with no internal zeros is preserved by
convolution (equivalently, by multiplying polynomials). This is used repeatedly
across this campaign and currently exists here only as a paper argument plus
finite numerical checks, so a formal version is a durable asset independent of
anything else in the queue.

The **application** is a root-location criterion: a positive-coefficient
polynomial all of whose roots lie within π/3 of the negative real axis has a
log-concave coefficient sequence. This project used its contrapositive on
2026-08-26 (`notes/matroid_intersection_number_lane_2026-08-26.md` §3) to
interpret a certified root computation, so it is currently load-bearing and
unproved here.

**Mathlib has no notion of log-concavity for sequences** (checked in the
2026-08-26 checkout: no `LogConcave` anywhere). The definitions below must
therefore be introduced, and getting them right is the single most important
thing in this packet: a correct Lean proof of a subtly wrong definition is worth
nothing. State them exactly as given, then prove.

## Definitions

Work with `f : ℕ → ℝ` of finite support, or with `Polynomial ℝ` and its
`coeff`, whichever is cleaner; say which you chose and keep it uniform.

```text
Nonneg f          :=  ∀ k, 0 ≤ f k
LogConcave f      :=  ∀ k, f (k+1) * f (k+1) ≥ f k * f (k+2)
NoInternalZeros f :=  ∀ i j k, i ≤ j → j ≤ k → f i ≠ 0 → f k ≠ 0 → f j ≠ 0
Conv f g          :=  fun n => ∑ i ∈ Finset.range (n+1), f i * g (n - i)
```

`NoInternalZeros` says the support is "interval-like": it is the hypothesis
without which the convolution theorem is **false**, so it must not be dropped or
weakened to `∀ k, f k ≠ 0`.

## G0 (required) — convolution closure

```text
theorem logConcave_conv
    (f g : ℕ → ℝ)
    (hf0 : Nonneg f) (hg0 : Nonneg g)
    (hf : LogConcave f) (hg : LogConcave g)
    (hfz : NoInternalZeros f) (hgz : NoInternalZeros g)
    (hfin : f.support.Finite) (hgin : g.support.Finite) :
    LogConcave (Conv f g) ∧ NoInternalZeros (Conv f g)
```

Equivalently, for `p q : Polynomial ℝ` with nonnegative, log-concave,
internal-zero-free coefficients, the coefficients of `p * q` are log-concave.
Deliver whichever form you prove; if you prove the sequence form, add the
polynomial corollary, since that is how the rest of the campaign will use it.

**Two proof routes; pick whichever Lean supports better.**

1. *Total positivity / Cauchy–Binet.* A nonnegative sequence with no internal
   zeros is log-concave iff its infinite Toeplitz matrix `A(i,j) = f (i - j)` is
   totally positive of order 2 (all 2×2 minors nonnegative). The Toeplitz matrix
   of a convolution is the product of the Toeplitz matrices, and TP2 is closed
   under products by the Cauchy–Binet formula. This is the conceptually clean
   route and the TP2 lemma is itself reusable.
2. *Direct.* Prove `c_k² - c_{k-1} c_{k+1} ≥ 0` for `c = Conv f g` by expanding
   both sides and pairing terms: the difference is a sum over index pairs of
   nonnegative contributions built from the two hypotheses. The bookkeeping is
   heavier but requires no matrix machinery.

Note the theorem genuinely needs `NoInternalZeros`; if you find a
counterexample to the statement as written, that is a valid and valuable outcome
— return it explicitly rather than silently adding hypotheses.

## G1 (required) — the two factor shapes

With `ρ > 0` and `φ` real:

```text
theorem logConcave_linear (w : ℝ) (hw : 0 < w) :
    LogConcave ![w, 1] ∧ Nonneg ![w, 1]

theorem logConcave_quadratic (ρ φ : ℝ) (hρ : 0 < ρ) :
    LogConcave ![ρ^2, 2*ρ*Real.cos φ, 1] ↔ Real.cos φ ^ 2 ≥ 1/4
```

and hence, for `|φ| ≤ π/3`, the quadratic `x² + 2ρcos(φ)x + ρ²` has log-concave
coefficients. The iff is the point: **π/3 is exactly the threshold**, since
`(2ρcos φ)² ≥ ρ²·1 ⟺ cos²φ ≥ 1/4 ⟺ |cos φ| ≥ 1/2`. Prove the equivalence,
not just the useful direction, so the sharpness is on the record.

## G2 (optional) — the sector criterion

```text
theorem logConcave_of_roots_in_sector (p : Polynomial ℝ)
    (hpos : ∀ k ≤ p.natDegree, 0 < p.coeff k)
    (hroots : ∀ z : ℂ, (p.map (algebraMap ℝ ℂ)).IsRoot z →
        |Complex.arg (-z)| ≤ Real.pi / 3) :
    LogConcave (fun k => p.coeff k)
```

Route: `p` has real coefficients, so its complex roots come in conjugate pairs;
factor `p` over `ℝ` into linear factors `(x + w)` with `w > 0` and quadratic
factors `x² + 2ρcos(φ)x + ρ²` with `|φ| ≤ π/3`; each has log-concave
coefficients by G1; conclude by G0 and induction on the number of factors.

**This is the part most likely to be expensive**, because it needs real-closed
factorization plumbing (`Mathlib/Analysis/Complex/Polynomial/Basic.lean` has the
fundamental theorem of algebra and conjugate-root material). If it turns out to
dominate the effort, **deliver G0 and G1 and stop** — they are the reusable
content and G2 is a corollary anyone can then assemble. Say plainly which
happened.

## Deliverable and grading

- Standalone current-Mathlib Lean project, source-labelled declarations, README
  traceability table.
- No `sorry`, `admit`, `axiom`, `implemented_by`, `native_decide`, hidden
  hypotheses, or a statement specialised to a fixed finite degree standing in for
  the general one.
- `lake build` clean, `#print axioms` for every named result.
- Grade **COMPLETE** only if G0 and G1 are both proved in the stated generality.
  G2 reached is a bonus, and say so explicitly. If G0 is false as stated, return
  the smallest explicit counterexample (two sequences, with their convolution and
  the failing index) and stop.

An honest PARTIAL with G0 proved beats a complete-looking submission that
weakens `NoInternalZeros`, assumes strict positivity everywhere to dodge the
support condition, or proves only the direction of G1 that this project happens
to want.
