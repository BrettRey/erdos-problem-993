# Independent proof/certificate audit of the Poisson-binomial manuscript
<!-- SUMMARY: Full independent audit of paper/poisson_binomial/main.tex: replay, fresh symbolic re-derivation, exact end-to-end tests, HJ source check, Lean build — all passed · status: audit passed, gates remaining noted · updated: 2026-07-17 -->

Auditor: Claude (Fable 5), 2026-07-17, at Brett's request ("do the audit").
Scope: the complete proof chain of Theorem 1.1, Proposition 1.2, Corollary 1.3,
and Proposition 3.1 in `paper/poisson_binomial/main.tex` (state at commit
260cfd7), plus the certificate supplement and the conditional Lean project.

## 1. Exact replay (per CERTIFICATE.md)

- Generator rerun: output byte-identical to both stored certificates
  (`cmp` clean against the 2026-07-10 summary and 2026-07-16 full files).
- Independent checker (`python3 -I -S`): `"status": "passed"`, 13 cells,
  275 coefficients, payload digest `64cacd6c…629a` as documented.
- `shasum -a 256 -c` on the supplement zip: OK.
- All six whole-file digests in CERTIFICATE.md reproduced exactly.

## 2. Fresh symbolic re-derivation (not reusing repo code)

`scripts/audit_pb_scalar_independent_20260717.py` re-implements the scalar
section from the manuscript's formulas alone (SymPy, exact rationals):

- L_r = b_r identity and the C_{r+1}/C_r ratio identity (eq 3.3) confirmed.
- eq (3.6): A(δ) − Q(H) = P(H)/(4H⁵(H+1)³) with the printed degree-10 P
  confirmed symbolically from the raw L_r/R_r definitions (K = 3 window).
- Cell [3,4]: 11 Bernstein coefficients, all positive, min 22272 (matches
  Table 1). Cells m = 4..15: integer polynomials of degree 2m+2, all
  Bernstein coefficients positive; min coefficients match Table 1 where
  compared (m = 4, 5, 6, 7, 15). Total coefficients: 275.
- Bernstein conversion formula (eq 3.9) validated against reconstruction.
- Analytic range: S₀/T₀ closed forms match direct sums; the J = 5 cell's five
  degree-4 Bernstein coefficients match the manuscript exactly (2360, 7500,
  25055/2, 254205/16, 31115/2); for general J = u+6 the five Bernstein
  coefficients equal (multiplier)·(printed u-polynomial)/2880 symbolically,
  and every printed u-polynomial has strictly positive coefficients.

## 3. Exact end-to-end tests (Fractions, no floats)

`scripts/audit_pb_endtoend_exact_20260717.py`:

- Proposition 3.1 tested directly from the K/L_r/R_r/A(δ) definitions
  (independent of the paper's cell decomposition) at 760 exact rationals,
  including ±1e-9 straddles of every K-transition δ = 1/(m+1) for
  m = 4..59 and triangular-cell boundaries: 0 failures. Smallest observed
  margin ≈ 0.358 near δ = 1/4 (consistent with the [3,4] cell being
  the tight one).
- Theorem 1.1 and Corollary 1.3 (ratio drop + geometric tail) verified on
  368 random Poisson-binomial laws with V ≥ 1 in exact rational
  arithmetic: 0 violations; worst random score 0.422.
- The variance-one binomial family reproduces δ_D = (n+1)/(3(n−1)) exactly
  (n = 5..89), decreasing toward 1/3, matching Proposition 1.2.

## 4. Hand-verified prose steps

- Well-definedness of D via f_n/f_{n−1} = (Σ(1−p_i)/p_i)⁻¹ and V < 1
  contradiction: checked.
- HJ cubics → normalized recurrence δ_{k±1}(1−δ_k) ≤ δ_k by the stated
  divisions: re-derived.
- δ_D > 0 propagation argument, iteration map x ↦ x/(1−x) giving
  δ_{D±r} ≤ δ/(1−rδ), endpoint exclusion, q_D ≥ a, and both telescoping
  mass bounds R_r, L_r: re-derived, including positivity ranges.
- Pairwise variance identity and the window relaxation V ≥ p²A(δ): checked.
- The inline smoothing proof of the maximal-mass bound V ≥ (p⁻²−1)/12:
  checked (self-contained; the BMM citation is corroborative).
- The final squeeze (V < 1/(4δ) ⇒ p² > δ/(3+δ) vs p² < 1/(4δA(δ))):
  checked.
- Proposition 1.2: q_k formula, p_n > 1/n, p_n < 2/(n+1) for n ≥ 5
  (n² − 4n − 1 > 0), δ₂ algebra: checked.
- ULC(n) and Johnson-bound counterexample computations in Section 1: checked.

## 5. External input and Lean status

- The only external mathematical input is Hillion–Johnson, Ann. Probab.
  44(1) (2016), Theorem A.2 (eq 78) and Corollary A.3 (eq 79). Verified
  against the local source text `notes/literature/arxiv_1303_3381.txt`:
  quoted term-for-term correctly, including the all-integer boundary
  convention. Their inductive proof of (78) was not re-verified here;
  it is peer-reviewed Annals of Probability material.
- Lean project `formalization/pb_effective_drop_aristotle`: `lake build`
  succeeds (8028 jobs, one unused-variable lint). No `sorry`, `admit`,
  `axiom`, or `implemented_by`. Scope matches the disclosure: recurrence
  propagation, endpoint exclusion, first-crossing bound, raw-drop
  corollary, all conditional on the HJ recurrence as hypothesis `hstep`;
  each Lean statement's mathematical content re-checked by hand.

## 6. Verdict and remaining gates

No mathematical, computational, or provenance defects found. This audit is
an AI audit; it satisfies the independent-replay requirement but Brett must
decide whether it discharges the "human proof/certificate audit" gate.
Still open per STATUS/venue record: authenticated MathSciNet + expert
novelty checks, durable certificate-archive DOI, prediction-ledger
forecast.
