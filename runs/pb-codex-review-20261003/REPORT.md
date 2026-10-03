# Codex review of the PB paper (2026-10-03)

Saved verbatim by the parent session (Claude Code, Opus 5.5) from the stdout of `tools/codex/bin/codex-readonly` (model reported in the session log: `gpt-6.1-sol`). Prompt: `prompt.md` in this folder.

---

**Verdict:** I found no mathematical error in the proof of Theorem 1.1 after reading §§2–4, independently recomputing all 32 compact-range endpoint checks and the requested Bernstein coefficients, and checking the reduction algebra and boundary indices. The Lean statements express the paper’s theorem and scalar proposition faithfully under \(0<p_i<1\); the inspected proof sources contain proofs of the required inputs rather than additional hypotheses assuming them. This is a source audit, **not an independent Lean replay**: I did not run `lake build` or obtain fresh axiom reports. The concrete defects are incorrect decimal-enclosure output, an overbroad description of the Newton result formalized, and incomplete release metadata.

**Findings**

1. **Location:** `scripts/verify_pb_cue_threshold.py`, lines 49, 56–57; its docstring’s claim of “rigorous enclosures.”  
   **Severity: error — numerical reporting, without affecting the threshold proof.**

   Both rational endpoints are converted to floating point and rounded to nearest. The program consequently prints singleton intervals such as
   \[
   V_{1996}\in[0.999965435231604,\;0.999965435231604].
   \]
   The printed decimal lies outside the rational enclosure computed internally. I checked this by converting the printed decimal back to an exact `Fraction`; the same defect occurs for \(V_{1998}\) and \(V_{10000}\).

   Valid outward-rounded enclosures are
   \[
   \begin{aligned}
   V_{1996}&\in[0.999965435231604,\;0.999965435231605],\\
   V_{1998}&\in[1.000066908642300,\;1.000066908642301],\\
   V_{10000}&\in[1.163238438868,\;1.163238438869].
   \end{aligned}
   \]
   The exact comparison `hi96 < 1 < lo98` passes, so the manuscript’s threshold remains supported.

   **Proposed fix:** Print exact rational endpoints or round lower endpoints downward and upper endpoints upward using integer arithmetic. Apply the same correction to the reciprocal-bound output. Replace the floating-point validation at line 60 with rational comparisons if the whole check is described as exact.

2. **Location:** Manuscript proof outline, lines 413–416; `PBDeduction/Statement.lean`, description of G3e; `LOCAL_REPLAY.md`, phrase “strict Newton inequalities.”  
   **Severity: imprecision — formalization scope.**

   `pb_strict_lc` proves
   \[
   f_{k-1}f_{k+1}<f_k^2.
   \]
   It does **not** state Newton’s stronger binomially normalized inequality, or the quantitative ULC bound (1.10). `PBLogConcave.lean` establishes strict log-concavity directly by convolution induction.

   Thus “including the results of … Newton … that it uses” is defensible only if it means the strict-log-concavity consequence used in the main proof. “Strict Newton inequalities” suggests a stronger formalized result than the source supplies. The supplement’s narrower wording, “strict log-concavity of the pmf,” is accurate.

   **Proposed fix:** Say explicitly that Lean proves the strict-log-concavity consequence needed for Theorem 1.1; avoid calling that theorem the formalization of Newton’s inequalities.

3. **Location:** Manuscript line 71 and the all-integer cubic inequalities at `eq:HJ-left`–`eq:HJ-right`.  
   **Severity: minor — zero-extension convention.**

   The initial convention assigns only \(f_{-1}=f_{n+1}=0\), whereas the cubic inequalities also use \(f_{-2}\), \(f_{n+2}\), and are stated for every integer. The later discussion mentions Hillion–Johnson’s full zero extension, and the Lean definition implements it correctly, so this is not a substantive proof gap.

   **Proposed fix:** Initially define \(f_k=0\) for every \(k\notin\{0,\ldots,n\}\).

4. **Location:** Supplement description, line 837, `ZENODO-DOI-PENDING`; `CERTIFICATE.md`, “Scope and archive status.”  
   **Severity: minor — incomplete release metadata.**

   The manuscript presents an archive URL that is still a placeholder. The manifest acknowledges that no immutable DOI exists. The local archive is present and internally consistent, but that does not establish public availability.

   **Proposed fix:** Supply the deposited supplement’s actual persistent link before publication, or describe the archive as accompanying the submission until one exists.

**What you checked and how**

- **Read-only scope.** I read the manuscript, manifest, three requested verification programs, certificate checker, archive builder, Lean definitions and statements, and the principal proof files. No repository files were modified.

- **§4 preamble and compact-range proof.** I checked the substitutions
  \[
  H=\delta^{-1}-1,\qquad K<H\le K+1,\qquad
  Q(H)=\frac{(3H+4)(H+1)}4,
  \]
  the product expression for \(L_r\), componentwise monotonicity, and the endpoint comparison’s direction.

  The integer endpoints are handled correctly: at \(H=m+1\), \(K=m\); at \(H=m\), \(K=m-1\). The proof uses the appropriate half-open cells. Evaluating \(A_m\) at the excluded endpoint \(\delta=1/(m+1)\) is legitimate for the monotonicity comparison, because the weights there remain nonnegative.

- **All 32 endpoint ratios independently recomputed.** I used a separate `fractions.Fraction` calculation, forming the weights directly and evaluating
  \[
  A=\Bigl(\sum_iw_i\Bigr)\Bigl(\sum_i i^2w_i\Bigr)
      -\Bigl(\sum_i iw_i\Bigr)^2.
  \]
  All 32 ratios were strictly greater than one. For \(K=4,\ldots,15\), their rounded-down values were exactly
  \[
  1.072,\ 1.439,\ 1.763,\ 2.050,\ 2.304,\ 2.531,\
  2.735,\ 2.920,\ 3.089,\ 3.243,\ 3.385,\ 3.516.
  \]
  Among the 20 \(K=3\) pieces, the minimum was on the last piece:
  \[
  \frac{65461}{64950}=1.007867\ldots,
  \]
  supporting the printed \(1.0078\).

- **§4.2 completely checked.** I checked the \(C_{r+1}/C_r\) identity, positivity of its denominator, symmetrization, \(A_{\rm sym}=ST\), triangular-cell coverage, \(J\le K\), the product lower bound, nonnegativity needed when multiplying the bounds, and the quartic expansion.

  Independently recovering the quartic’s coefficients from exact rational evaluations gave the five printed \(J=5\) Bernstein coefficients:
  \[
  2360,\quad7500,\quad25055/2,\quad254205/16,\quad31115/2.
  \]
  Separate exact checks at \(J=6,7,10,100\) reproduced all five identities
  \(\beta_i=\mu_i\pi_i(J-6)/2880\). The supplied SymPy checker additionally verified those identities symbolically for arbitrary \(u\).

  The claimed failure of symmetrization at \(H=7/2\) was reproduced:
  \[
  ST-Q=-\frac{3237589}{1882384}<0.
  \]

- **§§2–3 deductions and boundaries.** I checked:
  - positive denominators in deriving recurrence (2.4);
  - both directions of reciprocal bound (2.5), including the missing endpoint directions;
  - Lemma 2.1’s induction and endpoint exclusion;
  - the shift from the deficit center \(D\) to the mass center \(c=D-1\);
  - the telescoping right and left ratio products;
  - support inclusion for the entire variance window;
  - the pairwise variance restriction;
  - the maximal-mass proof by uniform smoothing;
  - both strict inequalities in the contradiction deducing Theorem 1.1.

  I found no off-by-one error or invalid inversion in these steps. The manuscript’s cubic inequalities agree with equations (78)–(79) of [Hillion–Johnson’s Appendix A](https://arxiv.org/pdf/1303.3381).

- **Additional exact tests.** A separate deterministic test generated 84 rational Poisson–binomial laws with \(n=4,5,8,12,20,40,80\), including symmetric binomials. It checked positivity, both cubic inequalities, the reciprocal bound, the pairwise variance identity and maximal-mass bound. Of these, 63 had \(V\ge1\), and 30 had \(\delta_D<1/4\); all 251 tested radii satisfied the propagation and mass bounds. All assertions passed. These tests supplement the algebraic audit; they do not prove universal claims.

- **Formal statement fidelity.** I read both projects’ `Defs.lean` and `Statement.lean`:
  - `pbPmf` is precisely the generating-polynomial coefficient for nonnegative indices. Its guard prevents negative integers from being mapped to coefficient zero through `toNat`.
  - `pb_zero_out` proves vanishing above \(n\).
  - `pbVar` is \(\sum_i p_i(1-p_i)\), and `pb_pairVar` supplies the variance identity.
  - `IsFirstDescent` specifies a strict descent and excludes every earlier strict descent. It is equivalent to the paper’s minimum definition.
  - `pb_first_descent_exists` separately proves existence under the theorem’s hypotheses.
  - Positivity on the support prevents division by zero in `deficit` at \(D\).
  - `PBScalar.scalar_inequality` uses the correct products, double sum, and characterization
    \((K+1)\delta<1\le(K+2)\delta\). Its hypotheses force \(K\ge3\).

  On this definition-and-signature audit, neither formal statement is weaker than its corresponding paper statement. The final theorem constructs the needed hypotheses from proved Poisson–binomial results.

- **Escape-hatch search.** Searches over all `.lean` files in both requested result projects, and the earlier conditional project, found no occurrences of `sorry`, `admit`, `axiom`, `native_decide`, `implemented_by` or `ofReduceBool`. An expanded search also found no `trustCompiler`, `unsafe`, `skipKernelTC`, external declarations or elaborator hooks. Compact checks use `decide +kernel`. The scalar certificate files contain recursion-depth settings, which do not disable kernel checking.

- **Computational and archive claims.** All four read-only checker runs passed. The second certificate checker reported 13 cells, 275 positive coefficients and the manifest’s payload digest. I also regenerated the full certificate **in memory** and obtained equality with the stored JSON. All nine whole-file hashes listed in the manifest matched. The existing ZIP contains 76 members; its internal hashes and external checksum passed, and all mapped archive members matched the current repository files. The deduction project’s scalar copy matched the separate scalar project. The older Bernstein-route comments are stale but explicitly identified as such in `CERTIFICATE.md`.

- **Other manuscript mathematics.** I read and checked the binomial \(1/3\) example, support corollary, tilted-law modal upper bound, ULC counterexample, balanced-family formulas, and CUE variance derivation. The combined tilt bound uses the indices appearing in [Pitman’s equation (20)](https://www.stat.berkeley.edu/~pitman/453.pdf). The CUE program certified \(V_{1996}<1<V_{1998}\), and its approximate values support the manuscript’s \(0.215\) versus \(0.0004\) comparison. I found no additional mathematical defect in these arguments on that reading.

**What you could not check**

- **I could not run `lake build`.** I also did not execute Lean’s axiom-printing files. Successful compilation and dependence on only `propext`, `Classical.choice` and `Quot.sound` therefore remain author-reported claims, although the inspected sources are consistent with them.
- I did not rebuild the ZIP or write generator outputs; I checked the existing archive and regenerated certificate data in memory.
- I did not independently establish novelty, exhaustively verify every bibliographic attribution, or authenticate the historical AI-review and provenance claims in the acknowledgements.
- I cannot certify the requested “different model family” independence. This report should not be represented as evidence of that provenance.
