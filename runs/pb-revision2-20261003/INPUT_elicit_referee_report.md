# Referee report: “A variance-scaled Turán inequality at the first descent of a Poisson–binomial mass function”

**Version reviewed:** the 12-page PDF supplied by the user, including its displayed supplementary-material section. The separate code/certificate archive referred to in the paper was not supplied; its DOI is still shown as `ZENODO-DOI-PENDING` in this version. [^1]

## Summary and contribution

The paper studies the normalized Turán deficit \(\delta_k=1-f_{k-1}f_{k+1}/f_k^2\) of a Poisson–binomial probability mass function and asks whether its variance controls the deficit near the mode. Its main result is the universal bound \(V\delta_D\ge 1/4\) at the first descent index \(D\), under \(V\ge1\). A variance-one binomial family gives an upper bound of 1/3 on the best possible universal constant, leaving its exact value open. The paper then propagates the first-descent estimate across the support using reciprocal-deficit inequalities and gives an application to eigenvalue counts for Haar unitary matrices. [^1]

The proof strategy is coherent and unusually explicit: Hillion–Johnson cubic inequalities constrain adjacent deficits; these yield local lower bounds on masses around the rightmost mode; the pairwise variance identity and an integer-valued maximal-mass inequality reduce the result to a scalar inequality. Section 4 establishes that inequality through Bernstein-basis positivity, with exact-rational calculations reported for bounded intervals and a symbolic family for the tail. [^1]

## Strengths

1. The theorem is cleanly stated with its nondegeneracy condition, and the paper supplies a concrete family showing that the constant cannot exceed 1/3. This appropriately distinguishes the proved lower bound from the unresolved sharp constant. [^1]
2. The argument’s dependencies are laid out in a readable sequence: a local recurrence, mass lower bounds, a variance comparison, then a scalar certificate. The paper also identifies three potentially lossy steps rather than implying the constant is sharp. [^1]
3. The paper addresses boundary cases in its standing reductions, including zero/one Bernoulli parameters, support endpoints, and the condition needed for the first-descent index to exist. [^1]
4. The exact-arithmetic approach is appropriate for the sign claims: the manuscript reports an independently written checker and makes clear which finite coefficient expansions establish the bounded range. [^1]

## Overall assessment

The central result appears to have a plausible and well-motivated proof architecture, and the manuscript gives enough detail to follow the reduction and the broad strategy. The principal publication-readiness concern is reproducibility of the computer-assisted scalar inequality: the version supplied points to a pending archive DOI, so the programs and full certificate data that the proof explicitly relies on cannot be independently checked from the reviewed materials. [^1] This is a material but readily repairable dissemination issue, not evidence that the theorem is false.

## Major comments

### 1. Make the computational certificate and code retrievable — major, Section 4 and Supplementary Material

The proof of Proposition 3.1 depends on positivity of 275 Bernstein coefficients over the bounded range and on symbolic coefficient formulae for the unbounded range. The text says that the exact coefficients, generating code, independent checker, and manifest are in an archive, but the DOI is `ZENODO-DOI-PENDING`; none of those files are part of the supplied PDF. [^1]

**Why it matters:** these computations are a load-bearing part of the proof, rather than an incidental illustration. Without the artifact, a reader cannot verify the coefficient-generation formulas, reproduce the expansions, or confirm the independent checker’s scope.

**Remedy:** replace the placeholder with a persistent public DOI before publication, or include the scripts, exact coefficient data, and run instructions as a journal supplement. Identify the exact archive version/hash used for the manuscript. The authors should also state whether the checker independently reconstructs the underlying polynomials as well as the Bernstein conversion, and document any shared assumptions or code.

### 2. Clarify the status and scope of the Lean formalization — major-to-moderate, Supplementary Material

The supplement calls the Lean formalization “conditional” on a one-sided form of recurrence (2.4), and lists several lemmas and deductions covered. [^1]

**Why it matters:** a reader could otherwise interpret the mention of Lean as formal verification of Theorem 1.1, whereas the stated condition indicates that at least part of the mathematical bridge remains outside the formalized result.

**Remedy:** provide the formalization and state precisely which propositions are machine-checked, which hypotheses are assumed, and which steps—including the scalar inequality and its computational certificates—are not formalized. Avoid wording that could be read as end-to-end formal verification unless that is what the code establishes.

## Minor comments

1. **Abstract and introduction:** make explicit that the best constant is not identified; the paper does state this in the body, but readers may otherwise overread the 1/4 theorem as sharp. [^1]
2. **Section 4.1:** the sentence reporting that “all 275” coefficients are positive is concise but not independently informative without the archive. After publication of the artifact, provide a short machine-readable index by interval/polynomial and coefficient count, or a compact checksum, to make the claim easier to audit. [^1]
3. **Example 1.6:** the numerical comparison at \(N=10^4\) is useful, but state directly that these are lower-bound values being compared, not the actual deficit \(\delta_D\). The surrounding prose mostly implies this, but the distinction matters for interpretation. [^1]
4. **Supplementary-material reference:** the placeholder should be searched globally in the final source and PDF; the same pending DOI appears in the supplementary-material description. [^1]

## Internal-consistency findings

The abstract’s main bound, the stated interval for the optimal constant, and the corollary’s reciprocal-deficit propagation are consistent with the theorem and proof outline as presented. The manuscript explicitly distinguishes the lower bound at the first descent from the upper-bound example and acknowledges that sharpness is unresolved. [^1]

I did not identify a definite algebraic contradiction in the displayed main proof. However, the computational proposition was not independently reproduced: the paper’s claimed certificate files were unavailable in the supplied version. Thus this is a check of the manuscript’s stated proof structure and displayed derivations, not independent validation of all exact-arithmetic coefficients.

## Recommendation

**Major revision (publication-readiness, not a demonstrated mathematical flaw).** The central issue is to make the load-bearing exact-arithmetic code and certificates persistently available and clarify precisely what the Lean formalization verifies. If the complete archive is supplied and matches the manuscript’s claims, the concerns raised here may be resolved without changing the theorem or proof strategy.

## Reviewer confidence and limits

**Medium.** The supplied PDF was inspected in full, including the main proof, displayed certificate formulas, supplement description, and references. Confidence is limited specifically by the unavailable separate code/certificate archive; no independent execution of the asserted computer checks was possible. The novelty claim is not independently assessed against a broad literature search.

[^1]: variance-scaled-turan-first-descent.pdf.