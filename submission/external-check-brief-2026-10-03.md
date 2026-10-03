# Brief for an external check of the Poisson–binomial paper
<!-- SUMMARY: Paste-ready brief for a cross-family check (ChatGPT Pro) of variance-local-log-concavity-poisson-binomial.pdf before ECP submission; lists the arguments that only Claude has checked by hand · status: ready to send · updated: 2026-10-03 -->

Attach `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.pdf` and, if the tool accepts it, `paper/poisson_binomial/poisson_binomial_certificate_supplement.zip`. Then paste the text below.

---

You are refereeing a short paper submitted to *Electronic Communications in Probability*. Please try to break it. A correct refutation of any claim is more useful than agreement.

The main result is Theorem 1.1: V·δ_D ≥ 1/4 for every Poisson–binomial law with variance V ≥ 1, where δ_k = 1 − f_{k−1}f_{k+1}/f_k² and D is the first index at which the pmf decreases. Exact-arithmetic scripts check every number in the paper. The written arguments below have been checked by hand only by models of one family, so they need an independent reader most.

1. **The reciprocal bound (2.5)**, |1/δ_{j+1} − 1/δ_j| ≤ 1 on the whole support. It is derived from the Hillion–Johnson cubic inequalities (2.2)–(2.4). Check the index ranges, including the two end pairs.
2. **Lemma 2.1** and its induction.
3. **Corollary 1.3.** Both inequalities in (1.5), which use Darroch's theorem; and (1.6).
4. **Proposition 1.4.** The upper bound δ_c < 2/V by exponential tilting. Check:
   - the tilted variance bound V(t) ≥ e^{−|t|}V;
   - that each tilted law has exactly two modes;
   - the strict use of Darroch's rule ("every mode differs from the mean by less than 1");
   - the integral bound.
5. **Proposition 1.5.** The ULC counterexample: strict log-concavity, the first-descent index m+1, the closed form of δ, and the variance limit 4.
6. **Example 1.6.** Eigenvalues of a Haar unitary matrix in a half circle. Check:
   - the Bernoulli representation (Meckes–Meckes Prop. 3(1); Hough–Krishnapur–Peres–Virág Thm. 7);
   - that the parameters are the eigenvalues of the Gram matrix Γ;
   - the exact formula for V_N, the O(N^{−2}) asymptotic, and the claim that V_N first reaches 1 at even N = 1998;
   - the symmetry argument giving D = N/2 + 1.
7. **The Pitman comparison** (the paragraph before Example 1.6): the derivation f_{k+1}/f_k < V/(V + k − EW).
8. **The main proof** (Sections 2–4). This means the reduction to the one-variable inequality A(δ) ≥ (3+δ)/(4δ²) and the proof of that inequality (Bernstein certificates for 3 < H ≤ 16; explicit polynomial identities for H ≥ 16). Spot-check whatever you can.
9. **Novelty.** Is any of the following already in the literature, in any form?
   - Theorem 1.1;
   - the support-wide bound δ_k ≥ 1/(4V + |k−D|);
   - the two-sided bound 1/(4V+1) ≤ δ_c < 2/V.

   Name the source if so.
10. **As an ECP referee**, would you recommend acceptance? What is the strongest objection?

For each problem, give the location (theorem, equation or paragraph), say whether it is an error or a presentation issue, and propose a fix. If you check a step and find it correct, say so; that is useful too.
