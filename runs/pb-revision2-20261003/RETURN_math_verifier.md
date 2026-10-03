<!-- Saved verbatim by the parent session from the subagent's final message (Claude Opus via responsibility-verifier, read-only, no shell; checked the passages added in the second revision of paper/poisson_binomial/variance-scaled-turan-first-descent.tex, 2026-10-03). Subagents cannot write report files. The parent then applied A1-A4 and B5-B10 (B9 in full: k->ell, Gamma for the Gram matrix, mu(t)/V(t) referenced to the proof of Proposition 1.4). -->

REPORT: independent check of the new and changed passages in variance-scaled-turan-first-descent.tex

No RESPONSIBILITY NOTICE. I checked passages (a) to (l) by hand against the four sources on disk, plus Johnson 2017 (`notes/literature/arxiv_1507_06268.txt`) and Dümbgen–Wellner (`notes/literature/arxiv_1910_03444.txt`). None of the new proofs contains a false mathematical step. I found four minor wrong statements, all at the wording level, and some attribution and presentation gaps. One of those gaps (item 5, Example 1.6) is likely to draw a referee comment.

Line numbers refer to `/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/paper/poisson_binomial/variance-scaled-turan-first-descent.tex`. Equation numbers follow the compiled PDF: (2.4) recurrence, (2.5) reciprocal bound, (2.8) deficit bound, (1.5)/(1.6) Corollary 1.3, (1.7) ULC bound.

## A. Wrong statements (all minor), with fixes

**1. Proof outline, lines 345–348 (passage l).** The text says δ ≤ 1/17 is handled "with coefficients that are polynomials with positive coefficients in an integer parameter". That is not true of the J = 5 interval, 16 ≤ H ≤ 21, which is 1/22 ≤ δ ≤ 1/17. That interval is proved with five positive rational Bernstein coefficients (lines 718–724). Only J ≥ 6 (H ≥ 21, δ ≤ 1/22) uses the polynomials π_i(u).
- Fix: "with positive rational coefficients for 1/22 ≤ δ < 1/4, and for δ ≤ 1/22 with coefficients that are polynomials with positive coefficients in an integer parameter."

**2. Introduction, lines 95–107 (outside the listed targets; possibly unchanged).** Dümbgen–Wellner Prop. 2 (source line 283) says (x+1)b(x+1)/b(x) strictly decreases. That gives (k+1)q_{k+1} < k q_k, so δ_k = 1 − q_{k+1}/q_k > 1/(k+1). So the D–W result does give a lower bound on δ_k. It depends on k alone, not on n or on the p_i, and it is weaker than (1.7). Two sentences therefore overstate:
- "do not control δ_k" (line 107);
- "The lower bounds on δ_k that we know of depend on n … or on the success probabilities" (lines 95–97).
- Fix: say they give no *variance-scaled* lower bound, and optionally note that D–W gives δ_k > 1/(k+1), which is implied by (1.7) and also degenerates at D in the family after Prop 1.5.

**3. Pitman paragraph, line 317 (passage k).** Pitman defines λ(x) only for ℓ < x < r, here 0 < x < n (source lines 347–350). So θ(n) is undefined.
- Fix: "for integers k with EW < k < n".

**4. Degeneracy paragraph, lines 278–280 (passage i).** The equation 2mε(1−ε) = 1 has no solution when m = 1 (it would need ε(1−ε) = 1/2 > 1/4).
- Fix: state m ≥ 2 and take ε_m the smaller root. At m = 2, ε = 1/2 and the law is Bin(4, 1/2).

## B. Presentation and attribution issues, with fixes

**5. Example 1.6, lines 297–301 (passage j): attribution gap.** Meckes–Meckes Prop 3(1) (source lines 155–167) only says that independent Bernoulli variables ξ_1, …, ξ_N exist with N_A equal in law to their sum. It does not identify the parameters. As written, the semicolon makes it look as if the cited result supplies the Gram-matrix eigenvalues. The identification is correct, but it needs one sentence:
- HKPV Theorem 7 (source lines 219–248) gives Bernoulli(λ_k) for a trace-class determinantal kernel. Apply it to the process restricted to A = (0, π), with kernel 1_A Π 1_A, where Π is the rank-N projection onto span{e^{ijθ}} in L²(dθ/2π).
- The operators 1_A Π 1_A and Π 1_A Π have the same nonzero spectrum (AB versus BA).
- The matrix of Π 1_A Π in the orthonormal basis e^{ijθ} is (2π)^{-1}∫_0^π e^{i(k−j)θ}dθ. That is the conjugate (equivalently the transpose) of the paper's K, so it has the same eigenvalues.

On the specific question whether this is the right matrix for the Meckes–Meckes kernel: yes. Their display sin(N(x−y)/2)/sin((x−y)/2) is the kernel with respect to dθ/2π, since K_N(x,x) = N and the trace must be N. It equals Σ_{j<N} e^{ij(x−y)} up to the unimodular factor e^{−i(N−1)(x−y)/2}. That factor is conjugation by a unitary multiplication operator, so it changes neither the determinantal process nor the spectrum.

Also add that all N eigenvalues lie strictly in (0,1). K and I − K are Gram matrices of the linearly independent functions e^{ijθ} on the upper and lower half-arcs respectively, so both are positive definite. This is what justifies n = N in the (1.7) comparison at line 313.

**6. Example 1.6, lines 310–311.**
- "The sequence V_N increases" has no proof. One line suffices: V_{N+1} − V_N = (2/π²) Σ_{d odd > N} d^{-2} > 0. The formula holds for every N.
- "First reaches 1 at N = 1998" holds among even N only. The same asymptotic formula gives V_1997 ≈ 1.000016 (my estimate from the asymptotic, not an exact computation). Say "among even N".
- "For larger even N": write "for even N ≥ 1998".

**7. Lines 207–211 (passage f).** The conclusion is correct, but the stated route is loose. Complementing a single summand does not preserve the deficits. The argument should be: for independent Bernoulli sums X and Z = Σ_{l≤m} Z_l, the sum X − Z + m = X + Σ(1 − Z_l) is a Bernoulli sum with variance Var X + Var Z. Its pmf is the pmf of X − Z translated by m, so D, δ_D, the δ_k − δ_D offsets and E(·) all shift together.

**8. Prop 1.5 proof, lines 269–272 (passage h).** The envelope j²2^{−|j|} does not dominate the normalizing sum: it vanishes at j = 0, where the weight is 1. Use (1+j²)2^{−|j|}, or monotone convergence, which the paper already sets up. Note too that the mean is m by symmetry, so Var = Σ j²w_j / Σ w_j.

**9. Notation.**
- k is used in two senses in lines 317–325: as Pitman's index and as the integer in Bin(n, p) with np = k+1−ε.
- μ(t) and V(t) are defined inside the proof of Prop 1.4 and then used outside it at lines 318–320.
- K denotes both the window size in §2 and the Gram matrix in Example 1.6.

**10. Optional, lines 424–427 (passage c).** The case split on δ_k < 1 versus δ_k = 1 is unnecessary. Since δ_k and δ_{k+1} are positive, dividing δ_{k+1}(1−δ_k) ≤ δ_k by δ_kδ_{k+1} gives 1/δ_{k+1} ≥ 1/δ_k − 1 directly.

## C. Passage-by-passage verdicts

**(a) Section 1 opening: nothing wrong.** exp(−C_k) = f_{k−1}f_{k+1}/f_k², so δ_k = 1 − e^{−C_k} wherever f_{k−1}f_{k+1} > 0 (that is, 1 ≤ k ≤ n−1). δ_0 = δ_n = 1 because f_0, f_n > 0 when 0 < p_i < 1. For a normal density, −(log φ)'' = 1/V, and the unit-step second difference of log φ is exactly 1/V. (The one exception is item 2 above, in the adjacent literature paragraph.)

**(b) Section 2.1, D ≥ 2: nothing wrong.** If f_1 < f_0, then D = 1 and 0 is the unique mode. Darroch gives EW < 1, and V = Σp_i(1−p_i) < Σp_i = EW < 1, contradicting V ≥ 1. The nonemptiness argument (lines 374–381) also checks out.

**(c) Reciprocal inequality (2.5): nothing wrong.**
- HJ (78) is g_{k−1}g_k² − 2g_{k−1}²g_{k+1} + g_kg_{k+1}g_{k−2} ≥ 0 and HJ (79) is g_k²g_{k+1} − 2g_{k+1}²g_{k−1} + g_{k+2}g_kg_{k−1} ≥ 0, both for all k ∈ ℤ (source lines 1435–1448). Expanding shows they are exactly (2.2) and (2.3).
- The divisions by f_{k−1}f_k² and f_{k+1}f_k² need positive masses, which hold for 1 ≤ k ≤ n−1. The identities at the edges use f_{−1} = f_{n+1} = 0, so δ_0 = δ_n = 1.
- Positivity: (1.7) follows from ULC(n). I rederived it: 1 − k(n−k)/((k+1)(n−k+1)) = (n+1)/((k+1)(n−k+1)). HJ (41) and the statement that Bernoulli sums are ULC(n) are at source lines 889–908.
- Index ranges: the second inequality of (2.4) gives the upward direction for j ∈ [1, n−1]. The first, applied at k = j+1, gives the downward direction for j ∈ [0, n−2].
- End pair (0,1): the nontrivial direction 1/δ_1 ≤ 2 comes from the first inequality at k = 1 (needs n ≥ 2). The other direction, 1/δ_1 ≥ 0 = 1/δ_0 − 1, is trivial.
- End pair (n−1, n): handled symmetrically.
- For n = 1 the bound is trivial; in any case V ≥ 1 forces n ≥ 4.

**(d) Lemma 2.1: nothing wrong.**
- Base case: D ≠ n because δ_n = 1 > δ.
- Step: all pairs on the path D, D±1, …, D±r lie in [0, n], so telescoping (2.5) gives 1/δ_{D±r} ≥ 1/δ − r. Since (r+1)δ < 1 is equivalent to 1/δ > r+1, this is greater than 1, which excludes the indices 0 and n.
- Inverting gives exactly δ/(1−rδ) in (2.8). The bound δ/(1−rδ) < 1 is again equivalent to (r+1)δ < 1.
- So D−K ≥ 1 and D+K ≤ n−1. Hence c−K ≥ 0 and c+K ≤ n−2, which covers every index later used: q_{D−K}, the deficits δ_{D−s} for s ≤ K, and q_{D+j} for j ≤ K−1.

**(e) Corollary 1.3: nothing wrong.**
- Telescoping (2.5) between any two indices in [0, n] gives 1/δ_k ≤ 1/δ_D + |k−D| ≤ 4V + |k−D|.
- Darroch gives |c − EW| < 1, so |D − EW| < 2 and |k−D| < |k−EW| + 2.
- For (1.6): η_− < 1 and η_+ ≥ 0 give η_+ ≤ η_+/η_−, which also covers D = n, where η_+ = 0.

**(f) Window claim: correct.** For V ≥ 1, 4V + A√V + 2 ≤ (6+A)V. The extension to differences of Bernoulli sums is correct; only the stated route needs fixing (item 7).

**(g) Proposition 1.4: nothing wrong in any sub-step.**
- The tilted parameters p_ie^t/(1−p_i+p_ie^t) are correct, and μ' = V(t).
- The bound e^t/(1−p+pe^t)² ≥ e^{−|t|} holds in both cases: for t ≥ 0 because 1−p+pe^t ≤ e^t, and for t ≤ 0 because 1−p+pe^t ≤ 1.
- c ≥ 1 gives q_c ≥ 1 by minimality of D, and D ≤ n makes q_{c+1} < 1 defined. Hence a ≤ 0 < b.
- Tilting preserves the deficits, and δ_k > 0 makes q_k strictly decreasing. So f^{(a)} has ratio 1 exactly at k = c (modes exactly c−1 and c), and f^{(b)} has ratio 1 exactly at k = c+1 (modes exactly c and c+1).
- Pitman's wording (source lines 290–294) is: "a PF sequence has either a unique index m or two consecutive indices m such that a_m = max … such a mode m differs from the mean μ by less than 1". Applied to each of the two modes, this gives the open intervals μ(a) ∈ (c−1, c) and μ(b) ∈ (c, c+1). Pitman's precise version (13) is consistent with this.
- The integral ∫_a^b e^{−|t|}dt equals 2 − e^a − e^{−b}.
- The final step is right: (1 − 1/q_c)(1 − q_{c+1}) ≥ 0, and 1 − q_{c+1}/q_c = δ_c by (2.1), since 1 ≤ c ≤ n−1.
- The strict "< 2/V" comes from Darroch's strict inequality.
- At line 242: V/(4V+1) ≥ 1/5 for V ≥ 1.

**(h) Proposition 1.5: nothing wrong apart from the envelope nit (item 8).**
- g_k/C(2m,k) ∝ 2^{−|k−m|} is log-concave. Multiplying by the strictly log-concave binomial row gives strict log-concavity.
- Symmetry plus strict log-concavity gives the unique mode m and first descent m+1.
- The closed form: the ratio is ((m−1)/(m+2))/(m/(m+1)), so δ_{m+1} = (2m+1)/(m(m+2)). This also holds at m = 1, where δ_2 = 1.
- Each factor (m−i+1)/(m+i) increases to 1 as m grows.
- The limit: Σ_ℤ 2^{−|j|} = 3 and Σ_ℤ j²2^{−|j|} = 12, so the variance tends to 4.
- The remark: δ_m = 1 − (m/(m+1))²/4 > 3/4, so 1/δ_m < 4/3. That is correct.

**(i) Degeneracy paragraph: correct for m ≥ 2 (item 4).**
- V = 1, and the law is symmetric about m: 2m − W has the same parameter multiset. So D = m+1, and the right-hand side of (1.7) at D equals (2m+1)/(m(m+2)).
- Johnson: Def. 1.2 is E(x) = (V(x)² − V(x−1)V(x+1))/(V(x)V(x+1)) (source lines 65–68). Lemma 5.1 and Ex. 5.2(2) give c = (Σ p_j/(1−p_j))^{−1} (source lines 530, 549–552). Rearranged, this is exactly the displayed bound, for k ≤ n−1.
- The odds sum is at least m(1−ε)/ε ≈ 2m², so its reciprocal tends to 0.

**(j) Example 1.6: the mathematics is correct. Attribution gaps are items 5 and 6.**
- For odd d, K_{jk} = i/(πd), so |K_{jk}|² = 1/(π²d²). K_{jk} = 0 for even d ≠ 0, and K_jj = 1/2.
- So V_N = N/2 − N/4 − (2/π²)Σ(N−d)/d², as stated.
- I rederived the asymptotic for even N. With Σ_{d odd ≥ N+1} d^{−2} = 1/(2N) − 1/(6N³) + … and Σ_{d odd < N} 1/d = ½log N + ½log 2 + γ/2 + 1/(12N²) + …, I get V_N = (log N + γ + log 2 + 1)/π² − 1/(6π²N²) + O(N^{−4}). So O(N^{−2}) is correct.
- Symmetry: U ↦ −U preserves Haar measure, and an eigenvalue at ±1 has probability 0. So X_N has the same law as N − X_N, and strict log-concavity gives c = N/2 and D = N/2+1.
- Threshold: the closed form crosses 1 at N ≈ 1996.7. It gives V_1996 ≈ 0.999965 and V_1998 ≈ 1.000067, consistent with the script's 0.99997 and 1.00007.
- At N = 10⁴: V ≈ 1.163238, so 1/(4V) ≈ 0.21492. The (1.7) value is 10001/25,010,000 ≈ 0.000400. Both match the paper.

**(k) Pitman paragraph: correct apart from the k-range (item 3).**
- Pitman (21), λ(k) < a_k/a_{k+1}, gives f_{k+1}/f_k < 1/θ(k).
- For t ≥ 0, (1−p+pe^t)² ≥ 1, so V(t) ≤ e^tV.
- Integrating, k − EW ≤ V(θ(k) − 1), so f_{k+1}/f_k < V/(V + k − EW).
- "Stronger far out": (1.5) implies f_{k+1}/f_k < (4V−1)/(4V+k−D) by a telescoping product. That is about 4V/k, against Pitman's about V/k.
- Binomial example: (n+1)p = k+1−ε+p lies in (k, k+1) when p < ε, so the mode is k, D = k+1, D − EW = ε, and Pitman's bound tends to 1. V ≥ 1 needs k ≥ 1.

**(l) Abstract and outline.** The abstract has nothing false. Every claim matches a proved result: the 1/4 bound, the 1/3 family, the support bound, 1/(4V+1) ≤ δ_c < 2/V, the ULC counterexample, and the CUE rates for large N. In the outline, the only problem is item 1. The heuristic (b_r ≈ e^{−r²/2H}, so ST ≍ H²) and the M² ≥ (1+12V)^{−1} rearrangement are correct.

## D. What I verified directly, what I inferred, and what would reverse the result

**Verified at source:**
- HJ (78)–(79), for all k ∈ ℤ: `arxiv_1303_3381.txt` lines 1435–1448.
- HJ Def. 3.11, (41), and Bernoulli sums being ULC(n): lines 889–908.
- HKPV Thm 7: `arxiv_math_0503110.txt` lines 219–248.
- Meckes–Meckes kernel display and Prop 3(1) wording: `arxiv_1612_08100.txt` lines 143–175.
- Johnson Def. 1.2, Lemma 5.1, Ex. 5.2(2): `arxiv_1507_06268.txt` lines 65–68, 530, 544–552.
- Dümbgen–Wellner Prop. 2: `arxiv_1910_03444.txt` line 283.
- Pitman's prose statement of Darroch's rule: `pitman_1997_coefficients_real_zeros.txt` lines 290–294.

**Inferred:** Pitman's displays (13) and (14)–(21) are garbled in the text extraction (lines 297–365). My reconstruction is: λ(k) < a_k/a_{k+1} < λ(k+1); λ(x) defined for ℓ < x < r; and a two-mode window [k + 1/(k+2), k+1 − 1/(n−k+1)] that lies strictly inside (k, k+1). It is consistent with the surrounding prose, but I did not read it verbatim. Check it against the PDF before relying on it.

**Not executed:** I could not run code. V_1996, V_1998 and the N = 10⁴ figures were checked only against the closed-form asymptotic (the O(N^{−2}) term is about 4×10^{−9}), not recomputed exactly.

**What would reverse my verdicts:**
- A restricted k-range in HJ (78)–(79). Ruled out: the source says "for all k ∈ ℤ".
- A non-strict Darroch bound in the two-mode case. Ruled out by Pitman's "less than 1" and by (13).

**Flag for the parent, not checkable locally:** supplement lines 767–772 describe the Lean formalization of Lemma 2.1. Confirm that the description still matches the lemma's current statement: membership in {1, …, n−1} plus the bound (2.8), now proved via (2.5).

**Sources read:** `/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/notes/literature/` (`pitman_1997_coefficients_real_zeros.txt`, `arxiv_1303_3381.txt`, `arxiv_math_0503110.txt`, `arxiv_1612_08100.txt`, `arxiv_1507_06268.txt`, `arxiv_1910_03444.txt`), plus the extracted PDF text at `/private/tmp/claude-502/-Users-brettreynolds-projects-LLM-CLI-projects-papers-queue-erdos-problem-993/eec3f615-7ddb-4430-8ddb-9f052413b82b/scratchpad/review2/fullpaper.txt`.

This check shares the parent's model family. Don't count it as an independent vote.
