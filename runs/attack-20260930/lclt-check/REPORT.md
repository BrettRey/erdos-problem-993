# lclt-check (resumed): route-lclt claims A–C, the hybrid N1, and the Lemma L gap
<!-- SUMMARY: Subagent report, saved verbatim by the parent session on 2026-09-30 (subagents cannot write .md reports). Grade: sound, marginal gain. Claims A-C correct (exact on 5,445 trees n<=14); hybrid removes Fang's R^-3 floor but best explicit N1 ~10^495 (10^1218 on provable constants); Lemma L (V, K, T) unproved · status: final · updated: 2026-09-30 -->

Saved verbatim from the agent's final message. The agent could not write this file itself because the harness blocks subagent `.md` reports. Only the whitespace was changed.

**Overall grade: sound, marginal gain.**

- Claims A, B and C are correct. I re-derived them and verified them exactly on all 5,445 trees with n ≤ 14 plus 7 larger trees.
- c_* = 0.26420 is right, but only as a limit constant. As a uniform bound it is refuted.
- The hybrid does remove the R^-3 local-step floor. But the mean-shift term in the centroid CLT is paid at exponent 2/(1−a). With a ≥ 0.8, N1 goes back above 10^100. With the certified constants the best N1 is about 10^495 (Fang sharpened: 10^2070).
- No RESPONSIBILITY NOTICE: nothing here is a material risk.

The first attempt's files were re-verified, not trusted. Reruns are in `rerun/`, and the census and exact_large outputs came back byte-identical. HEAD is now d15121a (it moved during the day). All files are in this directory, and `manifest.yaml` records the model, HEAD and SHA-256 for 9 inputs and every output.

## 1. Claims A, B, C: re-derived and tested exactly. Grade: sound

**My derivation.**
- **Bijection.** Take T bipartite with parts L and R, both independent. Map I to (J, K) with J = I∩R and K = I∩L.
  - J can be any subset of R.
  - K can be any subset of A_J = {v ∈ L : N(v) ∩ J = ∅}.
  - So I_T(x) = Σ_{J⊆R} x^{|J|}(1+x)^{B_J}, with B_J = |A_J|. Call this A0. It is λ-free and holds for any finite bipartite graph.
- **(A)** Under the hard-core measure, J has weight λ^{|J|}(1+λ)^{B_J}/Z, and |K| given J is Bin(B_J, p) with p = λ/(1+λ).
  - This is Fang's conditioning: Lemma 6.1, eq. (6.2), p. 10 ("Given I ∩ R, let B be the number of available vertices of L").
  - The second difference is linear, so D_k := 2p_k − p_{k−1} − p_{k+1} = E_J h_{B_J}(k − |J|), where h_B = −Δ²b_B.
- **Fixed kernel.** On E0 = {B_J ≥ B0}, split off a J-measurable block of B0 available vertices. Then X = M′ + Y0, with Y0 ~ Bin(B0, p) independent of M′ = |J| + Bin(B_J − B0, p). This gives the exact decomposition
  - D_k = P(E0)·E[h_{B0}(k − M′) | E0] + (2a_k − a_{k−1} − a_{k+1}), with a_j := P(X = j, E0^c).
- **(B)** Abel summation: Σ_m h(k−m)(P(m) − Q(m)) = Σ_m [h(k−m) − h(k−m−1)](F_P(m) − F_Q(m)). The boundary terms vanish. So the absolute value is at most TV(h)·d_K(P, Q).
- **(C)** For every lattice law G:
  - D_k ≥ P(E0)[E h_{B0}(k−G) − TV(h_{B0})·d_K(M′|E0, G)] − (a_{k−1} + a_{k+1}).
  - Since a_{k−1} + a_{k+1} ≤ P(E0^c), route-lclt's −2P(E0^c) is valid but loose by at least a factor of 2.
  - The inequality holds for any G. The discretised normal only decides whether the bound is positive, which is also why the exact form is clean for Lean.
  - D_k > 0 gives p_k > (p_{k−1} + p_{k+1})/2 ≥ √(p_{k−1}p_{k+1}), hence i_k² > i_{k−1}i_{k+1}; λ cancels.

**Exact tests (Fractions and ints).** Every assertion passed: A0, A, A′, the exact decomposition, (B), and D ≥ cert_sharp ≥ cert_old.

| Tree set | Trees | A0 tree-sides | A (tree, side, k) | A′/B/C cases (every B0 = 1..max B_J) |
|---|---|---|---|---|
| 3 ≤ n ≤ 12 (rerun, identical to first attempt) | 985 | 1,970 | 5,466 | 30,992 |
| n = 13 (new) | 1,301 | 2,602 | 6,176 | 40,144 |
| n = 14 (new) | 3,159 | 6,318 | 17,436 | 122,052 |
| **Total n ≤ 14** | **5,445** | **10,890** | **29,078** | **193,188** |

- Tree counts per n match A000055.
- The 2-D (|J|, B_J) DP matched brute force in every one of the 10,890 tree-side cases.

Larger trees (A0, A, A′, B, C exact; three window values of k; B0 = ⌊(1−η)E B_J⌋ with η ∈ {0.3, 0.5}):
- Rerun, identical to the first attempt: H(9,2) n = 28, path n = 40, double star (14,14) n = 30, spider 5×6 n = 31.
- New: complete binary tree n = 63, caterpillar 20×2 n = 60, random recursive tree n = 60 (seed 993). These three use the DP only.

**How to read the certificate counts.** Every exact test ran at the tie fugacity λ = i_{k−1}/i_k. There, D_k = p_k(1 − i_{k−1}i_{k+1}/i_k²), so D_k > 0 is exactly log-concavity.
- So the identity tests are sharp at that λ, and A0 covers every λ.
- The `cert_old_pos` / `cert_sharp_pos` tallies (for example 14,078/30,992 at n ≤ 12) are **not** tests of the certificate as the route uses it (mean-matching λ, asymptotic). They are not failure rates.

**Sharp-tail diagnostic (float64; route-lclt's own KT5 code with only the tail term replaced).**
- Replacing −2P(E0^c) by −(a_{k−1} + a_{k+1}) removes every route-lclt KT5 failure:
  - n = 8: 6/53 → 0
  - n = 12: 65/1,682 → 0
  - n = 14: 378/8,717 → 0
  - n = 16 (new): 969/71,100 → 0; worst normalised margin 0.338
- The old counts reproduce route-lclt's exactly. Because it reuses route-lclt code, this is a consistency check, not an independent one.

**Kernel constants.**
- φ‴(x) = (3x − x³)φ(x), with zeros at 0 and ±√3. So ∫|φ‴| = 2(φ(0) + 4φ(√3)).
  - Arb (python-flint 0.9.0, 200 bits): [1.510013000130477132617541 ± 1.0e−25].
  - An independent `acb.integral` over the sign-definite pieces gives the same ball.
- c_* = φ(0)/∫|φ‴| = 1/(2 + 8e^{−3/2}) = [0.2641979111219313291032677 ± 1.8e−26]. Confirmed. The first attempt's `kernel_arb.txt` agrees and is superseded; its `kernel_constants.py` cannot run under the venv, which lacks mpmath.
- **Refuted as a uniform bound: TV(h_B)·v^{3/2} ≤ ∫|φ‴|.** In exact rationals at p = 12/13 (λ = 12):
  - B = 256 gives 1.51364, above the rational 1.51002 > ∫|φ‴|.
  - B = 22 gives a value in (1.579, 1.581).
  - Float scan over λ ∈ [1/3, 12] and B ≤ 1200: the supremum is 1.57996, at λ = 12, B = 22 (4.6% above the limit).
- So c_* is a limit constant only. An effective certificate needs a uniform C_TV (about 1.58 on the scanned range, unproved for B > 1200), giving c_eff = φ(0)/C_TV ≈ 0.2525. That is harmless for orders of magnitude, but it kills the naive inequality as a formal target.

## 2. Best explicit N1 for the hybrid. Grade: sound, marginal gain

**Constants and sources.** `n1_check.py` does not import audit-b code, and it reproduces the first attempt's `n1_hybrid.json` to within 0.1.

- **Certificate tolerance.** ε_K = c_*((1−η)ρ0)^{3/2} = 0.0459 at η = 0.3 and ρ0 = 0.4448.
  - ρ0 = 0.4448 is route-lclt's KT2 census minimum at n = 14 and is unproved.
  - Half of ε_K goes to the Kolmogorov chain, split over 4 terms, so each term must be ≤ t = 0.00574. "Fixed accuracy about 0.05" means ε_K; the per-term tolerance is about 0.006.
- **c_σ** = 12/(2·13⁴) = 2.10e−4: Fang (5.3) at λ = 12, via audit-b step 3.
- **Raič constant** 58 (d = 1) or 65.95 (d = 2): arXiv:1802.06475v4 eq. (1.2), via the audits. Switching changes log10 N1 by 0.1.
- **K_M** = C_δ/(1 − 2^{a−1}), with C_δ = 12·sup_y e^{−y}(1 + C0·y)^{2/p}, C0 = 1/(1−ρ), a = 2 − 2/p. Source: audit-b steps 1–2, on Fang Lemma 4.2 and Prop. 4.1. I recomputed it in closed form (y* = 2/p − 1/C0):
  - p = 19/10, ρ ≤ 199/200: C_δ = 1174.3, K_M = 32,780 (matches audit-b's C_up = 3.3e4).
  - Paper as written: K_M = 3.76e7 (matches audit-b's C_up = 3.8e7).
- **K_S** comes from audit-b's C_γ and is subdominant everywhere. For a ≤ 1/2 its formula breaks, so I set a = 0.55 in that term only; this is labelled in the script.
- **The chain.** This is Fang (5.2)–(5.7) as audit-b wrote it out, transferred to M′ with normaliser Var M′ ≥ η·c_σ·n. Its terms are:
  - Berry–Esseen: Raič·b·√(2/(η c_σ n))
  - mean shift: √(K_M b^{a−1}/(π η c_σ))
  - variance shift: (2 + 1/√(2πe))·(K_M b^{a−1} + K_S n^{a−1})/(η c_σ)

**Results: log10 N1, fixed-b route / paper route (b = n^{1/4}).**

| Scenario | 1−a | Hybrid, ρ0 = 0.4448 | Hybrid, ρ0 provable = 1/(104 C_up) | Fang N0 (audits) |
|---|---|---|---|---|
| a = 2/3, C = 1 (like-for-like floor) | 0.333 | **66** / 107 | 116 | **240** / 537 (audit-b) |
| a = 0, C = 1 | 1 | 29.5 / 34 | 48 | 103 (audit-a; its c_σ may differ) |
| certified sharpened, p = 19/10 | 0.0526 | **495** / 965 | 1,218 / 2,375 | **2,070** / 5,375 |
| p = 193/100 | 0.0363 | 686 / 1,346 | 1,643 / 3,227 | 2,828 / 7,430 |
| as written (PDF best) | 0.0015 | 21,062 / 42,099 | n/a | 2.44e5 (paper route) |
| a = 0.88, C = 1 (audit-a slope at λ = 12) | 0.12 | 167 / 310 | 314 | n/a |
| a = 0.80, C = 1 (audit-b slope) | 0.20 | 103 / 182 | 188 | n/a |
| λ_max = 4 (unproved), a = 0.69, C = 1 | 0.31 | 61 / 100 | n/a | n/a |

**Dominant term.** In every scenario the L² mean-shift term sets b, and Berry–Esseen then sets n ≈ b². In formula form:
- log10 N1 ≈ (2/(1−a))·F + 2·log10(Raič/t) + log10(2/(η c_σ)), where F = log10(K_M/(π η c_σ t²)).
- For the certified case F = 12.70 digits. Of that, K_M contributes 4.52, 1/c_σ contributes 3.68, 1/(πη) contributes 0.03 and 1/t² contributes 4.48. Multiplying by 38 gives 483, and the remaining terms add 13, for 495.

**Answer to the key question.**
- **Yes, the R^-3 floor goes.** The hybrid never uses the characteristic function, so the tail cut R, ε1 ~ R^-3 and the 10^−40 to 10^−263 tolerances disappear.
  - Like for like (a = 2/3, C = 1, same c_σ): Fang's floor is 10^240 and the hybrid's is 10^66.
  - At the certified constants: 10^2070 falls to 10^495.
  - Roughly 41 of the ~54 digits per 1/(1−a) factor vanish.
- **But the exponent a brings the floor back.** The mean shift is still paid at exponent 2/(1−a) on the remaining ~12.7 digits.
  - If the root-moment exponent really is ≥ 0.8 (audit-b's slope 0.82) to 0.88 (audit-a, λ = 12), then N1 ≥ 10^103 to 10^167 even with unit constants. Those slopes are float heuristics that only bound a from below.
  - With the constants actually certified, N1 ≈ 10^495.
  - Removing the 1/t² cost entirely (a CLT that treats M − μ as Gaussian rather than as an L² error) would still leave about 10^80 to 10^90 at a = 0.88. That figure is hand arithmetic read off the decomposition ((12.70 − 4.48 → 4.80) × 16.7), not a script output.
- **Using only a ρ0 the audits can prove** (≥ 1/(104 C_up) = 2.9e−7), the tolerance collapses to about 10^−11 and log10 N1 = 1,218.
- **The only scenario near enumeration** is route-lclt's hypothetical: if (K) held as d_K ≤ K·n^{−1/2}, then N1 = (K/ε_K)² = 475K² at η = 0.3, or 1,900K² with the half budget. That is exactly the unproved content.
- **So the best explicit N1** is 10^495 if (V) is granted empirically, or 10^1218 using only what the audits can prove. Both are also conditional on (T) and on the two transfer assumptions in section 3. Neither is a theorem.

## 3. What remains unproved in Lemma L

**Setup.**
- T is a tree on n ≥ N_* vertices, k ∈ [⌈n/4⌉, min(q, α−1)], and λ ∈ [1/3, 12] satisfies E_λ X = k.
- L is the side maximising E|I∩L|, J = I∩R, B0 = ⌈(1−η)·E B_J⌉, E0 = {B_J ≥ B0}, and M′ = |J| + Bin(B_J − B0, p) on E0.
- Everything needs explicit constants.

**The three conditions.**
- **(V)** p(1−p)·E B_J ≥ ρ0·σ², equivalently Var E[X|J] ≤ (1−ρ0)σ².
  - The empirical minimum is 0.4448.
  - The provable bound is only ρ0 ≥ 1/(104 C_up): E|I∩L| ≥ k/2 ≥ n/8, λ ≤ 12, and σ² ≤ C_up·n.
- **(K)** d_K(M′|E0, G_d) ≤ θ·c_eff·((1−η)ρ0)^{3/2}, where G_d is the discretised normal with the conditional mean and variance of M′.
  - Any explicit o(1) rate would do. Route-lclt's K·n^{−1/2} is stronger than needed.
- **(T)** a_{k−1} + a_{k+1} ≤ θ″·P(E0)·E h_{B0}(k − G_d), with a_j = P(X = j, B_J < B0).
  - Route-lclt's P(B_J < B0) ≤ K_T·n^{−2} suffices but is stronger than necessary.
  - My derivation (untested; whether the constants work is unchecked):
    - Bound a_j ≤ P(B_J < B0/2) + P(B_J < B0)·max_{B ≥ B0/2} max_i b_{B,p}(i).
    - Then (T) reduces to a deep tail P(B_J < B0/2) = o(n^{−3/2}), for which a fourth-moment bound E(B_J − E B_J)⁴ = O(n²) suffices.
    - It also needs P(B_J < B0) ≤ κ/n with κ small enough: a constant-level inequality at the n^{−3/2} scale.
    - Chebyshev with Var B_J = O(n) gives O(1/n), so the requirement is borderline, not hopeless.

**Assumptions used in section 2 that are part of (K).**
- (i) *Transfer of Fang §5 to M′.* Given the centroid reveal, M′ is not a sum of independent summands. A linearisation through |J| + p·B_J is plausible, but revealing L-centroids disturbs the product-Bernoulli structure of K given J. I priced this at Raič d = 2 as a hypothesis; it changes N1 by 0.1 in log10.
- (ii) *Conditioning on the global event E0.* E0 breaks independence between components. The workaround M″ = |J| + Bin((B_J − B0)⁺, p) costs P(E0^c) in Kolmogorov distance.

**Elementary but unwritten.**
- A uniform explicit C_TV. It is at least 1.579 (exact); 1.580 is the float supremum over B ≤ 1200.
- A lower bound on E h_{B0}(k − G_d) with an explicit discretisation (Euler–Maclaurin) error.
- Control of the offset between E[X|E0] and k.

**No proof mechanism for (V), (K) or (T) was identified.** As section 2 shows, a Fang-type centroid chain cannot deliver (K) at a threshold anywhere near the enumeration record (n ≤ 32).

## 4. A candidate for Aristotle

**Packet: the exact, finite form of Lemma C.** Low mathematical value, low risk, cheap. It would give a formally verified reduction but does not advance Lemma L.

**Setting.**
- G is a finite simple graph with V = L ⊔ R and both L and R independent.
- A_J = {v ∈ L : ∀w ∈ J, ¬Adj v w}, and B_J = |A_J|.
- λ > 0, p = λ/(1+λ), Z = Σ_{I indep} λ^{|I|}, and p_j = #{I indep : |I| = j}·λ^j/Z.
- B0 ≥ 1, b_B(i) = C(B,i)p^i(1−p)^{B−i}, and h(i) = 2b_{B0}(i) − b_{B0}(i−1) − b_{B0}(i+1).
- w_J = λ^{|J|}(1+λ)^{B_J}/Z.
- μ′(m) = Σ_{B_J ≥ B0} w_J·b_{B_J−B0}(m − |J|), a sub-probability law of total mass P0.
- a_j = Σ_{B_J < B0} w_J·b_{B_J}(j − |J|).

**Top target.** For every k ∈ ℤ and every finitely supported g ≥ 0 on ℤ with Σg = P0:
- 2p_k − p_{k−1} − p_{k+1} ≥ Σ_m g(m)h(k−m) − TV(h)·sup_m |Σ_{i≤m}(μ′(i) − g(i))| − (a_{k−1} + a_{k+1}), where TV(h) = Σ_i |h(i+1) − h(i)|.

**Graded sub-targets.**
1. A0 in `Polynomial ℕ`, for any finite bipartite graph: Σ_{I indep} X^{|I|} = Σ_{J ⊆ R} X^{|J|}(1+X)^{B_J}.
2. The Abel summation inequality as a standalone lemma on finitely supported functions ℤ → ℝ.

Refutation is not expected; the statement is exact-checked on all trees with n ≤ 14 and 7 larger ones.

**Do not send** TV(h_B) ≤ ∫|φ‴|/v^{3/2}: it is exactly refuted. A uniform "TV·v^{3/2} ≤ 8/5" version is definition-complete but not settled, and it needs Stirling or Edgeworth analysis. It is not a bounded packet.

## Not done
- The hub variant.
- Any proof attempt on Lemma L.
- A re-audit of Fang §4–5. I recomputed C_δ from audit-b's closed form and imported C_γ.
- Beyond B ≤ 1200, the TV supremum is float-only.
- The KT5 n = 16 rerun reuses route-lclt code.

## Files
- **New scripts:** `kernel_check_exact.py`, `tv_exact_point.py`, `tv_sup_scan.py`, `n1_check.py`, `census_range.py`, `exact_large2.py`, `make_manifest2.py`
- **Re-verified first-attempt files:** `exact_AC.py`, `exact_large.py`, `kt5_sharp.py`, `n1_hybrid.py`
- **Outputs:** `rerun/` (census_n3_12, census_n13, census_n14, exact_large2, kernel_check_exact, tv_exact_point, tv_sup_scan, n1_check, kt5_sharp_n16)
- **Manifest:** `manifest.yaml`
