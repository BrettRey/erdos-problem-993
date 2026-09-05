## Claim

For a finite forest with independence number α ≥ 5 in which every independent set of size α−2 extends to a maximum independent set (b₂ = 0), the non‑extendable counts satisfy 5b₃ ≤ e₃ and 5b₄ ≤ e₄.

## Findings (by severity)

**1. No error in the main line.** I re‑derived, symbolically and independently, every load‑bearing identity: (5)–(7) (I_{T⁺}, E_{T⁺}, B_{T⁺}, reversed profiles), (9), (10), (11). All are correct as displayed, including b_{T,3} = b_{G,3}+f₀ and b_{T,4} = b_{G,4}+3b_{G,3}+f₁ in the three‑leaf case. Lemma 3 (blocked (α−1)-set ⟹ 3(α−1) ≤ n ≤ 2α), Lemma 4 (only one blocked factor through degree 4; the images (4) and (3,1) are distinct terms of 𝓔_F), Lemma 5, the α_i ≤ 4 classification (blocked singleton ⟹ universal vertex ⟹ K_{1,4}), the m = 2 impossibility (ι₂ = [x^{a−1}]I_{G−w} ≥ 1), and the m ≥ 4 edge‑addition argument (maximum sets, hence all e_d, unchanged; b_d can only drop) all check out.

**2. Repair needed — Lemma 5, the a = 1 line.** "e₂ = 0, e₁ = 1, e₀ ≥ 1 … also suffice" does not give 4e₂+6e₁−e₀ = 6−e₀ ≥ 0. It needs e₀ ≤ 6, which follows from the unconditional incidence bound e₀ ≤ 2e₁ proved two lines earlier (here e₀ ≤ 2). Insert that citation; no consequence downstream.

**3. Repair needed — (11) is stated without its derivation.** The step uses ι₃ = f₀+b₀(F) and ι₄ = f₁+b₁(F), i.e. it consumes b₀(F) = b₁(F) = 0. §11 proves b₁(F) = 0 immediately before but never says where it is used; without it the displayed margins are wrong. One sentence fixes this.

**4. Vacuous branch (harmless).** In the m = 1 case, "a = 4 and R = K_{1,4}" cannot occur: then G is K_{1,4} with a pendant w at a leaf, and {centre, w} is a maximal independent 2‑set, so b₂(G) = 1 ≠ 0. The branch's inequality (5h₃ = 5 ≤ e_{G,2}+q₃, e_{G,2} ≥ 6) is nevertheless true; it can simply be deleted.

**5. Hypotheses to state explicitly.** §11's support lemma needs n ≥ 3 plus "a tree with one non‑leaf is a star", and the root must be a non‑leaf so that w is a non‑leaf. §10's "non‑star components" must be read as "components ≠ K_{1,4}". The m ≥ 4 reduction silently uses b₂(K_{1,m}) = 0, true precisely because m ≥ 4.

**Independent corroboration at the tight end.** For α = 5, e₄+b₄ = n ≤ 2α = 10, so the theorem asserts b₄ ≤ 1. I verified this separately: under b₂ = 0 a blocked vertex u satisfies α(T−N[u]) ≤ 1, so deg u ≥ α−1; two blocked vertices must be non‑adjacent with {u,v} maximal, and the tree‑edge count forces exactly one common neighbour, giving α = deg u + deg v − 1 ≥ 7. Also checked: E_T = (1+2x)(1+x)⁴, B_T = x+x² for the K_{1,4}-centre extension ((e₃,e₄,b₃,b₄) = (14,6,1,1)); the nine‑vertex m = 3 case (35,35,2,6), margins (25,5); the two‑star base (56,70,2,8); and Lemma 5 on K_{1,4}⊔K₁⊔K_{1,4} → (84,126,2,10), margins (74,76), matching (9). Outside results: bipartite α ≥ n/2 (checked directly, not recalled); Woodroofe and Adiprasito–Björner–Goodarzi are only recalled, and §12 does not use them. I confirmed §5's split graph G_m (α = 6, b₂ = 0, e₄ = 15, b₄ = m), so forest structure (n ≤ 2α, N(u) independent) is genuinely load‑bearing and is in fact used.

## Verdict

**VALID WITH EXPLICIT REPAIR** (items 2 and 3; each is one line, using facts already proved in the artifact).

Decisive reasoning: the induction is on vertex count; the split m = 1/2/3/≥4 on a deepest support vertex is exhaustive for non‑stars; m = 2 is killed by an exact coefficient; m ≥ 4 reduces to a disjoint union with identical maximum sets; m = 3 and m = 1 reduce through the exact margin identities, whose added terms are non‑negative by the two incidence bounds; and the sole base violating 5b₃ ≤ e₃, namely K_{1,4}, occurs in exactly three places, each settled by explicit computation (margins (9,1), (25,5), and Lemma 5). No step applies the theorem to a tree of the same order, so the dependency audit is sound.
